#ifdef HEADEARPATH
#include "interface/hdsdp_schur.h"
#include "interface/hdsdp_utils.h"
#include "interface/hdsdp_conic.h"
#else
#include "hdsdp_schur.h"
#include "hdsdp_utils.h"
#include "hdsdp_conic.h"
#endif

static hdsdp_retcode HKKTIAllocDenseKKT( hdsdp_kkt *HKKT ) {
    
    hdsdp_retcode retcode = HDSDP_RETCODE_OK;
    
    /* This is one of the most memory intensive part of HDSDP */
    HDSDP_INIT(HKKT->kktMatElem, double, HKKT->nRow * HKKT->nRow);
    HDSDP_MEMCHECK(HKKT->kktMatElem);
    
    HDSDP_CALL(HFpLinsysCreate(&HKKT->kktM, HKKT->nRow, HDSDP_LINSYS_DENSE_ITERATIVE));
    
    double dCGAcc = KKT_ACCURACY;
    int nCGIteration = -1;
    
    if ( HKKT->nRow > 20000 ) {
        dCGAcc = KKT_ACCURACY * 100.0;
        nCGIteration = 500;
    } else if ( HKKT->nRow > 15000 ) {
        dCGAcc = KKT_ACCURACY * 50.0;
        nCGIteration = 450;
    } else if ( HKKT->nRow > 5000 ) {
        dCGAcc = KKT_ACCURACY * 5.0;
        nCGIteration = 120;
    }
    
    HFpLinsysSetParam(HKKT->kktM, 5.0 * dCGAcc, dCGAcc, -1, nCGIteration, -1);
    
    /* Mark the diagonal entries of the KKT solver */
    for ( int iCol = 0; iCol < HKKT->nRow; ++iCol ) {
        HKKT->kktDiag[iCol] = &HKKT->kktMatElem[iCol + iCol * HKKT->nRow];
    }
    
exit_cleanup:
    return retcode;
}

static hdsdp_retcode HKKTAllocSparseKKT( hdsdp_kkt *HKKT ) {
    
    /* Till here the aggregated sparsity pattern has not been tested
       Hence HDSDP tries to build up a sparse Schur complement and will switch back to
       dense Schur complement if aggregation of sparsity pattern shows that the KKT system becomes dense
       after accumulation.
     
       To improve efficieny, we only loop once and use realloc to extend memory when #nnz runs out
     */
    
    hdsdp_retcode retcode = HDSDP_RETCODE_OK;
    
    int nAllocatedNnz = HKKT->nRow * 3;
    int iNz = 0;
    int64_t nDenseKKTNnz = HDSDP_SPARSE_SCHUR_THRESHOLD * HKKT->nRow * HKKT->nRow;
    
    int *colBuffer = NULL;
    
    /* Allocate initial memory */
    HDSDP_INIT(HKKT->kktMatBeg, int, HKKT->nRow + 1);
    HDSDP_MEMCHECK(HKKT->kktMatBeg);
    
    HDSDP_INIT(HKKT->kktMatIdx, int, nAllocatedNnz);
    HDSDP_MEMCHECK(HKKT->kktMatIdx);
    
    /* Allocate buffer to collect column sparsity */
    HDSDP_INIT(colBuffer, int, HKKT->nRow);
    
    for ( int iCol = 0; iCol < HKKT->nRow; ++iCol ) {
        /* Zero out the buffer */
        HDSDP_ZERO(colBuffer, int, HKKT->nRow);
        for ( int iCone = 0; iCone < HKKT->nCones; ++iCone ) {
            HConeAddSymNz(HKKT->cones[iCone], iCol, colBuffer);
        }
        /* Collect the sparsity pattern in the lower triangular part of the buffer */
        for ( int iRow = iCol; iRow < HKKT->nRow; ++iRow ) {
            if ( colBuffer[iRow] ) {
                colBuffer[iRow] = iNz;
                HKKT->kktMatIdx[iNz] = iRow;
                iNz += 1;
            }
        }
        
        /* Retrieve sparsity mapping */
        for ( int iCone = 0; iCone < HKKT->nCones; ++iCone ) {
            HConeGetSymMapping(HKKT->cones[iCone], iCol, colBuffer);
        }
        
        /* Update index */
        HKKT->kktMatBeg[iCol + 1] = iNz;
        
        /* If there is risk of running out of memory in the next round, we do re-allocation */
        if ( nAllocatedNnz - iNz <= HKKT->nRow - iCol ) {
            nAllocatedNnz += HKKT->nRow * 2;
            HDSDP_REALLOC(HKKT->kktMatIdx, int, nAllocatedNnz);
            HDSDP_MEMCHECK(HKKT->kktMatIdx);
        }
        
        /* The aggregated sparsity pattern becomes dense. We switch back to dense matrix */
        if ( iNz >= nDenseKKTNnz ) {
            HKKT->isKKTSparse = 0;
            break;
        }
    }
    
    if ( HKKT->isKKTSparse ) {
        /* The Schur complement is sparse. And we use sparse direct solver */
        HDSDP_REALLOC(HKKT->kktMatIdx, int, iNz);
        HDSDP_INIT(HKKT->kktMatElem, double, iNz);
        
        /* Mark the diagonal elements. We assert that the start of each column must be its
           diagonal. Otherwise the column is empty and problem data is invalid.
         */
        for ( int iCol = 0; iCol < HKKT->nRow; ++iCol ) {
            if ( HKKT->kktMatIdx[HKKT->kktMatBeg[iCol]] != iCol ) {
                printf("KKT solver detects an empty column.\n");
                retcode = HDSDP_RETCODE_FAILED;
                goto exit_cleanup;
            }
            HKKT->kktDiag[iCol] = &HKKT->kktMatElem[HKKT->kktMatBeg[iCol]];
        }
        
        HDSDP_CALL(HFpLinsysCreate(&HKKT->kktM, HKKT->nRow, HDSDP_LINSYS_SPARSE_DIRECT));
        HDSDP_CALL(HFpLinsysSymbolic(HKKT->kktM, HKKT->kktMatBeg, HKKT->kktMatIdx))
    } else {
        HDSDP_FREE(HKKT->kktMatBeg);
        HDSDP_FREE(HKKT->kktMatIdx);
        HDSDP_CALL(HKKTIAllocDenseKKT(HKKT));
    }
    
exit_cleanup:
    HDSDP_FREE(colBuffer);
    return retcode;
}

static void HKKTClean( hdsdp_kkt *HKKT, int typeKKT ) {
    
    /* Clean up for the next KKT solve */
    HDSDP_ZERO(HKKT->dASinvVec, double, HKKT->nRow);
    HDSDP_ZERO(HKKT->dASinvRdSinvVec, double, HKKT->nRow);
    
    if ( typeKKT == KKT_TYPE_HOMOGENEOUS ) {
        HDSDP_ZERO(HKKT->dASinvCSinvVec, double, HKKT->nRow);
        HKKT->dCSinv = 0.0;
        HKKT->dCSinvCSinv = 0.0;
        HKKT->dCSinvRdSinv = 0.0;
    }
    
    HKKT->dTraceSinv = 0.0;
    
    if ( typeKKT == KKT_TYPE_INFEASIBLE || typeKKT == KKT_TYPE_HOMOGENEOUS || typeKKT == KKT_TYPE_PRIMAL ) {
        if ( HKKT->isKKTSparse ) {
            HDSDP_ZERO(HKKT->kktMatElem, double, HKKT->kktMatBeg[HKKT->nRow]);
        } else {
            HDSDP_ZERO(HKKT->kktMatElem, double, HKKT->nRow * HKKT->nRow);
        }
    }

    return;
}

extern hdsdp_retcode HKKTCreate( hdsdp_kkt **pHKKT ) {
    
    hdsdp_retcode retcode = HDSDP_RETCODE_OK;
    HDSDP_NULLCHECK(pHKKT);
    hdsdp_kkt *HKKT = NULL;
    HDSDP_INIT(HKKT, hdsdp_kkt, 1);
    HDSDP_MEMCHECK(HKKT);
    HDSDP_ZERO(HKKT, hdsdp_kkt, 1);
    *pHKKT = HKKT;
    
exit_cleanup:
    return retcode;
}

extern hdsdp_retcode HKKTInit( hdsdp_kkt *HKKT, int nRow, int nCones, hdsdp_cone **cones ) {
    
    hdsdp_retcode retcode = HDSDP_RETCODE_OK;
    
    HKKT->nRow = nRow;
    HKKT->nCones = nCones;
    HKKT->cones = cones;
    
    int coneDim = 0;
    int maxConeDim = 0;
    
    /* Get the maximum conic dimension */
    for ( int iCone = 0; iCone < nCones; ++iCone ) {
        coneDim = HConeGetDim(cones[iCone]);
        maxConeDim = HDSDP_MAX(coneDim, maxConeDim);
    }
    
    HKKT->maxConeDim = maxConeDim;
    
    /* Prepare the buffer to export S^-1 */
    HDSDP_INIT(HKKT->invBuffer, double, maxConeDim * maxConeDim);
    HDSDP_MEMCHECK(HKKT->invBuffer);
    
    /* Prepare the buffer to store internal computation results */
    HDSDP_INIT(HKKT->kktBuffer, double, maxConeDim * maxConeDim);
    HDSDP_MEMCHECK(HKKT->kktBuffer);
    HDSDP_INIT(HKKT->kktBuffer2, double, maxConeDim * maxConeDim);
    HDSDP_MEMCHECK(HKKT->kktBuffer2);
    
    /* Allocate vector space */
    HDSDP_INIT(HKKT->dASinvVec, double, nRow);
    HDSDP_MEMCHECK(HKKT->dASinvVec);
    
    HDSDP_INIT(HKKT->dASinvCSinvVec, double, nRow);
    HDSDP_MEMCHECK(HKKT->dASinvCSinvVec);
    
    HDSDP_INIT(HKKT->dASinvRdSinvVec, double, nRow);
    HDSDP_MEMCHECK(HKKT->dASinvRdSinvVec);
    
    HDSDP_INIT(HKKT->kktDiag, double *, nRow);
    HDSDP_MEMCHECK(HKKT->kktDiag);
    
    /* Choose the sparsity pattern of the Schur complement */
    /* Initial test by investigating the separate sparsity patterns */
    int64_t coneNnzs = 0;
    int64_t maxConicNnzs = 0;
    int64_t nDenseKKTNnz = HDSDP_SPARSE_SCHUR_THRESHOLD * nRow * nRow;
    
    HKKT->isKKTSparse = 1;
    
    for ( int iCone = 0; iCone < nCones; ++iCone ) {
        coneNnzs = HConeGetSymNnz(cones[iCone]);
        maxConicNnzs = HDSDP_MAX(coneNnzs, maxConicNnzs);
        if ( maxConicNnzs >= nDenseKKTNnz ) {
            HKKT->isKKTSparse = 0;
            break;
        }
    }
    
    if ( !HKKT->isKKTSparse ) {
        /* Use dense Schur complement */
        HDSDP_CALL(HKKTIAllocDenseKKT(HKKT));
        printf("    Using dense Schur complement\n");
    } else {
        /* Try using sparse Schur complement */
        HDSDP_CALL(HKKTAllocSparseKKT(HKKT));
        printf("    Using sparse Schur complement (%d nnzs)\n", HKKT->kktMatBeg[HKKT->nRow]);
    }
    
    HKKT->dPrimalX = NULL;
    
exit_cleanup:
    return retcode;
}

extern hdsdp_retcode HKKTBuildUp( hdsdp_kkt *HKKT, int typeKKT ) {
    
    hdsdp_retcode retcode = HDSDP_RETCODE_OK;
    
    HKKTClean(HKKT, typeKKT);
    
    for ( int iCone = 0; iCone < HKKT->nCones; ++iCone ) {
        HDSDP_CALL(HConeBuildSchurComplement(HKKT->cones[iCone], HKKT, typeKKT));
    }
    
exit_cleanup:
    return retcode;
}

extern hdsdp_retcode HKKTBuildUpExtraCone( hdsdp_kkt *HKKT, hdsdp_cone *cone, int typeKKT ) {
    
    hdsdp_retcode retcode = HDSDP_RETCODE_OK;
    HDSDP_CALL(HConeBuildSchurComplement(cone, HKKT, typeKKT));
    
exit_cleanup:
    return retcode;
}

extern hdsdp_retcode HKKTBuildUpFixed( hdsdp_kkt *HKKT, int typeKKT, int kktStrategy ) {
    
    hdsdp_retcode retcode = HDSDP_RETCODE_OK;
    
    HKKTClean(HKKT, typeKKT);
    
    for ( int iCone = 0; iCone < HKKT->nCones; ++iCone ) {
        HDSDP_CALL(HConeBuildSchurComplementFixed(HKKT->cones[iCone], HKKT, typeKKT, kktStrategy));
    }
    
exit_cleanup:
    return retcode;
}

extern void HKKTExport( hdsdp_kkt *HKKT, double *dKKTASinvVec, double *dKKTASinvRdSinvVec, double *dKKTASinvCSinvVec,
                       double *dCSinvCSinv, double *dCSinv, double *dCSinvRdCSinv, double *dTraceSinv ) {
    
    /* Export the vectors and information from KKT solver */
    if ( dKKTASinvVec ) {
        HDSDP_MEMCPY(dKKTASinvVec, HKKT->dASinvVec, double, HKKT->nRow);
    }
    
    if ( dKKTASinvRdSinvVec ) {
        HDSDP_MEMCPY(dKKTASinvRdSinvVec, HKKT->dASinvRdSinvVec, double, HKKT->nRow);
    }
    
    if ( dKKTASinvCSinvVec ) {
        HDSDP_MEMCPY(dKKTASinvCSinvVec, HKKT->dASinvCSinvVec, double, HKKT->nRow);
    }
    
    if ( dCSinvCSinv ) {
        *dCSinvCSinv = HKKT->dCSinvCSinv;
    }
    
    if ( dCSinv ) {
        *dCSinv = HKKT->dCSinv;
    }
    
    if ( dCSinvRdCSinv ) {
        *dCSinvRdCSinv = HKKT->dCSinvRdSinv;
    }
    
    if ( dTraceSinv ) {
        *dTraceSinv = HKKT->dTraceSinv;
    }

    return;
}

extern hdsdp_retcode HKKTFactorize( hdsdp_kkt *HKKT ) {
    
    hdsdp_retcode retcode = HDSDP_RETCODE_OK;
    
    HDSDP_CALL(HFpLinsysNumeric(HKKT->kktM, HKKT->kktMatBeg, HKKT->kktMatIdx, HKKT->kktMatElem));
    
exit_cleanup:
    return retcode;
}

extern hdsdp_retcode HKKTSolve( hdsdp_kkt *HKKT, double *dRhsVec, double *dLhsVec ) {
    
    /* Apply KKT solver */
    hdsdp_retcode retcode = HDSDP_RETCODE_OK;
    HDSDP_CALL(HFpLinsysSolve(HKKT->kktM, 1, dRhsVec, dLhsVec));
    
exit_cleanup:
    return retcode;
}

extern void HKKTRegularize( hdsdp_kkt *HKKT, double dKKTReg ) {
    
    /* Regularize the diagonal of the Schur complement */
    double dKKTMinDiag = HDSDP_INFINITY;
    
    for ( int iCol = 0; iCol < HKKT->nRow; ++iCol ) {
        dKKTMinDiag = HDSDP_MIN(*HKKT->kktDiag[iCol], dKKTMinDiag);
    }
    
    dKKTReg = dKKTReg * dKKTMinDiag;
    dKKTReg = HDSDP_MIN(dKKTReg, 1e-05);
    
    if ( dKKTReg < 1e-14 ) {
        dKKTReg = 0.0;
    }
    
#ifdef HDSDP_KKT_DEBUG
    hdsdp_printf("Regularizing KKT system with %5.2e\n", dKKTReg);
#endif
    
    for ( int iCol = 0; iCol < HKKT->nRow; ++iCol ) {
        *HKKT->kktDiag[iCol] += dKKTReg;
    }
    
    return;
}

extern void HKKTRegisterPSDP( hdsdp_kkt *HKKT, double **dPrimalX ) {
    
    HKKT->dPrimalX = dPrimalX;
    
    return;
}

extern void HKKTClear( hdsdp_kkt *HKKT ) {
    
    if ( !HKKT ) {
        return;
    }
    
    HDSDP_FREE(HKKT->dASinvVec);
    HDSDP_FREE(HKKT->dASinvCSinvVec);
    HDSDP_FREE(HKKT->dASinvRdSinvVec);
    HDSDP_FREE(HKKT->invBuffer);
    HDSDP_FREE(HKKT->kktBuffer);
    HDSDP_FREE(HKKT->kktBuffer2);
    
    HDSDP_FREE(HKKT->kktMatBeg);
    HDSDP_FREE(HKKT->kktMatIdx);
    HDSDP_FREE(HKKT->kktMatElem);
    
    HDSDP_FREE(HKKT->kktDiag);
    
    HFpLinsysDestroy(&HKKT->kktM);
    HDSDP_ZERO(HKKT, hdsdp_kkt, 1);
 
    return;
}

extern void HKKTDestroy( hdsdp_kkt **pHKKT ) {
    
    if ( !pHKKT ) {
        return;
    }
    
    HKKTClear(*pHKKT);
    HDSDP_FREE(*pHKKT);
    
    return;
}
