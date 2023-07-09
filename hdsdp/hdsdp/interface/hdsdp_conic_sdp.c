#ifdef HEADERPATH
#include "interface/hdsdp_conic_sdp.h"
#include "interface/def_hdsdp_user_data.h"
#include "interface/hdsdp_utils.h"
#include "linalg/hdsdp_sdpdata.h"
#include "linalg/sparse_opts.h"
#include "linalg/hdsdp_linsolver.h"
#include "external/hdsdp_cs.h"
#else
#include "hdsdp_conic_sdp.h"
#include "def_hdsdp_user_data.h"
#include "hdsdp_utils.h"
#include "hdsdp_sdpdata.h"
#include "sparse_opts.h"
#include "hdsdp_linsolver.h"
#include "hdsdp_cs.h"
#endif

#include <math.h>

#ifndef SMALL_DUAL_THRESHOLD
#define SMALL_DUAL_THRESHOLD (0)
#endif

#ifndef SPARSE_DUAL_THRESHOLD
#define SPARSE_DUAL_THRESHOLD (0.4)
#endif
static hdsdp_retcode sdpDenseConeIAllocDualMat( hdsdp_cone_sdp_dense *cone ) {
    /* This routine allocates the dual information for a dense SDP cone
       According to the aggregated sparsity pattern of the matrices in the cone,
       the SDP dual matrix will be either represented using sparse csc or dense packed format.
     
       If sparse data structure is employed, then dualMatBeg, dualMatIdx will be filled by the sparsity pattern
       of the dual matrix; dualPosToElemMap will be an n * (n + 1) / 2 integer array that maps an element
       in the lower triangular part to the position of dualMatElem. dualMatElem and dualCheckeerElem are
       two concurrent representations of the dual matrix and will be factorized.
       dualStep only serves as a buffer to store the SDP step matrix dS and will not be factorized.
     
       If dense data structure is employed. dualMatBeg, dualMatIdx and dualPosToElemMap will not be allocated
       and dualMatElem, dualCheckerElem and dualStep will adopt dense packed format.
     
       To build the conic structure, we first check if there is any dense coefficient matrix in the cone.
       Whenever there is a dense coefficient, the dual matrix will use dense data structure.
     
       Otherwise, we will allocate dualPosToElemMap and iterate through the coefficient matrices for their
       sparsity pattern, column-wise. Once the number of nonzeros exceeds a pre-defined threshold, the dual
       matrix will use a dense structure.

     */
    hdsdp_retcode retcode = HDSDP_RETCODE_OK;
    
    /* Initialize dual matrix as sparse */
    cone->isDualSparse = 1;
    int nDualNonzeros = 0;
    
    /* If there is a dense matrix, the dual matrix will be dense */
    if ( cone->sdpConeStats[SDP_COEFF_DENSE] > 0 ||
         cone->sdpConeStats[SDP_COEFF_DSR1] > 0  ||
         cone->nCol < SMALL_DUAL_THRESHOLD ) {
        cone->isDualSparse = 0;
    }
    
    /* Allocate the dual sparsity pattern */
    if ( cone->isDualSparse ) {
        HDSDP_INIT(cone->dualPosToElemMap, int, PACK_NNZ(cone->nCol));
        HDSDP_MEMCHECK(cone->dualPosToElemMap);
        
        /* Accumulate the sparsity pattern */
        for ( int iRow = 0; iRow < cone->nRow; ++iRow ) {
            sdpDataMatGetMatNz(cone->sdpRow[iRow], cone->dualPosToElemMap);
        }
        sdpDataMatGetMatNz(cone->sdpObj, cone->dualPosToElemMap);
        
        /* Count statistic of nonzeros and determine matrix type */
        for ( int iElem = 0; iElem < PACK_NNZ(cone->nCol); ++iElem ) {
            nDualNonzeros += cone->dualPosToElemMap[iElem];
        }
        
        if ( nDualNonzeros >= SPARSE_DUAL_THRESHOLD * cone->nCol * cone->nCol ) {
            /* We decide to use dense representation */
            cone->isDualSparse = 0;
        }
    }
    
    if ( !cone->isDualSparse ) {
        /* We use dense representation of the dual matrix. No*/
        HDSDP_FREE(cone->dualPosToElemMap);
        HDSDP_INIT(cone->dualMatElem, double, cone->nCol * cone->nCol);
        HDSDP_MEMCHECK(cone->dualMatElem);
        HDSDP_INIT(cone->dualCheckerElem, double, cone->nCol * cone->nCol);
        HDSDP_MEMCHECK(cone->dualCheckerElem);
        HDSDP_INIT(cone->dualStep, double, cone->nCol * cone->nCol);
        HDSDP_MEMCHECK(cone->dualStep);
        HDSDP_CALL(HFpLinsysCreate(&cone->dualFactor, cone->nCol, HDSDP_LINSYS_DENSE_DIRECT));
        /* Done */
    } else {
        /* We are using the sparse dual representation */
        assert( nDualNonzeros > 0 );
        HDSDP_INIT(cone->dualMatBeg, int, cone->nCol + 1);
        HDSDP_MEMCHECK(cone->dualMatBeg);
        HDSDP_INIT(cone->dualMatIdx, int, nDualNonzeros);
        HDSDP_MEMCHECK(cone->dualMatIdx);
        HDSDP_INIT(cone->dualMatElem, double, nDualNonzeros);
        HDSDP_MEMCHECK(cone->dualMatElem);
        HDSDP_INIT(cone->dualCheckerElem, double, nDualNonzeros);
        HDSDP_MEMCHECK(cone->dualCheckerElem);
        HDSDP_INIT(cone->dualStep, double, nDualNonzeros);
        HDSDP_MEMCHECK(cone->dualStep);
        
        /* Now start enumerating over the mapping and construct cone data structure */
        int iNz = 0;
        int *colPtr = cone->dualPosToElemMap;
        for ( int iCol = 0; iCol < cone->nCol; ++iCol ) {
            for ( int iRow = 0; iRow < cone->nCol - iCol; ++iRow ) {
                if ( colPtr[iRow] ) {
                    colPtr[iRow] = iNz;
                    cone->dualMatIdx[iNz] = iRow + iCol;
                    iNz += 1;
                }
            }
            cone->dualMatBeg[iCol + 1] = iNz;
            colPtr += (cone->nCol - iCol);
        }
        
        /* Symbolic factorization is left in the presolve routine */
        HDSDP_CALL(HFpLinsysCreate(&cone->dualFactor, cone->nCol, HDSDP_LINSYS_SPARSE_DIRECT));
        HFpLinsysSetParam(cone->dualFactor, -1.0, -1.0, 4, -1, -1);
        /* Done */
    }
    
exit_cleanup:
    return retcode;
}

static hdsdp_retcode sdpSparseConeIAllocDualMat( hdsdp_cone_sdp_sparse *cone ) {
    
    hdsdp_retcode retcode = HDSDP_RETCODE_OK;
    
    /* Build up the dual matrix for SDP sparse cone. The logic is exactly the same as for
       the dense cone
     */
    cone->isDualSparse = 1;
    int nDualNonzeros = 0;
    
    /* If there is a dense matrix, the dual matrix will be dense */
    if ( cone->sdpConeStats[SDP_COEFF_DENSE] > 0 ||
         cone->sdpConeStats[SDP_COEFF_DSR1] > 0  ||
         cone->nCol < SMALL_DUAL_THRESHOLD ) {
        cone->isDualSparse = 0;
    }
    
    /* Allocate the dual sparsity pattern */
    if ( cone->isDualSparse ) {
        HDSDP_INIT(cone->dualPosToElemMap, int, PACK_NNZ(cone->nCol));
        HDSDP_MEMCHECK(cone->dualPosToElemMap);
        
        /* Accumulate the sparsity pattern */
        for ( int iElem = 0; iElem < cone->nRow; ++iElem ) {
            sdpDataMatGetMatNz(cone->sdpRow[iElem], cone->dualPosToElemMap);
        }
        sdpDataMatGetMatNz(cone->sdpObj, cone->dualPosToElemMap);
        
        /* Count statistic of nonzeros and determine matrix type */
        for ( int iElem = 0; iElem < PACK_NNZ(cone->nCol); ++iElem ) {
            nDualNonzeros += cone->dualPosToElemMap[iElem];
        }
        
        if ( nDualNonzeros >= SPARSE_DUAL_THRESHOLD * cone->nCol * cone->nCol ) {
            /* We decide to use dense representation */
            cone->isDualSparse = 0;
        }
    }
    
    /* Exactly the same as the dense cone */
    if ( !cone->isDualSparse ) {
        /* We use dense representation of the dual matrix. No*/
        HDSDP_FREE(cone->dualPosToElemMap);
        HDSDP_INIT(cone->dualMatElem, double, cone->nCol * cone->nCol);
        HDSDP_MEMCHECK(cone->dualMatElem);
        HDSDP_INIT(cone->dualCheckerElem, double, cone->nCol * cone->nCol);
        HDSDP_MEMCHECK(cone->dualCheckerElem);
        HDSDP_INIT(cone->dualStep, double, cone->nCol * cone->nCol);
        HDSDP_MEMCHECK(cone->dualStep);
        HDSDP_CALL(HFpLinsysCreate(&cone->dualFactor, cone->nCol, HDSDP_LINSYS_DENSE_DIRECT));
        /* Done */
    } else {
        /* We are using the sparse dual representation */
        assert( nDualNonzeros > 0 );
        HDSDP_INIT(cone->dualMatBeg, int, cone->nCol + 1);
        HDSDP_MEMCHECK(cone->dualMatBeg);
        HDSDP_INIT(cone->dualMatIdx, int, nDualNonzeros);
        HDSDP_MEMCHECK(cone->dualMatIdx);
        HDSDP_INIT(cone->dualMatElem, double, nDualNonzeros);
        HDSDP_MEMCHECK(cone->dualMatElem);
        HDSDP_INIT(cone->dualCheckerElem, double, nDualNonzeros);
        HDSDP_MEMCHECK(cone->dualCheckerElem);
        HDSDP_INIT(cone->dualStep, double, nDualNonzeros);
        HDSDP_MEMCHECK(cone->dualStep);
        
        /* Now start enumerating over the mapping and construct cone data structure */
        int iNz = 0;
        int *colPtr = cone->dualPosToElemMap;
        for ( int iCol = 0; iCol < cone->nCol; ++iCol ) {
            for ( int iRow = 0; iRow < cone->nCol - iCol; ++iRow ) {
                if ( colPtr[iRow] ) {
                    colPtr[iRow] = iNz;
                    cone->dualMatIdx[iNz] = iRow + iCol;
                    iNz += 1;
                }
            }
            cone->dualMatBeg[iCol + 1] = iNz;
            colPtr += (cone->nCol - iCol);
        }
        
        HDSDP_CALL(HFpLinsysCreate(&cone->dualFactor, cone->nCol, HDSDP_LINSYS_SPARSE_DIRECT));
        HFpLinsysSetParam(cone->dualFactor, -1.0, -1.0, 4, -1, -1);
        
        /* Done */
    }
    
exit_cleanup:
    return retcode;
}

static void sdpDenseConeIFreeDualMat( hdsdp_cone_sdp_dense *cone ) {
    
    if ( !cone ) {
        return;
    }
    
    if ( !cone->isDualSparse ) {
        HDSDP_FREE(cone->dualMatElem);
        HDSDP_FREE(cone->dualCheckerElem);
        HDSDP_FREE(cone->dualStep);
    } else {
        HDSDP_FREE(cone->dualMatBeg);
        HDSDP_FREE(cone->dualMatIdx);
        HDSDP_FREE(cone->dualMatElem);
        HDSDP_FREE(cone->dualCheckerElem);
        HDSDP_FREE(cone->dualStep);
        HDSDP_FREE(cone->dualPosToElemMap);
    }
    
    HFpLinsysDestroy(&cone->dualFactor);
    
    return;
}

static void sdpSparseConeIFreeDualMat( hdsdp_cone_sdp_sparse *cone ) {
    
    if ( !cone ) {
        return;
    }
    
    if ( !cone->isDualSparse ) {
        HDSDP_FREE(cone->dualMatElem);
        HDSDP_FREE(cone->dualCheckerElem);
        HDSDP_FREE(cone->dualStep);
    } else {
        HDSDP_FREE(cone->dualMatBeg);
        HDSDP_FREE(cone->dualMatIdx);
        HDSDP_FREE(cone->dualMatElem);
        HDSDP_FREE(cone->dualCheckerElem);
        HDSDP_FREE(cone->dualStep);
        HDSDP_FREE(cone->dualPosToElemMap);
    }
    
    HFpLinsysDestroy(&cone->dualFactor);
    
    return;
}

/** @brief Create a dense sdp cone
 *
 */
extern hdsdp_retcode sdpDenseConeCreateImpl( hdsdp_cone_sdp_dense **pCone ) {
    
    hdsdp_retcode retcode = HDSDP_RETCODE_OK;
    
    HDSDP_NULLCHECK(pCone);
    hdsdp_cone_sdp_dense *cone = NULL;
    
    HDSDP_INIT(cone, hdsdp_cone_sdp_dense, 1);
    HDSDP_MEMCHECK(cone);
    HDSDP_ZERO(cone, hdsdp_cone_sdp_dense, 1);
    *pCone = cone;
    
exit_cleanup:
    
    return retcode;
}

extern hdsdp_retcode sdpSparseConeCreateImpl( hdsdp_cone_sdp_sparse **pCone ) {
    
    hdsdp_retcode retcode = HDSDP_RETCODE_OK;
    HDSDP_NULLCHECK(pCone);
    hdsdp_cone_sdp_sparse *cone = NULL;
    HDSDP_INIT(cone, hdsdp_cone_sdp_sparse, 1);
    HDSDP_MEMCHECK(cone);
    HDSDP_ZERO(cone, hdsdp_cone_sdp_sparse, 1);
    *pCone = cone;
    
exit_cleanup:
    
    return retcode;
}

/** @brief Set data into a dense cone
 *
 */
extern hdsdp_retcode sdpDenseConeProcDataImpl( hdsdp_cone_sdp_dense *cone, int nRow, int nCol,
                                               int *coneMatBeg, int *coneMatIdx, double *coneMatElem ) {
    
    hdsdp_retcode retcode = HDSDP_RETCODE_OK;
    
    /* Create the cone */
    cone->nRow = nRow;
    cone->nCol = nCol;
    
    HDSDP_INIT(cone->sdpRow, sdp_coeff *, nRow);
    HDSDP_MEMCHECK(cone->sdpRow);
    
    HDSDP_INIT(cone->sdpConePerm, int, nRow);
    HDSDP_MEMCHECK(cone->sdpConePerm);
    
    HDSDP_ZERO(cone->sdpConeStats, int, 5);
    
    /* Primal objective */
    HDSDP_CALL(sdpDataMatCreate(&cone->sdpObj));
    HDSDP_CALL(sdpDataMatSetData(cone->sdpObj, nCol, coneMatBeg[1], coneMatIdx, coneMatElem));
    cone->sdpConeStats[sdpDataMatGetType(cone->sdpObj)] += 1;
    
    /* Constraint data */
    for ( int iRow = 0; iRow < nRow; ++iRow ) {
        HDSDP_CALL(sdpDataMatCreate(&cone->sdpRow[iRow]));
        HDSDP_CALL(sdpDataMatSetData(cone->sdpRow[iRow], nCol, coneMatBeg[iRow + 2] - coneMatBeg[iRow + 1],
                                     coneMatIdx + coneMatBeg[iRow + 1], coneMatElem + coneMatBeg[iRow + 1]));
        sdp_coeff_type coeffType = sdpDataMatGetType(cone->sdpRow[iRow]);
        cone->sdpConeStats[coeffType] += 1;
    }
    
    /* Allocate the dual matrix */
    HDSDP_CALL(sdpDenseConeIAllocDualMat(cone));
        
exit_cleanup:
    
    return retcode;
}

extern hdsdp_retcode sdpSparseConeProcDataImpl( hdsdp_cone_sdp_sparse *cone, int nRow, int nCol,
                                               int *coneMatBeg, int *coneMatIdx, double *coneMatElem ) {
    
    hdsdp_retcode retcode = HDSDP_RETCODE_OK;
    
    /* Create the cone */
    cone->nRow = nRow;
    cone->nCol = nCol;
    
    HDSDP_ZERO(cone->sdpConeStats, int, 5);
    
    /* Primal objective */
    HDSDP_CALL(sdpDataMatCreate(&cone->sdpObj));
    HDSDP_CALL(sdpDataMatSetData(cone->sdpObj, nCol, coneMatBeg[1], coneMatIdx, coneMatElem));
    cone->sdpConeStats[sdpDataMatGetType(cone->sdpObj)] += 1;
    
    /* Sparse constraint data */
    int nRowElem = 0;
    
    /* Count #Nonzeros */
    nRowElem = csp_nnz_cols(nRow, &coneMatBeg[1]);
    
    HDSDP_INIT(cone->sdpRow, sdp_coeff *, nRowElem);
    HDSDP_INIT(cone->rowIdx, int, nRowElem);
    HDSDP_MEMCHECK(cone->sdpRow);
    HDSDP_MEMCHECK(cone->rowIdx);
    
    int nRowElemTmp = 0;
    for ( int iRow = 0; iRow < nRow; ++iRow ) {
        if ( coneMatBeg[iRow + 2] - coneMatBeg[iRow + 1] > 0 ) {
            HDSDP_CALL(sdpDataMatCreate(&cone->sdpRow[nRowElemTmp]));
            HDSDP_CALL(sdpDataMatSetData(cone->sdpRow[nRowElemTmp], nCol, coneMatBeg[iRow + 2] - coneMatBeg[iRow + 1],
                                         coneMatIdx + coneMatBeg[iRow + 1], coneMatElem + coneMatBeg[iRow + 1]));
            sdp_coeff_type coeffType = sdpDataMatGetType(cone->sdpRow[nRowElemTmp]);
#ifdef HDSDP_CONIC_DEBUG
            assert( coeffType != SDP_COEFF_ZERO );
#endif
            cone->sdpConeStats[coeffType] += 1;
            cone->rowIdx[nRowElemTmp] = iRow;
            nRowElemTmp += 1;
        } else {
            cone->sdpConeStats[SDP_COEFF_ZERO] += 1;
        }
    }
    
    cone->nRowElem = nRowElem;
    
#ifdef HDSDP_CONIC_DEBUG
    assert( nRowElem == nRowElemTmp );
#endif
    
    /* Allocate the dual matrix */
    HDSDP_CALL(sdpSparseConeIAllocDualMat(cone));
    
exit_cleanup:
    
    return retcode;
}

extern hdsdp_retcode sdpDenseConePresolveImpl( hdsdp_cone_sdp_dense *cone ) {
    
    hdsdp_retcode retcode = HDSDP_RETCODE_OK;
    
    double *preConeAuxi = NULL;
    /* cone->nCol * cone->nCol for a Full matrix */
    HDSDP_INIT(preConeAuxi, double, cone->nCol);
    HDSDP_MEMCHECK(preConeAuxi);
    
    HDSDP_CALL(sdpDataMatBuildUpEigs(cone->sdpObj, preConeAuxi));
    for ( int iRow = 0; iRow < cone->nRow; ++iRow ) {
        /* Update conic statistics */
        cone->sdpConeStats[sdpDataMatGetType(cone->sdpRow[iRow])] -= 1;
        HDSDP_CALL(sdpDataMatBuildUpEigs(cone->sdpRow[iRow], preConeAuxi));
        cone->sdpConeStats[sdpDataMatGetType(cone->sdpRow[iRow])] += 1;
    }
    
    /* Symbolically factorize the dual matrix */
    HDSDP_CALL(HFpLinsysSymbolic(cone->dualFactor, cone->dualMatBeg, cone->dualMatIdx));
    
exit_cleanup:
    
    HDSDP_FREE(preConeAuxi);
    return retcode;
}

extern hdsdp_retcode sdpSparseConePresolveImpl( hdsdp_cone_sdp_sparse *cone ) {
    
    hdsdp_retcode retcode = HDSDP_RETCODE_OK;
    
    double *preConeAuxi = NULL;
    /* cone->nCol * cone->nCol for a Full matrix */
    HDSDP_INIT(preConeAuxi, double, cone->nCol);
    HDSDP_MEMCHECK(preConeAuxi);
    
    HDSDP_CALL(sdpDataMatBuildUpEigs(cone->sdpObj, preConeAuxi));
    for ( int iElem = 0; iElem < cone->nRowElem; ++iElem ) {
        cone->sdpConeStats[sdpDataMatGetType(cone->sdpRow[iElem])] -= 1;
        HDSDP_CALL(sdpDataMatBuildUpEigs(cone->sdpRow[iElem], preConeAuxi));
        cone->sdpConeStats[sdpDataMatGetType(cone->sdpRow[iElem])] += 1;
    }
    
    /* Symbolically factorize the dual matrix */
    HDSDP_CALL(HFpLinsysSymbolic(cone->dualFactor, cone->dualMatBeg, cone->dualMatIdx));
    
exit_cleanup:
    
    HDSDP_FREE(preConeAuxi);
    return retcode;
}

extern void sdpDenseConeSetStartImpl( hdsdp_cone_sdp_dense *cone, double rResi ) {
    
    cone->dualResidual = rResi;
    return;
}

extern void sdpSparseConeSetStartImpl( hdsdp_cone_sdp_sparse *cone, double rResi ) {
    
    cone->dualResidual = rResi;
    return;
}

extern double sdpDenseConeGetObjNorm( hdsdp_cone_sdp_dense *cone, int whichNorm ) {
    
    return sdpDataMatNorm(cone->sdpObj, whichNorm);
}

extern double sdpSparseConeGetObjNorm( hdsdp_cone_sdp_sparse *cone, int whichNorm ) {
    
    return sdpDataMatNorm(cone->sdpObj, whichNorm);
}

extern double sdpDenseConeGetCoeffNorm( hdsdp_cone_sdp_dense *cone, int whichNorm ) {
    
    double norm = 0.0;
    if ( whichNorm == ABS_NORM ) {
        for ( int iRow = 0; iRow < cone->nRow; ++iRow ) {
            norm += sdpDataMatNorm(cone->sdpRow[iRow], ABS_NORM);
        }
    } else {
        for ( int iRow = 0; iRow < cone->nRow; ++iRow ) {
            double tmpNorm = sdpDataMatNorm(cone->sdpRow[iRow], FRO_NORM);
            norm += tmpNorm * tmpNorm;
        }
        norm = sqrt(norm);
    }
    
    return norm;
}

extern double sdpSparseConeGetCoeffNorm( hdsdp_cone_sdp_sparse *cone, int whichNorm ) {
    
    double norm = 0.0;
    if ( whichNorm == ABS_NORM ) {
        for ( int iElem = 0; iElem < cone->nRowElem; ++iElem ) {
            norm += sdpDataMatNorm(cone->sdpRow[iElem], ABS_NORM);
        }
    } else {
        for ( int iElem = 0; iElem < cone->nRowElem; ++iElem ) {
            double tmpNorm = sdpDataMatNorm(cone->sdpRow[iElem], FRO_NORM);
            norm += tmpNorm * tmpNorm;
        }
        norm = sqrt(norm);
    }
    
    return norm;
}

extern void sdpDenseConeUpdateImpl( hdsdp_cone_sdp_dense *cone, double barHsdTau, double *rowDual ) {
    
    return;
}

extern void sdpSparseConeUpdateImpl( hdsdp_cone_sdp_sparse *cone, double barHsdTau, double *rowDual ) {
    
    return;
}

extern double sdpDenseConeRatioTestImpl( hdsdp_cone_sdp_dense *cone, double barHsdTauStep, double *rowDualStep ) {
    
    return 0.0;
}

extern double sdpSparseConeRatioTestImpl( hdsdp_cone_sdp_sparse *cone, double barHsdTauStep, double *rowDualStep ) {
    
    return 0.0;
}

extern int64_t sdpDenseConeGetSymNnzImpl( hdsdp_cone_sdp_dense *cone ) {
    
    
    return 0;
}

extern int64_t sdpSparseConeGetSymNnzImpl( hdsdp_cone_sdp_sparse *cone ) {
    
    
    return 0;
}


extern void sdpDenseConeAddSymNnzImpl( hdsdp_cone_sdp_dense *cone, int *schurMatCol ) {
    
    
    return;
}

extern void sdpSparseConeAddSymNnzImpl( hdsdp_cone_sdp_sparse *cone, int *schurMatCol ) {
    
    
    return;
}

extern void sdpDenseConeClearImpl( hdsdp_cone_sdp_dense *cone ) {
    
    if ( !cone ) {
        return;
    }
    
    HDSDP_FREE(cone->sdpConePerm);
    sdpDataMatDestroy(&cone->sdpObj);
    
    for ( int iRow = 0; iRow < cone->nRow; ++iRow ) {
        sdpDataMatDestroy(&cone->sdpRow[iRow]);
    }
    
    HDSDP_FREE(cone->sdpRow);
    sdpDenseConeIFreeDualMat(cone);
    
    HDSDP_ZERO(cone, hdsdp_cone_sdp_dense, 1);
    
    return;
}

extern void sdpSparseConeClearImpl( hdsdp_cone_sdp_sparse *cone ) {
    
    if ( !cone ) {
        return;
    }
    
    HDSDP_FREE(cone->rowIdx);
    sdpDataMatDestroy(&cone->sdpObj);
    
    for ( int iRow = 0; iRow < cone->nRowElem; ++iRow ) {
        sdpDataMatDestroy(&cone->sdpRow[iRow]);
    }
    
    HDSDP_FREE(cone->sdpRow);
    sdpSparseConeIFreeDualMat(cone);
    
    HDSDP_ZERO(cone, hdsdp_cone_sdp_sparse, 1);
    
    return;
}

extern void sdpDenseConeDestroyImpl( hdsdp_cone_sdp_dense **pCone ) {
    
    if ( !pCone ) {
        return;
    }
    
    sdpDenseConeClearImpl(*pCone);
    HDSDP_FREE(*pCone);
    
    return;
}

extern void sdpSparseConeDestroyImpl( hdsdp_cone_sdp_sparse **pCone ) {
    
    if ( !pCone ) {
        return;
    }
    
    sdpSparseConeClearImpl(*pCone);
    HDSDP_FREE(*pCone);
    
    return;
}

extern void sdpDenseConeViewImpl( hdsdp_cone_sdp_dense *cone ) {
    
    printf("- Dense SDP cone of %d rows. \n", cone->nRow);
    printf("- Objective: \n");
    sdpDataMatView(cone->sdpObj);
    
    printf("- Constraint: \n");
    for ( int iRow = 0; iRow < cone->nRow; ++iRow ) {
        sdpDataMatView(cone->sdpRow[iRow]);
    }
    
    printf("- Conic statistics: Zero %d Sp %d Ds %d SpR1 %d DsR1 %d \n\n", cone->sdpConeStats[SDP_COEFF_ZERO],
           cone->sdpConeStats[SDP_COEFF_SPARSE], cone->sdpConeStats[SDP_COEFF_DENSE],
           cone->sdpConeStats[SDP_COEFF_SPR1], cone->sdpConeStats[SDP_COEFF_DSR1]);
    
    printf("- Dual sparsity: \n ");
    
    if ( cone->isDualSparse ) {
        dcs A;
        A.nz = -1;
        A.m = cone->nCol;
        A.n = cone->nCol;
        A.p = cone->dualMatBeg;
        A.i = cone->dualMatIdx;
        A.x = cone->dualMatElem;
        
        dcs_print(&A, 1);
    } else {
        printf("Using dense dual matrix. \n");
    }
    
    return;
}

extern void sdpSparseConeViewImpl( hdsdp_cone_sdp_sparse *cone ) {
    
    printf("- Sparse SDP cone of %d rows and %d nonzeros. \n", cone->nRow, cone->nRowElem);
    printf("- Objective: \n");
    sdpDataMatView(cone->sdpObj);
    
    printf("- Constraint: \n");
    for ( int iRow = 0; iRow < cone->nRowElem; ++iRow ) {
        printf("%d: ", cone->rowIdx[iRow]);
        sdpDataMatView(cone->sdpRow[iRow]);
    }
    
    printf("Conic statistics: Zero %d Sp %d Ds %d SpR1 %d DsR1 %d \n\n", cone->sdpConeStats[SDP_COEFF_ZERO],
           cone->sdpConeStats[SDP_COEFF_SPARSE], cone->sdpConeStats[SDP_COEFF_DENSE],
           cone->sdpConeStats[SDP_COEFF_SPR1], cone->sdpConeStats[SDP_COEFF_DSR1]);
    
    printf("- Dual sparsity: \n ");
    
    if ( cone->isDualSparse ) {
        dcs A;
        A.nz = -1;
        A.m = cone->nCol;
        A.n = cone->nCol;
        A.p = cone->dualMatBeg;
        A.i = cone->dualMatIdx;
        A.x = cone->dualMatElem;
        
        dcs_print(&A, 1);
    } else {
        printf("Using dense dual matrix. \n");
    }
    
    return;
}

