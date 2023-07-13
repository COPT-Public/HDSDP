#ifdef HEADERPATH
#include "interface/hdsdp_conic_sdp.h"
#include "interface/def_hdsdp_user_data.h"
#include "interface/hdsdp_utils.h"
#include "interface/def_hdsdp_schur.h"
#include "linalg/hdsdp_sdpdata.h"
#include "linalg/sparse_opts.h"
#include "linalg/dense_opts.h"
#include "linalg/vec_opts.h"
#include "linalg/hdsdp_linsolver.h"
#include "external/hdsdp_cs.h"
#else
#include "hdsdp_conic_sdp.h"
#include "def_hdsdp_user_data.h"
#include "hdsdp_utils.h"
#include "def_hdsdp_schur.h"
#include "hdsdp_sdpdata.h"
#include "sparse_opts.h"
#include "dense_opts.h"
#include "vec_opts.h"
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
        /* First the diagonal must be filled */
        int *iPos = cone->dualPosToElemMap;
        for ( int iCol = 0; iCol < cone->nCol; ++iCol ) {
            iPos[0] = 1;
            iPos += (cone->nCol - iCol);
        }
        
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
        /* First accumulate the diagonal matrix */
        int *iPos = cone->dualPosToElemMap;
        for ( int iCol = 0; iCol < cone->nCol; ++iCol ) {
            iPos[0] = 1;
            iPos += (cone->nCol - iCol);
        }
        for ( int iElem = 0; iElem < cone->nRowElem; ++iElem ) {
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

#define BUFFER_DUALVAR   (0)
#define BUFFER_DUALCHECK (1)
#define BUFFER_DUALSTEP  (2)
static void sdpDenseConeIZeroBuffer( hdsdp_cone_sdp_dense *cone, int whichBuffer ) {
    
    double *target = NULL;
    if ( whichBuffer == BUFFER_DUALVAR ) {
        target = cone->dualMatElem;
    } else if ( whichBuffer == BUFFER_DUALCHECK ) {
        target = cone->dualCheckerElem;
    } else {
        target = cone->dualStep;
    }
    
    if ( cone->isDualSparse ) {
        HDSDP_ZERO(target, double, cone->dualMatBeg[cone->nCol]);
    } else {
        HDSDP_ZERO(target, double, cone->nCol * cone->nCol);
    }
    
    return;
}

static void sdpSparseConeIZeroBuffer( hdsdp_cone_sdp_sparse *cone, int whichBuffer ) {
    
    double *target = NULL;
    if ( whichBuffer == BUFFER_DUALVAR ) {
        target = cone->dualMatElem;
    } else if ( whichBuffer == BUFFER_DUALCHECK ) {
        target = cone->dualCheckerElem;
    } else {
        target = cone->dualStep;
    }
    
    if ( cone->isDualSparse ) {
        HDSDP_ZERO(target, double, cone->dualMatBeg[cone->nCol]);
    } else {
        HDSDP_ZERO(target, double, cone->nCol * cone->nCol);
    }
    
    return;
}

static inline void sdpDenseConeIUpdateBuffer( hdsdp_cone_sdp_dense *cone, double dCCoef, double dACoefScal,
                                       double *dACoef, double dEyeCoef, int whichBuffer ) {
    /* Assemble the buffer by taking a free linear combination between coefficients
     
       B <- dResiCoef * I + dACoefScal * A' * y + C * dCCoef
     */
    
    double *target = NULL;
    
    switch (whichBuffer) {
        case BUFFER_DUALVAR:
            target = cone->dualMatElem;
            break;
        case BUFFER_DUALCHECK:
            target = cone->dualCheckerElem;
            break;
        case BUFFER_DUALSTEP:
            target = cone->dualStep;
            break;
        default:
            break;
    }
    
    /* Zero out buffer space*/
    if ( cone->isDualSparse ) {
        HDSDP_ZERO(target, double, cone->dualMatBeg[cone->nCol]);
    } else {
        HDSDP_ZERO(target, double, cone->nCol * cone->nCol);
    }
    
    /* Aggregate -A' * y */
    for ( int iRow = 0; iRow < cone->nRow; ++iRow ) {
        sdpDataMatAddToBuffer(cone->sdpRow[iRow], dACoefScal * dACoef[iRow], cone->dualPosToElemMap, target);
    }
    /* Add c * tau */
    sdpDataMatAddToBuffer(cone->sdpObj, dCCoef, cone->dualPosToElemMap, target);
    
    /* Add residual */
    if ( dEyeCoef != 0.0 ) {
        if ( cone->isDualSparse ) {
            int *iPos = cone->dualPosToElemMap;
            for ( int iCol = 0; iCol < cone->nCol; ++iCol ) {
                target[iPos[0]] += dEyeCoef;
                iPos += (cone->nCol - iCol);
            }
        } else {
            for ( int iCol = 0; iCol < cone->nCol; ++iCol ) {
                cone->dualMatElem[iCol + iCol * cone->nCol] += dEyeCoef;
            }
        }
    }
    
    return;
}

static inline void sdpSparseConeIUpdateBuffer( hdsdp_cone_sdp_sparse *cone, double dCCoef, double dACoefScal,
                                        double *dACoef, double dEyeCoef, int whichBuffer ) {
    
    double *target = NULL;
    
    switch (whichBuffer) {
        case BUFFER_DUALVAR:
            target = cone->dualMatElem;
            break;
        case BUFFER_DUALCHECK:
            target = cone->dualCheckerElem;
            break;
        case BUFFER_DUALSTEP:
            target = cone->dualStep;
            break;
        default:
            break;
    }
    
    /* Zero out buffer space*/
    if ( cone->isDualSparse ) {
        HDSDP_ZERO(target, double, cone->dualMatBeg[cone->nCol]);
    } else {
        HDSDP_ZERO(target, double, cone->nCol * cone->nCol);
    }
    
    /* Aggregate -A' * y */
    for ( int iElem = 0; iElem < cone->nRowElem; ++iElem ) {
        sdpDataMatAddToBuffer(cone->sdpRow[iElem], dACoefScal * dACoef[cone->rowIdx[iElem]],
                              cone->dualPosToElemMap, target);
    }
    /* Add c * tau */
    sdpDataMatAddToBuffer(cone->sdpObj, dCCoef, cone->dualPosToElemMap, target);
    
    /* Add residual */
    if ( dEyeCoef != 0.0 ) {
        if ( cone->isDualSparse ) {
            int *iPos = cone->dualPosToElemMap;
            for ( int iCol = 0; iCol < cone->nCol; ++iCol ) {
                target[iPos[0]] += dEyeCoef;
                iPos += (cone->nCol - iCol);
            }
        } else {
            for ( int iCol = 0; iCol < cone->nCol; ++iCol ) {
                cone->dualMatElem[iCol + iCol * cone->nCol] += dEyeCoef;
            }
        }
    }
    
    return;
}

static void sdpDenseConeILanczosMultiply( void *cone, double *dLhsVec, double *dRhsVec ) {
    /*
     Implement the Matrix vector multiplication for the Lanczos solver.
     Recall that when computing the dual stepsize in HDSDP, we need to evaluate alpha such that
     
       S + alpha * dS >= 0 and S = L * L', which is equivalent to finding alpha such that
     
       I + alpha * L^-1 * dS * L^-T >= 0 and each matrix-vector multiplication involves
     
     1. a backward solve L^T * x = w
     2. a general multiplication y = - dS * x
     3. a forward solve L * z = y
     
     Depending on the conic dual sparsity pattern, the solve and multiplication procedures will be different.
     */
    
    hdsdp_cone_sdp_dense *dsCone = (hdsdp_cone_sdp_dense *) cone;
    
    /* First step: backward solve */
    HFpLinsysBSolve(dsCone->dualFactor, 1, dLhsVec, dRhsVec);
    
    /* Second step: multiplication. Multiply dRhsVec by -dS */
    if ( dsCone->isDualSparse ) {
        HDSDP_ZERO(dsCone->dVecBuffer, double, dsCone->nCol);
        /* Compute -dS * x. Note that only the lower triangular part is stored in dS. */
        for ( int iCol = 0; iCol < dsCone->nCol; ++iCol ) {
            /* The first element of each column must be on the diagonal due to the identity dual residual */
            dsCone->dVecBuffer[dsCone->dualMatIdx[dsCone->dualMatBeg[iCol]]] -= \
            dsCone->dualStep[iCol] * dRhsVec[dsCone->dualMatIdx[dsCone->dualMatBeg[iCol]]];
            /* For each of element not in the diagonal, we have to map it to its symmetric position */
            for ( int iRow = dsCone->dualMatBeg[iCol] + 1; iRow < dsCone->dualMatBeg[iCol + 1]; ++iRow ) {
                dsCone->dVecBuffer[dsCone->dualMatIdx[iRow]] -= \
                dsCone->dualStep[iCol] * dRhsVec[dsCone->dualMatIdx[iRow]];
                dsCone->dVecBuffer[dsCone->dualMatIdx[iCol]] -= \
                dsCone->dualStep[iRow] * dRhsVec[dsCone->dualMatIdx[iCol]];
            }
        }
    } else {
        fds_symv(dsCone->nCol, -1.0, dsCone->dualStep, dRhsVec, 0.0, dsCone->dVecBuffer);
    }
    
    /* Last step: forward solve */
    HFpLinsysFSolve(dsCone->dualFactor, 1, dsCone->dVecBuffer, dRhsVec);
    
    return;
}

static void sdpSparseConeILanczosMultiply( void *cone, double *dLhsVec, double *dRhsVec ) {
    
    hdsdp_cone_sdp_sparse *spCone = (hdsdp_cone_sdp_sparse *) cone;
    
    /* First step: backward solve */
    HFpLinsysBSolve(spCone->dualFactor, 1, dLhsVec, dRhsVec);
    
    /* Second step: multiplication. Multiply dRhsVec by -dS */
    if ( spCone->isDualSparse ) {
        HDSDP_ZERO(spCone->dVecBuffer, double, spCone->nCol);
        /* Compute -dS * x. Note that only the lower triangular part is stored in dS. */
        for ( int iCol = 0; iCol < spCone->nCol; ++iCol ) {
            /* The first element of each column must be on the diagonal due to the identity dual residual */
            spCone->dVecBuffer[spCone->dualMatIdx[spCone->dualMatBeg[iCol]]] -= \
            spCone->dualStep[iCol] * dRhsVec[spCone->dualMatIdx[spCone->dualMatBeg[iCol]]];
            /* For each of element not in the diagonal, we have to map it to its symmetric position */
            for ( int iRow = spCone->dualMatBeg[iCol] + 1; iRow < spCone->dualMatBeg[iCol + 1]; ++iRow ) {
                spCone->dVecBuffer[spCone->dualMatIdx[iRow]] -= \
                spCone->dualStep[iCol] * dRhsVec[spCone->dualMatIdx[iRow]];
                spCone->dVecBuffer[spCone->dualMatIdx[iCol]] -= \
                spCone->dualStep[iRow] * dRhsVec[spCone->dualMatIdx[iCol]];
            }
        }
    } else {
        fds_symv(spCone->nCol, -1.0, spCone->dualStep, dRhsVec, 0.0, spCone->dVecBuffer);
    }
    
    /* Last step: forward solve */
    HFpLinsysFSolve(spCone->dualFactor, 1, spCone->dVecBuffer, dRhsVec);
    
    return;
}

#define SPARSE_EFFICIENCY (1.5)
static int sdpDenseConeIChooseKKTStrategy( int *rowRanks, int *rowSparsity, int *rowPerm, int nRow, int nCol, int iPerm ) {
    
    int bestKKTStrategy = KKT_M1;
    double bestKKTScore = HDSDP_INFINITY;
    
    double KKTScore1 = HDSDP_INFINITY;
    double KKTScore2 = 0.0;
    double KKTScore3 = 0.0;
    double KKTScore4 = 0.0;
    double KKTScore5 = 0.0;
    
    /* Get the rank information of the current matrix */
    int rowIdx = rowPerm[iPerm];
    int rowRank = rowRanks[rowIdx];
    double nCubed = (double) nCol * nCol * nCol;
    
    /* Compute the cumulative sparsity after iPerm */
    double nNzAfter = 0.0;
    
    for ( int iRow = iPerm; iRow < nRow; ++iRow ) {
        nNzAfter += rowSparsity[iRow];
    }
    
    /* Compute the score function */
    KKTScore2 = rowRank * ((double) rowSparsity[iPerm] * nCol + 3 * SPARSE_EFFICIENCY * nNzAfter);
    KKTScore3 = (double) nCol * SPARSE_EFFICIENCY * rowSparsity[rowIdx] + nCubed + SPARSE_EFFICIENCY * nNzAfter;
    KKTScore4 = (double) nCol * SPARSE_EFFICIENCY * rowSparsity[rowIdx] + SPARSE_EFFICIENCY * (nCol + 1) * nNzAfter;
    KKTScore5 = SPARSE_EFFICIENCY * (2.0 * SPARSE_EFFICIENCY * rowSparsity[rowIdx] + 1) * nNzAfter;
        
    /* Choose the best KKT strategy */
    if ( KKTScore2 <= bestKKTScore ) {
        bestKKTStrategy = KKT_M2;
        bestKKTScore = KKTScore2;
    }
    
    if ( KKTScore3 < bestKKTScore ) {
        bestKKTStrategy = KKT_M3;
        bestKKTScore = KKTScore3;
    }
    
    if ( KKTScore4 < bestKKTScore ) {
        bestKKTStrategy = KKT_M4;
        bestKKTScore = KKTScore4;
    }
    
    if ( KKTScore5 < bestKKTScore ) {
        bestKKTStrategy = KKT_M5;
        bestKKTScore = KKTScore5;
    }
    
#if 0
    printf("Row %d M1: %e M2: %e M3: %e M4: %e M5: %e \n",
           rowPerm[iPerm], KKTScore1, KKTScore2, KKTScore3, KKTScore4, KKTScore5);
#endif
    
    return bestKKTStrategy;
}

static hdsdp_retcode sdpDenseConeIGetKKTOrdering( hdsdp_cone_sdp_dense *cone ) {
    
    /* Compute the KKT ordering of the SDP cone
     
     In HDSDP, when setting up the Schur complement matrix, we adopt 4 strategies that are called
     M1 to M4 in the solver.
      
     --------------------
     Rank based
     --------------------
     M1: Not in use
     --
     M2: Make use of both low-rank and sparsity of the coefficient matrices
     --
     
     --------------------
     
     --------------------
     Sparsity based
     --------------------
     M2:
     M3:
     M4:
     
     
     
     Note: different from written in the paper, the original M1 strategy is abandoned since it's not competitive
           for problems where eigen-decomposition is expensive. Moreover, the three recently introduced sparsity-based
           strategies are more competitive
     */
    hdsdp_retcode retcode = HDSDP_RETCODE_OK;
    
    /* Prepare rank information */
    int *rowRanks = NULL;
    int *rowSparsity = NULL;
    
    HDSDP_INIT(rowRanks, int, cone->nRow);
    HDSDP_MEMCHECK(rowRanks);
    
    /* Prepare sparsity information */
    HDSDP_INIT(rowSparsity, int, cone->nRow);
    HDSDP_MEMCHECK(rowSparsity);
    
    /* Enumerate all the rows to collect all the needed information */
    for ( int iRow = 0; iRow < cone->nRow; ++iRow ) {
        cone->sdpConePerm[iRow] = iRow;
        rowRanks[iRow] = sdpDataMatGetRank(cone->sdpRow[iRow]);
        rowSparsity[iRow] = sdpDataMatGetNnz(cone->sdpRow[iRow]);
    }
    
    /* Generate permutation based on the sparsity pattern */
    HUtilSortIntByInt(cone->sdpConePerm, rowSparsity, 0, cone->nRow - 1);
    
    /* Choose KKT strategies */
    for ( int iPerm = 0; iPerm < cone->nRow; ++iPerm ) {
        cone->KKTStrategies[iPerm] = \
        sdpDenseConeIChooseKKTStrategy(rowRanks, rowSparsity,
                                       cone->sdpConePerm, cone->nRow, cone->nCol, iPerm);
    }
    
exit_cleanup:
        
    HDSDP_FREE(rowRanks);
    HDSDP_FREE(rowSparsity);
        
    return retcode;
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
    
    /* Allocate the dual diagonal (of factorization to compute barrier/potential) */
    HDSDP_INIT(cone->dualDiag, double, cone->nCol);
    HDSDP_MEMCHECK(cone->dualDiag);
    
    /* Allocate the Lanczos solver and set the method */
    HDSDP_CALL(HLanczosCreate(&cone->Lanczos));
    HDSDP_CALL(HLanczosInit(cone->Lanczos, cone->nCol, 30));
    HLanczosSetData(cone->Lanczos, cone, sdpDenseConeILanczosMultiply);
    
    /* Allocate the dual vector buffer */
    HDSDP_INIT(cone->dVecBuffer, double, cone->nCol);
    HDSDP_MEMCHECK(cone->dVecBuffer);
    
    /* Allocate the KKT permutation and strategies */
    HDSDP_INIT(cone->sdpConePerm, int, cone->nCol);
    HDSDP_INIT(cone->KKTStrategies, int, cone->nCol);
        
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
    
    /* Allocate the dual diagonal (of factorization to compute barrier/potential) */
    HDSDP_INIT(cone->dualDiag, double, cone->nCol);
    HDSDP_MEMCHECK(cone->dualDiag);
    
    /* Allocate the Lanczos solver and set the method */
    HDSDP_CALL(HLanczosCreate(&cone->Lanczos));
    HDSDP_CALL(HLanczosInit(cone->Lanczos, cone->nCol, 30));
    HLanczosSetData(cone->Lanczos, cone, sdpSparseConeILanczosMultiply);
    
    /* Allocate the dual vector buffer */
    HDSDP_INIT(cone->dVecBuffer, double, cone->nCol);
    HDSDP_MEMCHECK(cone->dVecBuffer);
    
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
    
    /* KKT re-ordering */
    HDSDP_CALL(sdpDenseConeIGetKKTOrdering(cone));
    
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
    /* Implement the numerical aggregation of conic coefficients
       When this routine is called, the dual matrix will be overwritten by
        
            S = -Rd - A' * y + c * tau
       
     With respect to the buffer utility, we take
     
     dCCoef = tau
     dACoef = y
     dAScal = -1.0
     dEyeCoef = - Rd

     */
    sdpDenseConeIUpdateBuffer(cone, barHsdTau, -1.0, rowDual, -cone->dualResidual, BUFFER_DUALVAR);
    
    return;
}

extern void sdpSparseConeUpdateImpl( hdsdp_cone_sdp_sparse *cone, double barHsdTau, double *rowDual ) {
    
    sdpSparseConeIUpdateBuffer(cone, barHsdTau, -1.0, rowDual, -cone->dualResidual, BUFFER_DUALVAR);
    
    return;
}

extern hdsdp_retcode sdpDenseConeRatioTestImpl( hdsdp_cone_sdp_dense *cone, double barHsdTauStep, double *rowDualStep, double dAdaRatio, double *maxStep ) {
    /* Given conic quantity dy, get maximum distance from the cone
       Given dTau and dy, we first evaluate
     
            dS = dGamma * R_d - A' * dy + c * dTau
       
     In terms of the buffer utility, we take
     
     dCCoef = dTau
     dACoef = dy
     dAScal = -1.0
     dEyeCoef = Rd * dGamma
     
     Then Lanczos solver will take over to compute the largest step we can take so that
     
        S + alpha * dS >= 0
     */
    
    hdsdp_retcode retcode = HDSDP_RETCODE_OK;
    sdpDenseConeIUpdateBuffer(cone, barHsdTauStep, -1.0, rowDualStep, dAdaRatio * cone->dualResidual, BUFFER_DUALSTEP);
    /* Finished with setting up dS. Now do Lanczos computation */
    HDSDP_CALL(HLanczosSolve(cone->Lanczos, NULL, maxStep));
    
exit_cleanup:
    return retcode;
}

extern hdsdp_retcode sdpSparseConeRatioTestImpl( hdsdp_cone_sdp_sparse *cone, double barHsdTauStep, double *rowDualStep, double dAdaRatio, double *maxStep ) {
    
    hdsdp_retcode retcode = HDSDP_RETCODE_OK;
    sdpSparseConeIUpdateBuffer(cone, barHsdTauStep, -1.0, rowDualStep, dAdaRatio * cone->dualResidual, BUFFER_DUALSTEP);
    /* Finished with setting up dS. Now do Lanczos computation */
    HDSDP_CALL(HLanczosSolve(cone->Lanczos, NULL, maxStep));
    
exit_cleanup:
    return retcode;
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

extern hdsdp_retcode sdpDenseConeGetBarrier( hdsdp_cone_sdp_dense *cone, double barHsdTau, double *rowDual, double *logdet ) {
    
    /* Compute the barrier function at current dual slack
       If not supplied the barrier of the most recently updated dual solution is used
     */
    hdsdp_retcode retcode = HDSDP_RETCODE_OK;
    double dLogDeterminant = 0.0;
    
    if ( rowDual ) {
        sdpDenseConeUpdateImpl(cone, barHsdTau, rowDual);
        HDSDP_CALL(HFpLinsysNumeric(cone->dualFactor, cone->dualMatBeg, cone->dualMatIdx, cone->dualMatElem));
    }
    
    HDSDP_CALL(HFpLinsysGetDiag(cone->dualFactor, cone->dualDiag));
    
    /* The log determinant is directly obtained from diagonal of the Cholesky factor
       since log det(S) = log det(L * L') = log (det(L)^2) = 2 log det(L)
     */
    
    for ( int i = 0; i < cone->nCol; ++i ) {
        dLogDeterminant += log(cone->dualDiag[i]);
    }
    
    *logdet = 2.0 * dLogDeterminant;
    
exit_cleanup:
    return retcode;
}

extern hdsdp_retcode sdpSparseConeGetBarrier( hdsdp_cone_sdp_sparse *cone, double barHsdTau, double *rowDual, double *logdet ) {
    
    /* Compute the barrier function at current dual slack
       If not supplied the barrier of the most recently updated dual solution is used
     */
    hdsdp_retcode retcode = HDSDP_RETCODE_OK;
    double dLogDeterminant = 0.0;
    
    if ( rowDual ) {
        sdpSparseConeUpdateImpl(cone, barHsdTau, rowDual);
        HDSDP_CALL(HFpLinsysNumeric(cone->dualFactor, cone->dualMatBeg, cone->dualMatIdx, cone->dualMatElem));
    }
    
    HDSDP_CALL(HFpLinsysGetDiag(cone->dualFactor, cone->dualDiag));
    
    /* The log determinant is directly obtained from diagonal of the Cholesky factor
       since log det(S) = log det(L * L') = log (det(L)^2) = 2 log det(L)
     */
    
    for ( int i = 0; i < cone->nCol; ++i ) {
        dLogDeterminant += log(cone->dualDiag[i]);
    }
    
    *logdet = 2.0 * dLogDeterminant;
    
exit_cleanup:
    return retcode;
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

    HDSDP_FREE(cone->dualDiag);
    HLanczosDestroy(&cone->Lanczos);
    HDSDP_FREE(cone->dVecBuffer);
    
    HDSDP_FREE(cone->sdpConePerm);
    HDSDP_FREE(cone->KKTStrategies);
    
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
    
    HDSDP_FREE(cone->dualDiag);
    HLanczosDestroy(&cone->Lanczos);
    HDSDP_FREE(cone->dVecBuffer);
    
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
    
    printf("- Conic statistics: Zero %d Sp %d Ds %d SpR1 %d DsR1 %d \n", cone->sdpConeStats[SDP_COEFF_ZERO],
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
        
        dcs_print(&A, 0);
    } else {
        printf("- Using dense dual matrix. \n");
        fds_print(cone->nCol, cone->dualMatElem);
    }
    
    printf("- KKT ordering: \n");
    
    int KKTStrategies[5] = {0};
    for ( int iRow = 0; iRow < cone->nRow; ++iRow ) {
        KKTStrategies[cone->KKTStrategies[iRow]] += 1;
    }
    
    printf("- KKT Statistics: M1: %d M2: %d M3: %d M4: %d M5: %d \n",
           KKTStrategies[KKT_M1], KKTStrategies[KKT_M2], KKTStrategies[KKT_M3],
           KKTStrategies[KKT_M4], KKTStrategies[KKT_M5]);
    
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
    
    printf("- Conic statistics: Zero %d Sp %d Ds %d SpR1 %d DsR1 %d \n", cone->sdpConeStats[SDP_COEFF_ZERO],
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
        
        dcs_print(&A, 0);
    } else {
        printf("- Using dense dual matrix. \n");
        fds_print(cone->nCol, cone->dualMatElem);
    }
    
    return;
}

