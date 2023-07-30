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
#include "external/def_hdsdp.h"
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
#include "def_hdsdp.h"
#endif

#include <math.h>

#ifndef SMALL_DUAL_THRESHOLD
#define SMALL_DUAL_THRESHOLD (0)
#endif

#ifndef SPARSE_DUAL_THRESHOLD
#define SPARSE_DUAL_THRESHOLD (0.6)
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
        HDSDP_CALL(HFpLinsysCreate(&cone->dualChecker, cone->nCol, HDSDP_LINSYS_DENSE_DIRECT));
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
        HDSDP_CALL(HFpLinsysCreate(&cone->dualChecker, cone->nCol, HDSDP_LINSYS_SPARSE_DIRECT));
        HFpLinsysSetParam(cone->dualChecker, -1.0, -1.0, 4, -1, -1);
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
        HDSDP_CALL(HFpLinsysCreate(&cone->dualChecker, cone->nCol, HDSDP_LINSYS_DENSE_DIRECT));
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
        HDSDP_CALL(HFpLinsysCreate(&cone->dualChecker, cone->nCol, HDSDP_LINSYS_SPARSE_DIRECT));
        HFpLinsysSetParam(cone->dualChecker, -1.0, -1.0, 4, -1, -1);
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
    HFpLinsysDestroy(&cone->dualChecker);
    
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
    HFpLinsysDestroy(&cone->dualChecker);
    
    return;
}

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
    
    /* Add cone perturbation */
    dEyeCoef += cone->dualPerturb;
    
    /* Add residual and perturbation */
    if ( dEyeCoef != 0.0 ) {
        if ( cone->isDualSparse ) {
            int *iPos = cone->dualPosToElemMap;
            for ( int iCol = 0; iCol < cone->nCol; ++iCol ) {
                target[iPos[0]] += dEyeCoef;
                iPos += (cone->nCol - iCol);
            }
        } else {
            for ( int iCol = 0; iCol < cone->nCol; ++iCol ) {
                target[iCol + iCol * cone->nCol] += dEyeCoef;
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
    
    /* Add cone perturbation */
    dEyeCoef += cone->dualPerturb;
    
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
                target[iCol + iCol * cone->nCol] += dEyeCoef;
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
    HFpLinsysBSolve(dsCone->LTarget, 1, dLhsVec, dRhsVec);
    
    /* Second step: multiplication. Multiply dRhsVec by -dS */
    if ( dsCone->isDualSparse ) {
        HDSDP_ZERO(dsCone->dVecBuffer, double, dsCone->nCol);
        /* Compute -dS * x. Note that only the lower triangular part is stored in dS. */
        for ( int iCol = 0, iRow = 0; iCol < dsCone->nCol; ++iCol ) {
            /* The first element of each column must be on the diagonal due to the identity dual residual */
            iRow = dsCone->dualMatBeg[iCol];
            dsCone->dVecBuffer[dsCone->dualMatIdx[iRow]] -= dsCone->dualStep[iRow] * dRhsVec[iCol];
            /* For each of element not in the diagonal, we have to map it to its symmetric position */
            for ( iRow = dsCone->dualMatBeg[iCol] + 1; iRow < dsCone->dualMatBeg[iCol + 1]; ++iRow ) {
                dsCone->dVecBuffer[dsCone->dualMatIdx[iRow]] -= dsCone->dualStep[iRow] * dRhsVec[iCol];
                dsCone->dVecBuffer[iCol] -= dsCone->dualStep[iRow] * dRhsVec[dsCone->dualMatIdx[iRow]];
            }
        }
    } else {
        fds_symv(dsCone->nCol, -1.0, dsCone->dualStep, dRhsVec, 0.0, dsCone->dVecBuffer);
    }
    
    /* Last step: forward solve */
    HFpLinsysFSolve(dsCone->LTarget, 1, dsCone->dVecBuffer, dRhsVec);
    
    return;
}

static void sdpSparseConeILanczosMultiply( void *cone, double *dLhsVec, double *dRhsVec ) {
    
    hdsdp_cone_sdp_sparse *spCone = (hdsdp_cone_sdp_sparse *) cone;
    
    /* First step: backward solve */
    HFpLinsysBSolve(spCone->LTarget, 1, dLhsVec, dRhsVec);
    
    /* Second step: multiplication. Multiply dRhsVec by -dS */
    if ( spCone->isDualSparse ) {
        HDSDP_ZERO(spCone->dVecBuffer, double, spCone->nCol);
        /* Compute -dS * x. Note that only the lower triangular part is stored in dS. */
        for ( int iCol = 0, iRow = 0; iCol < spCone->nCol; ++iCol ) {
            /* The first element of each column must be on the diagonal due to the identity dual residual */
            iRow = spCone->dualMatBeg[iCol];
            spCone->dVecBuffer[spCone->dualMatIdx[iRow]] -= spCone->dualStep[iRow] * dRhsVec[iCol];
            /* For each of element not in the diagonal, we have to map it to its symmetric position */
            for ( iRow = spCone->dualMatBeg[iCol] + 1; iRow < spCone->dualMatBeg[iCol + 1]; ++iRow ) {
                spCone->dVecBuffer[spCone->dualMatIdx[iRow]] -= spCone->dualStep[iRow] * dRhsVec[iCol];
                spCone->dVecBuffer[iCol] -= spCone->dualStep[iRow] * dRhsVec[spCone->dualMatIdx[iRow]];
            }
        }
    } else {
        fds_symv(spCone->nCol, -1.0, spCone->dualStep, dRhsVec, 0.0, spCone->dVecBuffer);
    }
    
    /* Last step: forward solve */
    HFpLinsysFSolve(spCone->LTarget, 1, spCone->dVecBuffer, dRhsVec);
    
    return;
}

#define SPARSE_EFFICIENCY (1.5)
static int sdpDenseConeIChooseKKTStrategy( int *rowRanks, int *rowSparsity, int *rowPerm, int nRow, int nCol, int iPerm ) {
    
    int bestKKTStrategy = KKT_M1;
    double bestKKTScore = HDSDP_INFINITY;
    
    double KKTScore1 = HDSDP_INFINITY;
    (void) KKTScore1;
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
    KKTScore3 = (double) nCol * SPARSE_EFFICIENCY * rowSparsity[iPerm] + nCubed + SPARSE_EFFICIENCY * nNzAfter + nCubed / nRow;
    KKTScore4 = (double) nCol * SPARSE_EFFICIENCY * rowSparsity[iPerm] + SPARSE_EFFICIENCY * (nCol + 1) * nNzAfter + nCubed / nRow;
    KKTScore5 = SPARSE_EFFICIENCY * (2.0 * SPARSE_EFFICIENCY * rowSparsity[iPerm] + 1) * nNzAfter + nCubed / nRow;
        
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
     M1: Make use of low-rank property. Requires eigen decomposition 
     --
     M2: Make use of both low-rank and sparsity of the coefficient matrices
     --
     
     --------------------
     
     --------------------
     Sparsity based
     --------------------
     M2: Make use of sparsity
     M3: Make use of sparsity
     M4: Make use of sparsity
     
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
    HUtilDescendSortIntByInt(cone->sdpConePerm, rowSparsity, 0, cone->nRow - 1);
    
    /* Choose KKT strategies */
    for ( int iPerm = 0; iPerm < cone->nRow; ++iPerm ) {
        cone->KKTStrategies[iPerm] = \
        sdpDenseConeIChooseKKTStrategy(rowRanks, rowSparsity,
                                       cone->sdpConePerm, cone->nRow, cone->nCol, iPerm);
    }
    
#if 1
    int KKTMethods[5] = {0};
    for ( int iRow = 0; iRow < cone->nRow; ++iRow ) {
        KKTMethods[cone->KKTStrategies[iRow]] += 1;
    }
    printf("    M1: %d M2: %d M3: %d M4: %d M5: %d \n",
           KKTMethods[0], KKTMethods[1], KKTMethods[2], KKTMethods[3], KKTMethods[4]);
    
#endif
    
exit_cleanup:
        
    HDSDP_FREE(rowRanks);
    HDSDP_FREE(rowSparsity);
        
    return retcode;
}

static hdsdp_retcode sdpDenseConeIGetKKTColumnByKKT1( hdsdp_cone_sdp_dense *cone, hdsdp_kkt *kkt, int iKKTCol, int typeKKT ) {
    
    hdsdp_retcode retcode = HDSDP_RETCODE_OK;
    /* Schur strategy M1 is disabled now */
    assert( 0 );
    
exit_cleanup:
    return retcode;
}

static hdsdp_retcode sdpDenseConeIGetKKTColumnByKKT2( hdsdp_cone_sdp_dense *cone, hdsdp_kkt *kkt, int iKKTCol, int typeKKT ) {
    
    hdsdp_retcode retcode = HDSDP_RETCODE_OK;
    assert ( typeKKT != KKT_TYPE_CORRECTOR );
        
    /* Apply Schur complement strategy M2 to build the Schur by exploiting the low rank structure.
    Recall that we need to setup a lower-triangular column of the permuted Schur complement matrix
    and we evaluate
     
     M_{i, j} = trace(A_i * S^-1 A_j * S^-1) for all j >= i (in the permuted index)
     
     As a bi-product, we simultaneously set up
     
     dASinvVec          KKT_TYPE_INFEASIBLE, KKT_TYPE_CORRECTOR
     dASinvRdSinvVec    KKT_TYPE_INFEASIBLE, KKT_TYPE_CORRECTOR
     
     dASinvCSinvVec     KKT_TYPE_HOMOGENEOUS
     dCSinvRdSinv       KKT_TYPE_HOMOGENEOUS
     dCSinvCSinv        KKT_TYPE_HOMOGENEOUS
     dCSinv             KKT_TYPE_HOMOGENEOUS
     
     at different stages of the algorithm.
     
     To implement M2 strategy, we
     
     For each row (i-th)
        Assert A_i = sign * a * a' is rank-one
        Setup v = S^-1 a and S^-1 A_i S^-1 = sign * v * v'
        (Always) Set up ASinvVec[i] = sign * v' * a
        (Always) Set up dASinvRdSinvVec[i] = sign * ||v||^2 * Rd
        (Optionally) Set up dASinvCSinvVec[i]
        For each row (j-th) >= i-th
            (Optionally) Set up M_{i, j} = trace(A_j * sign * v * v') = sign * v' * A_j * v
        End for
     End for
     
     The implementation requires
     1. a buffer vector v to store the result of solve and we use kkt->kktBuffer
     2. an auxiliary vector to compute quadratic form and we use kkt->kktBuffer + nCol
     */
    
    /* Assert that we are dealing with rank-one matrices */
    int iPermKKTCol = cone->sdpConePerm[iKKTCol];
    sdp_coeff *sdpTargetMatrix = cone->sdpRow[iPermKKTCol];
    assert( sdpDataMatGetRank(sdpTargetMatrix) == 1 );
    
    /* Get buffers and auxiliary information ready */
    double dSinAVecSign = 0.0;
    double *dSinvAVecBuffer = kkt->kktBuffer;
    double *dAuxiQuadFormVec = kkt->kktBuffer + cone->nCol;
    
    /* Do we compute the self-dual components */
    int doHSDComputation = 0;
    
    if ( typeKKT == KKT_TYPE_HOMOGENEOUS ) {
        doHSDComputation = 1;
    }
    
    /* Now start setting up the column of the KKT matrix */
    /* Set up the v vector */
    sdpDataMatKKT2SolveRankOne(sdpTargetMatrix, cone->dualFactor, kkt->invBuffer,
                               dSinvAVecBuffer, &dSinAVecSign);
    
    /* Set up ASinvVec[i] = sign * v' * a */
    kkt->dASinvVec[iPermKKTCol] += sdpDataMatKKT2TraceASinv(sdpTargetMatrix, dSinvAVecBuffer);
    /* Set up dASinvRdSinvVec[i] = sign * ||v||^2 * Rd */
    if ( cone->dualResidual ) {
        double dVnorm = nrm2(&cone->nCol, dSinvAVecBuffer, &HIntConstantOne);
        kkt->dASinvRdSinvVec[iPermKKTCol] += dSinAVecSign * cone->dualResidual * dVnorm * dVnorm;
    }
    
    if ( doHSDComputation ) {
        /* Set up dASinvCSinvVec[i] */
        double dASinvCSinvVal = dSinAVecSign * sdpDataMatKKT2QuadForm(cone->sdpObj, dSinvAVecBuffer, dAuxiQuadFormVec);
        kkt->dASinvCSinvVec[iPermKKTCol] += dASinvCSinvVal;
    }
    
    for ( int iRow = iKKTCol; iRow < kkt->nRow; ++iRow ) {
        /* Set up the KKT system */
        int iPermKKTRow = cone->sdpConePerm[iRow];
        
        if ( iPermKKTRow >= iPermKKTCol ) {
            FULL_ENTRY(kkt->kktMatElem, kkt->nRow, iPermKKTRow, iPermKKTCol) += \
            dSinAVecSign * sdpDataMatKKT2QuadForm(cone->sdpRow[iPermKKTRow], dSinvAVecBuffer, dAuxiQuadFormVec);
        } else {
            FULL_ENTRY(kkt->kktMatElem, kkt->nRow, iPermKKTCol, iPermKKTRow) += \
            dSinAVecSign * sdpDataMatKKT2QuadForm(cone->sdpRow[iPermKKTRow], dSinvAVecBuffer, dAuxiQuadFormVec);
        }
    }
        
exit_cleanup:
    return retcode;
}

static hdsdp_retcode sdpDenseConeIGetKKTColumnByKKT3( hdsdp_cone_sdp_dense *cone, hdsdp_kkt *kkt, int iKKTCol, int typeKKT ) {
    
    hdsdp_retcode retcode = HDSDP_RETCODE_OK;
    assert ( typeKKT != KKT_TYPE_CORRECTOR );
    
    /* Apply Schur complement strategy M3 to build the Schur by exploiting the coefficient sparsity.
        
     To implement M3 strategy, we
     
     For each row (i-th)
        Set up B = S^-1 A_i S^-1 explicitly
        (Always) Set up ASinvVec[i] when forming B
        (Always) Set up dASinvRdSinvVec[i] = trace(B) * Rd
        (Optionally) Set up dASinvCSinvVec[i] = trace(B * C)
        For each row (j-th) >= i-th
            (Optionally) Set up M_{i, j} = trace(A_j * B)
        End for
     End for
     
     The implementation requires
     1. a full buffer matrix B to store the dense intermediate result S^-1 A_i S^-1
     2. a full buffer matrix SinvA to store the dense intermediate result S^-1 A_i
     */
    
    int iPermKKTCol = cone->sdpConePerm[iKKTCol];
    sdp_coeff *sdpTargetMatrix = cone->sdpRow[iPermKKTCol];
    
    /* Get buffers and auxiliary information ready */
    double *dSinvASinvBuffer = kkt->kktBuffer;
    double *dAuxiMat = kkt->kktBuffer2;
    
    /* Do we compute the self-dual components */
    int doHSDComputation = 0;
    
    if ( typeKKT == KKT_TYPE_HOMOGENEOUS ) {
        doHSDComputation = 1;
    }
    
    /* Now start setting up the column of the KKT matrix */
    /* Compute the buffer and simultaneously compute trace(A * S^-1) */
    double dASinvVal = sdpDataMatKKT3ComputeSinvASinv(sdpTargetMatrix, cone->dualFactor,
                                                      kkt->invBuffer, dAuxiMat, dSinvASinvBuffer);
    kkt->dASinvVec[iPermKKTCol] += dASinvVal;
    
    /* Compute dASinvRdSinvVec by the diagonal of B */
    if ( cone->dualResidual ) {
        double dTraceBuffer = 0.0;
        for ( int iCol = 0; iCol < cone->nCol; ++iCol ) {
            dTraceBuffer += dSinvASinvBuffer[iCol + iCol * cone->nCol];
        }
        kkt->dASinvRdSinvVec[iPermKKTCol] += dTraceBuffer * cone->dualResidual;
    }
    
    if ( doHSDComputation ) {
        /* Set up dASinvCSinvVec[i] */
        kkt->dASinvCSinvVec[iPermKKTCol] += \
        sdpDataMatKKT3TraceABuffer(cone->sdpObj, dSinvASinvBuffer, dAuxiMat);
    }
    
    for ( int iRow = iKKTCol; iRow < kkt->nRow; ++iRow ) {
        int iPermKKTRow = cone->sdpConePerm[iRow];
        if ( iPermKKTRow >= iPermKKTCol ) {
            FULL_ENTRY(kkt->kktMatElem, kkt->nRow, iPermKKTRow, iPermKKTCol) += \
            sdpDataMatKKT3TraceABuffer(cone->sdpRow[iPermKKTRow], dSinvASinvBuffer, dAuxiMat);
        } else {
            FULL_ENTRY(kkt->kktMatElem, kkt->nRow, iPermKKTCol, iPermKKTRow) += \
            sdpDataMatKKT3TraceABuffer(cone->sdpRow[iPermKKTRow], dSinvASinvBuffer, dAuxiMat);
        }
    }
    
exit_cleanup:
    return retcode;
}

static hdsdp_retcode sdpDenseConeIGetKKTColumnByKKT4( hdsdp_cone_sdp_dense *cone, hdsdp_kkt *kkt, int iKKTCol, int typeKKT ) {
    
    hdsdp_retcode retcode = HDSDP_RETCODE_OK;
    assert ( typeKKT != KKT_TYPE_CORRECTOR );
    
    /* Apply Schur complement strategy M4 to build the Schur by exploiting sparsity.
     
     To implement M4 strategy, we
     
     For each row (i-th)
        Set up B = A_i * S^-1 explicitly
        (Always) Set up dASinvVec[i] = trace(B)
        (Always) Setup dASinvRdSinvVec[i] = Rd * trace(S^-1 * B)
        (Optionally) Set up dASinvCSinvVec[i] = trace(C * S^-1 * B)
        For each row (j-th) >= i-th
            (Optionally) Set up M_{i, j} = trace(A_j * S^-1 * B)
        End for
     End for
     
     The implementation requires
     1. a full buffer matrix B to store the dense intermediate result A_i * S^-1
     2. a buffer vector to help intermediate computation
     */
    
    int iPermKKTCol = cone->sdpConePerm[iKKTCol];
    sdp_coeff *sdpTargetMatrix = cone->sdpRow[iPermKKTCol];
    
    double *dASinvBuffer = kkt->kktBuffer;
    double *dAuxiMat = kkt->kktBuffer2;
    
    /* Do we compute the self-dual components */
    int doHSDComputation = 0;
    
    if ( typeKKT == KKT_TYPE_HOMOGENEOUS ) {
        doHSDComputation = 1;
    }
    
    /* Compute the buffer and simultaneously compute trace(S^-1 * A * S^-1)
       Now that trace(S^-1 * A * S^-1) is computed during the internal operations, we need to skip it when
       there is no dual residual. */
    kkt->dASinvRdSinvVec[iPermKKTCol] += cone->dualResidual * \
    sdpDataMatKKT4ComputeASinv(sdpTargetMatrix, cone->dualFactor, kkt->invBuffer, dAuxiMat, cone->dualResidual, dASinvBuffer);
    
    /* Compute trace(A * S^-1) */
    double dTraceBuffer = 0.0;
    for ( int iCol = 0; iCol < cone->nCol; ++iCol ) {
        dTraceBuffer += dASinvBuffer[iCol + cone->nCol * iCol];
    }
    kkt->dASinvVec[iPermKKTCol] += dTraceBuffer;
    
    if ( doHSDComputation ) {
        /* Set up dASinvCSinvVec[i] */
        kkt->dASinvCSinvVec[iPermKKTCol] += \
        sdpDataMatKKT4TraceASinvBuffer(cone->sdpObj, cone->dualFactor, kkt->invBuffer, dASinvBuffer, dAuxiMat);
    }
    
    for ( int iRow = iKKTCol; iRow < kkt->nRow; ++iRow ) {
        int iPermKKTRow = cone->sdpConePerm[iRow];
        if ( iPermKKTRow >= iPermKKTCol ) {
            FULL_ENTRY(kkt->kktMatElem, kkt->nRow, iPermKKTRow, iPermKKTCol) += \
            sdpDataMatKKT4TraceASinvBuffer(cone->sdpRow[iPermKKTRow], cone->dualFactor, kkt->invBuffer, dASinvBuffer, dAuxiMat);
        } else {
            FULL_ENTRY(kkt->kktMatElem, kkt->nRow, iPermKKTCol, iPermKKTRow) += \
            sdpDataMatKKT4TraceASinvBuffer(cone->sdpRow[iPermKKTRow], cone->dualFactor, kkt->invBuffer, dASinvBuffer, dAuxiMat);
        }
    }
    
exit_cleanup:
    return retcode;
}

static hdsdp_retcode sdpDenseConeIGetKKTColumnByKKT5( hdsdp_cone_sdp_dense *cone, hdsdp_kkt *kkt, int iKKTCol, int typeKKT ) {
    
    hdsdp_retcode retcode = HDSDP_RETCODE_OK;
    assert ( typeKKT != KKT_TYPE_CORRECTOR );
    
    /* Apply Schur complement strategy M5 to build the Schur by exploiting sparsity.
     
     To implement M5 strategy, we
     
     For each row (i-th)
        (Always) Set up dASinvVec[i] = trace(A_i * S^-1)
        (Always) Set up dASinvRdSinvVec[i] = Rd * trace(S^-1 * A_i * S^-1)
        (Optionally) Set up dASinvCSinvVec[i] = trace(C * S^-1 * A_i * S^-1)
        For each row (j-th) >= i-th
            (Optionally) Set up M_{i, j} = trace(A_i * S^-1 * A_j * S^-1)
        End for
     End for
     
     The implementation requires
     1. a full buffer matrix to help intermediate computation
    */
    
    int iPermKKTCol = cone->sdpConePerm[iKKTCol];
    sdp_coeff *sdpTargetMatrix = cone->sdpRow[iPermKKTCol];
    
    double *dAuxiMat = kkt->kktBuffer2;
    
    /* Do we compute the self-dual components */
    int doHSDComputation = 0;
    
    if ( typeKKT == KKT_TYPE_HOMOGENEOUS ) {
        doHSDComputation = 1;
    }
    
    /* Reuse KKT3 routine */
    double dASinvVal = sdpDataMatKKT3TraceABuffer(sdpTargetMatrix, kkt->invBuffer, dAuxiMat);
    kkt->dASinvVec[iPermKKTCol] += dASinvVal;
    
    if ( cone->dualResidual ) {
        double dASinvRdSinvVal = \
        sdpDataMatKKT5SinvADotSinv(sdpTargetMatrix, cone->dualFactor, kkt->invBuffer, dAuxiMat) * cone->dualResidual;
        kkt->dASinvRdSinvVec[iPermKKTCol] += dASinvRdSinvVal;
    }
    
    if ( doHSDComputation ) {
        /* Set up dASinvCSinvVec[i] */
        double dASinvCSinvVal = sdpDataMatKKT5TraceASinvBSinv(sdpTargetMatrix, cone->sdpObj, kkt->invBuffer, dAuxiMat);
        kkt->dASinvCSinvVec[iPermKKTCol] += dASinvCSinvVal;
    }
    
    for ( int iRow = iKKTCol; iRow < kkt->nRow; ++iRow ) {
        int iPermKKTRow = cone->sdpConePerm[iRow];
        if ( iPermKKTRow >= iPermKKTCol ) {
            FULL_ENTRY(kkt->kktMatElem, kkt->nRow, iPermKKTRow, iPermKKTCol) += \
            sdpDataMatKKT5TraceASinvBSinv(sdpTargetMatrix, cone->sdpRow[iPermKKTRow], kkt->invBuffer, dAuxiMat);
        } else {
            FULL_ENTRY(kkt->kktMatElem, kkt->nRow, iPermKKTCol, iPermKKTRow) += \
            sdpDataMatKKT5TraceASinvBSinv(sdpTargetMatrix, cone->sdpRow[iPermKKTRow], kkt->invBuffer, dAuxiMat);
        }
    }
    
exit_cleanup:
    return retcode;
}

static hdsdp_retcode sdpDenseConeIGetHSDComponents( hdsdp_cone_sdp_dense *cone, hdsdp_kkt *kkt ) {
    
    hdsdp_retcode retcode = HDSDP_RETCODE_OK;
    
    /* Special routine for the HSD components. In this routine we set up
     
     kkt->dCSinv
     kkt->dCSinvRdSinv
     kkt->dCSinvCSinv
     
     using M5 strategy
    */
    
    sdp_coeff *sdpTargetMatrix = cone->sdpObj;
    sdp_coeff_type dType = sdpDataMatGetType(sdpTargetMatrix);
    
    if ( dType == SDP_COEFF_ZERO ) {
        goto exit_cleanup;
    } else if ( dType == SDP_COEFF_SPR1 || dType == SDP_COEFF_SPR1 ) {
        /* Apply KKT M5 routines if C is sparse */
        kkt->dCSinvCSinv += sdpDataMatKKT5TraceASinvBSinv(sdpTargetMatrix, sdpTargetMatrix, kkt->invBuffer, kkt->kktBuffer);
        kkt->dCSinv += sdpDataMatKKT3TraceABuffer(sdpTargetMatrix, kkt->invBuffer, kkt->kktBuffer);
        
        if ( cone->dualResidual ) {
            kkt->dCSinvRdSinv += cone->dualResidual * \
            sdpDataMatKKT5SinvADotSinv(sdpTargetMatrix, cone->dualFactor, kkt->invBuffer, kkt->kktBuffer);
        }
        
    } else {
        /* Apply KKT M3 routine is C is dense */
        double *dSinvASinvBuffer = kkt->kktBuffer;
        double *dAuxiMat = kkt->kktBuffer2;
        kkt->dCSinv += sdpDataMatKKT3ComputeSinvASinv(sdpTargetMatrix, cone->dualFactor, kkt->invBuffer, dAuxiMat, dSinvASinvBuffer);
        kkt->dCSinvCSinv += sdpDataMatKKT3TraceABuffer(sdpTargetMatrix, dSinvASinvBuffer, dAuxiMat);
        
        if ( cone->dualResidual ) {
            double dTraceBuffer = 0.0;
            for ( int iCol = 0; iCol < cone->nCol; ++iCol ) {
                dTraceBuffer += dSinvASinvBuffer[iCol + iCol * cone->nCol];
            }
            kkt->dCSinvRdSinv += dTraceBuffer * cone->dualResidual;
        }
    }
    
exit_cleanup:
    return retcode;
}

static hdsdp_retcode sdpDenseConeIGetKKTCorrectorComponents( hdsdp_cone_sdp_dense *cone, hdsdp_kkt *kkt ) {
    
    hdsdp_retcode retcode = HDSDP_RETCODE_OK;
    
    /* Special routine for the corrector. In this routine we set up
       kkt->dASinvVec
       kkt->dASinvRdSinvVec
     */
    
    for ( int iRow = 0; iRow < kkt->nRow; ++iRow ) {
        kkt->dASinvVec[iRow] += sdpDataMatKKT3TraceABuffer(cone->sdpRow[iRow], kkt->invBuffer, kkt->kktBuffer);
    }
    
    if ( cone->dualResidual ) {
        for ( int iRow = 0; iRow < kkt->nRow; ++iRow ) {
            kkt->dASinvRdSinvVec[iRow] += cone->dualResidual * \
            sdpDataMatKKT5SinvADotSinv(cone->sdpRow[iRow], cone->dualFactor, kkt->invBuffer, kkt->kktBuffer);
        }
    }
 
exit_cleanup:
    return retcode;
}

static hdsdp_retcode sdpSparseConeIGetKKTColumnByKKT2( hdsdp_cone_sdp_sparse *cone, hdsdp_kkt *kkt, int iKKTNzCol, int typeKKT ) {
    
    hdsdp_retcode retcode = HDSDP_RETCODE_OK;
    assert( typeKKT != KKT_TYPE_CORRECTOR );
    
    /* Implement the sparse version of KKT strategy M2. There is no permutation for sparse KKT system
       Refer to the details in function sdpDenseConeIGetKKTColumnByKKT2 */
    
    int iKKTCol = cone->rowIdx[iKKTNzCol];
    sdp_coeff *sdpTargetMatrix = cone->sdpRow[iKKTNzCol];
    assert( sdpDataMatGetRank(sdpTargetMatrix) == 1 );
    
    double dSinAVecSign = 0.0;
    double *dSinvAVecBuffer = kkt->kktBuffer;
    double *dAuxiQuadFormVec = kkt->kktBuffer + cone->nCol;
    
    int doHSDComputation = 0;
    
    if ( typeKKT == KKT_TYPE_HOMOGENEOUS ) {
        doHSDComputation = 1;
    }
    
    sdpDataMatKKT2SolveRankOne(sdpTargetMatrix, cone->dualFactor, kkt->invBuffer,
                               dSinvAVecBuffer, &dSinAVecSign);
    
    kkt->dASinvVec[iKKTCol] += sdpDataMatKKT2TraceASinv(sdpTargetMatrix, dSinvAVecBuffer);
    
    if ( cone->dualResidual ) {
        double dVnorm = nrm2(&cone->nCol, dSinvAVecBuffer, &HIntConstantOne);
        kkt->dASinvRdSinvVec[iKKTCol] += dSinAVecSign * cone->dualResidual * dVnorm * dVnorm;
    }
    
    if ( doHSDComputation ) {
        /* Set up dASinvCSinvVec[i] */
        kkt->dASinvCSinvVec[iKKTCol] += \
        dSinAVecSign * sdpDataMatKKT2QuadForm(cone->sdpObj, dSinvAVecBuffer, dAuxiQuadFormVec);
    }
    
    int iKKTPosition = 0;
    for ( int iRowElem = iKKTNzCol; iRowElem < cone->nRowElem; ++iRowElem ) {
        if ( kkt->isKKTSparse ) {
            iKKTPosition = PACK_ENTRY(cone->kktMapping, cone->nRowElem, iRowElem, iKKTNzCol);
        } else {
            iKKTPosition = cone->rowIdx[iRowElem] + iKKTCol * kkt->nRow;
        }
        kkt->kktMatElem[iKKTPosition] += \
        dSinAVecSign * sdpDataMatKKT2QuadForm(cone->sdpRow[iRowElem], dSinvAVecBuffer, dAuxiQuadFormVec);
    }
    
exit_cleanup:
    return retcode;
}

static hdsdp_retcode sdpSparseConeIGetKKTColumnByKKT3( hdsdp_cone_sdp_sparse *cone, hdsdp_kkt *kkt, int iKKTNzCol, int typeKKT ) {
    
    hdsdp_retcode retcode = HDSDP_RETCODE_OK;
    assert( typeKKT != KKT_TYPE_CORRECTOR );
    
    /* Implement the sparse version of KKT strategy M3. There is no permutation for sparse KKT system
       Refer to the details in function sdpDenseConeIGetKKTColumnByKKT3 */
    
    int iKKTCol = cone->rowIdx[iKKTNzCol];
    sdp_coeff *sdpTargetMatrix = cone->sdpRow[iKKTNzCol];
    
    double *dSinvASinvBuffer = kkt->kktBuffer;
    double *dAuxiMat = kkt->kktBuffer2;
    
    int doHSDComputation = 0;
    
    if ( typeKKT == KKT_TYPE_HOMOGENEOUS ) {
        doHSDComputation = 1;
    }
    
    kkt->dASinvVec[iKKTCol] += \
    sdpDataMatKKT3ComputeSinvASinv(sdpTargetMatrix, cone->dualFactor,
                                   kkt->invBuffer, dAuxiMat, dSinvASinvBuffer);
    
    if ( cone->dualResidual ) {
        double dTraceBuffer = 0.0;
        for ( int iCol = 0; iCol < cone->nCol; ++iCol ) {
            dTraceBuffer += dSinvASinvBuffer[iCol + iCol * cone->nCol];
        }
        kkt->dASinvRdSinvVec[iKKTCol] += dTraceBuffer * cone->dualResidual;
    }
    
    if ( doHSDComputation ) {
        /* Set up dASinvCSinvVec[i] */
        kkt->dASinvCSinvVec[iKKTCol] += \
        sdpDataMatKKT3TraceABuffer(cone->sdpObj, dSinvASinvBuffer, dAuxiMat);
    }
    
    int iKKTPosition = 0;
    for ( int iRowElem = iKKTNzCol; iRowElem < cone->nRowElem; ++iRowElem ) {
        if ( kkt->isKKTSparse ) {
            iKKTPosition = PACK_ENTRY(cone->kktMapping, cone->nRowElem, iRowElem, iKKTNzCol);
        } else {
            iKKTPosition = cone->rowIdx[iRowElem] + iKKTCol * kkt->nRow;
        }
        kkt->kktMatElem[iKKTPosition] += \
        sdpDataMatKKT3TraceABuffer(cone->sdpRow[iRowElem], dSinvASinvBuffer, dAuxiMat);
    }
    
exit_cleanup:
    return retcode;
}

static hdsdp_retcode sdpSparseConeIGetKKTColumnByKKT4( hdsdp_cone_sdp_sparse *cone, hdsdp_kkt *kkt, int iKKTNzCol, int typeKKT ) {
    
    hdsdp_retcode retcode = HDSDP_RETCODE_OK;
    assert ( typeKKT != KKT_TYPE_CORRECTOR );
    
    /* In general this method should not be invoked. But we still implement it for cross-validation  */
    
    int iKKTCol = cone->rowIdx[iKKTNzCol];
    sdp_coeff *sdpTargetMatrix = cone->sdpRow[iKKTNzCol];
    
    double *dASinvBuffer = kkt->kktBuffer;
    double *dAuxiMat = kkt->kktBuffer2;
    
    /* Do we compute the self-dual components */
    int doHSDComputation = 0;
    
    if ( typeKKT == KKT_TYPE_HOMOGENEOUS ) {
        doHSDComputation = 1;
    }
    
    kkt->dASinvRdSinvVec[iKKTCol] += cone->dualResidual * \
    sdpDataMatKKT4ComputeASinv(sdpTargetMatrix, cone->dualFactor, kkt->invBuffer, dAuxiMat, cone->dualResidual, dASinvBuffer);
    
    /* Compute trace(A * S^-1) */
    double dTraceBuffer = 0.0;
    for ( int iCol = 0; iCol < cone->nCol; ++iCol ) {
        dTraceBuffer += dASinvBuffer[iCol + cone->nCol * iCol];
    }
    kkt->dASinvVec[iKKTCol] += dTraceBuffer;
    
    if ( doHSDComputation ) {
        /* Set up dASinvCSinvVec[i] */
        kkt->dASinvCSinvVec[iKKTCol] += \
        sdpDataMatKKT4TraceASinvBuffer(cone->sdpObj, cone->dualFactor, kkt->invBuffer, dASinvBuffer, dAuxiMat);
    }
    
    int iKKTPosition = 0;
    for ( int iRowElem = iKKTNzCol; iRowElem < cone->nRowElem ; ++iRowElem ) {
        if ( kkt->isKKTSparse ) {
            iKKTPosition = PACK_ENTRY(cone->kktMapping, cone->nRowElem, iRowElem, iKKTNzCol);
        } else {
            iKKTPosition = cone->rowIdx[iRowElem] + iKKTCol * kkt->nRow;
        }
        kkt->kktMatElem[iKKTPosition] += \
        sdpDataMatKKT4TraceASinvBuffer(cone->sdpRow[iRowElem], cone->dualFactor, kkt->invBuffer, dASinvBuffer, dAuxiMat);
    }
    
exit_cleanup:
    return retcode;
}

static hdsdp_retcode sdpSparseConeIGetKKTColumnByKKT5( hdsdp_cone_sdp_sparse *cone, hdsdp_kkt *kkt, int iKKTNzCol, int typeKKT ) {
    
    hdsdp_retcode retcode = HDSDP_RETCODE_OK;
    assert( typeKKT != KKT_TYPE_CORRECTOR );
    
    /* Implement the sparse version of KKT strategy M5. There is no permutation for sparse KKT system
       Refer to the details in function sdpDenseConeIGetKKTColumnByKKT5 */
    
    int iKKTCol = cone->rowIdx[iKKTNzCol];
    sdp_coeff *sdpTargetMatrix = cone->sdpRow[iKKTNzCol];
    
    double *dAuxiMat = kkt->kktBuffer2;
    
    int doHSDComputation = 0;
    
    if ( typeKKT == KKT_TYPE_HOMOGENEOUS ) {
        doHSDComputation = 1;
    }
    
    kkt->dASinvVec[iKKTCol] += \
    sdpDataMatKKT3TraceABuffer(sdpTargetMatrix, kkt->invBuffer, dAuxiMat);
    
    if ( cone->dualResidual ) {
        kkt->dASinvRdSinvVec[iKKTCol] += \
        sdpDataMatKKT5SinvADotSinv(sdpTargetMatrix, cone->dualFactor, kkt->invBuffer, dAuxiMat) * cone->dualResidual;
    }
    
    if ( doHSDComputation ) {
        /* Set up dASinvCSinvVec[i] */
        kkt->dASinvCSinvVec[iKKTCol] += \
        sdpDataMatKKT5TraceASinvBSinv(sdpTargetMatrix, cone->sdpObj, kkt->invBuffer, dAuxiMat);
    }
    
    int iKKTPosition = 0;
    for ( int iRowElem = iKKTNzCol; iRowElem < cone->nRowElem; ++iRowElem ) {
        
        if ( kkt->isKKTSparse ) {
            iKKTPosition = PACK_ENTRY(cone->kktMapping, cone->nRowElem, iRowElem, iKKTNzCol);
        } else {
            iKKTPosition = cone->rowIdx[iRowElem] + iKKTCol * kkt->nRow;
        }
        
        kkt->kktMatElem[iKKTPosition] += \
        sdpDataMatKKT5TraceASinvBSinv(sdpTargetMatrix, cone->sdpRow[iRowElem], kkt->invBuffer, dAuxiMat);
    }
    
exit_cleanup:
    return retcode;
}

static hdsdp_retcode sdpSparseConeIGetHSDComponents( hdsdp_cone_sdp_sparse *cone, hdsdp_kkt *kkt ) {
    
    hdsdp_retcode retcode = HDSDP_RETCODE_OK;
    /* Special routine for the HSD components */
    sdp_coeff *sdpTargetMatrix = cone->sdpObj;
    sdp_coeff_type dType = sdpDataMatGetType(sdpTargetMatrix);
    
    if ( dType == SDP_COEFF_ZERO ) {
        goto exit_cleanup;
    } else if ( dType == SDP_COEFF_SPR1 || dType == SDP_COEFF_SPR1 ) {
        /* Apply KKT M5 routines if C is sparse */
        kkt->dCSinvCSinv += sdpDataMatKKT5TraceASinvBSinv(sdpTargetMatrix, sdpTargetMatrix, kkt->invBuffer, kkt->kktBuffer);
        kkt->dCSinv += sdpDataMatKKT3TraceABuffer(sdpTargetMatrix, kkt->invBuffer, kkt->kktBuffer);
        
        if ( cone->dualResidual ) {
            kkt->dCSinvRdSinv += cone->dualResidual * \
            sdpDataMatKKT5SinvADotSinv(sdpTargetMatrix, cone->dualFactor, kkt->invBuffer, kkt->kktBuffer);
        }
        
    } else {
        /* Apply KKT M3 routine is C is dense */
        double *dSinvASinvBuffer = kkt->kktBuffer;
        double *dAuxiMat = kkt->kktBuffer2;
        kkt->dCSinv += sdpDataMatKKT3ComputeSinvASinv(sdpTargetMatrix, cone->dualFactor, kkt->invBuffer, dAuxiMat, dSinvASinvBuffer);
        kkt->dCSinvCSinv += sdpDataMatKKT3TraceABuffer(sdpTargetMatrix, dSinvASinvBuffer, dAuxiMat);
        
        if ( cone->dualResidual ) {
            double dTraceBuffer = 0.0;
            for ( int iCol = 0; iCol < cone->nCol; ++iCol ) {
                dTraceBuffer += dSinvASinvBuffer[iCol + iCol * cone->nCol];
            }
            kkt->dCSinvRdSinv += dTraceBuffer * cone->dualResidual;
        }
    }
 
exit_cleanup:
    return retcode;
}

static hdsdp_retcode sdpSparseConeIGetKKTCorrectorComponents( hdsdp_cone_sdp_sparse *cone, hdsdp_kkt *kkt ) {
    
    hdsdp_retcode retcode = HDSDP_RETCODE_OK;
    /* Special routine for the corrector */
    
    for ( int iRowElem = 0; iRowElem < cone->nRowElem; ++iRowElem ) {
        kkt->dASinvVec[cone->rowIdx[iRowElem]] += \
        sdpDataMatKKT3TraceABuffer(cone->sdpRow[iRowElem], kkt->invBuffer, kkt->kktBuffer);
    }
    
    if ( cone->dualResidual ) {
        for ( int iRowElem = 0; iRowElem < cone->nRowElem; ++iRowElem ) {
            kkt->dASinvRdSinvVec[cone->rowIdx[iRowElem]] += cone->dualResidual * \
            sdpDataMatKKT5SinvADotSinv(cone->sdpRow[iRowElem], cone->dualFactor, kkt->invBuffer, kkt->kktBuffer);
        }
    }
 
exit_cleanup:
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
    HDSDP_INIT(cone->sdpConePerm, int, cone->nRow);
    HDSDP_INIT(cone->KKTStrategies, int, cone->nRow);
    
    cone->coneKKTNnz = cone->nRow;
    cone->coneKKTNnz = cone->coneKKTNnz * cone->coneKKTNnz;
        
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
    
    /* Compute nonzeros in the Schur complement */
    cone->coneKKTNnz = (int64_t) cone->nRowElem;
    cone->coneKKTNnz = cone->coneKKTNnz * cone->coneKKTNnz;
    
    HDSDP_INIT(cone->kktMapping, int, PACK_NNZ(cone->nRowElem));
    
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
    HDSDP_CALL(HFpLinsysSymbolic(cone->dualChecker, cone->dualMatBeg, cone->dualMatIdx));
    
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
    HDSDP_CALL(HFpLinsysSymbolic(cone->dualChecker, cone->dualMatBeg, cone->dualMatIdx));
    
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

extern void sdpDenseConeScal( hdsdp_cone_sdp_dense *cone, double dScal ) {
    
    sdpDataMatScal(cone->sdpObj, dScal);
    return;
}

extern void sdpSparseConeScal( hdsdp_cone_sdp_sparse *cone, double dScal ) {
    
    sdpDataMatScal(cone->sdpObj, dScal);
    return;
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

extern hdsdp_retcode sdpDenseConeRatioTestImpl( hdsdp_cone_sdp_dense *cone, double barHsdTauStep, double *rowDualStep, double dAdaRatio,
                                                int whichBuffer, double *maxStep ) {
    /* Given conic quantity dy, get maximum distance from the cone
       Given dTau and dy, we first evaluate
     
            B <- dResiCoef * I + dACoefScal * A' * y + C * dCCoef
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
    
    /* Choose Lanczos target buffer */
    if ( whichBuffer == BUFFER_DUALVAR ) {
        cone->LTarget = cone->dualFactor;
    } else {
        cone->LTarget = cone->dualChecker;
    }
    
    /* Finished with setting up dS. Now do Lanczos computation */
    if ( cone->nCol == 1 ) {
        if ( whichBuffer == BUFFER_DUALVAR ) {
            *maxStep = (cone->dualStep[0] > 0.0) ? HDSDP_INFINITY : (-cone->dualMatElem[0] / cone->dualStep[0]);
        } else {
            *maxStep = (cone->dualStep[0] > 0.0) ? HDSDP_INFINITY : (-cone->dualCheckerElem[0] / cone->dualStep[0]);
        }
        
    }  else {
        HDSDP_CALL(HLanczosSolve(cone->Lanczos, NULL, maxStep));
    }
    
exit_cleanup:
    return retcode;
}

extern hdsdp_retcode sdpSparseConeRatioTestImpl( hdsdp_cone_sdp_sparse *cone, double barHsdTauStep, double *rowDualStep, double dAdaRatio,
                                                 int whichBuffer, double *maxStep ) {
    
    hdsdp_retcode retcode = HDSDP_RETCODE_OK;
    sdpSparseConeIUpdateBuffer(cone, barHsdTauStep, -1.0, rowDualStep, dAdaRatio * cone->dualResidual, BUFFER_DUALSTEP);
    
    /* Choose Lanczos target buffer */
    if ( whichBuffer == BUFFER_DUALVAR ) {
        cone->LTarget = cone->dualFactor;
    } else {
        cone->LTarget = cone->dualChecker;
    }
    
    if ( cone->nCol == 1 ) {
        if ( whichBuffer == BUFFER_DUALVAR ) {
            *maxStep = (cone->dualStep[0] > 0.0) ? HDSDP_INFINITY : (-cone->dualMatElem[0] / cone->dualStep[0]);
        } else {
            *maxStep = (cone->dualStep[0] > 0.0) ? HDSDP_INFINITY : (-cone->dualCheckerElem[0] / cone->dualStep[0]);
        }
        
    } else {
        HDSDP_CALL(HLanczosSolve(cone->Lanczos, NULL, maxStep));
    }
    
exit_cleanup:
    return retcode;
}

extern int sdpDenseConeGetDim( hdsdp_cone_sdp_dense *cone ) {
    
    return cone->nCol;
}

extern int sdpSparseConeGetDim( hdsdp_cone_sdp_sparse *cone ) {
    
    return cone->nCol;
}

extern hdsdp_retcode sdpDenseConeGetKKT( hdsdp_cone_sdp_dense *cone, void *kkt, int typeKKT ) {
    /*
     Set up the KKT system for a dense cone.
     Now that the cone is dense, the Schur complement is dense.
     
     We employ five (acturally four) pre-defined KKT strategies to set up rows (columns) of the KKT system
     
     with M_{i, j} = Tr(A_i, S^-1 A_j S^-1)
     
     and we simultaneously set up the RHS and LHS of the KKT system
     
     [ mu * A S^-1 A            -b - mu * A S^-1 C S^-1          ] [ dy ] = [ b * tau - mu * A S^-1 + mu * A S^-1 Rd S^-1              ]
     [ b - mu * A S^-1 C S^-1   mu * C S^-1 C S^-1 + mu * tau^-2 ] [dtau] = [ -b' * y + mu * tau^-1 + mu * C S^-1 - mu * C S^-1 Rd S^-1]
     
     */
    
    hdsdp_retcode retcode = HDSDP_RETCODE_OK;
    hdsdp_kkt *Hkkt = (hdsdp_kkt *) kkt;
    
    /* Prepare inverse */
    HFpLinsysInvert(cone->dualFactor, Hkkt->invBuffer, Hkkt->kktBuffer);
    
    if ( typeKKT == KKT_TYPE_CORRECTOR ) {
        HDSDP_CALL(sdpDenseConeIGetKKTCorrectorComponents(cone, kkt));
        goto exit_cleanup;
    }
    
    /* Get quantity for estimating primal objective */
    if ( cone->dualResidual ) {
        for ( int iCol = 0; iCol < cone->nCol; ++iCol ) {
            Hkkt->dTraceSinv += Hkkt->invBuffer[iCol + iCol * cone->nCol];
        }
    }
    
    for ( int iKKTCol = 0; iKKTCol < cone->nRow; ++iKKTCol ) {
        
        if ( sdpDataMatGetType(cone->sdpRow[cone->sdpConePerm[iKKTCol]]) == SDP_COEFF_ZERO ) {
            /* Continue if the target matrix is zero */
            continue;
        }
        
#if 0
        if ( cone->sdpConePerm[iKKTCol] == 0 ) {
            printf("Here. \n");
        }
#endif
        
        switch (cone->KKTStrategies[iKKTCol]) {
            case KKT_M1:
                HDSDP_CALL(sdpDenseConeIGetKKTColumnByKKT1(cone, Hkkt, iKKTCol, typeKKT));
                break;
            case KKT_M2:
                HDSDP_CALL(sdpDenseConeIGetKKTColumnByKKT2(cone, Hkkt, iKKTCol, typeKKT));
                break;
            case KKT_M3:
                HDSDP_CALL(sdpDenseConeIGetKKTColumnByKKT3(cone, Hkkt, iKKTCol, typeKKT));
                break;
            case KKT_M4:
                HDSDP_CALL(sdpDenseConeIGetKKTColumnByKKT4(cone, Hkkt, iKKTCol, typeKKT));
                break;
            case KKT_M5:
                HDSDP_CALL(sdpDenseConeIGetKKTColumnByKKT5(cone, Hkkt, iKKTCol, typeKKT));
                break;
            default:
                printf("Invalid KKT strategy. \n");
                retcode = HDSDP_RETCODE_FAILED;
                break;
        }
    }
    
    if ( typeKKT == KKT_TYPE_HOMOGENEOUS ) {
        HDSDP_CALL(sdpDenseConeIGetHSDComponents(cone, Hkkt));
    }
    
exit_cleanup:
    return retcode;
}

extern hdsdp_retcode sdpDenseConeGetKKTByFixedStrategy( hdsdp_cone_sdp_dense *cone, void *kkt, int typeKKT, int kktStrategy ) {
    /* Set up the Schur complement with a fixed strategy */
    hdsdp_retcode retcode = HDSDP_RETCODE_OK;
    hdsdp_kkt *Hkkt = (hdsdp_kkt *) kkt;
    
    /* Prepare inverse */
    HFpLinsysInvert(cone->dualFactor, Hkkt->invBuffer, Hkkt->kktBuffer);
    
    if ( typeKKT == KKT_TYPE_CORRECTOR ) {
        HDSDP_CALL(sdpDenseConeIGetKKTCorrectorComponents(cone, kkt));
        goto exit_cleanup;
    }
    
    /* Get quantity for estimating primal objective */
    if ( cone->dualResidual ) {
        for ( int iCol = 0; iCol < cone->nCol; ++iCol ) {
            Hkkt->dTraceSinv += Hkkt->invBuffer[iCol + iCol * cone->nCol];
        }
    }
    
    for ( int iKKTCol = 0; iKKTCol < cone->nRow; ++iKKTCol ) {
        
        if ( sdpDataMatGetType(cone->sdpRow[cone->sdpConePerm[iKKTCol]]) == SDP_COEFF_ZERO ) {
            continue;
        }
        
        switch (kktStrategy) {
            case KKT_M1:
                HDSDP_CALL(sdpDenseConeIGetKKTColumnByKKT1(cone, Hkkt, iKKTCol, typeKKT));
                break;
            case KKT_M2:
                HDSDP_CALL(sdpDenseConeIGetKKTColumnByKKT2(cone, Hkkt, iKKTCol, typeKKT));
                break;
            case KKT_M3:
                HDSDP_CALL(sdpDenseConeIGetKKTColumnByKKT3(cone, Hkkt, iKKTCol, typeKKT));
                break;
            case KKT_M4:
                HDSDP_CALL(sdpDenseConeIGetKKTColumnByKKT4(cone, Hkkt, iKKTCol, typeKKT));
                break;
            case KKT_M5:
                HDSDP_CALL(sdpDenseConeIGetKKTColumnByKKT5(cone, Hkkt, iKKTCol, typeKKT));
                break;
            default:
                printf("Invalid KKT strategy. \n");
                retcode = HDSDP_RETCODE_FAILED;
                break;
        }
    }
    
    if ( typeKKT == KKT_TYPE_HOMOGENEOUS ) {
        HDSDP_CALL(sdpDenseConeIGetHSDComponents(cone, Hkkt));
    }
    
exit_cleanup:
    return retcode;
}

extern hdsdp_retcode sdpSparseConeGetKKT( hdsdp_cone_sdp_sparse *cone, void *kkt, int typeKKT ) {
    
    /* When setting up the KKT system for a sparse cone,
     There is no pre-determined strategy and we choose from M2, M3 and M5 based on the sparsity pattern of the current matrix
     */
    
    hdsdp_retcode retcode = HDSDP_RETCODE_OK;
    hdsdp_kkt *Hkkt = (hdsdp_kkt *) kkt;
    
    /* Prepare inverse */
    HFpLinsysInvert(cone->dualFactor, Hkkt->invBuffer, Hkkt->kktBuffer);
    
    if ( typeKKT == KKT_TYPE_CORRECTOR ) {
        HDSDP_CALL(sdpSparseConeIGetKKTCorrectorComponents(cone, kkt));
        goto exit_cleanup;
    }
    
    /* Get quantity for estimating primal objective */
    if ( cone->dualResidual ) {
        for ( int iCol = 0; iCol < cone->nCol; ++iCol ) {
            Hkkt->dTraceSinv += Hkkt->invBuffer[iCol + iCol * cone->nCol];
        }
    }
    
    for ( int iKKTNzCol = 0; iKKTNzCol < cone->nRowElem; ++iKKTNzCol ) {
        
        int KKTStrategy = KKT_M1;
        (void) KKTStrategy;
        
#if 0
        if ( cone->rowIdx[iKKTNzCol] == 3 ) {
            printf("Here. \n");
        }
#endif
        
        sdp_coeff_type dataType = sdpDataMatGetType(cone->sdpRow[iKKTNzCol]);
        
        switch (dataType) {
            case SDP_COEFF_SPARSE:
                /* Use M5 */
                HDSDP_CALL(sdpSparseConeIGetKKTColumnByKKT5(cone, Hkkt, iKKTNzCol, typeKKT));
                break;
            case SDP_COEFF_DENSE:
                HDSDP_CALL(sdpSparseConeIGetKKTColumnByKKT3(cone, Hkkt, iKKTNzCol, typeKKT));
                break;
            case SDP_COEFF_DSR1:
                HDSDP_CALL(sdpSparseConeIGetKKTColumnByKKT2(cone, Hkkt, iKKTNzCol, typeKKT));
                break;
            case SDP_COEFF_SPR1:
                HDSDP_CALL(sdpSparseConeIGetKKTColumnByKKT5(cone, Hkkt, iKKTNzCol, typeKKT));
                break;
            default:
                assert( 0 );
                break;
        }
    }
    
    if ( typeKKT == KKT_TYPE_HOMOGENEOUS ) {
        HDSDP_CALL(sdpSparseConeIGetHSDComponents(cone, kkt));
    }
    
exit_cleanup:
    return retcode;
}

extern hdsdp_retcode sdpSparseConeGetKKTByFixedStrategy( hdsdp_cone_sdp_sparse *cone, void *kkt, int typeKKT, int kktStrategy ) {
    
    /* When setting up the KKT system for a sparse cone,
     There is no pre-determined strategy and we choose from M2, M3 and M5 based on the sparsity pattern of the current matrix
     */
    
    hdsdp_retcode retcode = HDSDP_RETCODE_OK;
    hdsdp_kkt *Hkkt = (hdsdp_kkt *) kkt;
    
    /* Prepare inverse */
    HFpLinsysInvert(cone->dualFactor, Hkkt->invBuffer, Hkkt->kktBuffer);
    
    /* Get quantity for estimating primal objective */
    if ( cone->dualResidual ) {
        for ( int iCol = 0; iCol < cone->nCol; ++iCol ) {
            Hkkt->dTraceSinv += Hkkt->invBuffer[iCol + iCol * cone->nCol];
        }
    }
    
    if ( typeKKT == KKT_TYPE_CORRECTOR ) {
        HDSDP_CALL(sdpSparseConeIGetKKTCorrectorComponents(cone, kkt));
        goto exit_cleanup;
    }
    
    for ( int iKKTNzCol = 0; iKKTNzCol < cone->nRowElem; ++iKKTNzCol ) {
                
        switch (kktStrategy) {
            case KKT_M1:
                assert( 0 );
                break;
            case KKT_M2:
                HDSDP_CALL(sdpSparseConeIGetKKTColumnByKKT2(cone, Hkkt, iKKTNzCol, typeKKT));
                break;
            case KKT_M3:
                HDSDP_CALL(sdpSparseConeIGetKKTColumnByKKT3(cone, Hkkt, iKKTNzCol, typeKKT));
                break;
            case KKT_M4:
                HDSDP_CALL(sdpSparseConeIGetKKTColumnByKKT4(cone, Hkkt, iKKTNzCol, typeKKT));
                break;
            case KKT_M5:
                HDSDP_CALL(sdpSparseConeIGetKKTColumnByKKT5(cone, Hkkt, iKKTNzCol, typeKKT));
                break;
            default:
                printf("Invalid KKT strategy. \n");
                retcode = HDSDP_RETCODE_FAILED;
                break;
        }
    }
    
    if ( typeKKT == KKT_TYPE_HOMOGENEOUS ) {
        HDSDP_CALL(sdpSparseConeIGetHSDComponents(cone, kkt));
    }
    
exit_cleanup:
    return retcode;
}

extern int64_t sdpDenseConeGetSymNnzImpl( hdsdp_cone_sdp_dense *cone ) {
    
    return cone->coneKKTNnz;
}

extern int64_t sdpSparseConeGetSymNnzImpl( hdsdp_cone_sdp_sparse *cone ) {
    
    return cone->coneKKTNnz;
}

extern void sdpDenseConeAddSymNnzImpl( hdsdp_cone_sdp_dense *cone, int iCol, int *schurMatCol ) {
    
    /* Accumulate the sparsity pattern in a column of the Schur complement
       For a dense cone the whole column is dense
     */
    for ( int iRow = 0; iRow < cone->nRow; ++iRow ) {
        schurMatCol[iRow] = 1;
    }
    
    return;
}

extern void sdpSparseConeAddSymNnzImpl( hdsdp_cone_sdp_sparse *cone, int iCol, int *schurMatCol ) {
    
    /* When counting the sparsity pattern for a sparse SDP cone, we allocate a column buffer
       and iterate through the columns of the Schur complement to extract the position of each nonzero
       of the current cone in the aggregated sparsity pattern of the Schur complement
     
       e.g., we have nnz of the cone as [0, 3], which corresponds to
     
            [ 1 0 0 0 ]
            [ 0 0 0 0 ]
            [ 0 0 0 0 ]
            [ 1 0 0 1 ]
     
       in the final Schur complement matrix (we consider only the lower-triangular part) .
     
       Assume that there is another cone of sparsity pattern [0, 1], then the aggregated sparsity pattern would by
            
            [ 1 0 0 0 ]  +  [ 1 0 0 0 ]  =  [ 1 0 0 0 ]
            [ 0 0 0 0 ]  +  [ 1 1 0 0 ]  =  [ 1 1 0 0 ]
            [ 0 0 0 0 ]  +  [ 0 0 0 0 ]  =  [ 0 0 0 0 ]
            [ 1 0 0 1 ]  +  [ 0 0 0 0 ]  =  [ 1 0 0 1 ]
     
       When the number of rows is large. We have to resort to a column-wise strategy.
       i.e., we create a buffer of size 4, and in each cone we record the number of elements
       that have been counted in the aggregated sparsity pattern (cone->iKKTCounted)
       
       For each column
            We zero-out the column buffer
            For each cone
                If column index == rowIdx[iKKTCounted]
                    We add the sparsity pattern of the cone to the column buffer
                End if
            End for
       End for
       
      After getting the sparsity of the column, we count the position of each nonzero element in the final sparse data structure.
      e.g., after the first iteration the buffer would be
       [ 1, 1, 0, 1 ] and we turn it into [ 0, 1, -1, 2 ]
     after the second iteration the buffer would be
       [ 0, 1, 0, 0 ] and we turn it into [ -1, 3, -1, -1 ]
     
     Then we let the cone retrieve their column sparsity from the transformed buffer
     */
    
    if ( cone->rowIdx[cone->iKKTCounted] == iCol ) {
        /* We get nonzero in this column */
        for ( int iElem = cone->iKKTCounted; iElem < cone->nRowElem; ++iElem ) {
            /* We only need lower-triangular part */
            schurMatCol[cone->rowIdx[iElem]] = 1;
        }
    }
    
    return;
}

extern void sdpDenseConeGetSymMapping( hdsdp_cone_sdp_dense *cone, int iCol, int *schurMatCol ) {
    
    /* Whenever there is a dense cone, this method will not be invoked */
    assert( 0 );
    return;
}

extern void sdpSparseConeGetSymMapping( hdsdp_cone_sdp_sparse *cone, int iCol, int *schurMatCol ) {
    
    /* Now we extract the sparsity pattern */
    int kktShift = ((2 * cone->nRowElem - cone->iKKTCounted + 1) * cone->iKKTCounted) / 2;
    int *kktStart = cone->kktMapping + kktShift;
    
    if ( cone->rowIdx[cone->iKKTCounted] == iCol ) {
        for ( int iElem = cone->iKKTCounted; iElem < cone->nRowElem; ++iElem ) {
            kktStart[iElem - cone->iKKTCounted] = schurMatCol[cone->rowIdx[iElem]];
        }
        cone->iKKTCounted += 1;
    }
    
    return;
}

extern hdsdp_retcode sdpDenseConeInteriorCheck( hdsdp_cone_sdp_dense *cone, double barHsdTau, double *rowDual, int *isInterior ) {
    
    hdsdp_retcode retcode = HDSDP_RETCODE_OK;
    sdpDenseConeUpdateImpl(cone, barHsdTau, rowDual);
    HDSDP_CALL(HFpLinsysPsdCheck(cone->dualFactor, cone->dualMatBeg, cone->dualMatIdx, cone->dualMatElem, isInterior));
    
exit_cleanup:
    return retcode;
}

extern hdsdp_retcode sdpSparseConeInteriorCheck( hdsdp_cone_sdp_sparse *cone, double barHsdTau, double *rowDual, int *isInterior ) {
    
    hdsdp_retcode retcode = HDSDP_RETCODE_OK;
    sdpSparseConeUpdateImpl(cone, barHsdTau, rowDual);
    HDSDP_CALL(HFpLinsysPsdCheck(cone->dualFactor, cone->dualMatBeg, cone->dualMatIdx, cone->dualMatElem, isInterior));
    
exit_cleanup:
    return retcode;
}

extern hdsdp_retcode sdpDenseConeInteriorCheckExpert( hdsdp_cone_sdp_dense *cone, double dCCoef, double dACoefScal, double *dACoef, double dEyeCoef,
                                                      int whichBuffer, int *isInterior ) {
    
    /* The expert routine for dense cone interior point check */
    hdsdp_retcode retcode = HDSDP_RETCODE_OK;
    sdpDenseConeIUpdateBuffer(cone, dCCoef, dACoefScal, dACoef, dEyeCoef, whichBuffer);
    
    if ( whichBuffer == BUFFER_DUALVAR ) {
        HDSDP_CALL(HFpLinsysPsdCheck(cone->dualFactor, cone->dualMatBeg, cone->dualMatIdx, cone->dualMatElem, isInterior));
    } else {
        HDSDP_CALL(HFpLinsysPsdCheck(cone->dualChecker, cone->dualMatBeg, cone->dualMatIdx, cone->dualCheckerElem, isInterior));
    }
    
exit_cleanup:
    return retcode;
}

extern hdsdp_retcode sdpSparseConeInteriorCheckExpert( hdsdp_cone_sdp_sparse *cone, double dCCoef, double dACoefScal, double *dACoef, double dEyeCoef,
                                                       int whichBuffer, int *isInterior ) {
    
    hdsdp_retcode retcode = HDSDP_RETCODE_OK;
    sdpSparseConeIUpdateBuffer(cone, dCCoef, dACoefScal, dACoef, dEyeCoef, whichBuffer);
    
    if ( whichBuffer == BUFFER_DUALVAR ) {
        HDSDP_CALL(HFpLinsysPsdCheck(cone->dualFactor, cone->dualMatBeg, cone->dualMatIdx, cone->dualMatElem, isInterior));
    } else {
        HDSDP_CALL(HFpLinsysPsdCheck(cone->dualChecker, cone->dualMatBeg, cone->dualMatIdx, cone->dualCheckerElem, isInterior));
    }
    
exit_cleanup:
    return retcode;
}

extern void sdpDenseConeReduceResidual( hdsdp_cone_sdp_dense *cone, double resiReduction ) {
    
    cone->dualResidual = resiReduction;
    return;
}

extern void sdpSparseConeReduceResidual( hdsdp_cone_sdp_sparse *cone, double resiReduction ) {
    
    cone->dualResidual = resiReduction;
    return;
}

extern void sdpDenseConeSetPerturb( hdsdp_cone_sdp_dense *cone, double dDualPerturb ) {
    
    assert( dDualPerturb >= 0.0 );
    cone->dualPerturb = dDualPerturb;
    
    return;
}

extern void sdpSparseConeSetPerturb( hdsdp_cone_sdp_sparse *cone, double dDualPerturb ) {
    
    assert( dDualPerturb >= 0.0 );
    cone->dualPerturb = dDualPerturb;
    
    return;
}

extern hdsdp_retcode sdpDenseConeGetBarrier( hdsdp_cone_sdp_dense *cone, double barHsdTau, double *rowDual, int whichBuffer, double *logdet ) {
    
    /* Compute the barrier function at current dual slack
       If not supplied the barrier of the most recently updated dual solution is used
     */
    hdsdp_retcode retcode = HDSDP_RETCODE_OK;
    double dLogDeterminant = 0.0;
    
    hdsdp_linsys_fp *sTarget = NULL;
    double *sElememt = NULL;
    
    if ( whichBuffer == BUFFER_DUALVAR ) {
        sTarget = cone->dualFactor;
        sElememt = cone->dualMatElem;
    } else {
        sTarget = cone->dualChecker;
        sElememt = cone->dualCheckerElem;
    }
    
    if ( rowDual ) {
        sdpDenseConeUpdateImpl(cone, barHsdTau, rowDual);
        HDSDP_CALL(HFpLinsysNumeric(sTarget, cone->dualMatBeg, cone->dualMatIdx, sElememt));
    }
    
    HDSDP_CALL(HFpLinsysGetDiag(sTarget, cone->dualDiag));
    
    /* The log determinant is directly obtained from diagonal of the Cholesky factor
       since log det(S) = log det(L * L') = log (det(L)^2) = 2 log det(L)
     */
    
    for ( int iCol = 0; iCol < cone->nCol; ++iCol ) {
        dLogDeterminant += log(cone->dualDiag[iCol]);
    }
    
    *logdet = 2.0 * dLogDeterminant;
    
exit_cleanup:
    return retcode;
}

extern hdsdp_retcode sdpSparseConeGetBarrier( hdsdp_cone_sdp_sparse *cone, double barHsdTau, double *rowDual, int whichBuffer, double *logdet ) {
    
    /* Compute the barrier function at current dual slack
       If not supplied the barrier of the most recently updated dual solution is used */
    hdsdp_retcode retcode = HDSDP_RETCODE_OK;
    double dLogDeterminant = 0.0;
    
    hdsdp_linsys_fp *sTarget = NULL;
    double *sElement = NULL;
    
    if ( whichBuffer == BUFFER_DUALVAR ) {
        sTarget = cone->dualFactor;
        sElement = cone->dualMatElem;
    } else {
        sTarget = cone->dualChecker;
        sElement = cone->dualCheckerElem;
    }
    
    if ( rowDual ) {
        sdpSparseConeUpdateImpl(cone, barHsdTau, rowDual);
        HDSDP_CALL(HFpLinsysNumeric(sTarget, cone->dualMatBeg, cone->dualMatIdx, sElement));
    }
    
    HDSDP_CALL(HFpLinsysGetDiag(sTarget, cone->dualDiag));
    
    /* The log determinant is directly obtained from diagonal of the Cholesky factor
       since log det(S) = log det(L * L') = log (det(L)^2) = 2 log det(L)
     */
    
    for ( int iCol = 0; iCol < cone->nCol; ++iCol ) {
        dLogDeterminant += log(cone->dualDiag[iCol]);
    }
    
    *logdet = 2.0 * dLogDeterminant;
    
exit_cleanup:
    return retcode;
}

extern hdsdp_retcode sdpDenseConeAddStepToBufferAndCheck( hdsdp_cone_sdp_dense *cone, double dStep, int whichBuffer, int *isInterior ) {
    
    /* This utility add dS to the conic buffer and directly factorize it  */
    hdsdp_retcode retcode = HDSDP_RETCODE_OK;
    
    hdsdp_linsys_fp *sTarget = NULL;
    double *sElement = NULL;
    int nDualNz = 0;
    
    if ( cone->isDualSparse ) {
        nDualNz = cone->dualMatBeg[cone->nCol];
    } else {
        nDualNz = cone->nCol * cone->nCol;
    }
    
    if ( whichBuffer == BUFFER_DUALVAR ) {
        sTarget = cone->dualFactor;
        sElement = cone->dualMatElem;
    } else {
        sTarget = cone->dualChecker;
        sElement = cone->dualCheckerElem;
        HDSDP_MEMCPY(cone->dualCheckerElem, cone->dualMatElem, double, nDualNz);
    }
    
    axpy(&nDualNz, &dStep, cone->dualStep, &HIntConstantOne, sElement, &HIntConstantOne);
    HDSDP_CALL(HFpLinsysPsdCheck(sTarget, cone->dualMatBeg, cone->dualMatIdx, sElement, isInterior));
    
exit_cleanup:
    return retcode;
}

extern hdsdp_retcode sdpSparseConeAddStepToBufferAndCheck( hdsdp_cone_sdp_sparse *cone, double dStep, int whichBuffer, int *isInterior ) {
    
    /* This utility add dS to the conic buffer and directly factorize it  */
    hdsdp_retcode retcode = HDSDP_RETCODE_OK;
    
    hdsdp_linsys_fp *sTarget = NULL;
    double *sElement = NULL;
    int nDualNz = 0;
    
    if ( cone->isDualSparse ) {
        nDualNz = cone->dualMatBeg[cone->nCol];
    } else {
        nDualNz = cone->nCol * cone->nCol;
    }
    
    if ( whichBuffer == BUFFER_DUALVAR ) {
        sTarget = cone->dualFactor;
        sElement = cone->dualMatElem;
    } else {
        sTarget = cone->dualChecker;
        sElement = cone->dualCheckerElem;
        HDSDP_MEMCPY(cone->dualCheckerElem, cone->dualMatElem, double, nDualNz);
    }
    
    axpy(&nDualNz, &dStep, cone->dualStep, &HIntConstantOne, sElement, &HIntConstantOne);
    HDSDP_CALL(HFpLinsysPsdCheck(sTarget, cone->dualMatBeg, cone->dualMatIdx, sElement, isInterior));
    
exit_cleanup:
    return retcode;
}

extern void sdpDenseConeClearImpl( hdsdp_cone_sdp_dense *cone ) {
    
    if ( !cone ) {
        return;
    }
    
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
    
    HDSDP_FREE(cone->kktMapping);
    
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

extern void sdpDenseConeFeatureDetectImpl( hdsdp_cone_sdp_dense *cone, double *rowRHS,
                                           int coneIntFeatures[20], double coneDblFeatures[20] ) {
    /* When there is a single SDP cone. This routine detects if the SDP has the following structures
     
      1. (Almost) no primal interior point
      2. (Almost) no dual interior point
      3. Implied trace bound tr(X) == Z
      4. Implied dual upperbound and lower bound l <= y <= u
      5. No objective
      6. Dense cone
     
     Also the number of different cones will be records as the conic features (including objective)
     */
    
    /* We detect the case of no primal interior by seeking constraints that look like
       trace(a * a' * X) = b \approx 0.0
     */
    
    /* Get statistics */
    coneIntFeatures[INT_FEATURE_N_ZEORMATS] = cone->sdpConeStats[SDP_COEFF_ZERO];
    coneIntFeatures[INT_FEATURE_N_DSMATS] = cone->sdpConeStats[SDP_COEFF_DENSE];
    coneIntFeatures[INT_FEATURE_N_SPMATS] = cone->sdpConeStats[SDP_COEFF_SPARSE];
    coneIntFeatures[INT_FEATURE_N_SPR1MATS] = cone->sdpConeStats[SDP_COEFF_SPR1];
    coneIntFeatures[INT_FEATURE_N_DSR1MATS] = cone->sdpConeStats[SDP_COEFF_DSR1];
    
    int isNoPrimalInterior = 0;
    int isImpliedTraceX = 0;
    double dImpliedTraceX = 0.0;
    int iUnitCol = 0;
    int *iUnitColIdx = NULL;
    
    /* Detect if there is no primal interior point */
    for ( int iRow = 0; iRow < cone->nRow; ++iRow ) {
        if ( sdpDataMatGetRank(cone->sdpRow[iRow]) == 1 &&
            fabs(rowRHS[iRow]) < 1e-03 * sdpDataMatNorm(cone->sdpRow[iRow], FRO_NORM) ) {
            isNoPrimalInterior = 1;
        }
    }
    
    if ( isNoPrimalInterior ) {
        coneIntFeatures[INT_FEATURE_I_NOPINTERIOR] = 1;
    }
    
    /* Detect implied trace bound
       We detect two cases of implied bound.
     
       The first case is the constraint trace(I * X) == a
       The second case is the diagonal constraint diag(X) = d */
    
    /* First detect is there is trace(I * X) == a */
    int isEyeMultiple = 0;
    for ( int iRow = 0; iRow < cone->nRow; ++iRow ) {
        isEyeMultiple = sdpDataMatIsEye(cone->sdpRow[iRow], &dImpliedTraceX);
        if ( isEyeMultiple && rowRHS[iRow] / dImpliedTraceX > 0.0 ) {
            isImpliedTraceX = 1;
            dImpliedTraceX = rowRHS[iRow] / dImpliedTraceX;
            break;
        }
    }
    
    /* Then detect if there is diag(X) constraint. These constraints are of rank-one sparse coeffients */
    HDSDP_INIT(iUnitColIdx, int, cone->nRow);
    dImpliedTraceX = 0.0;
    if ( iUnitColIdx && !isImpliedTraceX ) {
        for ( int iRow = 0; iRow < cone->nRow; ++iRow ) {
            if ( sdpDataMatIsUnitCol(cone->sdpRow[iRow], &iUnitCol) && !iUnitColIdx[iUnitCol] ) {
                iUnitColIdx[iUnitCol] = 1;
                dImpliedTraceX += rowRHS[iRow];
            }
        }
    }
    
    int nSumUnit = 0;
    for ( int iRow = 0; iRow < cone->nRow; ++iRow ) {
        nSumUnit += iUnitColIdx[iRow];
    }
    
    if ( nSumUnit == cone->nCol ) {
        isImpliedTraceX = 1;
    }
    
    HDSDP_FREE(iUnitColIdx);
    
    if ( isImpliedTraceX ) {
        coneIntFeatures[INT_FEATURE_I_IMPTRACE] = 1;
        coneDblFeatures[DBL_FEATURE_IMPTRACEX] = dImpliedTraceX;
    }
    
    /* Detect if there is no objective */
    if ( sdpDataMatGetType(cone->sdpObj) == SDP_COEFF_ZERO ) {
        coneIntFeatures[INT_FEATURE_I_NULLOBJ] = 1;
    }
    
    /* Detect dense cone */
    if ( cone->sdpConeStats[SDP_COEFF_DENSE] >= 0.7 * cone->nRow ) {
        coneIntFeatures[INT_FEATURE_I_VERYDENSE] = 1;
    }
    
    return;
}

extern void sdpSparseConeFeatureDetectImpl( hdsdp_cone_sdp_dense *cone, double *rowRHS,
                                            int coneIntFeatures[20], double coneDblFeatures[20] ) {
    
    coneIntFeatures[INT_FEATURE_N_ZEORMATS] = cone->sdpConeStats[SDP_COEFF_ZERO];
    coneIntFeatures[INT_FEATURE_N_DSMATS] = cone->sdpConeStats[SDP_COEFF_DENSE];
    coneIntFeatures[INT_FEATURE_N_SPMATS] = cone->sdpConeStats[SDP_COEFF_SPARSE];
    coneIntFeatures[INT_FEATURE_N_SPR1MATS] = cone->sdpConeStats[SDP_COEFF_SPR1];
    coneIntFeatures[INT_FEATURE_N_DSR1MATS] = cone->sdpConeStats[SDP_COEFF_DSR1];
    
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
