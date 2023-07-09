#ifdef HEADERPATH
#include "interface/hdsdp_conic_sdp.h"
#include "interface/def_hdsdp_user_data.h"
#include "interface/hdsdp_utils.h"
#include "linalg/hdsdp_sdpdata.h"
#include "linalg/sparse_opts.h"
#else
#include "hdsdp_conic_sdp.h"
#include "def_hdsdp_user_data.h"
#include "hdsdp_utils.h"
#include "hdsdp_sdpdata.h"
#include "sparse_opts.h"
#endif


static hdsdp_retcode sdpDenseConeAllocDualMat( hdsdp_cone_sdp_dense *cone ) {
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
    
    /* If there is a dense matrix, the dual matrix will be dense */
    if ( cone->sdpConeStats[SDP_COEFF_DENSE] > 0 || cone->sdpConeStats[SDP_COEFF_DSR1] > 0 ) {
        cone->isDualSparse = 0;
    }
    
    /* Allocate the dual sparsity pattern */
    if ( cone->isDualSparse ) {
        HDSDP_INIT(cone->dualPosToElemMap, int, PACK_NNZ(cone->nCol));
        HDSDP_MEMCHECK(cone->dualPosToElemMap);
    }
    
    /* Accumulate the sparsity pattern */
    
    
    
    
exit_cleanup:
    return retcode;
}

static hdsdp_retcode sdpSparseConeAllocDualMat( hdsdp_cone_sdp_sparse *cone ) {
    
    hdsdp_retcode retcode = HDSDP_RETCODE_OK;
    
    
    
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
    
    HDSDP_INIT(cone->sdpConePerm, int, nRow);
    HDSDP_MEMCHECK(cone->sdpConePerm);
    
    HDSDP_ZERO(cone->sdpConeStats, int, 5);
    
    /* Primal objective */
    HDSDP_CALL(sdpDataMatCreate(&cone->sdpObj));
    HDSDP_CALL(sdpDataMatSetData(cone->sdpObj, nCol, coneMatBeg[1], coneMatIdx, coneMatElem));
    
    /* Constraint data */
    for ( int iRow = 0; iRow < nRow; ++iRow ) {
        HDSDP_CALL(sdpDataMatCreate(&cone->sdpRow[iRow]));
        HDSDP_CALL(sdpDataMatSetData(cone->sdpRow[iRow], nCol, coneMatBeg[iRow + 2] - coneMatBeg[iRow + 1],
                                     coneMatIdx + coneMatBeg[iRow + 1], coneMatElem + coneMatBeg[iRow + 1]));
        sdp_coeff_type coeffType = sdpDataMatGetType(cone->sdpRow[iRow]);
        cone->sdpConeStats[coeffType] += 1;
    }
    
    /* TODO: Allocate dual variables */
        
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
    
    /* TODO: Allocate dual variables */
    
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
    for ( int iRow = 0; iRow < cone->nRowElem; ++iRow ) {
        cone->sdpConeStats[sdpDataMatGetType(cone->sdpRow[iRow])] -= 1;
        HDSDP_CALL(sdpDataMatBuildUpEigs(cone->sdpRow[iRow], preConeAuxi));
        cone->sdpConeStats[sdpDataMatGetType(cone->sdpRow[iRow])] += 1;
    }
    
exit_cleanup:
    
    HDSDP_FREE(preConeAuxi);
    return retcode;
}

extern void sdpDenseConeSetStartImpl( hdsdp_cone_sdp_dense *cone, double rResi ) {
    
    return;
}

extern void sdpSparseConeSetStartImpl( hdsdp_cone_sdp_sparse *cone, double rResi ) {
    
    return;
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
    
    /* TODO: Free memory of the dual iterations */
    
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
    
    /* TODO: Free memory of the dual iterations */
    
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
    
    return;
}

