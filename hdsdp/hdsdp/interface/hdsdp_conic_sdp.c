#include "interface/hdsdp_conic_sdp.h"
#include "interface/def_hdsdp_user_data.h"
#include "interface/hdsdp_utils.h"

#include "linalg/hdsdp_sdpdata.h"

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
    
    HDSDP_ZERO(cone->sdpConeStats, int, 5);
    
    /* Primal objective */
    HDSDP_CALL(sdpDataMatCreate(&cone->sdpObj));
    HDSDP_CALL(sdpDataMatSetData(cone->sdpObj, nCol, coneMatBeg[1], coneMatIdx, coneMatElem));
    
    /* Sparse constraint data */
    int nRowElem = 0;
    
    /* Count #Nonzeros */
    for ( int iRow = 0; iRow < nRow; ++iRow ) {
        if ( coneMatBeg[iRow + 2] - coneMatBeg[iRow + 1] > 0 ) {
            nRowElem += 1;
        }
    }
    
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
            nRowElemTmp += 1;
        } else {
            cone->sdpConeStats[SDP_COEFF_ZERO] += 1;
        }
    }
    
#ifdef HDSDP_CONIC_DEBUG
    assert( nRowElem == nRowElemTmp );
#endif
    
    /* TODO: Allocate dual variables */
    
exit_cleanup:
    
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
    
    printf("Dense SDP cone of %d rows. \n", cone->nRow);
    printf("Objective: \n");
    sdpDataMatView(cone->sdpObj);
    
    printf("Constraint: \n");
    for ( int iRow = 0; iRow < cone->nRow; ++iRow ) {
        sdpDataMatView(cone->sdpRow[iRow]);
    }
    
    printf("\nConic statistics: Zero %d Sp %d Ds %d SpR1%d DsR1 %d \n", cone->sdpConeStats[SDP_COEFF_ZERO],
           cone->sdpConeStats[SDP_COEFF_SPARSE], cone->sdpConeStats[SDP_COEFF_DENSE],
           cone->sdpConeStats[SDP_COEFF_SPR1], cone->sdpConeStats[SDP_COEFF_DSR1]);
    
    return;
}

extern void sdpSparseConeViewImpl( hdsdp_cone_sdp_sparse *cone ) {
    
    printf("Sparse SDP cone of %d rows and %d nonzeros. \n", cone->nRow, cone->nRowElem);
    printf("Objective: \n");
    sdpDataMatView(cone->sdpObj);
    
    printf("Constraint: \n");
    for ( int iRow = 0; iRow < cone->nRowElem; ++iRow ) {
        printf("\nRow index %d \n", cone->rowIdx[iRow]);
        sdpDataMatView(cone->sdpRow[iRow]);
    }
    
    printf("\nConic statistics: Zero %d Sp %d Ds %d SpR1%d DsR1 %d \n", cone->sdpConeStats[SDP_COEFF_ZERO],
           cone->sdpConeStats[SDP_COEFF_SPARSE], cone->sdpConeStats[SDP_COEFF_DENSE],
           cone->sdpConeStats[SDP_COEFF_SPR1], cone->sdpConeStats[SDP_COEFF_DSR1]);
    
    return;
}

