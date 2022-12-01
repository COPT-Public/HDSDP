/** @file hdsdp\_user\_data.c
 *
 */
#include "hdsdp_user_data.h"
#include "hdsdp_utils.h"
#include "hdsdp_conic.h"
#include "sparse_opts.h"

/* Interface of user data */

/** @struct hdsdp\_user\_data
 *  @brief HDSDP user conic data for SDP and LP
 *
 */
struct hdsdp_user_data {
    
    cone_type cone;
    
    int     nConicRow;
    int     nConicCol;
    int    *coneMatBeg;
    int    *coneMatIdx;
    double *coneMatElem;
    
};

/** @brief Check if LP data implies bound constraint on y
 *
 */
static int HUserDataICheckLpBound( int nCol, int nRow, int *colMatBeg, int *colMatIdx, double *colMatElem ) {
    
    
    
    return 0;
}

extern hdsdp_retcode HUserDataCreate( user_data **pHdata ) {
    
    hdsdp_retcode retcode = HDSDP_RETCODE_OK;
    
    if ( !pHdata ) {
        retcode = HDSDP_RETCODE_FAILED;
        goto exit_cleanup;
    }
    
    user_data *Hdata = NULL;
    HDSDP_INIT(Hdata, user_data, 1);
    HDSDP_MEMCHECK(Hdata);
    
    HDSDP_ZERO(Hdata, user_data, 1);
    *pHdata = Hdata;
    
exit_cleanup:
    
    return retcode;
}

extern void HUserDataSetConeInfo( user_data *Hdata, cone_type cone, int nRow, int nCol ) {
    
    Hdata->cone = cone;
    Hdata->nConicRow = nRow;
    Hdata->nConicCol = nCol;
    
    return;
}

extern void HUserDataSetConeData( user_data *Hdata, int *coneMatBeg, int *coneMatIdx, double *coneMatElem ) {
    
    Hdata->coneMatBeg = coneMatBeg;
    Hdata->coneMatIdx = coneMatIdx;
    Hdata->coneMatElem = coneMatElem;
    
    return;
}

extern cone_type HUserDataChooseCone( user_data *Hdata ) {
    
    assert( Hdata->cone != HDSDP_CONETYPE_SPARSE_SDP );
    if ( Hdata->cone == HDSDP_CONETYPE_SOCP || Hdata->cone == HDSDP_CONETYPE_BOUND ||
         Hdata->cone == HDSDP_CONETYPE_SPARSE_SDP ) {
        return Hdata->cone;
    } else if ( Hdata->cone == HDSDP_CONETYPE_DENSE_SDP ) {
        int nzCoeffs = csp_nnz_cols(Hdata->nConicRow + 1, Hdata->coneMatBeg);
        return ( nzCoeffs > 0.3 * Hdata->nConicRow ) ? \
                HDSDP_CONETYPE_DENSE_SDP : HDSDP_CONETYPE_SPARSE_SDP;
    } else if ( Hdata->cone == HDSDP_CONETYPE_LP ) {
        
    }
    
    return HDSDP_CONETYPE_UNKNOWN;
}

extern void HUserDataClear( user_data *Hdata ) {
    
    if ( !Hdata ) {
        return;
    }
    
    HDSDP_ZERO(Hdata, user_data, 1);
    
    return;
}

extern void HUserDataDestroy( user_data **pHdata ) {
    
    if ( !pHdata ) {
        return;
    }
    
    HUserDataClear(*pHdata);
    HDSDP_FREE(*pHdata);
    
    return;
}
