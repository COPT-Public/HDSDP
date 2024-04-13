/** @file hdsdp\_user\_data.c
 *
 */

#ifdef HEADERPATH
#include "interface/def_hdsdp_user_data.h"
#include "interface/hdsdp_user_data.h"
#include "interface/hdsdp_utils.h"
#include "interface/hdsdp_conic.h"
#include "linalg/sparse_opts.h"
#else
#include "def_hdsdp_user_data.h"
#include "hdsdp_user_data.h"
#include "hdsdp_utils.h"
#include "hdsdp_conic.h"
#include "sparse_opts.h"
#endif


/** @brief Check if LP data implies bound constraint on y
 *
 */
static int HUserDataICheckLpBound( user_data *Hdata ) {
    
    if ( Hdata->cone != HDSDP_CONETYPE_LP ) {
        return 0;
    }
    
    /* If each column holds at most one variable, then it is bound */
    for ( int i = 0; i < Hdata->nConicRow; ++i ) {
        if ( Hdata->coneMatBeg[i + 1] - Hdata->coneMatBeg[i] >= 2 ) {
            return 0;
        }
    }
    
    return 1;
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

extern void HUserDataSetConeData( user_data *Hdata, cone_type cone, int nRow, int nCol,
                                  int *coneMatBeg, int *coneMatIdx, double *coneMatElem ) {
    
    Hdata->cone = cone;
    Hdata->nConicRow = nRow;
    Hdata->nConicCol = nCol;
    Hdata->coneMatBeg = coneMatBeg;
    Hdata->coneMatIdx = coneMatIdx;
    Hdata->coneMatElem = coneMatElem;
    
    return;
}

extern cone_type HUserDataChooseCone( user_data *Hdata ) {
        
    /* Automatic choice between different cone types*/
    if ( Hdata->cone == HDSDP_CONETYPE_SOCP || Hdata->cone == HDSDP_CONETYPE_BOUND || Hdata->cone == HDSDP_CONETYPE_LP ||
         Hdata->cone == HDSDP_CONETYPE_SPARSE_SDP || Hdata->cone == HDSDP_CONETYPE_SCALAR_BOUND ) {
        
        return Hdata->cone;
        
    } else if ( Hdata->cone == HDSDP_CONETYPE_DENSE_SDP ) {
        
        int nzSDPCoeffs = csp_nnz_cols(Hdata->nConicRow, &Hdata->coneMatBeg[1]);
        return ( nzSDPCoeffs > HDSDP_SPARSE_CONE_THRESHOLD * Hdata->nConicRow ) ? \
                HDSDP_CONETYPE_DENSE_SDP : HDSDP_CONETYPE_SPARSE_SDP;
        
    } else if ( Hdata->cone == HDSDP_CONETYPE_LP ) {
        
        if ( HUserDataICheckLpBound(Hdata) ) {
            return HDSDP_CONETYPE_BOUND;
        } else {
            return HDSDP_CONETYPE_LP;
        }
        
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
