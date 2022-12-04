#include "hdsdp_conic_sdp.h"
#include "def_hdsdp_user_data.h"
#include "hdsdp_user_data.h"
#include "hdsdp_utils.h"

/** @brief Choose type of an SDP ceofficient matrix
 *
 */
static sdp_coeff_type sdpCoeffChooseType( int nnz, int *coeffIdx, double *coeffElem ) {
    
    
    
    return SDP_COEFF_ZERO;
}

/** @brief Create a dense sdp cone
 *
 */
extern hdsdp_retcode sdpDenseConeCreate( hdsdp_cone_sdp_dense **pSDPCone ) {
    
    hdsdp_retcode retcode = HDSDP_RETCODE_OK;
    
    if ( !pSDPCone ) {
        retcode = HDSDP_RETCODE_FAILED;
        goto exit_cleanup;
    }
    
    hdsdp_cone_sdp_dense *sdpCone = NULL;
    HDSDP_INIT(sdpCone, hdsdp_cone_sdp_dense, 1);
    HDSDP_MEMCHECK(sdpCone);
    
    HDSDP_ZERO(sdpCone, hdsdp_cone_sdp_dense, 1);
    *pSDPCone = sdpCone;
    
exit_cleanup:
    
    return retcode;
}

/** @brief Set data into a dense cone
 *
 */
extern hdsdp_retcode sdpDenseConeSetData( hdsdp_cone_sdp_dense *sdpCone, user_data *sdpInput ) {
    
    hdsdp_retcode retcode = HDSDP_RETCODE_OK;
    
    int *sdpMatBeg = sdpInput->coneMatBeg;
    int *sdpMatIdx = sdpInput->coneMatIdx;
    double *sdpMatElem = sdpInput->coneMatElem;
    
    /* Prepare conic statistics */
    int *sdpConeStats = sdpCone->sdpConeStats;
    
    
exit_cleanup:
    
    return retcode;
}
