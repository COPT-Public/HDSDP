#include "hdsdp_conic_sdp.h"
#include "hdsdp_utils.h"

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

extern hdsdp_retcode sdpDenseConeSetData( hdsdp_cone_sdp_dense *sdpCone ) {
    
    hdsdp_retcode retcode = HDSDP_RETCODE_OK;
    
    
    
exit_cleanup:
    
    return retcode;
}
