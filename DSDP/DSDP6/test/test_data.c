#include <stdio.h>
#include "SDPdata.h"
#include "dsdphsd.h"
#include "dsdpsolver.h"
/* Test the data interface of DSDP */

DSDP_INT test_data(void) {
    
    DSDP_INT retcode = DSDP_RETCODE_OK;
    
    Solver *hsdSolver = NULL;
    Solver **phsdSolver = &hsdSolver;
    
    retcode = DSDPCreate(phsdSolver);
    retcode = DSDPSetDim(hsdSolver, cdim1 + cdim2 + cdim3, ncones, m, 0);
    
    DSDP_INT typehint[11] = {
        MAT_TYPE_RANK1, MAT_TYPE_UNKNOWN, MAT_TYPE_UNKNOWN, MAT_TYPE_UNKNOWN, MAT_TYPE_UNKNOWN,
        MAT_TYPE_UNKNOWN, MAT_TYPE_UNKNOWN, MAT_TYPE_UNKNOWN, MAT_TYPE_UNKNOWN, MAT_TYPE_UNKNOWN,
        MAT_TYPE_UNKNOWN};
    
    retcode = DSDPSetSDPConeData(hsdSolver, 0, cdim1, typehint, Ap_1, Ai_1, Ax_1);
    retcode = DSDPSetSDPConeData(hsdSolver, 1, cdim2, NULL, Ap_2, Ai_2, Ax_2);
    retcode = DSDPSetSDPConeData(hsdSolver, 2, cdim3, NULL, Ap_3, Ai_3, Ax_3);
    
cleanup:
    
    retcode = DSDPDestroy(hsdSolver);
    
    return retcode;
}
