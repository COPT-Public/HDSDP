#include <stdio.h>
#include "dsdphsd.h"
#include "dsdpsolver.h"
#include "data.h"


DSDP_INT test_dsdp(void) {
    
    DSDP_INT retcode = DSDP_RETCODE_OK;
    
    
    DSDP_INT m = coneM;
    DSDP_INT n = coneN;
    
    
    Solver *hsdSolver = NULL;
    Solver **phsdSolver = &hsdSolver;
    
    retcode = DSDPCreate(phsdSolver);
    retcode = DSDPSetDim(hsdSolver, n, 1, m, 0);
    retcode = DSDPSetSDPConeData(hsdSolver, 0, n, NULL, coneAp, coneAi, coneAx);
    
    retcode = DSDPSetObj(hsdSolver, dObj);
    retcode = DSDPOptimize(hsdSolver);
    
    // mwTrace("End Profiling. \n");
    
cleanup:
    retcode = DSDPDestroy(hsdSolver);
    
    return retcode;
}
