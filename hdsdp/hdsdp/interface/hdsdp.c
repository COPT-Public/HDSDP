/** @file hdsdp.c
 *  @brief A homogeneous dual-scaling interior point solver
 */

#include "hdsdp.h"
#include "hdsdp_conic.h"
#include "hdsdp_utils.h"
#include "hdsdp_user_data.h"

struct hdsdp_solver_internal {
    
    char *coneModelName[100];
    
    /* Logging */
    
    
    /* User data */
    
    
    /* Schur complement */
    
    
    /* Step and vector */
    
    
    /* Monitor */
    
    
};


extern hdsdp_retcode HDSDPSolverCreate( hdsdp **pHSolver ) {
    
    hdsdp_retcode retcode = HDSDP_RETCODE_OK;
    
    if ( !pHSolver ) {
        retcode = HDSDP_RETCODE_FAILED;
        return retcode;
    }
    
    hdsdp *HSolver = NULL;
    HDSDP_INIT(HSolver, hdsdp, 1);
    HDSDP_MEMCHECK(HSolver);
    
    *pHSolver = HSolver;
    
exit_cleanup:
    
    return retcode;
}

extern void HDSDPSolverClear( hdsdp *HSolver ) {
    
    if ( !HSolver ) {
        return;
    }
    
    return;
}

extern void HDSDPSolverDestroy( hdsdp **pHSolver ) {
    
    if ( !pHSolver ) {
        return;
    }
    
    HDSDPSolverClear(*pHSolver);
    HDSDP_FREE(*pHSolver);
    
    return;
}
