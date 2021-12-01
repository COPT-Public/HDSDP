#include "hsd.h"

static char etype[] = "Dual scaling algorithm iteration";

// Implement the main frame of the dual-scaling algorithm for DSDP-HSD

/* Utility function to check iteration progress */
extern DSDP_INT checkIterProgress( HSDSolver *dsdpSolver, DSDP_INT iter ) {
    
    DSDP_INT retcode = DSDP_RETCODE_OK;
    DSDP_INT npassed = 0;
    
    for (DSDP_INT i = 0; i < iter; ++i) {
        npassed += dsdpSolver->iterProgress[i];
    }
    
    if (npassed <= iter) {
        error(etype, "Pre-requisite iterations are not completed. \n");
    }
    
    return retcode;
}

