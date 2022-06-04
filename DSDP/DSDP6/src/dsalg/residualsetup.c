#include "residualsetup.h"

// Setup all the residuals for the dual scaling algorithm to work
static char etype[] = "Residual setup";

static DSDP_INT getRy( HSDSolver *dsdpSolver ) {
    // Check the dual infeasibility
    dsdpSolver->Ry *= (1 - dsdpSolver->drate * dsdpSolver->alpha);
    return DSDP_RETCODE_OK;
}

extern DSDP_INT setupRes( HSDSolver *dsdpSolver ) {
    // Setup residuals used for LP and SDP
    DSDP_INT retcode = DSDP_RETCODE_OK;
    getRy(dsdpSolver);
    dsdpSolver->iterProgress[ITER_RESIDUAL] = TRUE;
    
    return retcode;
}
