#include "residualsetup.h"

// Setup all the residuals for the dual scaling algorithm to work
static char etype[] = "Residual setup";

static DSDP_INT getRkappaTau( HSDSolver *dsdpSolver ) {
    // Setup the complementarity residual rtk
    DSDP_INT retcode = DSDP_RETCODE_OK;
    dsdpSolver->rtk = dsdpSolver->tau * dsdpSolver->kappa - dsdpSolver->mu;
    return retcode;
}

static DSDP_INT getRy( HSDSolver *dsdpSolver ) {
    // Check the dual infeasibility
    dsdpSolver->Ry *= (1 - dsdpSolver->alpha);
    return DSDP_RETCODE_OK;
}

static DSDP_INT getLPResidualry( HSDSolver *dsdpSolver, vec *ry ) {
    // Compute ry = - A' * y + c * tau - s
    DSDP_INT retcode = DSDP_RETCODE_OK;
    
    lpMat *Adata = dsdpSolver->lpData;
    DSDP_INT m   = Adata->dimy;
    
    vec *y = dsdpSolver->y;
    vec *s = dsdpSolver->s;
    vec *c = dsdpSolver->lpObj;
    double tau = dsdpSolver->tau;
    double *ydata = ry->x;

    assert( ry->dim = m );
    
    // Setup - A' * y
    double alpha = - 1.0;
    retcode = lpMataATy(alpha, Adata, y, ydata);
    // Get - A' * y + c * tau
    retcode = vec_axpy(tau, c, ry);
    // Get - A' * y + c * tau - s
    retcode = vec_axpy(alpha, s, ry);
    
    return retcode;
}

extern DSDP_INT setupRes( HSDSolver *dsdpSolver ) {
    // Setup residuals used for LP and SDP
    DSDP_INT retcode = DSDP_RETCODE_OK;
    
    checkIterProgress(dsdpSolver, ITER_RESIDUAL);
    retcode = getRy(dsdpSolver); checkCode;
    retcode = getRkappaTau(dsdpSolver); checkCode;
    
    return retcode;
}
