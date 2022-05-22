#include "residualsetup.h"

// Setup all the residuals for the dual scaling algorithm to work
static char etype[] = "Residual setup";

static DSDP_INT getRy( HSDSolver *dsdpSolver ) {
    // Check the dual infeasibility
    dsdpSolver->Ry *= (1 - dsdpSolver->drate * dsdpSolver->alpha);
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
    lpMataATy(alpha, Adata, y, ydata);
    // Get - A' * y + c * tau
    vec_axpy(tau, c, ry);
    // Get - A' * y + c * tau - s
    vec_axpy(alpha, s, ry);
    
    return retcode;
}

extern DSDP_INT setupRes( HSDSolver *dsdpSolver ) {
    // Setup residuals used for LP and SDP
    DSDP_INT retcode = DSDP_RETCODE_OK;
    getRy(dsdpSolver);
    dsdpSolver->iterProgress[ITER_RESIDUAL] = TRUE;
    
    return retcode;
}
