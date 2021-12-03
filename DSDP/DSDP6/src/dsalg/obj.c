#include "obj.h"
#include "dsdpsolver.h"

/* Objective value computer */
extern DSDP_INT getDualObj( HSDSolver *dsdpSolver ) {
    // Compute b' * y
    // Iteration prerequisite: normal
    // Event prerequisite    : none
    DSDP_INT retcode = DSDP_RETCODE_OK;
    DSDP_INT m = dsdpSolver->m;
    DSDP_INT incx = 1;
    double *ydata = dsdpSolver->y->x;
    double *bdata = dsdpSolver->dObj->x;
    dsdpSolver->dObjVal = dot(&m, bdata, &incx, ydata, &incx);
    
    return retcode;
}

extern DSDP_INT getSDPPrimalObj( HSDSolver *dsdpSolver ) {
    // Compute primal objective given primal feasible projection
    // Iteration prerequisite: normal
    // Event prerequisite: EVENT_SDP_NO_RY, EVENT_PFEAS_FOUND
    
    DSDP_INT retcode = DSDP_RETCODE_OK;
    
    assert( dsdpSolver->eventMonitor[EVENT_SDP_NO_RY] &&
            dsdpSolver->eventMonitor[EVENT_PFEAS_FOUND] );
    
    double pObjVal  = 0.0;
    DSDP_INT nblock = dsdpSolver->nBlock;
    DSDP_INT m  = dsdpSolver->m;
    double n    = 0.0;
    double mu   = dsdpSolver->mu;
    DSDP_INT incx = 1;
    
    vec *asinv = dsdpSolver->asinv;

    double *x1 = asinv->x;
    double *d2 = dsdpSolver->d2->x;
        
    // TODO: Save \sum n_k and avoid the repeating summation
    for (DSDP_INT i = 0; i < nblock; ++i) {
        n = (double) dsdpSolver->S[i]->dim;
        pObjVal += n;
    }
    pObjVal += dot(&m, x1, &incx, d2, &incx);
    dsdpSolver->pObjVal = pObjVal * mu + dsdpSolver->dObjVal;
    
    return retcode;
}
