#include "dsdputils.h"
#include "dsdpsolver.h"
static char etype[] = "DSDP Utility";

/* Coefficient norm computer */
extern DSDP_INT getMatnrm( HSDSolver *dsdpSolver, DSDP_INT blockid, DSDP_INT constrid, double *nrm ) {
    
    DSDP_INT retcode = DSDP_RETCODE_OK;
    void *data = dsdpSolver->sdpData[blockid]->sdpData[constrid];
    
    switch (dsdpSolver->sdpData[blockid]->types[constrid]) {
        case MAT_TYPE_ZERO:
            *nrm = 0.0;
            break;
        case MAT_TYPE_DENSE:
            retcode = denseMatFnorm(data, nrm);
            checkCode;
            break;
        case MAT_TYPE_SPARSE:
            retcode = spsMatFnorm(data, nrm);
            checkCode;
            break;
        case MAT_TYPE_RANK1:
            retcode = r1MatFnorm(data, nrm);
            checkCode;
            break;
        default:
            error(etype, "Unknown matrix type. \n");
            break;
    }
    
    return retcode;
}

/* Compute S + alpha * some matrix */
extern DSDP_INT addMattoS( HSDSolver *dsdpSolver, DSDP_INT blockid, DSDP_INT constrid, double alpha ) {
    
    DSDP_INT retcode = DSDP_RETCODE_OK;
    void *data = dsdpSolver->sdpData[blockid]->sdpData[constrid];
    
    switch (dsdpSolver->sdpData[blockid]->types[constrid]) {
        case MAT_TYPE_ZERO:
            break;
        case MAT_TYPE_DENSE:
            retcode = spsMatAddds(dsdpSolver->S[blockid], alpha, data);
            checkCode;
            break;
        case MAT_TYPE_SPARSE:
            retcode = spsMataXpbY(alpha, data, 1.0, dsdpSolver->S[blockid]);
            checkCode;
            break;
        case MAT_TYPE_RANK1:
            retcode = spsMatAddr1(dsdpSolver->S[blockid], alpha, data);
            checkCode;
            break;
        default:
            error(etype, "Unknown matrix type. \n");
            break;
    }
    
    return retcode;
}

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

extern DSDP_INT getSDPPrimalObjB( HSDSolver *dsdpSolver ) {
    // Compute primal objective given primal feasible projection
    // Iteration prerequisite: normal
    // Event prerequisite: EVENT_SDP_NO_RY, EVENT_PFEAS_FOUND
    
    DSDP_INT retcode = DSDP_RETCODE_OK;
    
    assert( dsdpSolver->eventMonitor[EVENT_SDP_NO_RY] &&
            dsdpSolver->eventMonitor[EVENT_PFEAS_FOUND] );
    
    double pObjVal  = 0.0;
    DSDP_INT m  = dsdpSolver->m;
    double mu   = dsdpSolver->mu;
    DSDP_INT incx = 1;
    
    vec *asinv = dsdpSolver->asinv;

    double *x1 = asinv->x;
    double *d2 = dsdpSolver->d2->x;
        
    pObjVal = (double) dsdpSolver->n + dot(&m, x1, &incx, d2, &incx);
    dsdpSolver->pObjVal = pObjVal * mu + dsdpSolver->dObjVal;
    
    return retcode;
}
