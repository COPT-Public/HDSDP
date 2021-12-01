#include "obj.h"
#include "dsdpsolver.h"

/* Objective value computer */
extern DSDP_INT getDualObj( HSDSolver *dsdpSolver ) {
    // Compute b' * y
    // Iteration prerequisite: normal
    // Event prerequisite    : none
    DSDP_INT retcode = DSDP_RETCODE_OK;
    DSDP_INT m = dsdpSolver->m;
    
    double dObj = 0.0;
    double *ydata = dsdpSolver->y->x;
    double *bdata = dsdpSolver->dObj->x;
    
    for (DSDP_INT i = 0; i < (DSDP_INT) m / 4; i+=4) {
        dObj += ydata[i    ] * bdata[i    ];
        dObj += ydata[i + 1] * bdata[i + 1];
        dObj += ydata[i + 2] * bdata[i + 2];
        dObj += ydata[i + 3] * bdata[i + 3];
    }
    
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
    
    vec **asinv = dsdpSolver->asinv;
    
    double *x1 = NULL;
    double *d2 = dsdpSolver->d2->x;
    
    for (DSDP_INT i = 0; i < nblock; ++i) {
        n = (double) dsdpSolver->S[i]->dim;
        pObjVal += n * mu;
        x1 = asinv[i]->x;
        for (DSDP_INT i = 0; i < (DSDP_INT) m / 4; i+=4) {
            pObjVal += d2[i    ] * x1[i    ];
            pObjVal += d2[i + 1] * x1[i + 1];
            pObjVal += d2[i + 2] * x1[i + 2];
            pObjVal += d2[i + 3] * x1[i + 3];
        }
    }
    
    dsdpSolver->pObjVal = pObjVal * mu + dsdpSolver->dObjVal;
    
    return retcode;
}
