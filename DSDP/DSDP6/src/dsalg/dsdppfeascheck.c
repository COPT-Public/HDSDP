#include "dsdppfeascheck.h"
#include "dsdputils.h"

/* Implement DSDP primal feasibility check */
extern DSDP_INT dsdpCheckPhaseAPfeas( HSDSolver *dsdpSolver, double dtaudelta,
                                     vec *dydelta, DSDP_INT *ispfeas ) {
    // Phase A feasibility check
    // - Ry + C * (tau - dtaudelta) - dsdpgetATy(A, y - dydelta)
    DSDP_INT retcode = DSDP_RETCODE_OK;
    
    vec_axpby(1.0, dsdpSolver->y, -1.0, dydelta);
    retcode = getPhaseACheckerS(dsdpSolver, dydelta->x, dsdpSolver->tau - dtaudelta);
    retcode = dsdpCheckerInCone(dsdpSolver, ispfeas);
    
    return retcode;
}

extern DSDP_INT dsdpCheckBackwardNewton( HSDSolver *dsdpSolver, DSDP_INT *ispfeas ) {
    
    // Check backward newton step C - dsdpgetATy(A, y - dymuprimal)
    DSDP_INT retcode = DSDP_RETCODE_OK;
    vec_zaxpby(dsdpSolver->b2, 1.0, dsdpSolver->y,
               -1.0, dsdpSolver->b1);
    
    retcode = getPhaseBCheckerS(dsdpSolver, dsdpSolver->b2->x);
    retcode = dsdpCheckerInCone(dsdpSolver, ispfeas);
    
    return retcode;
}

extern DSDP_INT dsdpCheckPrimalInfeas( HSDSolver *dsdpSolver ) {
    
    // Check dual unboundedness (primal infeasibility) through iterations
    DSDP_INT retcode = DSDP_RETCODE_OK;
    
    if ((!dsdpSolver->eventMonitor[EVENT_INVALID_GAP]) &&
        (!dsdpSolver->eventMonitor[EVENT_LARGE_DOBJ])) {
        return retcode;
    }
    
    double bTdy = 0.0;
    DSDP_INT incone = FALSE;
    vec_dot(dsdpSolver->dObj, dsdpSolver->dy, &bTdy);
    
    if (bTdy <= 0.0) {
        return retcode;
    }
    
    for (DSDP_INT i = 0; i < dsdpSolver->nBlock; ++i) {
        memcpy(dsdpSolver->Scker[i]->x, dsdpSolver->dS[i]->x,
               sizeof(double) * dsdpSolver->Scker[i]->nnz);
        spsMatIspd(dsdpSolver->S[i], &incone);
        if (!incone) {
            break;
        }
    }
    
    if (incone) {
        dsdpSolver->eventMonitor[EVENT_PINFEAS_DETECTED] = TRUE;
    }
    
    return retcode;
}
