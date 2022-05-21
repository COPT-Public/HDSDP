#include "dsdppfeascheck.h"
#include "dsdputils.h"

/* Implement DSDP primal feasibility check */
extern void dsdpCheckPhaseAPfeas( HSDSolver *dsdpSolver, double dtaudelta,
                                     vec *dydelta, DSDP_INT *ispfeas ) {
    // Phase A feasibility check
    // - Ry + C * (tau - dtaudelta) - dsdpgetATy(A, y - dydelta)
    vec_axpby(1.0, dsdpSolver->y, -1.0, dydelta);
    DSDPConic( COPS_CONSTR_EXPR )(dsdpSolver, CHECKER, -1.0,
                                  dydelta, dsdpSolver->tau - dtaudelta, -1.0);
    *ispfeas = DSDPConic( COPS_CHECK_INCONE )(dsdpSolver, CHECKER);
//     getPhaseACheckerS(dsdpSolver, dydelta, dsdpSolver->tau - dtaudelta);
//     dsdpCheckerInCone(dsdpSolver, ispfeas);
}

extern void dsdpCheckBackwardNewton( HSDSolver *dsdpSolver, DSDP_INT *ispfeas ) {
    
    // Check backward newton step C - dsdpgetATy(A, y - dymuprimal)
    vec_zaxpby(dsdpSolver->b2, 1.0, dsdpSolver->y,
               -1.0, dsdpSolver->b1);
    
    // Check bound cone
    double *bwnewton = dsdpSolver->b2->x, bound = dsdpSolver->ybound;
    
    for (DSDP_INT i = 0; i < dsdpSolver->m; ++i) {
        if (bound < fabs(bwnewton[i])) {
            *ispfeas = FALSE; return;
        }
    }
    
    getPhaseBCheckerS(dsdpSolver, dsdpSolver->b2);
    dsdpCheckerInCone(dsdpSolver, ispfeas);
}

extern void dsdpCheckPrimalInfeas( HSDSolver *dsdpSolver ) {
    
    // Check dual unboundedness (primal infeasibility) through iterations
    if (dsdpSolver->pObjVal <= dsdpSolver->dObjVal) {
        dsdpSolver->eventMonitor[EVENT_INVALID_GAP] = TRUE;
        dsdpSolver->pObjVal = dsdpSolver->dObjVal + 1e+05;
        DSDPStatUpdate(&dsdpSolver->dsdpStats, STAT_GAP_BROKEN, 1.0);
    } else {
        dsdpSolver->eventMonitor[EVENT_INVALID_GAP] = FALSE;
    }
    
    if (dsdpSolver->dObjVal >= 1e+08) {
        dsdpSolver->eventMonitor[EVENT_LARGE_DOBJ] = TRUE;
    } else {
        dsdpSolver->eventMonitor[EVENT_LARGE_DOBJ] = FALSE;
    }
    
    if ((!dsdpSolver->eventMonitor[EVENT_INVALID_GAP]) &&
        (!dsdpSolver->eventMonitor[EVENT_LARGE_DOBJ])) {
        return;
    }
    
    double bTdy = 0.0;
    DSDP_INT incone = FALSE;
    vec_dot(dsdpSolver->dObj, dsdpSolver->dy, &bTdy);
    
    if (bTdy <= 0.0) {
        return;
    }
    
    for (DSDP_INT i = 0; i < dsdpSolver->nBlock; ++i) {
        getPhaseBS(dsdpSolver, dsdpSolver->dy);
        memcpy(dsdpSolver->Scker[i]->x, dsdpSolver->dS[i]->x,
               sizeof(double) * dsdpSolver->Scker[i]->nnz);
        spsMatIspd(dsdpSolver->Scker[i], &incone);
        if (!incone) {
            break;
        }
    }
    
    if (incone) {
        dsdpSolver->eventMonitor[EVENT_PINFEAS_DETECTED] = TRUE;
    }
        
    return;
}
