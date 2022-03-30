#include "dsdpproxmeasure.h"
#include "dsdppfeascheck.h"
#include "dsdputils.h"

extern DSDP_INT dsdpgetPhaseAProxMeasure( HSDSolver *dsdpSolver, double newmu ) {
    
    // Check primal feasibility and proximity measure given a new mu parameter
    DSDP_INT retcode = DSDP_RETCODE_OK;
    DSDP_INT ispfeas = FALSE;
    
    // Use b1 as auxiliary
    vec *dydelta = dsdpSolver->b1;
    vec_zaxpby(dydelta, dsdpSolver->tau / dsdpSolver->mu,
               dsdpSolver->d2, -1.0, dsdpSolver->d3);
    
    // Compute Pnrm
    dsdpSolver->Pnrm = denseMatxTAx(dsdpSolver->Msdp, dsdpSolver->M->schurAux, dydelta->x);
    dsdpSolver->Pnrm = sqrt(dsdpSolver->Pnrm);
    
    // Check feasibility
    retcode = dsdpCheckPhaseAPfeas(dsdpSolver, 0.0, dydelta, &ispfeas);
    
    // Update primal objective
    if (ispfeas) {
        dsdpSolver->eventMonitor[EVENT_PFEAS_FOUND] = TRUE;
        double pObj = dsdpSolver->dObjVal, tmp1, tmp2;
        vec_dot(dsdpSolver->u, dydelta, &tmp1);
        vec_dot(dsdpSolver->asinv, dydelta, &tmp2);
        pObj += dsdpSolver->mu / dsdpSolver->tau * \
                (dsdpSolver->rysinv + tmp1 + tmp2 + dsdpSolver->n + dsdpSolver->m * 2);
        pObj = pObj / dsdpSolver->tau;
    }
    
    return retcode;
}

extern DSDP_INT dsdpgetPhaseBProxMeasure( HSDSolver *dsdpSolver, double *muub, double *mulb ) {
    
    // Compute the proximal measure of the current iterate
    DSDP_INT retcode = DSDP_RETCODE_OK;
    DSDP_INT ispfeas = FALSE;
    
    double gap = 0.0, rho, bound;
    
    DSDPGetDblParam(dsdpSolver, DBL_PARAM_PRLX_PENTALTY, &bound);
    DSDPGetDblParam(dsdpSolver, DBL_PARAM_RHON, &rho);
    
    vec_zaxpby(dsdpSolver->b1, 1 / dsdpSolver->mu, dsdpSolver->d1, -1.0, dsdpSolver->d2);
    
    dsdpSolver->Pnrm = denseMatxTAx(dsdpSolver->Msdp, dsdpSolver->M->schurAux, dsdpSolver->b1->x);
    dsdpSolver->Pnrm = sqrt(MAX(dsdpSolver->Pnrm, 0.0));
    retcode = dsdpCheckBackwardNewton(dsdpSolver, &ispfeas);
    
    if (ispfeas) {
        dsdpSolver->eventMonitor[EVENT_PFEAS_FOUND] = TRUE;
        vec_dot(dsdpSolver->b1, dsdpSolver->asinv, &gap);
        gap = (gap + dsdpSolver->n + dsdpSolver->m * 2) * dsdpSolver->mu;
        dsdpSolver->pObjVal = dsdpSolver->dObjVal + gap;
        *muub = gap / (dsdpSolver->n + dsdpSolver->m * 2);
        
        double pinfeas = 0.0, tmp = 0.0, s;
        vec *sl = dsdpSolver->sl, *su = dsdpSolver->su, *dy = dsdpSolver->b1;
        
        // Get primal infeasibility
        for (DSDP_INT i = 0; i < dsdpSolver->m; ++i) {
            tmp = dy->x[i];
            s = sl->x[i]; pinfeas = MAX(pinfeas, fabs(1.0 / s + tmp / (s * s)));
            s = su->x[i]; pinfeas = MAX(pinfeas, fabs(1.0 / s - tmp / (s * s)));
        }
        
        dsdpSolver->pinfeas = pinfeas * dsdpSolver->mu;
        
        // Compute the dnewtonub and dnewtonlb
        vec_lslack(dsdpSolver->b2, sl, -bound); vec_uslack(dsdpSolver->b2, su, bound);
        // Save primal objective
        if (gap < 0.1 * (fabs(dsdpSolver->pObjVal) + fabs(dsdpSolver->dObjVal) + 1) &&
            gap >= 1e-06 * (fabs(dsdpSolver->pObjVal) + fabs(dsdpSolver->dObjVal) + 1)) {
            dsdpSolver->mumaker = dsdpSolver->mu;
            vec_copy(dsdpSolver->y, dsdpSolver->ymaker);
            vec_copy(dsdpSolver->b1, dsdpSolver->dymaker);
        }
        
    } else {
        *muub = (dsdpSolver->pObjVal - dsdpSolver->dObjVal) / (dsdpSolver->n + 2 *  dsdpSolver->m);
    }
    
    *mulb = (*muub) / rho;
    
    return retcode;
}
