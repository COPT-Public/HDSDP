#include "dsdpproxmeasure.h"
#include "dsdppfeascheck.h"
#include "dsdputils.h"

extern DSDP_INT dsdpgetPhaseAProxMeasure( HSDSolver *dsdpSolver, double newmu ) {
    
    // Check primal feasibility and proximity measure given a new mu parameter
    DSDP_INT retcode = DSDP_RETCODE_OK;
    
    DSDP_INT ispfeas = FALSE, ignoredobj;
    DSDPGetIntParam(dsdpSolver, INT_PARAM_DOBJ_IGNORE, &ignoredobj);
    double taudenom = 0.0, dtaudelta = 0.0, tau = dsdpSolver->tau,
           dObj = (ignoredobj) ? 0.0 : dsdpSolver->dObjVal, utd2 = 0.0;
    
    // Use b1 as auxiliary
    vec *dydelta = dsdpSolver->b1;
    vec_dot(dsdpSolver->u, dsdpSolver->d12, &taudenom);
    
    taudenom = - taudenom + dsdpSolver->csinvcsinv + 1 / (tau * tau);
    
    if (fabs(taudenom) < 1e-15) {
        dtaudelta = 0.0;
    } else {
        // dtaudelta = - dObj + mu / tau + mu * csinv + tau * u' * d2 - mu * u' * d3;
        vec_dot(dsdpSolver->u, dsdpSolver->d2, &utd2);
        dtaudelta = - dObj + newmu * (1 / tau + dsdpSolver->csinv) + tau * utd2;
        vec_dot(dsdpSolver->u, dsdpSolver->d3, &utd2);
        dtaudelta -= utd2 * newmu;
        
        // dtaudelta = (1 / mu) *  dtaudelta / taudenom;
        dtaudelta = dtaudelta / (newmu * taudenom);
    }
    
    // dydelta = d12 * dtaudelta + d2 * tau / mu - d3;
    vec_zaxpby(dydelta, dtaudelta, dsdpSolver->d12, tau / newmu, dsdpSolver->d2);
    vec_axpy(-1.0, dsdpSolver->d3, dydelta);
    vec_dot(dsdpSolver->u, dydelta, &utd2);
    
    dsdpSolver->Pnrm = (dsdpSolver->csinvcsinv + 1 / (tau * tau)) * (dtaudelta * dtaudelta);
    dsdpSolver->Pnrm -= 2 * utd2 * dtaudelta;
    utd2 = denseMatxTAx(dsdpSolver->Msdp, dsdpSolver->M->schurAux, dydelta->x);
    dsdpSolver->Pnrm = sqrt(dsdpSolver->Pnrm + dsdpSolver->Mscaler * utd2);
    
    retcode = dsdpCheckPhaseAPfeas(dsdpSolver, dtaudelta, dydelta, &ispfeas);
    
    if (ispfeas) {
        dsdpSolver->eventMonitor[EVENT_PFEAS_FOUND] = TRUE;
        utd2 = newmu / (tau * tau) * (tau - dtaudelta);
        dsdpSolver->pObjVal = dObj - utd2;
        dsdpSolver->pObjVal = dsdpSolver->pObjVal / dsdpSolver->tau;
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
