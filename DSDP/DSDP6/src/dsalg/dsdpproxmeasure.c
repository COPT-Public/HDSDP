#include "dsdpproxmeasure.h"
#include "dsdppfeascheck.h"
#include "dsdputils.h"

extern DSDP_INT dsdpgetPhaseAProxMeasure( HSDSolver *dsdpSolver, double newmu ) {
    
    // Check primal feasibility and proximity measure given a new mu parameter
    DSDP_INT retcode = DSDP_RETCODE_OK;
    
    double taudenom = 0.0;
    double dtaudelta = 0.0;
    double tau = dsdpSolver->tau;
    double dObj = dsdpSolver->dObjVal;
    double utd2 = 0.0;
    DSDP_INT ispfeas = FALSE;
    
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
    denseMatxTAx(dsdpSolver->Msdp, dydelta->x, &utd2);
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
    
    double gap = 0.0, rho;
    retcode = DSDPGetDblParam(dsdpSolver, DBL_PARAM_RHO, &rho);
    
    vec_zaxpby(dsdpSolver->b1, 1 / dsdpSolver->mu, dsdpSolver->d1, -1.0, dsdpSolver->d2);
    
    retcode = denseMatxTAx(dsdpSolver->Msdp, dsdpSolver->b1->x, &dsdpSolver->Pnrm);
    dsdpSolver->Pnrm = MAX(dsdpSolver->Pnrm * sqrt(dsdpSolver->Mscaler), 1e-08);
    retcode = dsdpCheckBackwardNewton(dsdpSolver, &ispfeas);
    
    if (ispfeas) {
        dsdpSolver->eventMonitor[EVENT_PFEAS_FOUND] = TRUE;
        vec_dot(dsdpSolver->b1, dsdpSolver->asinv, &gap);
        gap = (gap + dsdpSolver->n) * dsdpSolver->mu;
        dsdpSolver->pObjVal = dsdpSolver->dObjVal + gap;
        *muub = gap / dsdpSolver->n;
        
        if (gap < 0.1 * (fabs(dsdpSolver->pObjVal) + fabs(dsdpSolver->dObjVal))) {
            dsdpSolver->mumaker = dsdpSolver->mu;
            vec_copy(dsdpSolver->y, dsdpSolver->ymaker);
            vec_copy(dsdpSolver->b1, dsdpSolver->dymaker);
        }
        
    } else {
        *muub = (dsdpSolver->pObjVal - dsdpSolver->dObjVal) / dsdpSolver->n;
    }
    
    *mulb = (*muub) / rho;
    
    return retcode;
}
