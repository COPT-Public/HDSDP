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
        dtaudelta = - dObj + newmu * (1 / tau + dsdpSolver->csinv) + utd2;
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
    denseMatxTAx(dsdpSolver->Msdp, dydelta, &utd2);
    dsdpSolver->Pnrm = sqrt(dsdpSolver->Pnrm + utd2);
    
    retcode = dsdpCheckPhaseAPfeas(dsdpSolver, dtaudelta, dydelta, &ispfeas);
    
    if (ispfeas) {
        dsdpSolver->eventMonitor[EVENT_PFEAS_FOUND] = TRUE;
        utd2 = newmu / (tau * tau) * (tau - dtaudelta);
        dsdpSolver->pObjVal = dObj - utd2;
    }
    
    return retcode;
}
