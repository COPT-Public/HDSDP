#include "dsdpproxmeasure.h"
#include "dsdppfeascheck.h"
#include "stepheur.h"
#include "dsdputils.h"

extern DSDP_INT dsdpgetPhaseAProxMeasure( HSDSolver *dsdpSolver, double newmu ) {
    
    // Check primal feasibility and proximity measure given a new mu parameter
    DSDP_INT retcode = DSDP_RETCODE_OK;
    DSDP_INT ispfeas = FALSE;
    
    // Use b1 as auxiliary
    vec *dydelta = dsdpSolver->b1, *vecaux = dsdpSolver->cgSolver->aux;
    vec_zaxpby(dydelta, dsdpSolver->tau / dsdpSolver->mu,
               dsdpSolver->d2, -1.0, dsdpSolver->d3);
    
    // Compute Pnrm
    vec_zaxpby(vecaux, 1 / dsdpSolver->mu, dsdpSolver->dObj, -1.0, dsdpSolver->asinv);
    vec_dot(vecaux, dydelta, &dsdpSolver->Pnrm);
    dsdpSolver->Pnrm = sqrt(dsdpSolver->Pnrm);
    
    // Check feasibility
    dsdpCheckPhaseAPfeas(dsdpSolver, 0.0, dydelta, &ispfeas);
    
    // Update primal objective
    if (ispfeas) {
        // Fix the code here about primal objective evaluation
        dsdpSolver->eventMonitor[EVENT_PFEAS_FOUND] = TRUE;
        double pObj = dsdpSolver->dObjVal, tmp1, tmp2;
        
        vec_zaxpby(dydelta, dsdpSolver->tau / dsdpSolver->mu,
                   dsdpSolver->d2, -1.0, dsdpSolver->d3);
        vec_dot(dsdpSolver->asinvrysinv, dydelta, &tmp1);
        vec_dot(dsdpSolver->asinv, dydelta, &tmp2);
        
        pObj += dsdpSolver->mu / dsdpSolver->tau * \
                (dsdpSolver->rysinv + tmp1 + tmp2 + dsdpSolver->nall);
        pObj = pObj / dsdpSolver->tau;
        tmp1 = dsdpSolver->pObjVal;
        if (tmp1 == 1e+10 || tmp1 == 1e+05 || tmp1 == 0.0) {
            dsdpSolver->pObjVal = pObj;
        } else {
            dsdpSolver->pObjVal = MIN(pObj, dsdpSolver->pObjVal);
        }
        
        if (dsdpSolver->mu < 1e-03) {
            dsdpSolver->mumaker = dsdpSolver->mu;
            vec_copy(dsdpSolver->y, dsdpSolver->ymaker);
            vec_copy(dsdpSolver->b1, dsdpSolver->dymaker);
        }
    }
    
    return retcode;
}

extern DSDP_INT dsdpgetPhaseBProxMeasure( HSDSolver *dsdpSolver, double *muub, double *mulb ) {
    
    // Compute the proximal measure of the current iterate
    DSDP_INT retcode = DSDP_RETCODE_OK;
    DSDP_INT ispfeas = FALSE, usegold;
    
    double gap = 0.0, rho, bound = dsdpSolver->ybound;
    DSDPGetDblParam(dsdpSolver, DBL_PARAM_RHON, &rho);
    DSDPGetIntParam(dsdpSolver, INT_PARAM_GOLDSEARCH, &usegold);
    vec_zaxpby(dsdpSolver->b1, 1 / dsdpSolver->mu, dsdpSolver->d1, -1.0, dsdpSolver->d2);
    vec_zaxpby(dsdpSolver->d4, 1 / dsdpSolver->mu, dsdpSolver->dObj, -1.0, dsdpSolver->asinv);
    vec_dot(dsdpSolver->d4, dsdpSolver->dy, &dsdpSolver->Pnrm);
    dsdpSolver->Pnrm = sqrt(MAX(dsdpSolver->Pnrm, 0.0));
    dsdpCheckBackwardNewton(dsdpSolver, &ispfeas);
    
    if (ispfeas) {
        dsdpSolver->eventMonitor[EVENT_PFEAS_FOUND] = TRUE;
        vec_dot(dsdpSolver->b1, dsdpSolver->asinv, &gap);
        
        gap = (gap + dsdpSolver->nall) * dsdpSolver->mu;
        
        // Bound tightening
        double pinfeas = 0.0, tmp = 0.0, s;
        vec *sl = dsdpSolver->sl, *su = dsdpSolver->su, *dy = dsdpSolver->b1;
        
        if (usegold && dsdpSolver->mu > 1e-04 && dsdpSolver->pObjVal != 1e+10 && dsdpSolver->pObjVal != 1e+05) {
            double approxpObj;
            searchpObj(dsdpSolver, &approxpObj);
            dsdpSolver->pObjVal = MIN(dsdpSolver->pObjVal, \
                                      0.5 * approxpObj + 0.5 * dsdpSolver->pObjVal);
            dsdpSolver->pObjVal = MIN(dsdpSolver->dObjVal + gap, dsdpSolver->pObjVal);
            gap = dsdpSolver->pObjVal - dsdpSolver->dObjVal;
        } else {
            dsdpSolver->pObjVal = dsdpSolver->dObjVal + gap;
        }
        
        *muub = gap / dsdpSolver->nall;
        
        // Get primal infeasibility
        for (DSDP_INT i = 0; i < dsdpSolver->m; ++i) {
            tmp = dy->x[i];
            s = sl->x[i]; pinfeas = MAX(pinfeas, fabs(1.0 / s + tmp / (s * s)));
            s = su->x[i]; pinfeas = MAX(pinfeas, fabs(1.0 / s - tmp / (s * s)));
        }
        
        dsdpSolver->pinfeas = pinfeas * dsdpSolver->mu;
        
        // Compute the dnewtonub and dnewtonlb
        vec_lslack(dsdpSolver->b2, sl, -bound); vec_uslack(dsdpSolver->b2, su, bound);
        double denominator = fabs(dsdpSolver->pObjVal) + fabs(dsdpSolver->dObjVal) + 1 / dsdpSolver->cScaler;
        // Save primal objective
        if (gap >= 5e-03 * denominator) {
            dsdpSolver->mumaker2 = dsdpSolver->mu;
            vec_copy(dsdpSolver->y, dsdpSolver->ymaker2);
            vec_copy(dsdpSolver->b1, dsdpSolver->dymaker2);
        } else if (gap >= 1e-06 * denominator){
            dsdpSolver->mumaker = dsdpSolver->mu;
            vec_copy(dsdpSolver->y, dsdpSolver->ymaker);
            vec_copy(dsdpSolver->b1, dsdpSolver->dymaker);
        }
    } else {
        if (usegold && dsdpSolver->mu > 1e-04) {
            double approxpObj;
            searchpObj(dsdpSolver, &approxpObj);
            dsdpSolver->pObjVal = MIN(dsdpSolver->pObjVal, \
                                      0.9 * approxpObj + 0.1 * dsdpSolver->pObjVal);
        }
        *muub = (dsdpSolver->pObjVal - dsdpSolver->dObjVal) / dsdpSolver->nall;
    }
    
    *mulb = (*muub) / rho;
    return retcode;
}
