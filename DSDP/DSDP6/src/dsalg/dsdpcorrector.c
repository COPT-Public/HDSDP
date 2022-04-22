#include "dsdpcorrector.h"
#include "dsdputils.h"

/*
 Implement the dual corrector for DSDP Phase B
*/

static char etype[] = "Dual corrector step";

static double getCurrentLogBarrier( HSDSolver *dsdpSolver, vec *y, DSDP_INT *inCone ) {
    
    double logbarrier = 0.0, logdet = 0.0, *aux = dsdpSolver->M->schurAux, sl, su, bound;
    bound = dsdpSolver->ybound;
    vec *b = dsdpSolver->dObj; DSDP_INT i;
    
    if (inCone) { // Not sure if y is in the cone
        DSDP_INT psd = FALSE;
        for (i = 0; i < dsdpSolver->m; ++i) {
            su = bound - y->x[i]; sl = y->x[i] + bound;
            logdet += (log(su) + log(sl));
        }
        getPhaseBS(dsdpSolver, y->x); dsdpInCone(dsdpSolver, &psd);
        if (!psd) { *inCone = FALSE; return 0.0; }
        *inCone = TRUE;
    } else { // No check
        for (i = 0; i < dsdpSolver->m; ++i) {
            su = bound - y->x[i]; sl = y->x[i] + bound;
            logdet += (log(su) + log(sl));
        }
        getPhaseBS(dsdpSolver, y->x);
        for (i = 0; i < dsdpSolver->nBlock; ++i) {
            spsMatFactorize(dsdpSolver->S[i]);
        }
    }
    
    for (DSDP_INT i = 0; i < dsdpSolver->nBlock; ++i) {
        logdet += spsMatGetlogdet(dsdpSolver->S[i], aux);
    }
    
    vec_dot(b, y, &logbarrier);
    logbarrier += dsdpSolver->mu * logdet;
    return (-logbarrier);
}

static double getCurrentLogdet( HSDSolver *dsdpSolver, vec *y, double tau, DSDP_INT pRelax, DSDP_INT *inCone ) {
    
    double logbarrier = 0.0, logdet = 0.0, *aux = dsdpSolver->M->schurAux, sl, su, bound = dsdpSolver->ybound;
    DSDP_INT i;
    
    if (inCone) { // Not sure if y is in the cone
        DSDP_INT psd = FALSE;
        if (pRelax) {
            for (i = 0; i < dsdpSolver->m; ++i) {
                su = bound - y->x[i]; sl = y->x[i] + bound;
                logdet += (log(su) + log(sl));
            }
        }
        getPhaseAS(dsdpSolver, y->x, tau); dsdpInCone(dsdpSolver, &psd);
        if (!psd) { *inCone = FALSE; return 0.0; }
        *inCone = TRUE;
    } else { // No check
        
        if (pRelax) {
            for (i = 0; i < dsdpSolver->m; ++i) {
                su = bound - y->x[i]; sl = y->x[i] + bound;
                logdet += (log(su) + log(sl));
            }
        }
        getPhaseAS(dsdpSolver, y->x, tau);
        for (i = 0; i < dsdpSolver->nBlock; ++i) {
            spsMatFactorize(dsdpSolver->S[i]);
        }
    }
    
    for (DSDP_INT i = 0; i < dsdpSolver->nBlock; ++i) {
        logdet += spsMatGetlogdet(dsdpSolver->S[i], aux);
    }
    
    logbarrier = dsdpSolver->mu * logdet;
    return (-logbarrier);
}

static DSDP_INT adjCorrectorStep( HSDSolver *dsdpSolver ) {
    // Automatically adjust the number of currector steps by problem structure
    
    DSDP_INT nusercorr;
    DSDPGetIntParam(dsdpSolver, INT_PARAM_BCORRECTOR, &nusercorr);
    
    if (nusercorr == 0) { return 0; }
    if (dsdpSolver->Pnrm < 0.1) { nusercorr = 0; }
    if (dsdpSolver->alpha < 0.1 && dsdpSolver->mu < 1e-05) { nusercorr = 0; }
    if (dsdpSolver->mu < 1e-06) { nusercorr = 0; }
    
    return nusercorr;
}

extern DSDP_INT dInfeasCorrectorStep( HSDSolver *dsdpSolver, DSDP_INT isfinal )  {
    // Corrector step in Phase A.
    DSDP_INT retcode = DSDP_RETCODE_OK;
    CGSolver *cg = dsdpSolver->cgSolver; dsdpCGSetMaxIter(cg, 300); dsdpCGSetTol(cg, 1e-04);
    
    DSDP_INT prelax = (dsdpSolver->ybound != DSDP_INFINITY);
    vec *d2 = dsdpSolver->M->asinv; // asinv is directed solved using CG
    vec *d1 = dsdpSolver->d1, *ynew = dsdpSolver->d4; // d4 is reused for storing ynew
    vec *dycorr = dsdpSolver->b1, *sl = dsdpSolver->sl, *su = dsdpSolver->su, *ynow = dsdpSolver->y;
    double oldbarrier, newbarrier = 0.0, step, tmp, bound, rhon;
    
    bound = dsdpSolver->ybound; DSDPGetDblParam(dsdpSolver, DBL_PARAM_RHON, &rhon);
    DSDP_INT ncorrector, inCone = FALSE; DSDPGetIntParam(dsdpSolver, INT_PARAM_ACORRECTOR, &ncorrector);
    
    
    for (DSDP_INT i = 0, j; i < ncorrector && dsdpSolver->Ry; ++i) {
        
        if (i == 0) {
            for (j = 0; j < dsdpSolver->nBlock; ++j) {
                spsMatFactorize(dsdpSolver->S[j]);
            }
        }
        
        asinvSetup(dsdpSolver->M);  // Compute asinv
        
        if (prelax) {
            vec_lslack(ynow, sl, -bound); vec_uslack(ynow, su, bound);
            for (j = 0; j < dsdpSolver->m; ++j) {
                tmp = su->x[j]; d2->x[j] += 1.0 / tmp;
                tmp = sl->x[j]; d2->x[j] -= 1.0 / tmp;
            }
        }
        
        dsdpCGSolve(cg, d2, NULL); // Compute corrector step
        vec_zaxpby(dycorr, 0.0, d1, -1.0, d2);
        oldbarrier = (i == 0) ? getCurrentLogdet(dsdpSolver, ynow, dsdpSolver->tau, prelax, NULL) : newbarrier;
        
        getPhaseBdS(dsdpSolver, 1.0, dycorr->x, 0.0);
                
        step = DSDP_INFINITY;
        for (j = 0; j < dsdpSolver->nBlock; ++j) {
            dsdpGetAlpha(dsdpSolver->lczSolver[j], dsdpSolver->S[j], dsdpSolver->dS[j],
                         dsdpSolver->spaux[j], &tmp);
            step = MIN(step, tmp);
        }
        
        if (prelax) {
            tmp = vec_step(sl, dycorr,  1.0); step = MIN(tmp, step);
            tmp = vec_step(su, dycorr, -1.0); step = MIN(tmp, step);
        }
        
        
        if (isfinal) {
            step = 0.9 * step; // Update Pnrm ?
        } else {
            step = MIN(0.95 * step, 1.0); // Update Pnrm ?
            step = MIN(step, rhon / dsdpSolver->Pnrm);
        }
        
        // Do line search
        for (j = 0; ; ++j) {
            vec_zaxpby(ynew, step, dycorr, 1.0, ynow);
            if (j == 0) {
                newbarrier = getCurrentLogdet(dsdpSolver, ynew, dsdpSolver->tau, prelax, &inCone);
                if (!inCone) {
                    step *= 0.5; --j;
                    if (step <= 1e-04) { break; }
                    continue;
                }
            } else {
                newbarrier = getCurrentLogdet(dsdpSolver, ynow, dsdpSolver->tau, prelax, NULL);
            }
            
            if (step <= 1e-04 || newbarrier <= oldbarrier || isfinal) {
                break;
            } else {
                step *= 0.5;
            }
        }
        
        if (step <= 1e-04) { break; }
        vec_copy(ynew, ynow); getPhaseAS(dsdpSolver, ynow->x, dsdpSolver->tau);
    }
    
    dsdpSolver->iterProgress[ITER_CORRECTOR] = TRUE;
    return retcode;
}

extern DSDP_INT dualCorrectorStep( HSDSolver *dsdpSolver ) {
    // Dual corrector step in Phase B. When into this routine,
    // S is not inverted and M is ready
    DSDP_INT retcode = DSDP_RETCODE_OK;
    retcode = checkIterProgress(dsdpSolver, ITER_CORRECTOR);
    
    // Prepare CG solver
    CGSolver *cg = dsdpSolver->cgSolver; dsdpCGSetMaxIter(cg, 300); dsdpCGSetTol(cg, 1e-04);
    
    DSDP_INT nall = dsdpSolver->nall;
    vec *d2 = dsdpSolver->M->asinv; // asinv is directed solved using CG
    vec *d1 = dsdpSolver->d1, *b = dsdpSolver->dObj, *ynew = dsdpSolver->d4; // d4 is reused for storing ynew
    vec *dycorr = dsdpSolver->b1, *sl = dsdpSolver->sl, *su = dsdpSolver->su, *ynow = dsdpSolver->y;
    double shrink = nall / (nall + sqrt(nall)), bTd1, bTd2, bTdycorr, aval;
    double oldbarrier, newbarrier, step, tmp, bound, rhon;
    
    bound = dsdpSolver->ybound; DSDPGetDblParam(dsdpSolver, DBL_PARAM_RHON, &rhon);
    DSDP_INT ncorrector = adjCorrectorStep(dsdpSolver), inCone = FALSE;
    vec_dot(b, d1, &bTd1);
    
    
    for (DSDP_INT i = 0, j; i < ncorrector; ++i) {
        
        if (dsdpSolver->mu < 1e-05) {
            break;
        }
        
        if (i == 0) {
            for (j = 0; j < dsdpSolver->nBlock; ++j) {
                spsMatFactorize(dsdpSolver->S[j]);
            }
        }
        
        asinvSetup(dsdpSolver->M);  // Compute asinv
        
        vec_lslack(ynow, sl, -bound); vec_uslack(ynow, su, bound);
        for (j = 0; j < dsdpSolver->m; ++j) {
            tmp = su->x[j]; d2->x[j] += 1.0 / tmp;
            tmp = sl->x[j]; d2->x[j] -= 1.0 / tmp;
        }
        
        dsdpCGSolve(cg, d2, NULL); // Compute corrector step
        
        vec_dot(b, d2, &bTd2);
        if (bTd1 > 0 && bTd2 > 0) {
            dsdpSolver->mu = MIN(dsdpSolver->mu, bTd1 / bTd2);
//             dsdpSolver->mu = MAX(dsdpSolver->schurmu / rhon, dsdpSolver->mu);
        }
        dsdpSolver->mu *= shrink;
        
        // dycorr = dy1 / mu - dy2;
        vec_zaxpby(dycorr, 1 / dsdpSolver->mu, d1, -1.0, d2);
        vec_dot(b, dycorr, &bTdycorr);
        oldbarrier = getCurrentLogBarrier(dsdpSolver, ynow, NULL);
        getPhaseBdS(dsdpSolver, 1.0, dycorr->x, 0.0);
                
        step = DSDP_INFINITY;
        for (j = 0; j < dsdpSolver->nBlock; ++j) {
            dsdpGetAlpha(dsdpSolver->lczSolver[j], dsdpSolver->S[j], dsdpSolver->dS[j],
                         dsdpSolver->spaux[j], &tmp);
            step = MIN(step, tmp);
        }
        
        tmp = vec_step(sl, dycorr,  1.0); step = MIN(tmp, step);
        tmp = vec_step(su, dycorr, -1.0); step = MIN(tmp, step);
        
        if (step >= 5.0) {
            vec_zaxpby(ynew, 1.0, dycorr, 1.0, ynow);
            vec_copy(ynew, ynow); getPhaseBS(dsdpSolver, ynow->x);
            continue;
        }
        
        step = MIN(0.95 * step, 1.0);
        step = MIN(step, rhon / dsdpSolver->Pnrm);
        
        // Do line search
        for (j = 0; ; ++j) {
            vec_zaxpby(ynew, step, dycorr, 1.0, ynow);
            if (j == 0) {
                newbarrier = getCurrentLogBarrier(dsdpSolver, ynew, &inCone);
                if (!inCone) {
                    step *= 0.5; --j;
                    if (step <= 1e-04) { break; }
                    continue;
                }
            } else {
                newbarrier = getCurrentLogBarrier(dsdpSolver, ynew, NULL);
            }
            
            if (step <= 1e-04 || newbarrier <= oldbarrier - fabs(0.05 * bTdycorr * step)) {
                break;
            } else {
                aval = 2 * (newbarrier - oldbarrier + bTdycorr * step) / (step * step);
                if ((bTdycorr / aval < step) && (bTdycorr / aval > 0)) {
                    step = bTdycorr / aval;
                } else {
                    step *= 0.5;
                }
            }
        }
        
        if (step <= 1e-04) { break; }
        vec_copy(ynew, ynow); getPhaseBS(dsdpSolver, ynow->x);
    }
    
    dsdpSolver->iterProgress[ITER_CORRECTOR] = TRUE;
    return retcode;
}
