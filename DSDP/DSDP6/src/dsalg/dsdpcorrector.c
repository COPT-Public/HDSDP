#include "dsdpcorrector.h"
#include "dsdputils.h"

/*
 Implement the dual corrector for DSDP Phase B
*/

static char etype[] = "Dual corrector step";

static double getCurrentLogBarrier( HSDSolver *dsdpSolver, vec *y, DSDP_INT *inCone ) {
    
    double logbarrier = 0.0, logdet = 0.0, *aux = dsdpSolver->M->schurAux, sl, su, bound;
    DSDPGetDblParam(dsdpSolver, DBL_PARAM_PRLX_PENTALTY, &bound);
    vec *b = dsdpSolver->dObj; DSDP_INT i;
    
    if (inCone) { // Not sure if y is in the cone
        DSDP_INT psd = FALSE;
        for (i = 0; i < dsdpSolver->m; ++i) {
            su = bound - y->x[i]; sl = y->x[i] + bound;
            if (sl <= 0 || su <= 0) {
                *inCone = FALSE; return 0.0;
            }
            logdet += (log(su) + log(sl));
        }
        getPhaseBS(dsdpSolver, y->x); dsdpInCone(dsdpSolver, &psd);
        if (!psd) { *inCone = FALSE; return 0.0; }
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

static DSDP_INT adjCorrectorStep( HSDSolver *dsdpSolver ) {
    // Automatically adjust the number of currector steps by problem structure
    
    DSDP_INT m = dsdpSolver->m, n = dsdpSolver->n, nusercorr;
    DSDPGetIntParam(dsdpSolver, INT_PARAM_BCORRECTOR, &nusercorr);
    
    if (nusercorr == 0) { return 0; }
    if (dsdpSolver->Pnrm < 0.1) { nusercorr = 0; }
    if (dsdpSolver->alpha < 0.1 && dsdpSolver->mu < 1e-05) { nusercorr = 0; }
    if (dsdpSolver->mu < 1e-06) { nusercorr = 0; }
    
    if (n >= m) {
        nusercorr = MIN(nusercorr, 4);
    } else if (n >= 0.5 * m) {
        nusercorr = MIN(nusercorr, 6);
    } else if (n >= 0.1 * m) {
        nusercorr = MIN(nusercorr, 8);
    }
    
    if (m > 5 * n) {
        nusercorr = MAX(nusercorr, 4);
    } else if (m > 10 * n) {
        nusercorr = MAX(nusercorr, 8);
    } else if (m > 20 * n) {
        nusercorr = MAX(nusercorr, 12);
    }
    
    nusercorr = MIN(nusercorr, 12);
    return nusercorr;
}

extern DSDP_INT dualCorrectorStep( HSDSolver *dsdpSolver ) {
    // Dual corrector step in Phase B. When into this routine,
    // S is not inverted and M is ready
    DSDP_INT retcode = DSDP_RETCODE_OK;
    retcode = checkIterProgress(dsdpSolver, ITER_CORRECTOR);
    
    // Prepare CG solver
    CGSolver *cg = dsdpSolver->cgSolver; dsdpCGSetMaxIter(cg, 300); dsdpCGSetTol(cg, 1e-04);
    
    DSDP_INT nall = dsdpSolver->n + dsdpSolver->m * 2;
    vec *d2 = dsdpSolver->M->asinv; // asinv is directed solved using CG
    vec *d1 = dsdpSolver->d1, *b = dsdpSolver->dObj, *ynew = dsdpSolver->d4; // d4 is reused for storing ynew
    vec *dycorr = dsdpSolver->b1, *sl = dsdpSolver->sl, *su = dsdpSolver->su, *ynow = dsdpSolver->y;
    double shrink = nall / (nall + sqrt(nall)), bTd1, bTd2, bTdycorr, aval;
    double oldbarrier, newbarrier, step, tmp, bound, rhon;
    
    DSDPGetDblParam(dsdpSolver, DBL_PARAM_PRLX_PENTALTY, &bound);
    DSDPGetDblParam(dsdpSolver, DBL_PARAM_RHON, &rhon);
    
    DSDP_INT ncorrector = adjCorrectorStep(dsdpSolver), inCone = FALSE;
    vec_dot(b, d1, &bTd1);
    
    for (DSDP_INT i = 0, j; i < ncorrector; ++i) {
        
        if (dsdpSolver->mu < 1e-05) {
            break;
        }
        
        if (i == 0) { dsdpInCone(dsdpSolver, &inCone); assert(inCone); }
        vec_lslack(ynow, sl, -bound); vec_uslack(ynow, su, bound);
        asinvSetup(dsdpSolver->M);  // Compute asinv
        
        for (j = 0; j < dsdpSolver->m; ++j) {
            tmp = su->x[j]; d2->x[j] += 1.0 / tmp;
            tmp = sl->x[j]; d2->x[j] -= 1.0 / tmp;
        }
        
        dsdpCGSolve(cg, d2, NULL); // Compute corrector step
        
        vec_dot(b, d2, &bTd2);
        if (bTd1 > 0 && bTd2 > 0) {
            dsdpSolver->mu = MIN(dsdpSolver->mu, bTd1 / bTd2);
        }
        dsdpSolver->mu *= shrink; // May be placed outside the loop
        
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
            
            if (step <= 1e-04 || newbarrier <= oldbarrier - fabs(0.1 * bTdycorr * step)) {
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
