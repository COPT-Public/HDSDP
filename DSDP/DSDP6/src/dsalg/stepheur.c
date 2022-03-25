#include "stepheur.h"
#include "dsdputils.h"
/*
   Implement the strategies to choose the proper stepsize and update the
   dual iteration variables
*/
#ifdef ROUGH_MU
#undef ROUGH_MU
#endif

static char etype[] = "Step size computation";

void spsMatExport( spsMat *A );

static DSDP_INT getKappaTauStep( HSDSolver *dsdpSolver, double *kappatauStep ) {
    // Compute the maximum step size to take at kappa and tau
    DSDP_INT retcode = DSDP_RETCODE_OK;
    double Aalpha, tau, kappa, dtau, dkappa, step;
    retcode = DSDPGetDblParam(dsdpSolver, DBL_PARAM_AALPHA, &Aalpha);
    tau = dsdpSolver->tau; kappa = dsdpSolver->kappa;
    dtau = dsdpSolver->dtau; dkappa = dsdpSolver->dkappa;
    step = MIN((dtau / tau), (dkappa / kappa));
    *kappatauStep = (step < 0.0) ? fabs(Aalpha / step) : DSDP_INFINITY;
    return retcode;
}

static DSDP_INT takeKappaTauStep( HSDSolver *dsdpSolver ) {
    // Take step in kappa and tau
    DSDP_INT retcode = DSDP_RETCODE_OK;
    double step = dsdpSolver->alpha;
    
    dsdpSolver->tau = dsdpSolver->tau + step * dsdpSolver->dtau;
    dsdpSolver->kappa = dsdpSolver->kappa + step * dsdpSolver->dkappa;
    
    return retcode;
}

static DSDP_INT takeyStep( HSDSolver *dsdpSolver ) {
    // Take step in y
    DSDP_INT retcode = DSDP_RETCODE_OK;
    
    double step = dsdpSolver->alpha;
    vec *y  = dsdpSolver->y, *dy = dsdpSolver->dy;
    vec_axpy(step, dy, y);
    
    return retcode;
}

static DSDP_INT takelpsStep( HSDSolver *dsdpSolver ) {
    // Take step in LP s
    double step = dsdpSolver->alpha;
    vec *s  = dsdpSolver->s, *ds = dsdpSolver->ds;
    vec_axpy(step, ds, s);
    return DSDP_RETCODE_OK;
}

static DSDP_INT takeSDPSStep( HSDSolver *dsdpSolver ) {
    // Take step in SDP S
    return getPhaseAS(dsdpSolver, dsdpSolver->y->x, dsdpSolver->tau);;
}

static DSDP_INT getLPsStep( HSDSolver *dsdpSolver, double *sStep ) {
    // Compute the maixmum step size to take at s for LP
    DSDP_INT retcode = DSDP_RETCODE_OK;
    
    vec *s  = dsdpSolver->s, *ds = dsdpSolver->ds;
    
    double step = 100.0;
    
    for (DSDP_INT i = 0; i < s->dim; ++i) {
        step = MIN((ds->x[i] / s->x[i]), step);
    }
    
    if (step < 0.0) {
        *sStep = fabs(0.995 / step);
    }
    
    return retcode;
}

static DSDP_INT getBlockSDPSStep( HSDSolver *dsdpSolver, DSDP_INT k, double *SkStep ) {
    // Compute the maixmum step size to take at S for some cone
    // Note that we compute the step utilizing the Cholesky factorization of S by
    // S + a * dS = L (I + a L^-1 dS L^-T) L^T
    DSDP_INT retcode = DSDP_RETCODE_OK;
    assert( k < dsdpSolver->nBlock );
    retcode = dsdpGetAlpha(dsdpSolver->lczSolver[k], dsdpSolver->S[k], dsdpSolver->dS[k],
                           dsdpSolver->spaux[k], SkStep);
    checkCode;
    
    return retcode;
}

static DSDP_INT getSDPSStep( HSDSolver *dsdpSolver, double *SStep ) {
    // Compute the maximum step size to take for the SDP cones
    DSDP_INT retcode = DSDP_RETCODE_OK;
    DSDP_INT nblock = dsdpSolver->nBlock;
    
    double res = 0.0;
    double step = DSDP_INFINITY;
    
    for (DSDP_INT i = 0; i < nblock; ++i) {
        retcode = getBlockSDPSStep(dsdpSolver, i, &res);
        step = MIN(step, res);
    }
    
    *SStep = step;
    return retcode;
}

static double getBoundyStep( HSDSolver *dsdpSolver ) {
    // Compute the maimum stepsize to take for the bound cones
    vec *y = dsdpSolver->y, *dy = dsdpSolver->dy;
    double bound, yi, dyi, step = DSDP_INFINITY;
    DSDPGetDblParam(dsdpSolver, DBL_PARAM_PRLX_PENTALTY, &bound);
    for (DSDP_INT i = 0; i < y->dim; ++i) {
        yi = y->x[i]; dyi = dy->x[i];
        if (dyi == 0.0) continue;
        step = (dyi > 0.0) ? MIN(step, (bound - yi) / dyi) : MIN(step, (- bound - yi) / dyi);
    }
    return step;
}

static DSDP_INT getCurrentyPotential( HSDSolver *dsdpSolver, vec *y,
                                      double rho, double *potential, DSDP_INT *inCone ) {
    // Compute the current potential
    // phi(y) = rho * log(pObj - dObj) - log det S - \sum log sl_i - \sum log su_i
    DSDP_INT retcode = DSDP_RETCODE_OK;
    
    double pval = 0.0, dObjVal = 0.0, *aux = dsdpSolver->M->schurAux;
    DSDP_INT i; double sl, su, bound;
    DSDPGetDblParam(dsdpSolver, DBL_PARAM_PRLX_PENTALTY, &bound);
    
    if (inCone) {
        // Not sure whether y is in the cone
        DSDP_INT psd = FALSE;
        // Check and get potential of the bound cone
        for (i = 0; i < dsdpSolver->m; ++i) {
            su = bound - y->x[i]; sl = y->x[i] + bound;
            if (sl <= 0 || su <= 0) {
                *inCone = FALSE; return retcode;
            }
            pval -= (log(su) + log(sl));
        }
        
        retcode = getPhaseBS(dsdpSolver, y->x); dsdpInCone(dsdpSolver, &psd);
        
        if (psd) {
            *inCone = TRUE; vec_dot(dsdpSolver->dObj, y, &dObjVal);
            pval += rho * log(dsdpSolver->pObjVal - dObjVal);
            for (i = 0; i < dsdpSolver->nBlock; ++i) {
                pval -= spsMatGetlogdet(dsdpSolver->S[i], aux);
            }
            dsdpSolver->dPotential = pval; *potential = pval;
        } else {
            *inCone = FALSE;
        }
        return retcode;
    }
    
    if (y) {
        // Get potential of a new y that is in the cone
        retcode = getPhaseBS(dsdpSolver, y->x);
        vec_dot(dsdpSolver->dObj, y, &dObjVal);
        pval = rho * log(dsdpSolver->pObjVal - dObjVal);
        
        // Bound cone
        for (i = 0; i < dsdpSolver->m; ++i) {
            su = bound - y->x[i]; sl = y->x[i] + bound;
            pval -= (log(su) + log(sl));
        }
        
        for (i = 0; i < dsdpSolver->nBlock; ++i) {
            spsMatFactorize(dsdpSolver->S[i]);
            pval -= spsMatGetlogdet(dsdpSolver->S[i], aux);
        }
        
    } else {
        
        vec *yold = dsdpSolver->y;
        // Get potential of new rho (old y)
        pval = rho * log(dsdpSolver->pObjVal - dsdpSolver->dObjVal);
        for (i = 0; i < dsdpSolver->m; ++i) {
            su = bound - yold->x[i]; sl = yold->x[i] + bound;
            pval -= (log(su) + log(sl));
        }
        
        for (i = 0; i < dsdpSolver->nBlock; ++i) {
            pval -= spsMatGetlogdet(dsdpSolver->S[i], aux);
        }
        
        dsdpSolver->dPotential = pval;
    }
    
    *potential = pval;
    return retcode;
}

extern DSDP_INT getMaxStep( HSDSolver *dsdpSolver ) {
    // Compute the maximum step size to take for one iteration
    DSDP_INT retcode = DSDP_RETCODE_OK;
    
    retcode = checkIterProgress(dsdpSolver, ITER_COMPUTE_STEP);
    
    assert( !dsdpSolver->iterProgress[ITER_COMPUTE_STEP] );
    if (dsdpSolver->iterProgress[ITER_COMPUTE_STEP]) {
        error(etype, "Stepsize has been computed. \n");
    }
    
    double stepkappatau = 0.0, steplps = 100.0, sdpS = 0.0, Aalpha;
    retcode = DSDPGetDblParam(dsdpSolver, DBL_PARAM_AALPHA, &Aalpha);
    retcode = getKappaTauStep(dsdpSolver, &stepkappatau);
    // retcode = getLPsStep(dsdpSolver, &steplps);
    retcode = getSDPSStep(dsdpSolver, &sdpS); checkCode;
    sdpS = MIN(sdpS, steplps);
    dsdpSolver->alpha = MIN(sdpS * Aalpha, stepkappatau);
    dsdpSolver->alpha = MIN(dsdpSolver->alpha, 1.0);
    
    dsdpSolver->iterProgress[ITER_COMPUTE_STEP] = TRUE;
    
    return retcode;
}

extern DSDP_INT takeStep( HSDSolver *dsdpSolver ) {
    // Take step towards next iterate
    DSDP_INT retcode = DSDP_RETCODE_OK;
    
    retcode = checkIterProgress(dsdpSolver, ITER_TAKE_STEP);
    assert( !dsdpSolver->iterProgress[ITER_TAKE_STEP] );
    
    if (dsdpSolver->iterProgress[ITER_TAKE_STEP]) {
        error(etype, "Step has been taken. \n");
    }
    
    retcode = takeKappaTauStep(dsdpSolver);
    retcode = takeyStep(dsdpSolver);
    // retcode = takelpsStep(dsdpSolver, step);
    retcode = takeSDPSStep(dsdpSolver); checkCode;
    dsdpSolver->iterProgress[ITER_TAKE_STEP] = TRUE;
    
    return retcode;
}

extern DSDP_INT selectMu( HSDSolver *dsdpSolver, double *newmu ) {
    // Choose the next barrier parameter
    /*
     The backward newton step is stored in b2 and Scker, dy1 is in d1; dy is in b1
     
     At this stage, if a new primal feasible solution if found, then sl and su are filled by
     backward newton steps. Otherwise sl and su are untouched
     
     */
    
    DSDP_INT retcode = DSDP_RETCODE_OK;
    
    double alpha = DSDP_INFINITY, alphap = 0.0, tmp = 1000;
    
    if (dsdpSolver->eventMonitor[EVENT_PFEAS_FOUND]) {
        retcode = getPhaseBdS(dsdpSolver, -1.0 / dsdpSolver->mu, dsdpSolver->d1->x, 0.0);
        for (DSDP_INT i = 0; i < dsdpSolver->nBlock; ++i) {
            dsdpGetAlpha(dsdpSolver->lczSolver[i], dsdpSolver->Scker[i],
                         dsdpSolver->dS[i], dsdpSolver->spaux[i], &tmp);
            alpha = MIN(alpha, tmp);
        }
        
        // Compute step of bound cone
        alpha = MIN(vec_step(dsdpSolver->su, dsdpSolver->d1, 1.0 / dsdpSolver->mu), alpha);
        alpha = MIN(vec_step(dsdpSolver->sl, dsdpSolver->d1, - 1.0 / dsdpSolver->mu), alpha);
        *newmu = dsdpSolver->mu / (1 + alpha * 0.95);
        
    } else {
        // dS = dsdpgetATy(A, dy);
        retcode = getPhaseBdS(dsdpSolver, -1.0, dsdpSolver->b1->x, 0.0);
        // alphap = dsdpgetalpha(S, dS, 0.95 / 1.0);
        for (DSDP_INT i = 0; i < dsdpSolver->nBlock; ++i) {
            dsdpGetAlpha(dsdpSolver->lczSolver[i], dsdpSolver->S[i],
                         dsdpSolver->dS[i], dsdpSolver->spaux[i], &tmp);
            alpha = MIN(alpha, tmp);
        }
        
        // sl and su are still untouched
        alpha = MIN(vec_step(dsdpSolver->su, dsdpSolver->b1,  1.0), alpha);
        alpha = MIN(vec_step(dsdpSolver->sl, dsdpSolver->b1, -1.0), alpha);
        
        alphap = alpha;
        // Shat = S + 0.95 * alphap * dS;
        for (DSDP_INT i = 0; i < dsdpSolver->nBlock; ++i) {
            memcpy(dsdpSolver->Scker[i]->x, dsdpSolver->S[i]->x,
                   sizeof(double) * dsdpSolver->S[i]->nnz);
            // This step sometimes fails due to inaccurate Lanczos
            spsMataXpbY(MIN(0.95 * alphap, 1.0), dsdpSolver->dS[i],
                        1.0, dsdpSolver->Scker[i], dsdpSolver->symS[i]);
        }
        
        // Overwrite sl and su by slhat and suhat
        vec_axpy(- MIN(0.95 * alphap, 1.0), dsdpSolver->b1, dsdpSolver->sl);
        vec_axpy(  MIN(0.95 * alphap, 1.0), dsdpSolver->b1, dsdpSolver->su);

        for (DSDP_INT i = 0; i < dsdpSolver->nBlock; ++i) {
            spsMatFactorize(dsdpSolver->Scker[i]);
        }
        
        // dS = - alphap * dsdpgetATy(A, dy1) / muk;
        getPhaseBdS(dsdpSolver, alphap / dsdpSolver->mu, dsdpSolver->d1->x, 0.0);
        
        tmp = DSDP_INFINITY;
        for (DSDP_INT i = 0; i < dsdpSolver->nBlock; ++i) {
            dsdpGetAlpha(dsdpSolver->lczSolver[i], dsdpSolver->Scker[i], dsdpSolver->dS[i],
                         dsdpSolver->spaux[i], &alpha);
            tmp = MIN(tmp, alpha);
        }
        
        tmp = MIN(vec_step(dsdpSolver->sl, dsdpSolver->d1,  alphap / dsdpSolver->mu), tmp);
        tmp = MIN(vec_step(dsdpSolver->su, dsdpSolver->d1, -alphap / dsdpSolver->mu), tmp);
        
//        tmp = MIN(tmp, 1000.0);
        
        *newmu = (alphap * dsdpSolver->mu) / (1 + tmp) + \
                  + (1 - alphap) * (dsdpSolver->pObjVal - dsdpSolver->dObjVal) / (dsdpSolver->n + dsdpSolver->m * 2);
        
        // tmp = MIN(1, 0.97 * tmp); tmp = dsdpSolver->mu / (1.0 + alphap); alphap *= 0.6;
        // *newmu = alphap * tmp + (1.0 - alphap) * dsdpSolver->mu;
    }
    
    return retcode;
}

extern DSDP_INT dualPotentialReduction( HSDSolver *dsdpSolver ) {
    
    // Implement the dual potential reduction method
    // dy is filled
    DSDP_INT retcode = DSDP_RETCODE_OK;
    retcode = checkIterProgress(dsdpSolver, ITER_COMPUTE_STEP);
    
    double rho, oldpotential = 0.0, maxstep = 0.0;
    DSDPSetDblParam(dsdpSolver, DBL_PARAM_RHO,
                    (dsdpSolver->pObjVal - dsdpSolver->dObjVal) / dsdpSolver->mu);
    getDblParam(dsdpSolver->param, DBL_PARAM_RHO, &rho);
    
    vec *ytarget = dsdpSolver->d4, *y = dsdpSolver->y, *dy = dsdpSolver->dy;
    getCurrentyPotential(dsdpSolver, NULL, rho, &oldpotential, NULL);
    getSDPSStep(dsdpSolver, &maxstep);
    maxstep = MIN(maxstep, getBoundyStep(dsdpSolver));
    maxstep = MIN(maxstep * 0.95, 1.0);
    
    double alpha = maxstep, newpotential = 0.0, bestpotential = oldpotential, bestalpha = alpha;
    DSDP_INT inCone = FALSE;
    
    // Start line search
    for (DSDP_INT i = 0; ; ++i) {
        // y = y + alpha * dy
        vec_zaxpby(ytarget, 1.0, y, alpha, dy);
        if (i == 0) {
            getCurrentyPotential(dsdpSolver, ytarget, rho, &newpotential, &inCone);
            if (!inCone) {
                alpha *= 0.8; bestalpha *= 0.8; --i;
                if (alpha <= 1e-03) {
                    break;
                }
                continue;
            }
        } else {
            getCurrentyPotential(dsdpSolver, ytarget, rho, &newpotential, NULL);
        }
        
        if (alpha <= 1e-02 || newpotential <= oldpotential - 0.05) {
            break;
        } else {
            if (newpotential < bestpotential) {
                bestalpha = alpha; bestpotential = newpotential;
            }
        }
        alpha *= 0.8;
    }
    
    if (alpha <= 1e-02) {
        alpha = bestalpha;
    }
    
    // Take step
    vec_axpy(alpha, dy, y); getPhaseBS(dsdpSolver, y->x);
    dsdpSolver->alpha = alpha;
    dsdpSolver->dPotential = bestpotential;
    
    dsdpSolver->iterProgress[ITER_COMPUTE_STEP] = TRUE;
    dsdpSolver->iterProgress[ITER_TAKE_STEP] = TRUE;
    
    return retcode;
}
