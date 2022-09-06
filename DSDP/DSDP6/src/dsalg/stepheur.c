#include "stepheur.h"
#include "dsdplapack.h"
#include "dsdplog.h"
#include "dsdputils.h"
#include "sparseopts.h"
#include "vec.h"

/*
   Implement the strategies to choose the proper stepsize and update the
   dual iteration variables
*/
#ifdef ROUGH_MU
#undef ROUGH_MU
#endif

static char etype[] = "Step size computation";

static DSDP_INT getKappaTauStep( HSDSolver *dsdpSolver, double *kappatauStep ) {
    // Compute the maximum step size to take at kappa and tau
    DSDP_INT retcode = DSDP_RETCODE_OK;
    
    if (!dsdpSolver->eventMonitor[EVENT_HSD_UPDATE]) {
        *kappatauStep = DSDP_INFINITY; return retcode;
    }
    
    double Aalpha = 0.9, tau, dtau, step;
    tau = dsdpSolver->tau; dtau = dsdpSolver->dtau;
#ifdef KAPPATAU
    double kappa, dkappa;
    kappa = dsdpSolver->kappa; dkappa = dsdpSolver->dkappa;
    step = MIN((dtau / tau), (dkappa / kappa));
#else
    step = dtau / tau;
#endif
    *kappatauStep = (step < 0.0) ? fabs(Aalpha / step) : DSDP_INFINITY;
    return retcode;
}

static DSDP_INT takeKappaTauStep( HSDSolver *dsdpSolver ) {
    // Take step in kappa and tau
    DSDP_INT retcode = DSDP_RETCODE_OK;
    
    if (!dsdpSolver->eventMonitor[EVENT_HSD_UPDATE]) {
        return retcode;
    }
    
    double step = dsdpSolver->alpha;
    dsdpSolver->tau = dsdpSolver->tau + step * dsdpSolver->dtau;
#ifdef KAPPATAU
    dsdpSolver->kappa = dsdpSolver->kappa + step * dsdpSolver->dkappa;
#else
    dsdpSolver->kappa = dsdpSolver->mu / dsdpSolver->tau;
#endif
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

static void takeSDPSStep( HSDSolver *dsdpSolver ) {
    // Take step in SDP S
    getPhaseAS(dsdpSolver, dsdpSolver->y, dsdpSolver->tau);;
}

static double getBoundyStep( HSDSolver *dsdpSolver ) {
    // Compute the maximum stepsize to take for the bound cones
    vec *y = dsdpSolver->y, *dy = dsdpSolver->dy;
    double bound, yi, dyi, tmp, step = DSDP_INFINITY;
    bound = dsdpSolver->ybound;
    
    if (dsdpSolver->eventMonitor[EVENT_IN_PHASE_A]) {
        double dtau = dsdpSolver->dtau, ds, tmpl = DSDP_INFINITY, tmpu = DSDP_INFINITY;
        double *sl = dsdpSolver->sl->x, *su = dsdpSolver->su->x;
        for (DSDP_INT i = 0; i < y->dim; ++i) {
            dyi = dy->x[i];
            ds =  dyi + bound * dtau; tmp = ds / sl[i]; tmpl = MIN(tmpl, tmp);
            ds = -dyi + bound * dtau; tmp = ds / su[i]; tmpu = MIN(tmpu, tmp);
        }
        step = (tmpl >= 0) ? step : MIN(step, - 1.0 / tmpl);
        step = (tmpu >= 0) ? step : MIN(step, - 1.0 / tmpu);
    } else {
        for (DSDP_INT i = 0; i < y->dim; ++i) {
            yi = y->x[i]; dyi = dy->x[i];
            if (dyi == 0.0) continue;
            if (dyi > 0.0) {
                tmp = (bound - yi) / dyi;
            } else{
                tmp = (- bound - yi) / dyi;
            }
            step = MIN(step, tmp);
        }
    }
    return step;
}

static DSDP_INT getCurrentyPotential( HSDSolver *dsdpSolver, vec *y,
                                      double rho, double *potential, DSDP_INT *inCone ) {
    // Compute the current potential
    // phi(y) = rho * log(pObj - dObj) - log det S - \sum log sl_i - \sum log su_i
    DSDP_INT retcode = DSDP_RETCODE_OK;
    
    double pval = 0.0, dObjVal = 0.0, *aux = dsdpSolver->M->schurAux;
    DSDP_INT i; double sl, su, bound = dsdpSolver->ybound, logdet;
    
    if (inCone) {
        logdet = DSDPConic( COPS_GET_LOGDET )(dsdpSolver, y, inCone);
        if (!*inCone) {
            *potential = 0.0; return retcode;
        }
        
        vec_dot(dsdpSolver->dObj, y, &dObjVal);
        pval = rho * log(dsdpSolver->pObjVal - dObjVal);
        pval -= logdet; *potential = pval;
        dsdpSolver->dPotential = pval;
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
        for (i = 0; i < dsdpSolver->lpDim; ++i) {
            pval -= log(dsdpSolver->s->x[i]);
        }
        dsdpSolver->dPotential = pval;
    }
    
    *potential = pval;
    return retcode;
}

extern void computeAdaptivedRate( HSDSolver *dsdpSolver ) {
    // Use heuristic to determine the rate for eliminating dual infeasibility
    double alphac = 0.0, alphainf = DSDP_INFINITY, tmp;
    // dSd3 = A^T * d3
    DSDPConic( COPS_CONSTR_EXPR )(dsdpSolver, DELTAS, 1.0, dsdpSolver->d3, 0.0, 0.0);
    alphac = DSDPConic( COPS_GET_MAXSTEP )(dsdpSolver, DUALVAR);
    alphac = MIN(alphac * 0.98, 1.0);
    DSDP_INT incone = FALSE, j;
    
    vec *aux = dsdpSolver->cg->aux, *y = dsdpSolver->y;
    for (j = 0; ; ++j) {
        vec_zaxpby(aux, -alphac, dsdpSolver->d3, 1.0, y);
        getPhaseACheckerS(dsdpSolver, aux, 1.0); getPhaseALpCheckers(dsdpSolver, aux, 1.0);
        dsdpCheckerInCone(dsdpSolver, &incone); dsdpLpCheckerInCone(dsdpSolver, &incone);
        if (incone) {
            break;
        } else {
            alphac *= 0.8;
        }
        if (alphac < 1e-02) {
            fatal_error_msg(etype);
        }
    }
    dsdpSolver->alpha = alphac;
    vec_copy(dsdpSolver->d4, dsdpSolver->d12);
    getPhaseAdS(dsdpSolver, 1.0, dsdpSolver->d12, 0.0);
    getPhaseALpds(dsdpSolver, 1.0, dsdpSolver->d12, 0.0);
    alphainf = getMaxSDPstep(dsdpSolver, CHECKER);
    tmp = getMaxLpstep(dsdpSolver, CHECKER);
    alphainf = MIN(alphainf, tmp);
    
    dsdpSolver->drate = MIN(1.0, 0.98 * alphainf / alphac);
    
    if (dsdpSolver->Pnrm < 1.0) {
        dsdpSolver->drate = MAX(0.9, dsdpSolver->drate);
    } else if (dsdpSolver->Pnrm < 10.0) {
        dsdpSolver->drate = MAX(0.3, dsdpSolver->drate);
    } else if (dsdpSolver->Pnrm < 50.0) {
        dsdpSolver->drate = MAX(0.1, dsdpSolver->drate);
    }
}

extern void getMaxStep( HSDSolver *dsdpSolver ) {
    // Compute the maximum step size to take for one iteration
    dsdpSolver->alpha = DSDPConic( COPS_GET_MAXSTEP )(dsdpSolver, DUALVAR);
    // MIN(dsdpSolver->alpha * XXX, 1.0): XXX is the most critical parameter
    dsdpSolver->alpha = MIN(dsdpSolver->alpha * 0.95, 1.0);
    dsdpSolver->iterProgress[ITER_COMPUTE_STEP] = TRUE;
}

extern void takeStep( HSDSolver *dsdpSolver ) {
    // Take step towards next iterate
    takeKappaTauStep(dsdpSolver);
    takeyStep(dsdpSolver);
    takeSDPSStep(dsdpSolver);
    DSDPConic( COPS_GET_SLACK )(dsdpSolver, DUALVAR);
    dsdpSolver->iterProgress[ITER_TAKE_STEP] = TRUE;
}

extern void selectMu( HSDSolver *dsdpSolver, double *newmu ) {
    // Choose the next barrier parameter
    /*
     The backward newton step is stored in b2 and Scker, dy1 is in d1; dy is in b1
     At this stage, if a new primal feasible solution if found, then sl and su are filled by
     backward newton steps. Otherwise sl and su are untouched
    */
        
    double alpha = DSDP_INFINITY, alphap = 0.0, tmp = 1000, tmp2;
    vec *aux = dsdpSolver->cg->aux;
    
    if (dsdpSolver->eventMonitor[EVENT_PFEAS_FOUND]) {
        // IMPORTANT: Is it correct ?
        getPhaseBdS(dsdpSolver, -1.0 / dsdpSolver->mu, dsdpSolver->d1, 0.0);
        getPhaseBLpds(dsdpSolver, -1.0 / dsdpSolver->mu, dsdpSolver->d1, 0.0);
        
        alpha = getMaxSDPstep(dsdpSolver, CHECKER);
        tmp2 = getMaxLpstep(dsdpSolver, CHECKER);
        alpha = MIN(alpha, tmp2);
        
        // Compute step of bound cone
        tmp2 = vec_step(dsdpSolver->su, dsdpSolver->d1, 1.0 / dsdpSolver->mu);
        alpha = MIN(tmp2, alpha);
        tmp2 = vec_step(dsdpSolver->sl, dsdpSolver->d1, - 1.0 / dsdpSolver->mu);
        alpha = MIN(tmp2, alpha);
        
        alpha = MIN(alpha * 0.97, 1000.0);
        *newmu = dsdpSolver->mu / (1.0 + alpha);
        
    } else {
        
        // dS = dsdpgetATy(A, dy);
        getPhaseBdS(dsdpSolver, -1.0, dsdpSolver->b1, 0.0);
        getPhaseBLpds(dsdpSolver, -1.0, dsdpSolver->b1, 0.0);
        
        alpha = getMaxSDPstep(dsdpSolver, DUALVAR);
        tmp2 = getMaxLpstep(dsdpSolver, DUALVAR);
        alpha = MIN(alpha, tmp2);

        // sl and su are still untouched
        tmp2 = vec_step(dsdpSolver->su, dsdpSolver->b1, 1.0);
        alpha = MIN(tmp2, alpha);
        tmp2 = vec_step(dsdpSolver->sl, dsdpSolver->b1, -1.0);
        alpha = MIN(tmp2, alpha);
        
        // Line-search
        alphap = (alpha < 1.0) ? MIN(1.0, 0.97 * alpha) : 1.0;
        vec_copy(dsdpSolver->sl, dsdpSolver->d12); // Backup
        vec_copy(dsdpSolver->su, dsdpSolver->d4);
        
        for (DSDP_INT i = 0, j, inCone = FALSE; ; ++i) {
            // Shat = S + alphap * dS = C - A' * y + alphap * A' * dy = C - A' * (y - alphap * dy)
            vec_zaxpby(aux, 1.0, dsdpSolver->y, -alphap, dsdpSolver->b1);
            getPhaseBCheckerS(dsdpSolver, aux);
            
            inCone = TRUE;
            if (dsdpSolver->isLPset) {
                vec_copy(dsdpSolver->s, dsdpSolver->scker);
                vec_axpy(alphap, dsdpSolver->ds, dsdpSolver->scker);
                inCone = vec_incone(dsdpSolver->scker);
            }
            
            if (inCone) {
                // Overwrite sl and su by slhat and suhat
                vec_axpy(- alphap, dsdpSolver->b1, dsdpSolver->sl);
                vec_axpy(  alphap, dsdpSolver->b1, dsdpSolver->su);
                inCone = vec_incone(dsdpSolver->sl);
                inCone = vec_incone(dsdpSolver->su);
            }
            
            if (inCone) {
                for (j = 0; j < dsdpSolver->nBlock; ++j) {
                    spsMatIspd(dsdpSolver->Scker[j], &inCone);
                    if (!inCone) { break; }
                }
            }

            if (!inCone) {
                alphap = (i > 2) ? 0.5 * alphap : 0.97 * alphap;
                vec_copy(dsdpSolver->d12, dsdpSolver->sl);
                vec_copy(dsdpSolver->d4, dsdpSolver->su);
            }
            
            if (inCone || alphap < 1e-08) { break; }
        }
        
        // IMPORTANT: Is it correct ?
        getPhaseBdS(dsdpSolver, -alphap / dsdpSolver->mu, dsdpSolver->d1, 0.0);
        getPhaseBLpds(dsdpSolver, -alphap / dsdpSolver->mu, dsdpSolver->d1, 0.0);
        
        tmp = DSDP_INFINITY;
        for (DSDP_INT i = 0; i < dsdpSolver->nBlock; ++i) {
            spsMatGetAlpha(dsdpSolver->lczs[i], dsdpSolver->Scker[i], dsdpSolver->dS[i], &alpha);
            tmp = MIN(tmp, alpha);
        }
        
        if (dsdpSolver->isLPset) {
            tmp2 = vec_step(dsdpSolver->scker, dsdpSolver->ds, 1.0);
            tmp = MIN(tmp, tmp2);
        }
        tmp2 = vec_step(dsdpSolver->sl, dsdpSolver->d1,  alphap / dsdpSolver->mu);
        tmp = MIN(tmp, tmp2);
        tmp2 = vec_step(dsdpSolver->su, dsdpSolver->d1, -alphap / dsdpSolver->mu);
        tmp = MIN(tmp, tmp2);
        
        tmp = MIN(0.97 * tmp, 1000.0);
        
        *newmu = (alphap * dsdpSolver->mu) / (1.0 + tmp) + \
                  + (1.0 - alphap) * (dsdpSolver->pObjVal - dsdpSolver->dObjVal) / dsdpSolver->nall;
    }
}

extern void dualPotentialReduction( HSDSolver *dsdpSolver ) {
    
    // Implement the dual potential reduction method
    // dy is filled
    double rho, alpha, oldpotential = 0.0, newpotential = 0.0, maxstep = 0.0, better = 0.0, tmp1;
    DSDP_INT inCone;
    getDblParam(dsdpSolver->param, DBL_PARAM_RHO, &rho);
    
    better = (dsdpSolver->Pnrm < 0.5) ? 0.0 : 0.05;
    vec *ytarget = dsdpSolver->d4, *y = dsdpSolver->y, *dy = dsdpSolver->dy;
    
    maxstep = getMaxSDPstep(dsdpSolver, DUALVAR);
    tmp1 = getBoundyStep(dsdpSolver);
    // getSDPSStep(dsdpSolver, &maxstep);
    maxstep = MIN(maxstep, tmp1);
    tmp1 = getMaxLpstep(dsdpSolver, DUALVAR);
    maxstep = MIN(maxstep, tmp1);
    
    if (maxstep > 1000.0 && FALSE) {
        alpha = 1.0;
    } else {
        
        inCone = TRUE;
        getCurrentyPotential(dsdpSolver, NULL, rho, &oldpotential, NULL);
        if (!inCone) {
            fatal_error_msg(etype);
        }
        maxstep = MIN(maxstep * 0.95, 1.0); alpha = maxstep;
        // Start line search
        for (DSDP_INT i = 0; ; ++i) {
            // y = y + alpha * dy
            vec_zaxpby(ytarget, 1.0, y, alpha, dy);
            if (i == 0) {
                getCurrentyPotential(dsdpSolver, ytarget, rho, &newpotential, &inCone);
                if (!inCone) {
                    alpha /= 3; --i;
                    if (alpha <= 1e-06) { break; assert( FALSE );}
                    continue;
                }
            } else {
                getCurrentyPotential(dsdpSolver, ytarget, rho, &newpotential, &inCone);
            }
            
            if (alpha <= 1e-04 || (newpotential <= oldpotential - better)
                || (alpha * dsdpSolver->Pnrm <= 0.001)) {
                break;
            }
            alpha *= 0.3;
        }
    }
    
    // Take step
    if (alpha <= 1e-04) {
        // vec_axpy(alpha, dy, y);
        getPhaseBS(dsdpSolver, y);
        dsdpInCone(dsdpSolver, &inCone); // TODO: Fix possible error in assert( inCone );
        getPhaseBLps(dsdpSolver, y);
        getBslack(dsdpSolver, y, DUALVAR);
    } else {
        vec_axpy(alpha, dy, y);
        getBslack(dsdpSolver, y, DUALVAR);
    }
    
    dsdpSolver->alpha = alpha;
    dsdpSolver->dPotential = newpotential;
    
    dsdpSolver->iterProgress[ITER_COMPUTE_STEP] = TRUE;
    dsdpSolver->iterProgress[ITER_TAKE_STEP] = TRUE;
}
