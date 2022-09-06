#include "dsdpcorrector.h"
#include "dsdplapack.h"
#include "dsdplog.h"
#include "dsdputils.h"
#include "adpcg.h"
#include "vec.h"

/*
 Implement the dual corrector for DSDP Phase B
*/

static char etype[] = "Dual corrector step";

static double getCurrentLogBarrier( HSDSolver *dsdpSolver, vec *y, DSDP_INT *inCone ) {
    
    double logbarrier = 0.0, logdet = 0.0;
    vec *b = dsdpSolver->dObj;
    logdet = DSDPConic( COPS_GET_LOGDET )(dsdpSolver, y, inCone);
    
    if (inCone) {
        if (!(*inCone)) { return 0.0; }
    }
    vec_dot(b, y, &logbarrier);
    logbarrier += dsdpSolver->mu * logdet;
    return (-logbarrier);
}

static DSDP_INT adjCorrectorStep( HSDSolver *dsdpSolver ) {
    // Automatically adjust the number of currector steps by problem structure
    DSDP_INT nusercorr;
    DSDPGetIntParam(dsdpSolver, INT_PARAM_BCORRECTOR, &nusercorr);
    
    if (nusercorr == 0) { return 0; }
    if (dsdpSolver->Pnrm < 0.1) { nusercorr = 0; }
    if ((dsdpSolver->alpha < 0.1 && dsdpSolver->mu < 1e-05) || (dsdpSolver->alpha < 1e-03)) {
        DSDPSetIntParam(dsdpSolver, INT_PARAM_BCORRECTOR, 0);
        nusercorr = 0;
    }
    if (dsdpSolver->mu < 1e-06) {
        nusercorr = 0;
        DSDPSetIntParam(dsdpSolver, INT_PARAM_BCORRECTOR, 0);
    }
    
    return nusercorr;
}

extern DSDP_INT dInfeasCorrectorStep( HSDSolver *dsdpSolver, DSDP_INT isfinal )  {
    // Corrector step in Phase A.
    DSDP_INT retcode = DSDP_RETCODE_OK;
    DSDP_INT nall = dsdpSolver->nall;
    
    adpcg *cg = dsdpSolver->cg;
    vec *d4 = dsdpSolver->M->asinvrysinv; // asinvrysinv is directly solved using CG
    vec *d2 = dsdpSolver->M->asinv; // asinv is directly solved using CG
    vec *d1 = dsdpSolver->d1, *ynew = dsdpSolver->d4; // d4 is reused for storing ynew
    vec *dycorr = dsdpSolver->dy, *sl = dsdpSolver->sl, *su = dsdpSolver->su, *ynow = dsdpSolver->y;
    double oldbarrier, newbarrier = 0.0, step, tmp, tmp2, tmp3, bound, rhon, dratemax, alphainf, drate, shrink = (double) nall / (nall + sqrt(nall));
    bound = dsdpSolver->ybound; DSDPGetDblParam(dsdpSolver, DBL_PARAM_RHON, &rhon);
    DSDP_INT ncorrector, inCone = FALSE;
    DSDPGetIntParam(dsdpSolver, INT_PARAM_ACORRECTOR, &ncorrector);
    dratemax = (isfinal) ? 0.0 : 0.8; // 0.6 37
    
    if (!dsdpSolver->Ry) {
        if (!DSDPConic( COPS_CHECK_INCONE )(dsdpSolver, DUALVAR)) {
            fatal_error_msg(etype);
        }
    }
    
    for (DSDP_INT i = 0, j; i < ncorrector && dsdpSolver->Ry; ++i) {
        
        if (i == 0) {
            if (!DSDPConic( COPS_CHECK_INCONE )(dsdpSolver, DUALVAR)) {
                fatal_error_msg(etype);
            }
        }
        
        DSDPConic( COPS_GET_SCHURVEC )(dsdpSolver, (DSDP_INT) (dratemax != 0));
        cg_solve(cg, d2, NULL);
        if (dratemax) {
            cg_solve(cg, d4, NULL);
        }
        vec_zaxpby(dycorr, 0.0, d1, -1.0, d2);

        if (i == 0) {
            oldbarrier = - dsdpSolver->mu * \
                         DSDPConic( COPS_GET_LOGDET )(dsdpSolver, ynow, NULL);
        } else {
            oldbarrier = newbarrier;
        }
         
        dsdpSolver->drate = 0.0;
        DSDPConic( COPS_STEPDIRECTION )(dsdpSolver);
        step = DSDPConic( COPS_GET_MAXSTEP )(dsdpSolver, DUALVAR);
        
        step = (isfinal) ? MIN(step, 1.0) : MIN(0.8 * step, 1.0); // 0.75 ? 0.8 ?
        // First line search
        for (j = 0; ; ++j) {
            vec_zaxpby(ynew, step, dycorr, 1.0, ynow);
            if (j == 0) {
                newbarrier = - dsdpSolver->mu * DSDPConic( COPS_GET_LOGDET )(dsdpSolver, ynew, &inCone);
                // newbarrier = getCurrentLogdet(dsdpSolver, ynew, dsdpSolver->tau, prelax, &inCone);
                if (!inCone) {
                    step *= 0.5; --j;
                    if (step <= 0.005) break;
                    continue;
                }
            } else { 
                newbarrier = - dsdpSolver->mu * DSDPConic( COPS_GET_LOGDET )(dsdpSolver, ynew, NULL); // getCurrentLogdet(dsdpSolver, ynew, dsdpSolver->tau, prelax, NULL);
                newbarrier = -1e+100;
            }
            
            if (step <= 0.005 || newbarrier <= oldbarrier || isfinal) {
                break;
            } else {
                step *= 0.5;
            }
        }
        
        if (step <= 0.005) {
            DSDPConic( COPS_GET_SLACK )(dsdpSolver, DUALVAR);
            if (!DSDPConic( COPS_CHECK_INCONE )(dsdpSolver, DUALVAR)) {
                fatal_error_msg(etype);
            }
            break;
        }
        
        if (dratemax) {
            alphainf = DSDP_INFINITY;
            // Bite off more infeasibility. Shat is now factorized
            vec_copy(d4, dycorr);
            
            // Refractor it
            getPhaseAdS(dsdpSolver, 1.0, dycorr, 0.0);
            getPhaseALpds(dsdpSolver, 1.0, dycorr, 0.0);
            tmp = getMaxSDPstep(dsdpSolver, DUALVAR);
            tmp2 = getMaxLpstep(dsdpSolver, DUALVAR);
            tmp = MIN(tmp, tmp2);
            alphainf = MIN(alphainf, tmp);
            
            // Check bound of y
            tmp = DSDP_INFINITY;
            for (j = 0; j < dsdpSolver->m; ++j) {
                tmp2 = dycorr->x[j];
                if (tmp2 == 0.0) { continue; }
                tmp3 = (tmp2 > 0.0) ? (bound - ynew->x[j]) / tmp2 : (- bound - ynew->x[j]) / tmp2;
                tmp = MIN(tmp, tmp3);
            }
            alphainf = MIN(alphainf, tmp);
            drate = MIN(1.0, alphainf * dratemax / step); tmp = dsdpSolver->Ry;
            for (inCone = FALSE; !inCone;) {
                dsdpSolver->Ry = tmp * (1 - step * drate);
                vec_copy(ynow, ynew);
                vec_zaxpby(dycorr, -1.0, d2, drate, d4);
                vec_axpy(step, dycorr, ynew);
                getPhaseAS(dsdpSolver, ynew, dsdpSolver->tau);
                getPhaseALps(dsdpSolver, ynew, dsdpSolver->tau);
                dsdpInCone(dsdpSolver, &inCone);
                if (inCone) {
                    dsdpLpInCone(dsdpSolver, &inCone);
                }
                if (inCone) { break; }
                else { drate *= 0.8; }
            }
            
            // printf("Reduce dual infeasibility %e \n", drate * step);
            if (drate * step < 0.1) {
                dratemax *= 0.9;
            }

            if (drate * step < 5e-04) {
                dratemax *= 0.0;
                ncorrector = MIN(ncorrector, i + 1);
            }
            
            if (drate * step > 0.3) {
                dsdpSolver->mu *= 0.95;
                dratemax = MIN(drate * 2.0, 0.8);
            }
            
            if (drate * step > 0.8) {
                dsdpSolver->mu *= 0.9;
                dratemax = MIN(drate * 2.0, 0.8);
            }
            
             dsdpSolver->mu *= shrink;
        }
        vec_copy(ynew, ynow);
        vec_lslack(ynow, sl, -bound);
        vec_uslack(ynow, su, bound);
    }

    cg_finish(cg);
    dsdpSolver->iterProgress[ITER_CORRECTOR] = TRUE;
    return retcode;
}

extern DSDP_INT dualCorrectorStep( HSDSolver *dsdpSolver ) {
    // Dual corrector step in Phase B. When into this routine,
    // S is not inverted and M is ready
    DSDP_INT retcode = DSDP_RETCODE_OK;
    retcode = checkIterProgress(dsdpSolver, ITER_CORRECTOR);
    
    // Prepare CG solver
    adpcg *cg = dsdpSolver->cg;
    DSDP_INT nall = dsdpSolver->nall;
    vec *d2 = dsdpSolver->M->asinv; // asinv is directed solved using CG
    vec *d1 = dsdpSolver->d1, *b = dsdpSolver->dObj, *ynew = dsdpSolver->d4; // d4 is reused for storing ynew
    vec *dycorr = dsdpSolver->b1, *sl = dsdpSolver->sl, *su = dsdpSolver->su, *ynow = dsdpSolver->y;
    double shrink = nall / (nall + sqrt(nall)), bTd1, bTd2, bTdycorr, aval;
    double oldbarrier, newbarrier, step, tmp, rhon;
    
    DSDPGetDblParam(dsdpSolver, DBL_PARAM_RHON, &rhon);
    DSDP_INT ncorrector = adjCorrectorStep(dsdpSolver), inCone = FALSE;
    vec_dot(b, d1, &bTd1);
    
    for (DSDP_INT i = 0, j; i < ncorrector; ++i) {
        
        if (dsdpSolver->mu < 1e-05) { break; }
        
        DSDPConic( COPS_GET_SCHURVEC )(dsdpSolver, FALSE);
        cg_solve(cg, d2, NULL); vec_dot(b, d2, &bTd2);
        if (bTd1 > 0 && bTd2 > 0) {
            dsdpSolver->mu = bTd1 / bTd2; // MIN(dsdpSolver->mu, bTd1 / bTd2);
        }
        dsdpSolver->mu *= shrink;
        
        // dycorr = dy1 / mu - dy2;
        vec_zaxpby(dycorr, 1 / dsdpSolver->mu, d1, -1.0, d2);
        vec_dot(b, dycorr, &bTdycorr);
        oldbarrier = getCurrentLogBarrier(dsdpSolver, ynow, NULL);
        getPhaseBdS(dsdpSolver, 1.0, dycorr, 0.0);
        getPhaseBLpds(dsdpSolver, 1.0, dycorr, 0.0);
        step = DSDP_INFINITY;
        
        tmp = getMaxSDPstep(dsdpSolver, DUALVAR);
        step = MIN(step, tmp);
        tmp = getMaxLpstep(dsdpSolver, DUALVAR);
        step = MIN(step, tmp);
        
        tmp = vec_step(sl, dycorr,  1.0); step = MIN(tmp, step);
        tmp = vec_step(su, dycorr, -1.0); step = MIN(tmp, step);
        
        if (step >= 1000.0) {
            vec_zaxpby(ynew, 1.0, dycorr, 1.0, ynow);
            vec_copy(ynew, ynow);
            getPhaseBS(dsdpSolver, ynow);
            getPhaseBLps(dsdpSolver, ynow);
            dsdpInCone(dsdpSolver, &inCone);
            if (!inCone) {
                fatal_error_msg(etype);
            }
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
                    if (step <= 1e-04) {
                        break;
                    }
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
        
        if (step <= 1e-04) {
            // infeasible or not sufficient descent
            if (inCone) {
                vec_copy(ynew, ynow);
            }
            getBslack(dsdpSolver, ynow, DUALVAR);
            getPhaseBS(dsdpSolver, ynow);
            dsdpInCone(dsdpSolver, &inCone);
            assert( inCone );
            getPhaseBLps(dsdpSolver, ynow);
            break;
        } else {
            // Feasible and sufficient descent
            vec_copy(ynew, ynow);
            getBslack(dsdpSolver, ynow, DUALVAR);
            getPhaseBLps(dsdpSolver, ynow);
            // getPhaseBS(dsdpSolver, ynow);
        }
    }
    
    cg_finish(cg);
    dsdpSolver->iterProgress[ITER_CORRECTOR] = TRUE;
    return retcode;
}
