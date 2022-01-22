#include "stepheur.h"
#include "dsdputils.h"
/*
   Implement the strategies to choose the proper stepsize and update the
   dual iteration variables
*/

static char etype[] = "Step size computation";

static DSDP_INT getKappaTauStep( HSDSolver *dsdpSolver, double *kappatauStep ) {
    // Compute the maximum step size to take at kappa and tau
    DSDP_INT retcode = DSDP_RETCODE_OK;
    
    double tau    = dsdpSolver->tau;
    double kappa  = dsdpSolver->kappa;
    double dtau   = dsdpSolver->dtau;
    double dkappa = dsdpSolver->dkappa;
    
    double step = MIN((dtau / tau), (dkappa / kappa));
    
    if (step < 0.0) {
        *kappatauStep = fabs(dsdpSolver->param->Aalpha / step);
    } else {
        *kappatauStep = DSDP_INFINITY;
    }
    
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
    vec *y  = dsdpSolver->y;
    vec *dy = dsdpSolver->dy;
    retcode = vec_axpy(step, dy, y);
    
    return retcode;
}

static DSDP_INT takelpsStep( HSDSolver *dsdpSolver ) {
    // Take step in LP s
    DSDP_INT retcode = DSDP_RETCODE_OK;
    
    double step = dsdpSolver->alpha;
    vec *s  = dsdpSolver->s;
    vec *ds = dsdpSolver->ds;
    retcode = vec_axpy(step, ds, s);
    
    return retcode;
}

static DSDP_INT takeSDPSStep( HSDSolver *dsdpSolver ) {
    // Take step in SDP S
    DSDP_INT retcode = DSDP_RETCODE_OK;
    retcode = getPhaseAS(dsdpSolver, dsdpSolver->y->x, dsdpSolver->tau);
    return retcode;
}

static DSDP_INT getLPsStep( HSDSolver *dsdpSolver, double *sStep ) {
    // Compute the maixmum step size to take at s for LP
    DSDP_INT retcode = DSDP_RETCODE_OK;
    
    vec *s  = dsdpSolver->s;
    vec *ds = dsdpSolver->ds;
    
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
    retcode = dsdpGetAlpha(dsdpSolver->S[k], dsdpSolver->dS[k],
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
        step = MIN(step, dsdpSolver->param->Aalpha * res);
    }
    
    *SStep = step;
    
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
    
    double stepkappatau = 0.0;
    double steplps = 100.0;
    double sdpS = 0.0;
    
    retcode = getKappaTauStep(dsdpSolver, &stepkappatau);
    // retcode = getLPsStep(dsdpSolver, &steplps);
    retcode = getSDPSStep(dsdpSolver, &sdpS); checkCode;
    sdpS = MIN(sdpS, steplps);
    dsdpSolver->alpha = MIN(sdpS, stepkappatau);
    dsdpSolver->alpha = MIN(dsdpSolver->param->Aalpha * dsdpSolver->alpha, 1.0);
    
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
    // The backward newton step is stored in b2 and Scker, dy1 is in d1; dy is in b1
    
    DSDP_INT retcode = DSDP_RETCODE_OK;
    
    double alpha = DSDP_INFINITY;
    double tmp = 1.0;
    
    if (dsdpSolver->eventMonitor[EVENT_PFEAS_FOUND]) {
        retcode = getPhaseBdS(dsdpSolver, 1.0, dsdpSolver->d1->x, 0.0);
        
        for (DSDP_INT i = 0; i < dsdpSolver->nBlock; ++i) {
            retcode = dsdpGetAlpha(dsdpSolver->Scker[i], dsdpSolver->dS[i],
                                   dsdpSolver->spaux[i], &tmp);
            checkCode;
            alpha = MIN(alpha, tmp);
        }
        *newmu = dsdpSolver->mu / (1 + 0.95 * alpha);
        
    } else {
        
        // dS = dsdpgetATy(A, dy);
        retcode = getPhaseBdS(dsdpSolver, -1.0, dsdpSolver->b1->x, 0.0);
        
        // alphap = dsdpgetalpha(S, dS);
        retcode = getSDPSStep(dsdpSolver, &tmp);
        assert( alpha != DSDP_INFINITY );
        
        // Shat = S + 0.95 * alphap * dS;
        for (DSDP_INT i = 0; i < dsdpSolver->nBlock; ++i) {
            retcode = spsMataXpbY(0.95 * tmp, dsdpSolver->dS[i],
                                  1.0, dsdpSolver->S[i], dsdpSolver->symS[i]);
        }
        
        
        
        
    }
    
    
    
    
    return retcode;
}
