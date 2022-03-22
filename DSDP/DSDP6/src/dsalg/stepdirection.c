#include "stepdirection.h"
#include "dsdputils.h"
/* Recover the steps to take after solving the Schur system */

static char etype[] = "Stepsize recovery";

static DSDP_INT assemblePhaseAArrs( HSDSolver *dsdpSolver ) {
    // Assemble auxiliary arrays
    /*
     b1 = b - mu * u;
     b2 = d2 * tau / mu - d3 + d4;
     d11 = d2 / mu;
     d1  = d11 + d12;
    */
    DSDP_INT retcode = DSDP_RETCODE_OK;
    
    // b1 = b - mu * u
    vec_zaxpby(dsdpSolver->b1, 1.0, dsdpSolver->dObj,
               - dsdpSolver->mu, dsdpSolver->u);
    // b2 = tau / mu * d2 - d3
    vec_zaxpby(dsdpSolver->b2, dsdpSolver->tau / dsdpSolver->mu,
               dsdpSolver->d2, -1.0, dsdpSolver->d3);
    // b2 = b2 + d4
    vec_axpy(1.0, dsdpSolver->d4, dsdpSolver->b2);
    // d1 = d2
    vec_copy(dsdpSolver->d2, dsdpSolver->d1);
    // d1 = d1 / mu
    vec_rscale(dsdpSolver->d1, dsdpSolver->mu);
    // d1 = d1 + d2
    vec_axpy(1.0, dsdpSolver->d12, dsdpSolver->d1);
    
    return retcode;
}

static DSDP_INT getdTau( HSDSolver *dsdpSolver ) {
    
    DSDP_INT retcode = DSDP_RETCODE_OK;
    
    vec *b1 = dsdpSolver->b1, *b2 = dsdpSolver->b2, *d1 = dsdpSolver->d1;
    double mu = dsdpSolver->mu, tau = dsdpSolver->tau;
    double kappa = dsdpSolver->kappa, taudenom = 0.0, csinv = dsdpSolver->csinv;
    double csinvcsinv = dsdpSolver->csinvcsinv, csinvrysinv = dsdpSolver->csinvrysinv;
    
    vec_dot(b1, d1, &taudenom);
    taudenom += mu * csinvcsinv + kappa / tau;
    if (fabs(taudenom) < 1e-15) {
        dsdpSolver->dtau = 0.0;
    } else {
        dsdpSolver->dtau = - dsdpSolver->dObjVal;
        vec_dot(b1, b2, &csinvcsinv);
        dsdpSolver->dtau += mu * (1 / tau + csinv - csinvrysinv) - csinvcsinv;
        dsdpSolver->dtau = dsdpSolver->dtau / taudenom;
    }
   
    return retcode;
}

static DSDP_INT getdy( HSDSolver *dsdpSolver ) {
    // dy = d1 * dtau + d2
    DSDP_INT retcode = DSDP_RETCODE_OK;
    
    vec *dy = dsdpSolver->dy;
    vec *d1 = dsdpSolver->d1;
    vec *b2 = dsdpSolver->b2;
    
    retcode = vec_zaxpby(dy, dsdpSolver->dtau, d1, 1.0, b2);
    
    return retcode;
}

static DSDP_INT getdyB( HSDSolver *dsdpSolver ) {
    // dy = dy1 / muprimal - dy2
    DSDP_INT retcode = DSDP_RETCODE_OK;
    
    vec_zaxpby(dsdpSolver->dy, 1 / dsdpSolver->mu,
               dsdpSolver->d1, -1.0, dsdpSolver->d2);
    
    return retcode;
}

static DSDP_INT getdS( HSDSolver *dsdpSolver ) {
    
    DSDP_INT retcode = DSDP_RETCODE_OK;
    double *dy = dsdpSolver->dy->x;
    double dtau = dsdpSolver->dtau;
    retcode = getPhaseAdS(dsdpSolver, dy, dtau); checkCode;
    
    return retcode;
}

static DSDP_INT getdSB( HSDSolver *dsdpSolver ) {
    
    DSDP_INT retcode = DSDP_RETCODE_OK;
    retcode = getPhaseBdS(dsdpSolver, 1.0, dsdpSolver->dy->x, 0.0);
    return retcode;
}

static DSDP_INT getdKappa( HSDSolver *dsdpSolver ) {
    
    DSDP_INT retcode = DSDP_RETCODE_OK;
    
    dsdpSolver->dkappa = - dsdpSolver->kappa
                         + dsdpSolver->mu / dsdpSolver->tau
                         - (dsdpSolver->kappa * dsdpSolver->dtau) / dsdpSolver->tau;
    return retcode;
}

static DSDP_INT getdsLP( HSDSolver *dsdpSolver ) {
    // Compute the dual direction for LP ry - A' * dy + c * dtau
    DSDP_INT retcode = DSDP_RETCODE_OK;

    vec *ds = dsdpSolver->ds;
    double *dydata = dsdpSolver->dy->x;
    double *dsdata = ds->x;
    DSDP_INT *Ap = dsdpSolver->lpData->lpdata->p;
    DSDP_INT *Ai = dsdpSolver->lpData->lpdata->i;
    double   *Ax = dsdpSolver->lpData->lpdata->x;
    DSDP_INT n = ds->dim;
    
    retcode = vec_copy(dsdpSolver->ry, ds);
    retcode = vec_axpy(dsdpSolver->dtau, dsdpSolver->lpObj, ds);
    
    double tmp = 0.0;
    
    for (DSDP_INT i = 0; i < n; ++i) {
        tmp = 0.0;
        for (DSDP_INT j = Ap[i]; j < Ap[i + 1]; ++i) {
            tmp += dydata[Ai[j]] * Ax[j];
        }
        dsdata[i] -= tmp;
    }

    return retcode;
}

static DSDP_INT getSDPDirs( HSDSolver *dsdpSolver ) {
    
    DSDP_INT retcode = DSDP_RETCODE_OK;
    
    if (dsdpSolver->eventMonitor[EVENT_IN_PHASE_A]) {
        retcode = getdTau(dsdpSolver); checkCode;
        retcode = getdy(dsdpSolver); checkCode;
        retcode = getdS(dsdpSolver); checkCode;
        retcode = getdKappa(dsdpSolver); checkCode;
    } else {
        retcode = getdyB(dsdpSolver); checkCode;
        retcode = getdSB(dsdpSolver); checkCode;
    }
    
    return retcode;
}

extern DSDP_INT getStepDirs( HSDSolver *dsdpSolver ) {
    
    DSDP_INT retcode = DSDP_RETCODE_OK;
    
    retcode = checkIterProgress(dsdpSolver, ITER_STEP_DIRECTION);
    assert( !dsdpSolver->iterProgress[ITER_STEP_DIRECTION] );
    
    if (dsdpSolver->iterProgress[ITER_STEP_DIRECTION]) {
        error(etype, "SDP directions have been set up. \n");
    }
    
    if (dsdpSolver->eventMonitor[EVENT_IN_PHASE_A]) {
        retcode = assemblePhaseAArrs(dsdpSolver);
    }
    retcode = getSDPDirs(dsdpSolver); checkCode;
    // retcode = getdsLP(dsdpSolver); checkCode;
    dsdpSolver->iterProgress[ITER_STEP_DIRECTION] = TRUE;
    
    return retcode;
}
