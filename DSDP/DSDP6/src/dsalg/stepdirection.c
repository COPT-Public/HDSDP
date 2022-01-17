#include "stepdirection.h"
#include "dsdputils.h"
/* Recover the steps to take after solving the Schur system */

static char etype[] = "Stepsize recovery";

static DSDP_INT assembleArrs( HSDSolver *dsdpSolver ) {
    // Assemble auxiliary arrays
    /*
     b1 = b - mu * u;
     b2 = d2 * tau / mu - d3 + d4;
     d11 = d2 / mu;
     d1  = d11 + d12;
    */
    DSDP_INT retcode = DSDP_RETCODE_OK;
    vec_zaxpby(dsdpSolver->b1, 1.0, dsdpSolver->dObj,
               dsdpSolver->mu, dsdpSolver->u);
    vec_zaxpby(dsdpSolver->b2, dsdpSolver->tau / dsdpSolver->mu,
               dsdpSolver->d2, -1.0, dsdpSolver->d3);
    vec_axpy(1.0, dsdpSolver->d4, dsdpSolver->b2);
    vec_copy(dsdpSolver->d2, dsdpSolver->d1);
    vec_rscale(dsdpSolver->d1, dsdpSolver->mu);
    vec_axpy(1.0, dsdpSolver->d12, dsdpSolver->d1);
    
    return retcode;
}

static DSDP_INT getdTau( HSDSolver *dsdpSolver ) {
    
    DSDP_INT retcode = DSDP_RETCODE_OK;
    
    vec *b1 = dsdpSolver->b1;
    vec *b2 = dsdpSolver->b2;
    vec *d1 = dsdpSolver->d1;
    
    double mu = dsdpSolver->mu;
    double tau = dsdpSolver->tau;
    double kappa = dsdpSolver->kappa;
    double taudenom = 0.0;
    double csinv = dsdpSolver->csinv;
    double csinvcsinv = dsdpSolver->csinvcsinv;
    double csinvrysinv = dsdpSolver->csinvrysinv;
    
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

static DSDP_INT getdS( HSDSolver *dsdpSolver ) {
    
    DSDP_INT retcode = DSDP_RETCODE_OK;
    double *dy = dsdpSolver->dy->x;
    double dtau = dsdpSolver->dtau;
    retcode = getPhaseAdS(dsdpSolver, dy, dtau); checkCode;
    
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
    
    retcode = checkIterProgress(dsdpSolver, ITER_RECOVER_LP_DIR);
    assert( !dsdpSolver->iterProgress[ITER_RECOVER_LP_DIR] );
    
    if (dsdpSolver->iterProgress[ITER_RECOVER_LP_DIR]) {
        error(etype, "LP directions have been set up. \n");
    }
    
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
    
    dsdpSolver->iterProgress[ITER_RECOVER_LP_DIR] = TRUE;
    
    return retcode;
}

static DSDP_INT getSDPDirs( HSDSolver *dsdpSolver ) {
    
    DSDP_INT retcode = DSDP_RETCODE_OK;
    retcode = checkIterProgress(dsdpSolver, ITER_RECOVER_SDP_DIR);
    assert( !dsdpSolver->iterProgress[ITER_RECOVER_SDP_DIR] );
    
    if (dsdpSolver->iterProgress[ITER_RECOVER_SDP_DIR]) {
        error(etype, "SDP directions have been set up. \n");
    }
    
    retcode = getdTau(dsdpSolver); checkCode;
    retcode = getdS(dsdpSolver); checkCode;
    retcode = getdKappa(dsdpSolver); checkCode;
    
    dsdpSolver->iterProgress[ITER_RECOVER_SDP_DIR] = TRUE;
    
    return retcode;
}

extern DSDP_INT getStepDirs( HSDSolver *dsdpSolver ) {
    
    DSDP_INT retcode = DSDP_RETCODE_OK;
    
    retcode = assembleArrs(dsdpSolver);
    retcode = getSDPDirs(dsdpSolver); checkCode;
    // retcode = getdsLP(dsdpSolver); checkCode;
    
    return retcode;
}
