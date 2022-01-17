#include "stepdirection.h"
/* Recover the steps to take after solving the Schur system */

static char etype[] = "Stepsize recovery";

static DSDP_INT getdTau( HSDSolver *dsdpSolver, double rM, vec *b2 ) {
    
    DSDP_INT retcode = DSDP_RETCODE_OK;
    
    vec *d1 = dsdpSolver->d1;
    vec *d2 = dsdpSolver->d2;
    
    DSDP_INT m = d1->dim;
    DSDP_INT one = 1;
    
    double mu = dsdpSolver->mu;
    double T = 0.0;
    double b2dotd1 = 0.0;
    double b2dotd2 = 0.0;
    
    b2dotd1 = ddot(&m, b2->x, &one, d1->x, &one);
    T = b2dotd1 + mu * dsdpSolver->csinvcsinv + dsdpSolver->kappa / dsdpSolver->tau;
    if (fabs(T) < 1e-12) {
        dsdpSolver->dtau = 0.0;
    } else {
        b2dotd2 = ddot(&m, b2->x, &one, d2->x, &one);
        dsdpSolver->dtau = (rM - b2dotd2) / T;
    }
    
    return retcode;
}

static DSDP_INT getdy( HSDSolver *dsdpSolver ) {
    // dy = d1 * dtau + d2
    DSDP_INT retcode = DSDP_RETCODE_OK;
    
    vec *dy = dsdpSolver->dy;
    vec *d1 = dsdpSolver->d1;
    vec *d2 = dsdpSolver->d2;
    
    retcode = vec_zaxpby(dy, dsdpSolver->dtau, d1, 1.0, d2);
    
    return retcode;
}

static DSDP_INT getBlockdS( HSDSolver *dsdpSolver, DSDP_INT blockid ) {
    // Compute dS for some block (dS = Ry - ATdy + C * dtau)
    // dS is computed exactly in the same way as Ry
    DSDP_INT retcode = DSDP_RETCODE_OK;
    
    return retcode;
}

static DSDP_INT getdS( HSDSolver *dsdpSolver ) {
    
    DSDP_INT retcode = DSDP_RETCODE_OK;
    DSDP_INT nblock = dsdpSolver->nBlock;
    
    for (DSDP_INT i = 0; i < nblock; ++i) {
        retcode = getBlockdS(dsdpSolver, i); checkCode;
    }
    
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

static DSDP_INT getSDPDirs( HSDSolver *dsdpSolver, double rM, vec *b2 ) {
    
    DSDP_INT retcode = DSDP_RETCODE_OK;
    retcode = checkIterProgress(dsdpSolver, ITER_RECOVER_SDP_DIR);
    assert( !dsdpSolver->iterProgress[ITER_RECOVER_SDP_DIR] );
    
    if (dsdpSolver->iterProgress[ITER_RECOVER_SDP_DIR]) {
        error(etype, "SDP directions have been set up. \n");
    }
    
    retcode = getdTau(dsdpSolver, rM, b2); checkCode;
    retcode = getdy(dsdpSolver); checkCode;
    retcode = getdS(dsdpSolver); checkCode;
    retcode = getdKappa(dsdpSolver); checkCode;
    
    dsdpSolver->iterProgress[ITER_RECOVER_SDP_DIR] = TRUE;
    
    return retcode;
}

extern DSDP_INT getStepDirs( HSDSolver *dsdpSolver, double rM, vec *b2 ) {
    
    DSDP_INT retcode = DSDP_RETCODE_OK;
    
    retcode = getSDPDirs(dsdpSolver, rM, b2); checkCode;
    retcode = getdsLP(dsdpSolver); checkCode;
    
    return retcode;
}
