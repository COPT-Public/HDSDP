#include "stepdirection.h"
#include "dsdputils.h"
/* Recover the steps to take after solving the Schur system */

static char etype[] = "Stepsize recovery";

static DSDP_INT assemblePhaseAArrs( HSDSolver *dsdpSolver ) {
    // Assemble auxiliary arrays
    /*
     b1 = b * pweight - mu * u;
     b2 = d2 * pweight * tau / mu - d3 + drate * d4;
     d11 = d2 / mu;
     d1  = d11 * pweight + d12;
    */
    DSDP_INT retcode = DSDP_RETCODE_OK;
    double pweight; DSDPGetDblParam(dsdpSolver, DBL_PARAM_OBJ_WEIGHT, &pweight);
    
    // b2 = d2 * pweight * tau / mu - d3
    vec_zaxpby(dsdpSolver->b2, pweight * dsdpSolver->tau / dsdpSolver->mu,
               dsdpSolver->d2, -1.0, dsdpSolver->d3);
    // b2 = b2 + d4
    vec_axpy(dsdpSolver->drate, dsdpSolver->d4, dsdpSolver->b2);
    // d1 = d2
    vec_copy(dsdpSolver->d2, dsdpSolver->d1);
    // d1 = d1 / mu * pweight
    vec_scale(dsdpSolver->d1, pweight / dsdpSolver->mu);
    
    if (dsdpSolver->eventMonitor[EVENT_HSD_UPDATE]) {
        // b1 = b * pweight - mu * u;
        vec_zaxpby(dsdpSolver->b1, pweight, dsdpSolver->dObj,
                   - dsdpSolver->mu, dsdpSolver->u);
        // d1 = d1 + d2
        vec_axpy(1.0, dsdpSolver->d12, dsdpSolver->d1);
    }
    
    return retcode;
}

static void getdTau( HSDSolver *dsdpSolver ) {
    
    if (!dsdpSolver->eventMonitor[EVENT_HSD_UPDATE]) {
        dsdpSolver->dtau = 0.0; return;
    }
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
}

static void getdy( HSDSolver *dsdpSolver ) {
    // dy = d1 * dtau + d2
    vec *dy = dsdpSolver->dy, *d1 = dsdpSolver->d1, *b2 = dsdpSolver->b2;
    vec_zaxpby(dy, dsdpSolver->dtau, d1, 1.0, b2);
}

static void getdyB( HSDSolver *dsdpSolver ) {
    // dy = dy1 / muprimal - dy2
    vec_zaxpby(dsdpSolver->dy, 1 / dsdpSolver->mu, dsdpSolver->d1, -1.0, dsdpSolver->d2);
    vec_zaxpby(dsdpSolver->d4, 1 / dsdpSolver->mu, dsdpSolver->dObj, -1.0, dsdpSolver->asinv);
    vec_dot(dsdpSolver->d4, dsdpSolver->dy, &dsdpSolver->Pnrm);
    dsdpSolver->Pnrm = sqrt(MAX(dsdpSolver->Pnrm, 0.0));
    if (dsdpSolver->Pnrm < 0.1) {
        dsdpSolver->mu *= 0.1;
        vec_zaxpby(dsdpSolver->dy, 1 / dsdpSolver->mu, dsdpSolver->d1, -1.0, dsdpSolver->d2);
        vec_zaxpby(dsdpSolver->d4, 1 / dsdpSolver->mu, dsdpSolver->dObj, -1.0, dsdpSolver->asinv);
        vec_dot(dsdpSolver->d4, dsdpSolver->dy, &dsdpSolver->Pnrm);
        dsdpSolver->Pnrm = sqrt(MAX(dsdpSolver->Pnrm, 0.0));
    }
}

static void getdS( HSDSolver *dsdpSolver ) {
    double *dy = dsdpSolver->dy->x, dtau = dsdpSolver->dtau;
    getPhaseAdS(dsdpSolver, dsdpSolver->drate, dy, dtau);
}

static void getdSB( HSDSolver *dsdpSolver ) {
    getPhaseBdS(dsdpSolver, 1.0, dsdpSolver->dy->x, 0.0);
}

static void getdKappa( HSDSolver *dsdpSolver ) {
#ifdef KAPPATAU
    dsdpSolver->dkappa = - dsdpSolver->kappa
                         + dsdpSolver->mu / dsdpSolver->tau
                         - (dsdpSolver->kappa * dsdpSolver->dtau) / dsdpSolver->tau;
#else
    dsdpSolver->dkappa = 0.0;
#endif
}

static void getSDPDirs( HSDSolver *dsdpSolver ) {
    if (dsdpSolver->eventMonitor[EVENT_IN_PHASE_A]) {
        // The four directions must be set up in order
        getdTau(dsdpSolver); getdy(dsdpSolver);
        getdS(dsdpSolver); getdKappa(dsdpSolver);
    } else {
        getdyB(dsdpSolver); getdSB(dsdpSolver);
    }
}

extern DSDP_INT getStepDirs( HSDSolver *dsdpSolver ) {
    
    DSDP_INT retcode = DSDP_RETCODE_OK;
    retcode = checkIterProgress(dsdpSolver, ITER_STEP_DIRECTION);
    if (dsdpSolver->iterProgress[ITER_STEP_DIRECTION]) {
        error(etype, "SDP directions have been set up. \n");
    }
    
    if (dsdpSolver->eventMonitor[EVENT_IN_PHASE_A]) {
        retcode = assemblePhaseAArrs(dsdpSolver);
    }
    getSDPDirs(dsdpSolver);
    dsdpSolver->iterProgress[ITER_STEP_DIRECTION] = TRUE;
    return retcode;
}
