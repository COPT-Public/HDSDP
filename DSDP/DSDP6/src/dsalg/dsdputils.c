#include "dsdputils.h"
#include "dsdpsolver.h"
#include "sparsemat.h"
#include "schurmat.h"
static char etype[] = "DSDP Conic Utility";

// #define CONIC

/*
  Utility routines that manage the operations between
  different types of matrices
*/

/*
 Conic operations: Scaling operations
 Compute one-norm for SDP, LP and other cones
 */

static double SDPConic( COPS_GET_A_ONE_NORM )
( HSDSolver *dsdpSolver ) {
    
    DSDP_INT k = dsdpSolver->nBlock, i;
    double onenorm = 0.0;
    for (i = 0; i < k; ++i) {
        onenorm += sdpMatGetAOneNorm(dsdpSolver->sdpData[i]);
    }
    return onenorm;
}

static double LPConic( COPS_GET_A_ONE_NORM )
( HSDSolver *dsdpSolver ) {
    
    if (!dsdpSolver->isLPset) { return 0.0; }
    double *Ax = dsdpSolver->lpData->Ax;
    DSDP_INT nnz = dsdpSolver->lpData->nnz, one = 1;
    return dasum(&nnz, Ax, &one);
}

static double BConic( COPS_GET_A_ONE_NORM )
( HSDSolver *dsdpSolver ) {
    
    return 0.0;
}

static double SDPConic( COPS_GET_C_FNORM )
( HSDSolver *dsdpSolver ) {
    
    DSDP_INT k = dsdpSolver->nBlock, i;
    double fnrm = 0.0, tmp = 0.0;
    for (i = 0; i < k; ++i) {
        tmp = sdpMatGetCFnorm(dsdpSolver->sdpData[i]);
        fnrm += tmp * tmp;
    }
    return sqrt(fnrm);
}

static double LPConic( COPS_GET_C_FNORM )
( HSDSolver *dsdpSolver ) {
    
    if (!dsdpSolver->isLPset) { return 0.0; }
    double *Ax = dsdpSolver->lpData->Ax;
    DSDP_INT nnz = dsdpSolver->lpData->nnz, one = 1;
    return dnrm2(&nnz, Ax, &one);
}

static double BConic( COPS_GET_C_FNORM )
( HSDSolver *dsdpSolver ) {
    
    return 0.0;
}

static double SDPConic( COPS_GET_C_ONE_NORM )
( HSDSolver *dsdpSolver ) {
    
    DSDP_INT k = dsdpSolver->nBlock, i;
    double onenorm = 0.0;
    for (i = 0; i < k; ++i) {
        onenorm += sdpMatGetCOneNorm(dsdpSolver->sdpData[i]);
    }
    return onenorm;
}

static double LPConic( COPS_GET_C_ONE_NORM )
( HSDSolver *dsdpSolver ) {
    
    if (!dsdpSolver->isLPset) { return 0.0; }
    return vec_onenorm(dsdpSolver->lpObj);
}

static double BConic( COPS_GET_C_ONE_NORM )
( HSDSolver *dsdpSolver ) {
    
    return 0.0;
}

static void SDPConic( COPS_DO_C_SCALE )
( HSDSolver *dsdpSolver ) {
    
    for (DSDP_INT k = 0; k < dsdpSolver->nBlock; ++k) {
        sdpMatRScaleC(dsdpSolver->sdpData[k], dsdpSolver->cScaler);
    }
    return;
}

static void LPConic( COPS_DO_C_SCALE )
( HSDSolver *dsdpSolver ) {
    
    if (!dsdpSolver->isLPset) { return; }
    vec_rscale(dsdpSolver->lpObj, dsdpSolver->cScaler);
}

static void BConic( COPS_DO_C_SCALE )
( HSDSolver *dsdsSpolver ) {
    
    return;
}

static void SDPConic( COPS_DO_VAR_SCALE )
( HSDSolver *dsdpSolver, double tau ) {
    
    for (DSDP_INT i = 0; i < dsdpSolver->nBlock; ++i) {
        spsMatRscale(dsdpSolver->S[i], tau);
    }
}

static void LPConic( COPS_DO_VAR_SCALE )
( HSDSolver *dsdpSolver, double tau ) {
    
    if (!dsdpSolver->isLPset) { return; }
    vec_rscale(dsdpSolver->s, dsdpSolver->tau);
}

static void BConic( COPS_DO_VAR_SCALE )
( HSDSolver *dsdpSolver, double tau ) {
    
    return;
}

static DSDP_INT SDPConic( COPS_CHECK_INCONE )
( HSDSolver *dsdpSolver, DSDP_INT type ) {
    
    if (type == DUALVAR) {
        for (DSDP_INT i = 0, incone = FALSE; i < dsdpSolver->nBlock; ++i) {
                    spsMatIspd(dsdpSolver->S[i], &incone);
                    if (!incone) { return FALSE; }
        }
    } else {
        for (DSDP_INT i = 0, incone = FALSE; i < dsdpSolver->nBlock; ++i) {
                    spsMatIspd(dsdpSolver->Scker[i], &incone);
                    if (!incone) { return FALSE; }
        }
    }
    return TRUE;
}

static DSDP_INT LPConic( COPS_CHECK_INCONE )
( HSDSolver *dsdpSolver, DSDP_INT type ) {
    
    if (!dsdpSolver->isLPset) { return TRUE; }
    return vec_incone(dsdpSolver->s);
}

static DSDP_INT BConic( COPS_CHECK_INCONE )
( HSDSolver *dsdpSolver, DSDP_INT type ) {
    
    if (type == DUALVAR) {
        return (vec_incone(dsdpSolver->sl) && vec_incone(dsdpSolver->su));
    } else {
        return (vec_incone(dsdpSolver->slcker) && vec_incone(dsdpSolver->sucker));
    }
}

static void SDPConic( COPS_SYMFAC )
( HSDSolver *dsdpSolver ) {
    
    for (DSDP_INT i = 0; i < dsdpSolver->nBlock; ++i) {
        spsMatSymbolic(dsdpSolver->S[i]);
        spsMatSymbolic(dsdpSolver->Scker[i]);
    }
}

static void LPConic( COPS_SYMFAC )
( HSDSolver *dsdpSolver ) {
    
    return;
}

static void BConic ( COPS_SYMFAC )
( HSDSolver *dsdpSolver ) {
    
    return;
}

static void SDPConic ( COPS_GET_SCHUR )
( HSDSolver *dsdpSolver ) {
    
    DSDPSchurSetup(dsdpSolver->M);
    dsdpSolver->iterProgress[ITER_SCHUR] = TRUE;
}

static void LPConic( COPS_GET_SCHUR )
( HSDSolver *dsdpSolver ) {
    // Compute A s^-2 A^T
    if (!dsdpSolver->isLPset) { return; }
    if (dsdpSolver->Msdp->stype != SCHUR_TYPE_DENSE) {
        printf("| Fatal Error in dsdputils.c -> Line 206 -> LPConic( COPS_GET_SCHUR ). Give up. \n");
        printf("---------------------------------------"
               "---------------------------------------"
               "--------------------\n");
        printf("| DSDP Ends by Fatal Error. No solution available. \n");
        printf("---------------------------------------"
               "---------------------------------------"
               "--------------------\n");
        exit(0);
    }
    double *M = dsdpSolver->Msdp->denseM->array;
    double *asinv = dsdpSolver->asinv->x;
    double *asinvrysinv = dsdpSolver->asinvrysinv->x;
    double *csinv = &dsdpSolver->csinv, *rysinv = &dsdpSolver->rysinv;
    
    DSDP_INT *Ap = dsdpSolver->lpData->Ap;
    DSDP_INT *Ai = dsdpSolver->lpData->Ai;
    double *Ax = dsdpSolver->lpData->Ax, *s = dsdpSolver->s->x;
    double Ry = dsdpSolver->Ry, *lpObj = dsdpSolver->lpObj->x;
    DSDP_INT i, j, k, m = dsdpSolver->m, n = dsdpSolver->lpDim;
    double Mij = 0.0, *a = dsdpSolver->x->x, tmp;
    
    if (Ry) {
        for (i = 0; i < m; ++i) {
            memset(a, 0, sizeof(double) * n);
            for (k = Ap[i]; k < Ap[i + 1]; ++k) {
                tmp = 1 / s[Ai[k]]; asinv[i] += Ax[k] * tmp;
                asinvrysinv[i] += Ax[k] * tmp * tmp * Ry;
                a[Ai[k]] = Ax[k] * (tmp * tmp);
            }
            for (j = 0; j <= i; ++j) {
                Mij = 0.0;
                for (k = Ap[j]; k < Ap[j + 1]; ++k) {
                    Mij += Ax[k] * a[Ai[k]];
                }
                fullIdx(M, m, i, j) += Mij;
            }
        }
    } else {
        for (i = 0; i < m; ++i) {
            memset(a, 0, sizeof(double) * n);
            for (k = Ap[i]; k < Ap[i + 1]; ++k) {
                tmp = s[Ai[k]]; asinv[i] += Ax[k] / tmp;
                a[Ai[k]] = Ax[k] / (tmp * tmp);
            }
            for (j = 0; j <= i; ++j) {
                Mij = 0.0;
                for (k = Ap[j]; k < Ap[j + 1]; ++k) {
                    Mij += Ax[k] * a[Ai[k]];
                }
                fullIdx(M, m, i, j) += Mij;
            }
        }
    }
    
    double tmp2 = 0.0; tmp = 0.0;
    for (i = 0; i < n; ++i) {
        tmp += 1 / s[i];
        tmp2 += lpObj[i] / s[i];
    }
    
    *csinv += tmp2;
    *rysinv += Ry * tmp;
    
    return;
}

static void BConic( COPS_GET_SCHUR )
( HSDSolver *dsdpSolver ) {
    
    DSDP_INT m = dsdpSolver->m, i;
    double **Mdiag = dsdpSolver->Msdp->diag, *asinv = dsdpSolver->asinv->x, s = 0.0;
    double *sl = dsdpSolver->sl->x, *su = dsdpSolver->su->x, bound = dsdpSolver->ybound;
    
    if (dsdpSolver->eventMonitor[EVENT_IN_PHASE_A]) {
        double *u = dsdpSolver->u->x;
        double sinvsqr, csinvsuml = 0.0, \
               csinvsumu = 0.0, cscssuml = 0.0, cscssumu = 0.0;
        
        for (i = 0; i < m; ++i) {
            // Upperbound
            s = su[i]; sinvsqr = 1.0 / (s * s);
            *Mdiag[i] += sinvsqr; asinv[i] += 1.0 / s;
            u[i] += bound * sinvsqr; csinvsumu += 1.0 / s;
            cscssumu += sinvsqr;
            // Lowerbound
            s = sl[i]; sinvsqr = 1.0 / (s * s);
            *Mdiag[i] += sinvsqr; asinv[i] -= 1.0 / s;
            u[i] += bound * sinvsqr; csinvsuml += 1.0 / s;
            cscssuml += sinvsqr;
        }
        
        dsdpSolver->csinv += bound * csinvsumu - bound * csinvsuml;
        dsdpSolver->csinvcsinv += bound * bound * (cscssumu + cscssuml);
        
    } else {
        for (i = 0; i < m; ++i) {
            s = su[i]; *Mdiag[i] += 1.0 / (s * s); asinv[i] += 1.0 / s;
            s = sl[i]; *Mdiag[i] += 1.0 / (s * s); asinv[i] -= 1.0 / s;
        }
    }
}

static void SDPConic( COPS_CONSTR_EXPR )
( HSDSolver *dsdpSolver, DSDP_INT type,
  double ycoef, vec *y, double tau, double r ) {
    // r * Ry + tau * C + ycoef * ATy
    spsMat **targets = (type == DUALVAR) ? \
    dsdpSolver->S : dsdpSolver->Scker;
    void (*f)(HSDSolver*, DSDP_INT, DSDP_INT, double) = \
    (type == DUALVAR) ? addMattoS : addMattoChecker;
    double dperturb = dsdpSolver->dperturb;
    if (type == DELTAS) {
        targets = dsdpSolver->dS; f = addMattodS; dperturb = 0.0;
    }
    dperturb += r * dsdpSolver->Ry;
    spsMat *target = NULL;
    for (DSDP_INT i = 0; i < dsdpSolver->nBlock; ++i) {
        target = targets[i]; spsMatReset(target);
        if (ycoef) { sdpMatATy(dsdpSolver->sdpData[i], ycoef, y,
                               tau, target, dsdpSolver->symS[i]); }
        if (dperturb) { spsMatAdddiag(target, dperturb, dsdpSolver->symS[i]); }
    }
    return;
}

static void LPConic( COPS_CONSTR_EXPR )
( HSDSolver *dsdpSolver, DSDP_INT type,
  double ycoef, vec *y, double tau, double r ) {
    // r * Ry + tau * C + ycoef * ATy
    if (!dsdpSolver->isLPset) { return; }
    vec *target = (type == DUALVAR) ? dsdpSolver->s : dsdpSolver->scker;
    double dperturb = dsdpSolver->dperturb;
    if (type == DELTAS) {
        target = dsdpSolver->ds; dperturb = 0.0;
    }
    vec_reset(target); double *slack = target->x, *c = dsdpSolver->lpObj->x;
    DSDP_INT i, j, *Ap = dsdpSolver->lpData->Ap, *Ai = dsdpSolver->lpData->Ai;
    double *Ax = dsdpSolver->lpData->Ax;
    for (i = 0; i < dsdpSolver->m; ++i) {
        for (j = Ap[i]; j < Ap[i + 1]; ++j) {
            slack[Ai[j]] += ycoef * Ax[j] * y->x[i];
        }
    }
    dperturb += r * dsdpSolver->Ry;
    if (tau || dperturb) {
        for (i = 0; i < dsdpSolver->lpDim; ++i) {
            slack[i] += tau * c[i]; slack[i] += dperturb;
        }
    }
    return;
}

static void BConic( COPS_CONSTR_EXPR )
( HSDSolver *dsdpSolver, DSDP_INT type,
  double ycoef, vec *y, double tau, double r ) {
    // r = 0 and S = tau * C + ycoef * ATy
    // -I * y + sl = -l and I * y + su =  u
    double bd = dsdpSolver->ybound;
    if (type == DELTAS) { // Take care of this !!
        double *dy = dsdpSolver->dy->x;
        for (DSDP_INT i = 0; i < dsdpSolver->m; ++i) {
            dy[i] = tau * bd - ycoef * y->x[i];
        }
    }
    double *sl, *su, *ydata = y->x;
    if (type == DUALVAR) {
        sl = dsdpSolver->sl->x;
        su = dsdpSolver->su->x;
    } else {
        sl = dsdpSolver->slcker->x;
        su = dsdpSolver->sucker->x;
    }
    for (DSDP_INT i = 0; i < dsdpSolver->m; ++i) {
        sl[i] = tau * bd - ycoef * ydata[i];
        su[i] = tau * bd + ycoef * ydata[i];
    }
    return;
}

static void SDPConic( COPS_GET_SLACK )
( HSDSolver *dsdpSolver, DSDP_INT type ) {
    
    if (type == DUALVAR) {
        if (dsdpSolver->eventMonitor[EVENT_IN_PHASE_A]) {
            getPhaseAS(dsdpSolver, dsdpSolver->y, dsdpSolver->tau);
        } else {
            getPhaseBS(dsdpSolver, dsdpSolver->y);
        }
    } else {
        if (dsdpSolver->eventMonitor[EVENT_IN_PHASE_A]) {
            getPhaseACheckerS(dsdpSolver, dsdpSolver->y, dsdpSolver->tau);
        } else {
            getPhaseBCheckerS(dsdpSolver, dsdpSolver->y);
        }
    }
}

static void LPConic( COPS_GET_SLACK )
( HSDSolver *dsdpSolver, DSDP_INT type ) {
    // c - ATy + Ry
    if (!dsdpSolver->isLPset) { return; }
    LPConic( COPS_CONSTR_EXPR )(dsdpSolver, type, -1.0,
                                dsdpSolver->y, dsdpSolver->tau, -1.0);
    return;
}

static void BConic( COPS_GET_SLACK )
( HSDSolver *dsdpSolver, DSDP_INT type ) {
    
    if (type == DUALVAR) {
        vec_lslack(dsdpSolver->y, dsdpSolver->sl,
                   -dsdpSolver->ybound * dsdpSolver->tau);
        vec_uslack(dsdpSolver->y, dsdpSolver->su,
                   dsdpSolver->ybound * dsdpSolver->tau);
    } else {
        vec_lslack(dsdpSolver->y, dsdpSolver->slcker,
                   -dsdpSolver->ybound * dsdpSolver->tau);
        vec_uslack(dsdpSolver->y, dsdpSolver->sucker,
                   dsdpSolver->ybound * dsdpSolver->tau);
    }
}

static double SDPConic( COPS_GET_MAXSTEP )
( HSDSolver *dsdpSolver, DSDP_INT type ) {
    
    double step = DSDP_INFINITY, tmp;
    spsMat **targets = (type == DUALVAR) ? dsdpSolver->S : dsdpSolver->Scker;
    for (DSDP_INT i = 0; i < dsdpSolver->nBlock; ++i) {
        dsdpGetAlpha(dsdpSolver->lczSolver[i], targets[i], dsdpSolver->dS[i], &tmp);
        step = MIN(step, tmp);
    }
    return step;
}

static double LPConic( COPS_GET_MAXSTEP )
( HSDSolver *dsdpSolver, DSDP_INT type ) {
    
    if (!dsdpSolver->isLPset) { return DSDP_INFINITY; };
    vec *target = (type == DUALVAR) ? dsdpSolver->s : dsdpSolver->scker;
    return vec_step(target, dsdpSolver->ds, 1.0);
}

static double BConic( COPS_GET_MAXSTEP )
( HSDSolver *dsdpSolver, DSDP_INT type ) {
    // dsl = dy, dsu = -dy;
    double step = DSDP_INFINITY, tmp;
    vec *target = (type == DUALVAR) ? dsdpSolver->sl : dsdpSolver->slcker;
    tmp = vec_step(target, dsdpSolver->dy, 1.0); step = MIN(step, tmp);
    target = (type == DUALVAR) ? dsdpSolver->su : dsdpSolver->sucker;
    tmp = vec_step(target, dsdpSolver->dy, -1.0); step = MIN(step, tmp);
    return step;
}

static void SDPConic( COPS_STEPDIRECTION )
( HSDSolver *dsdpSolver ) {
    
    if (dsdpSolver->eventMonitor[EVENT_IN_PHASE_A]) {
        getPhaseAdS(dsdpSolver, dsdpSolver->drate,
                    dsdpSolver->dy, dsdpSolver->dtau);
    } else {
        getPhaseBdS(dsdpSolver, 1.0, dsdpSolver->dy, 0.0);
    }
}

static void LPConic( COPS_STEPDIRECTION )
( HSDSolver *dsdpSolver ) {
    
    if (!dsdpSolver->isLPset) { return; }
    if (dsdpSolver->eventMonitor[EVENT_IN_PHASE_A]) {
        getPhaseALpds(dsdpSolver, dsdpSolver->drate,
                      dsdpSolver->dy, dsdpSolver->dtau);
    } else {
        getPhaseALpds(dsdpSolver, 0.0, dsdpSolver->dy, 0.0);
    }
}

static void BConic( COPS_STEPDIRECTION )
( HSDSolver *dsdpSolver ) {
    
    return;
}

static void SDPConic( COPS_GET_SCHURVEC )
( HSDSolver *dsdpSolver, DSDP_INT dInfeas ) {
    
    if (dInfeas) {
        arysinvSetup(dsdpSolver->M);
    } else {
        asinvSetup(dsdpSolver->M);
    }
}

static void LPConic( COPS_GET_SCHURVEC )
( HSDSolver *dsdpSolver, DSDP_INT dInfeas ) {
    
    if (!dsdpSolver->isLPset) { return; }
    
    double *a = dsdpSolver->x->x, tmp;
    double *s = dsdpSolver->s->x;
    DSDP_INT i, k, m = dsdpSolver->m, n = dsdpSolver->lpDim;
    DSDP_INT *Ap = dsdpSolver->lpData->Ap;
    DSDP_INT *Ai = dsdpSolver->lpData->Ai;
    double *Ax = dsdpSolver->lpData->Ax;
    double *asinv = dsdpSolver->asinv->x;
    
    if (dInfeas) {
        double Ry = dsdpSolver->Ry, *asinvrysinv = dsdpSolver->asinvrysinv->x;
        for (i = 0; i < m; ++i) {
            memset(a, 0, sizeof(double) * n);
            for (k = Ap[i]; k < Ap[i + 1]; ++k) {
                tmp = 1 / s[Ai[k]]; asinv[i] += Ax[k] * tmp;
                asinvrysinv[i] += Ax[k] * tmp * tmp * Ry;
            }
        }
    } else {
        for (i = 0; i < m; ++i) {
            memset(a, 0, sizeof(double) * n);
            for (k = Ap[i]; k < Ap[i + 1]; ++k) {
                asinv[i] += Ax[k] / s[Ai[k]];
            }
        }
    }
    return;
}

static void BConic( COPS_GET_SCHURVEC )
( HSDSolver *dsdpSolver, DSDP_INT dInfeas ) {
    
    double *asinv = dsdpSolver->asinv->x, \
           *sl = dsdpSolver->sl->x, *su = dsdpSolver->su->x;
    DSDP_INT j;
    vec_lslack(dsdpSolver->y, dsdpSolver->sl, -dsdpSolver->ybound);
    vec_uslack(dsdpSolver->y, dsdpSolver->su,  dsdpSolver->ybound);
    for (j = 0; j < dsdpSolver->m; ++j) {
        asinv[j] += 1.0 / su[j]; asinv[j] -= 1.0 / sl[j];
    }
}

static double SDPConic( COPS_GET_LOGDET )
( HSDSolver *dsdpSolver, vec *y, DSDP_INT *inCone ) {
    
    double logdet = 0.0;
    if (inCone) { if (!(*inCone)) return 0.0; }
    if (dsdpSolver->eventMonitor[EVENT_IN_PHASE_A]) {
        getPhaseAS(dsdpSolver, y, dsdpSolver->tau);
    } else {
        getPhaseBS(dsdpSolver, y);
    }
    if (inCone) {
        *inCone = SDPConic( COPS_CHECK_INCONE )(dsdpSolver, DUALVAR);
        if (!(*inCone)) return 0.0;
    } else {
        if (!SDPConic( COPS_CHECK_INCONE )(dsdpSolver, DUALVAR)) {
            printf("| Fatal Error in dsdputils.c -> Line 560 -> SDPConic( COPS_CHECK_INCONE ). Give up. \n");
            printf("---------------------------------------"
                   "---------------------------------------"
                   "--------------------\n");
            printf("| DSDP Ends by Fatal Error. No solution available. \n");
            printf("---------------------------------------"
                   "---------------------------------------"
                   "--------------------\n");
            exit(0);
        }
    }
    
    for (DSDP_INT i = 0; i < dsdpSolver->nBlock; ++i) {
        logdet += spsMatGetlogdet(dsdpSolver->S[i],
                                  dsdpSolver->M->schurAux);
    }
    return logdet;
}

static double LPConic( COPS_GET_LOGDET )
( HSDSolver *dsdpSolver, vec *y, DSDP_INT *inCone ) {
    
    if (!dsdpSolver->isLPset) { return 0.0; }
    double logdet = 0.0;
    if (inCone) { if (!(*inCone)) return 0.0; }
    
    if (dsdpSolver->eventMonitor[EVENT_IN_PHASE_A]) {
        getPhaseALps(dsdpSolver, y, dsdpSolver->tau);
    } else {
        getPhaseBLps(dsdpSolver, y);
    }
    
    for (DSDP_INT i = 0; i < dsdpSolver->lpDim; ++i) {
        logdet += log(dsdpSolver->s->x[i]);
    }
    
    if (inCone) {
        *inCone = (logdet == logdet) ? TRUE : FALSE;
        if (!(*inCone)) return 0.0;
    } else {
        if (logdet != logdet) {
            printf("| Fatal Error in dsdputils.c -> Line 601 -> LPConic( COPS_GET_LOGDET ). Give up. \n");
            printf("---------------------------------------"
                   "---------------------------------------"
                   "--------------------\n");
            printf("| DSDP Ends by Fatal Error. No solution available. \n");
            printf("---------------------------------------"
                   "---------------------------------------"
                   "--------------------\n");
            exit(0);
        }
    }
    
    return logdet;
}

static double BConic( COPS_GET_LOGDET )
( HSDSolver *dsdpSolver, vec *y, DSDP_INT *inCone ) {
    
    if (inCone) { if (!(*inCone)) return 0.0; }
    double logdet = 0.0, bound = dsdpSolver->ybound;
    for (DSDP_INT i = 0; i < dsdpSolver->m; ++i) {
        logdet += log(bound - y->x[i]) + log(y->x[i] + bound);
    }
    if (logdet != logdet && inCone) {*inCone = FALSE; logdet = 0.0;}
    return logdet;
}

static void SDPConic( COPS_GET_AX )
( HSDSolver *dsdpSolver, vec *AX ) {
    
    for (DSDP_INT i = 0; i < dsdpSolver->nBlock; ++i) {
        sdpMatAX(dsdpSolver->sdpData[i], dsdpSolver->dsaux[i], AX);
    }
    return;
}

static void LPConic( COPS_GET_AX )
( HSDSolver *dsdpSolver, vec *Ax ) {
    
    if (!dsdpSolver->isLPset) { return; }
    lpMatAx(dsdpSolver->lpData, dsdpSolver->x, Ax);
    return;
}

static void BConic( COPS_GET_AX )
( HSDSolver *dsdpSolver, vec *AX ) {
    
    return;
}

static double SDPConic( COPS_GET_CX )
( HSDSolver *dsdpSolver ) {
    
    double CX = 0.0; DSDP_INT i;
    for (i = 0; i < dsdpSolver->nBlock; ++i) {
        CX += sdpMatCX(dsdpSolver->sdpData[i], dsdpSolver->dsaux[i]);
    }
    return CX;
}

static double LPConic( COPS_GET_CX )
( HSDSolver *dsdpSolver ) {
    
    if (!dsdpSolver->isLPset) { return 0.0; }
    return lpMatcx(dsdpSolver->lpObj, dsdpSolver->x);
}

static double BConic( COPS_GET_CX )
( HSDSolver *dsdpSolver ) {
    
    return 0.0;
}

extern double DSDPConic( COPS_GET_CX )
( HSDSolver *dsdpSolver ) {
    
    return SDPConic( COPS_GET_CX )(dsdpSolver) + \
           LPConic ( COPS_GET_CX )(dsdpSolver) + \
           BConic  ( COPS_GET_CX )(dsdpSolver);
}

extern void DSDPConic( COPS_GET_AX )
( HSDSolver *dsdpSolver, vec *AX ) {
    
    vec_reset(AX);
    SDPConic( COPS_GET_AX )(dsdpSolver, AX);
    LPConic ( COPS_GET_AX )(dsdpSolver, AX);
    BConic  ( COPS_GET_AX )(dsdpSolver, AX);
}

extern double DSDPConic( COPS_GET_C_FNORM )
( HSDSolver *dsdpSolver ) {
    
    double cnrm = 0.0, tmp;
    tmp = SDPConic( COPS_GET_C_FNORM )(dsdpSolver); cnrm += tmp * tmp;
    tmp = LPConic ( COPS_GET_C_FNORM )(dsdpSolver); cnrm += tmp * tmp;
    tmp = BConic  ( COPS_GET_C_FNORM )(dsdpSolver); cnrm += tmp * tmp;
    return sqrt(cnrm);
}

extern void DSDPConic( COPS_DO_VAR_SCALE )
( HSDSolver *dsdpSolver, double tau ) {
    
    SDPConic( COPS_DO_VAR_SCALE )(dsdpSolver, tau);
    LPConic ( COPS_DO_VAR_SCALE )(dsdpSolver, tau);
    BConic  ( COPS_DO_VAR_SCALE )(dsdpSolver, tau);
}

extern double DSDPConic( COPS_GET_LOGDET )
( HSDSolver *dsdpSolver, vec *y, DSDP_INT *inCone ) {
    
    double logdet = 0.0; if (inCone) *inCone = TRUE;
    logdet += SDPConic( COPS_GET_LOGDET )(dsdpSolver, y, inCone);
    logdet += LPConic ( COPS_GET_LOGDET )(dsdpSolver, y, inCone);
    logdet += BConic  ( COPS_GET_LOGDET )(dsdpSolver, y, inCone);
    return logdet;
}

extern void DSDPConic( COPS_GET_SCHURVEC )
( HSDSolver *dsdpSolver, DSDP_INT dInfeas ) {
    
    SDPConic( COPS_GET_SCHURVEC ) (dsdpSolver, dInfeas);
    LPConic ( COPS_GET_SCHURVEC ) (dsdpSolver, dInfeas);
    BConic  ( COPS_GET_SCHURVEC ) (dsdpSolver, dInfeas);
}

extern void DSDPConic( COPS_STEPDIRECTION )
( HSDSolver *dsdpSolver ) {
    
    SDPConic( COPS_STEPDIRECTION )(dsdpSolver);
    LPConic ( COPS_STEPDIRECTION )(dsdpSolver);
    BConic  ( COPS_STEPDIRECTION )(dsdpSolver);
}

extern double DSDPConic( COPS_GET_MAXSTEP )
( HSDSolver *dsdpSolver, DSDP_INT type ) {
    
    double step = DSDP_INFINITY, tmp;
    tmp = SDPConic( COPS_GET_MAXSTEP )(dsdpSolver, type); step = MIN(step, tmp);
    tmp = LPConic ( COPS_GET_MAXSTEP )(dsdpSolver, type); step = MIN(step, tmp);
    tmp = BConic  ( COPS_GET_MAXSTEP )(dsdpSolver, type); step = MIN(step, tmp);
    return step;
}

extern void DSDPConic( COPS_CONSTR_EXPR )
( HSDSolver *dsdpSolver, DSDP_INT type,
 
    double ycoef, vec *y, double tau, double r ) {
    SDPConic( COPS_CONSTR_EXPR )(dsdpSolver, type, ycoef, y, tau, r);
    LPConic ( COPS_CONSTR_EXPR )(dsdpSolver, type, ycoef, y, tau, r);
    BConic  ( COPS_CONSTR_EXPR )(dsdpSolver, type, ycoef, y, tau, r);
}

extern void DSDPConic( COPS_SYMFAC )
( HSDSolver *dsdpSolver ) {
    
    SDPConic( COPS_SYMFAC )(dsdpSolver);
    LPConic ( COPS_SYMFAC )(dsdpSolver);
    BConic  ( COPS_SYMFAC )(dsdpSolver);
}

extern void DSDPConic ( COPS_GET_SLACK )
( HSDSolver *dsdpSolver, DSDP_INT type ) {
    
    SDPConic( COPS_GET_SLACK )(dsdpSolver, type);
    LPConic ( COPS_GET_SLACK )(dsdpSolver, type);
    BConic  ( COPS_GET_SLACK )(dsdpSolver, type);
}

extern void DSDPConic( COPS_GET_SCHUR )
( HSDSolver *dsdpSolver ) {
    
    SDPConic( COPS_GET_SCHUR )(dsdpSolver);
    LPConic ( COPS_GET_SCHUR )(dsdpSolver);
    // Do not add bound cone for Golden linesearch
    vec_copy(dsdpSolver->asinv, dsdpSolver->d12);
    BConic  ( COPS_GET_SCHUR )(dsdpSolver);
}

extern DSDP_INT DSDPConic( COPS_CHECK_INCONE )
( HSDSolver *dsdpSolver, DSDP_INT type ) {
    
    if (!BConic  ( COPS_CHECK_INCONE )(dsdpSolver, type)) { return FALSE; }
    if (!LPConic ( COPS_CHECK_INCONE )(dsdpSolver, type)) { return FALSE; }
    if (!SDPConic( COPS_CHECK_INCONE )(dsdpSolver, type)) { return FALSE; }
    return TRUE;
}

extern double DSDPConic( COPS_GET_A_ONE_NORM )
( HSDSolver *dsdpSolver ) {
    
    return SDPConic( COPS_GET_A_ONE_NORM )(dsdpSolver) + \
           LPConic ( COPS_GET_A_ONE_NORM )(dsdpSolver) + \
           BConic  ( COPS_GET_A_ONE_NORM )(dsdpSolver);
}

extern double DSDPConic( COPS_GET_C_ONE_NORM )
( HSDSolver *dsdpSolver ) {
    
    return SDPConic( COPS_GET_C_ONE_NORM )(dsdpSolver) + \
           LPConic ( COPS_GET_C_ONE_NORM )(dsdpSolver) + \
           BConic  ( COPS_GET_C_ONE_NORM )(dsdpSolver);
}

extern void DSDPConic( COPS_DO_C_SCALE )
( HSDSolver *dsdpSolver ) {
    
    SDPConic( COPS_DO_C_SCALE )(dsdpSolver);
    LPConic ( COPS_DO_C_SCALE )(dsdpSolver);
    BConic  ( COPS_DO_C_SCALE )(dsdpSolver);
}

extern void DSDPConic( COPS_GET_DOBJ )
( HSDSolver *dsdpSolver ) {
    // Compute b' * y
    checkIterProgress(dsdpSolver, ITER_DUAL_OBJ);
    vec_dot(dsdpSolver->dObj, dsdpSolver->y, &dsdpSolver->dObjVal);
    dsdpSolver->iterProgress[ITER_DUAL_OBJ] = TRUE;
}

// Special operations for Bound cone
extern void getBslack( HSDSolver *dsdpSolver, vec *y, DSDP_INT type ) {
    
    if (type == DUALVAR) {
        vec_lslack(y, dsdpSolver->sl, -dsdpSolver->ybound);
        vec_uslack(y, dsdpSolver->su,  dsdpSolver->ybound);
    } else {
        vec_lslack(y, dsdpSolver->slcker, -dsdpSolver->ybound);
        vec_uslack(y, dsdpSolver->sucker,  dsdpSolver->ybound);
    }
}

// Special operations for LP Cone
extern void getPhaseALps( HSDSolver *dsdpSolver, vec *y, double tau ) {
    
    LPConic( COPS_CONSTR_EXPR ) (dsdpSolver, DUALVAR, -1.0, y, tau, -1.0);
}

extern void getPhaseBLpCheckers( HSDSolver *dsdpSolver, vec *y ) {
    // C - ATy
    LPConic( COPS_CONSTR_EXPR ) (dsdpSolver, CHECKER, -1.0, y, 1.0, 0.0);
}

/* Retrieve S = C - ATy */
extern void getPhaseBLps( HSDSolver *dsdpSolver, vec *y ) {
    
    LPConic( COPS_CONSTR_EXPR ) (dsdpSolver, DUALVAR, -1.0, y, 1.0, 0.0);
}

extern void getPhaseBLpds( HSDSolver *dsdpSolver, double alpha, vec *dy, double beta ) {
    
    LPConic( COPS_CONSTR_EXPR ) (dsdpSolver, DELTAS, -alpha, dy, beta, 0.0);
}

extern void getPhaseALpCheckers( HSDSolver *dsdpSolver, vec *y, double tau ) {
    
    LPConic( COPS_CONSTR_EXPR ) (dsdpSolver, CHECKER, -1.0, y, tau, -1.0);
}

extern void getPhaseALpds( HSDSolver *dsdpSolver, double drate, vec *dy, double dtau ) {
    
    LPConic( COPS_CONSTR_EXPR ) (dsdpSolver, DELTAS, -1.0, dy, dtau, drate);
}

extern void dsdpLpCheckerInCone( HSDSolver *dsdpSolver, DSDP_INT *ispsd ) {
    
    *ispsd = LPConic( COPS_CHECK_INCONE )(dsdpSolver, CHECKER);
}

extern void dsdpLpInCone( HSDSolver *dsdpSolver, DSDP_INT *ispsd ) {
    
    *ispsd = LPConic( COPS_CHECK_INCONE )(dsdpSolver, DUALVAR);
}

extern double getMaxLpstep( HSDSolver *dsdpSolver, DSDP_INT type ) {
    
    return LPConic( COPS_GET_MAXSTEP )(dsdpSolver, type);
}

// Special operations for SDP Cone
/* Retrieve S = - Ry + C * tau - ATy across all the blocks */
extern void getPhaseAS( HSDSolver *dsdpSolver, vec *y, double tau ) {
    
    SDPConic( COPS_CONSTR_EXPR ) (dsdpSolver, DUALVAR, -1.0, y, tau, -1.0);
}

/* Retrieve S = C - ATy */
extern void getPhaseBS( HSDSolver *dsdpSolver, vec *y ) {
    
    SDPConic( COPS_CONSTR_EXPR ) (dsdpSolver, DUALVAR, -1.0, y, 1.0, 0.0);
}

/* Retrieve S for verifying positive definiteness */
extern void getPhaseACheckerS( HSDSolver *dsdpSolver, vec *y, double tau ) {
    
    SDPConic( COPS_CONSTR_EXPR ) (dsdpSolver, CHECKER, -1.0, y, tau, -1.0);
}

extern void getPhaseBCheckerS( HSDSolver *dsdpSolver, vec *y ) {
    // C - ATy
    SDPConic( COPS_CONSTR_EXPR ) (dsdpSolver, CHECKER, -1.0, y, 1.0, 0.0);
}

/* Retrieve dS = drate * Ry + C * dtau - ATdy across all the blocks */
extern void getPhaseAdS( HSDSolver *dsdpSolver, double drate, vec *dy, double dtau ) {
    
    SDPConic( COPS_CONSTR_EXPR ) (dsdpSolver, DELTAS, -1.0, dy, dtau, drate);
}

/* Retrieve dS = C * beta - ATdy * alpha across all the blocks */
extern void getPhaseBdS( HSDSolver *dsdpSolver, double alpha, vec *dy, double beta ) {
    
    SDPConic( COPS_CONSTR_EXPR ) (dsdpSolver, DELTAS, -alpha, dy, beta, 0.0);
}

/* DSDP routine for checking positive definite-ness of matrix */
extern void dsdpCheckerInCone( HSDSolver *dsdpSolver, DSDP_INT *ispsd ) {
    
    *ispsd = SDPConic( COPS_CHECK_INCONE )(dsdpSolver, CHECKER);
}

extern void dsdpInCone( HSDSolver *dsdpSolver, DSDP_INT *ispsd ) {
    
    *ispsd = SDPConic( COPS_CHECK_INCONE )(dsdpSolver, DUALVAR);
}

extern double getMaxSDPstep( HSDSolver *dsdpSolver, DSDP_INT type ) {
    
    return SDPConic( COPS_GET_MAXSTEP )(dsdpSolver, type);
}

/* Coefficient norm computer */
extern void getMatFnorm( HSDSolver *dsdpSolver, DSDP_INT blockid,
                         DSDP_INT constrid, double *nrm ) {
    
    void *data = dsdpSolver->sdpData[blockid]->sdpData[constrid];
    switch (dsdpSolver->sdpData[blockid]->types[constrid]) {
        case MAT_TYPE_ZERO  : *nrm = 0.0; break;
        case MAT_TYPE_DENSE : denseMatFnorm(data, nrm); break;
        case MAT_TYPE_SPARSE: spsMatFnorm(data, nrm); break;
        case MAT_TYPE_RANKK : rkMatFnorm(data, nrm); break;
        default             : break;
    }
}

extern double getMatOneNorm( HSDSolver *dsdpSolver, DSDP_INT blockid, DSDP_INT constrid ) {
    
    void *data = dsdpSolver->sdpData[blockid]->sdpData[constrid];
    switch (dsdpSolver->sdpData[blockid]->types[constrid]) {
        case MAT_TYPE_ZERO  : return 0.0;
        case MAT_TYPE_DENSE : return denseMatOneNorm(data);
        case MAT_TYPE_SPARSE: return spsMatOneNorm(data);
        case MAT_TYPE_RANKK : return r1MatOneNorm(((rkMat *) data)->data[0]);
        default             : assert( FALSE );
    }
}

/* Matrix Scaler */
extern void matRScale( HSDSolver *dsdpSolver, DSDP_INT blockid,
                       DSDP_INT constrid, double scaler) {
    
    if (scaler == 1.0) { return; }
    void *data = dsdpSolver->sdpData[blockid]->sdpData[constrid];
    switch (dsdpSolver->sdpData[blockid]->types[constrid]) {
        case MAT_TYPE_ZERO  : break;
        case MAT_TYPE_DENSE : denseMatRscale((dsMat *) data, scaler); break;
        case MAT_TYPE_SPARSE: spsMatRscale((spsMat *) data, scaler); break;
        case MAT_TYPE_RANKK : rkMatRscale((rkMat *) data, scaler); break;
        default             : break;
    }
}

/* Compute S + alpha * some matrix */
extern void addMattoIter( void *data, DSDP_INT mattype, double alpha,
                             spsMat *iterS, DSDP_INT *sumHash ) {
    // Add a data matrix to iteration S, dS or Scker
    switch (mattype) {
        case MAT_TYPE_SPARSE: spsMataXpbY(alpha, data, 1.0, iterS, sumHash); break;
        case MAT_TYPE_DENSE : spsMatAddds(iterS, alpha, data); break;
        case MAT_TYPE_RANKK : spsMatAddrk(iterS, alpha, data, sumHash); break;
        default             : break;
    }
}

extern void addMattoS( HSDSolver *dsdpSolver,
                       DSDP_INT   blockid,
                       DSDP_INT   constrid,
                       double     alpha ) {
    
    if (dsdpSolver->sdpData[blockid]->types[constrid] == MAT_TYPE_ZERO
        || alpha == 0.0) { return; }
    addMattoIter(dsdpSolver->sdpData[blockid]->sdpData[constrid],
                 dsdpSolver->sdpData[blockid]->types[constrid],
                 alpha, dsdpSolver->S[blockid], dsdpSolver->symS[blockid]);
}

/* Compute checkerS + alpha * some matrix */
extern void addMattoChecker( HSDSolver *dsdpSolver,
                             DSDP_INT   blockid,
                             DSDP_INT   constrid,
                             double     alpha ) {
    
    if (dsdpSolver->sdpData[blockid]->types[constrid] == MAT_TYPE_ZERO
        || alpha == 0.0) { return; }
    addMattoIter(dsdpSolver->sdpData[blockid]->sdpData[constrid],
                 dsdpSolver->sdpData[blockid]->types[constrid],
                 alpha, dsdpSolver->Scker[blockid], dsdpSolver->symS[blockid]);
}

extern void addMattodS( HSDSolver *dsdpSolver,
                        DSDP_INT   blockid,
                        DSDP_INT   constrid,
                        double     alpha ) {
    
    if (dsdpSolver->sdpData[blockid]->types[constrid] == MAT_TYPE_ZERO
        || alpha == 0.0) { return; }
    addMattoIter(dsdpSolver->sdpData[blockid]->sdpData[constrid],
                 dsdpSolver->sdpData[blockid]->types[constrid],
                 alpha, dsdpSolver->dS[blockid], dsdpSolver->symS[blockid]);
}
