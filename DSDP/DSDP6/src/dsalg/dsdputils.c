#include "dsdputils.h"
#include "dsdpsolver.h"
static char etype[] = "DSDP Utility";

/*
  Utility routines that manage the operations between
  different types of matrices
*/

extern void invertDualVars( HSDSolver *dsdpSolver ) {
    // Compute inverse of the dual matrix when M3, M4 or M5 techniques are used
    // Also used if corrector is being computed
    DSDPSchur *M = dsdpSolver->M;
    for (DSDP_INT i = 0; i < M->nblock; ++i) {
        spsMatInverse(M->S[i], M->Sinv[i], M->schurAux);
    }
}

/* Retrieve S = - Ry + C * tau - ATy across all the blocks */
extern DSDP_INT getPhaseAS( HSDSolver *dsdpSolver, double *y, double tau ) {
    
    DSDP_INT retcode = DSDP_RETCODE_OK;
    spsMat *S = NULL; DSDP_INT m = dsdpSolver->m;
    for (DSDP_INT i = 0; i < dsdpSolver->nBlock; ++i) {
        S = dsdpSolver->S[i]; retcode = spsMatReset(S);
        for (DSDP_INT j = 0; j < m; ++j) {
            retcode = addMattoS(dsdpSolver, i, j, - y[j]);
        }
        addMattoS(dsdpSolver, i, m, tau);
        spsMatAdddiag(S, - dsdpSolver->Ry, dsdpSolver->symS[i]);
    }
    return retcode;
}

/* Retrieve S = C - ATy */
extern DSDP_INT getPhaseBS( HSDSolver *dsdpSolver, double *y ) {
    
    DSDP_INT retcode = DSDP_RETCODE_OK;
    spsMat *S = NULL; DSDP_INT m = dsdpSolver->m;
    for (DSDP_INT i = 0; i < dsdpSolver->nBlock; ++i) {
        S = dsdpSolver->S[i]; retcode = spsMatReset(S);
        for (DSDP_INT j = 0; j < m; ++j) {
            retcode = addMattoS(dsdpSolver, i, j, - y[j]);
        }
        addMattoS(dsdpSolver, i, m, 1.0);
        spsMatAdddiag(S, dsdpSolver->dperturb, dsdpSolver->symS[i]);
    }
    return retcode;
}

/* Retrieve S for verifying positive definiteness */
extern DSDP_INT getPhaseACheckerS( HSDSolver *dsdpSolver, double *y, double tau ) {
    
    DSDP_INT retcode = DSDP_RETCODE_OK;
    spsMat *checker = NULL; DSDP_INT m = dsdpSolver->m;
    for (DSDP_INT i = 0; i < dsdpSolver->nBlock; ++i) {
        checker = dsdpSolver->Scker[i];
        retcode = spsMatReset(checker);
        for (DSDP_INT j = 0; j < m; ++j) {
            retcode = addMattoChecker(dsdpSolver, i, j, - y[j]);
        }
        retcode = addMattoChecker(dsdpSolver, i, m, tau);
        retcode = spsMatAdddiag(checker, - dsdpSolver->Ry, dsdpSolver->symS[i]);
    }
    return retcode;
}

extern DSDP_INT getPhaseBCheckerS( HSDSolver *dsdpSolver, double *y ) {
    // C - ATy
    DSDP_INT retcode = DSDP_RETCODE_OK;
    spsMat *checker = NULL; DSDP_INT m = dsdpSolver->m;
    for (DSDP_INT i = 0; i < dsdpSolver->nBlock; ++i) {
        checker = dsdpSolver->Scker[i];
        retcode = spsMatReset(checker);
        for (DSDP_INT j = 0; j < m; ++j) {
            retcode = addMattoChecker(dsdpSolver, i, j, - y[j]);
        }
        retcode = addMattoChecker(dsdpSolver, i, m, 1.0);
        retcode = spsMatAdddiag(dsdpSolver->Scker[i], dsdpSolver->dperturb, dsdpSolver->symS[i]);
    }
    return retcode;
}

/* Retrieve dS = drate * Ry + C * dtau - ATdy across all the blocks */
extern DSDP_INT getPhaseAdS( HSDSolver *dsdpSolver, double drate, double *dy, double dtau ) {
    
    DSDP_INT retcode = DSDP_RETCODE_OK;
    spsMat *dS = NULL; DSDP_INT m = dsdpSolver->m;
    for (DSDP_INT i = 0; i < dsdpSolver->nBlock; ++i) {
        dS = dsdpSolver->dS[i]; retcode = spsMatReset(dS);
        for (DSDP_INT j = 0; j < m; ++j) {
            retcode = addMattodS(dsdpSolver, i, j, - dy[j]);
        }
        addMattodS(dsdpSolver, i, m, dtau);
        spsMatAdddiag(dS, drate * dsdpSolver->Ry, dsdpSolver->symS[i]);
    }
    return retcode;
}

/* Retrieve dS = C * beta - ATdy * alpha across all the blocks */
extern DSDP_INT getPhaseBdS( HSDSolver *dsdpSolver, double alpha,
                             double *dy, double beta ) {
    
    DSDP_INT retcode = DSDP_RETCODE_OK;
    spsMat *dS = NULL; DSDP_INT m = dsdpSolver->m;
    for (DSDP_INT i = 0; i < dsdpSolver->nBlock; ++i) {
        dS = dsdpSolver->dS[i]; retcode = spsMatReset(dS);
        if (alpha) {
            for (DSDP_INT j = 0; j < m; ++j) {
                retcode = addMattodS(dsdpSolver, i, j, - alpha * dy[j]);
            }
        }
        if (beta) { retcode = addMattodS(dsdpSolver, i, m, beta); }
    }
    return retcode;
}

/* DSDP routine for checking positive definite-ness of matrix */
extern DSDP_INT dsdpCheckerInCone( HSDSolver *dsdpSolver, DSDP_INT *ispsd ) {
    // Determine whether the current checking auxiliary variable lies in the cone
    DSDP_INT retcode = DSDP_RETCODE_OK;
    DSDP_INT incone = FALSE;
    for (DSDP_INT i = 0; i < dsdpSolver->nBlock; ++i) {
        spsMatIspd(dsdpSolver->Scker[i], &incone);
        if (!incone) { break; }
    }
    *ispsd = incone; return retcode;
}

extern DSDP_INT dsdpInCone( HSDSolver *dsdpSolver, DSDP_INT *ispsd ) {
    // Determine whether the current dual variable lies in the cone
    DSDP_INT retcode = DSDP_RETCODE_OK;
    DSDP_INT incone = FALSE;
    for (DSDP_INT i = 0; i < dsdpSolver->nBlock; ++i) {
        spsMatIspd(dsdpSolver->S[i], &incone);
        if (!incone) { break; }
    }
    *ispsd = incone; return retcode;
}

/* Coefficient norm computer */
extern DSDP_INT getMatFnorm( HSDSolver *dsdpSolver, DSDP_INT blockid, DSDP_INT constrid, double *nrm ) {
    
    DSDP_INT retcode = DSDP_RETCODE_OK;
    void *data = dsdpSolver->sdpData[blockid]->sdpData[constrid];
    switch (dsdpSolver->sdpData[blockid]->types[constrid]) {
        case MAT_TYPE_ZERO  : *nrm = 0.0; break;
        case MAT_TYPE_DENSE : denseMatFnorm(data, nrm); break;
        case MAT_TYPE_SPARSE: spsMatFnorm(data, nrm); break;
        case MAT_TYPE_RANKK : rkMatFnorm(data, nrm); break;
        default             : error(etype, "Unknown matrix type. \n"); break;
    }
    
    return retcode;
}

extern DSDP_INT getMatOneNorm( HSDSolver *dsdpSolver, DSDP_INT blockid, DSDP_INT constrid, double *nrm ) {
    
    DSDP_INT retcode = DSDP_RETCODE_OK;
    void *data = dsdpSolver->sdpData[blockid]->sdpData[constrid];
    switch (dsdpSolver->sdpData[blockid]->types[constrid]) {
        case MAT_TYPE_ZERO  : *nrm = 0.0; break;
        case MAT_TYPE_DENSE : denseMatOneNorm(data, nrm); break;
        case MAT_TYPE_SPARSE: spsMatOneNorm(data, nrm); break;
        case MAT_TYPE_RANKK : r1MatOneNorm(data, nrm); break;
        default             : error(etype, "Unknown matrix type. \n"); break;
    }
    return retcode;
}

/* Matrix Scaler */
extern DSDP_INT matRScale( HSDSolver *dsdpSolver, DSDP_INT blockid, DSDP_INT constrid, double scaler) {
    
    DSDP_INT retcode = DSDP_RETCODE_OK;
    if (scaler == 1.0) { return retcode; }
    void *data = dsdpSolver->sdpData[blockid]->sdpData[constrid];
    switch (dsdpSolver->sdpData[blockid]->types[constrid]) {
        case MAT_TYPE_ZERO  : break;
        case MAT_TYPE_DENSE : retcode = denseMatRscale((dsMat *) data, scaler); break;
        case MAT_TYPE_SPARSE: retcode = spsMatRscale((spsMat *) data, scaler); break;
        case MAT_TYPE_RANKK : retcode = rkMatRscale((rkMat *) data, scaler); break;
        default             : error(etype, "Unknown matrix type. \n"); break;
    }
    return retcode;
}

/* Compute S + alpha * some matrix */
extern DSDP_INT addMattoIter( void *data, DSDP_INT mattype, double alpha, spsMat *iterS, DSDP_INT *sumHash ) {
    // Add a data matrix to iteration S, dS or Scker
    DSDP_INT retcode = DSDP_RETCODE_OK;
    if (alpha == 0.0) { return retcode; }
    switch (mattype) {
        case MAT_TYPE_SPARSE: spsMataXpbY(alpha, data, 1.0, iterS, sumHash); break;
        case MAT_TYPE_DENSE : spsMatAddds(iterS, alpha, data); break;
        case MAT_TYPE_RANKK : spsMatAddrk(iterS, alpha, data, sumHash); break;
        default             : error(etype, "Invalid matrix type. \n"); break;
    }
    return retcode;
}

extern DSDP_INT addMattoS( HSDSolver *dsdpSolver, DSDP_INT blockid, DSDP_INT constrid, double alpha ) {
    
    DSDP_INT retcode = DSDP_RETCODE_OK;
    if (dsdpSolver->sdpData[blockid]->types[constrid] == MAT_TYPE_ZERO ||
        alpha == 0.0) { return retcode; }
    void *data = dsdpSolver->sdpData[blockid]->sdpData[constrid];
    retcode = addMattoIter(data, dsdpSolver->sdpData[blockid]->types[constrid],
                           alpha, dsdpSolver->S[blockid], dsdpSolver->symS[blockid]);
    return retcode;
}

/* Compute checkerS + alpha * some matrix */
extern DSDP_INT addMattoChecker( HSDSolver *dsdpSolver, DSDP_INT blockid, DSDP_INT constrid, double alpha ) {
    
    DSDP_INT retcode = DSDP_RETCODE_OK;
    if (dsdpSolver->sdpData[blockid]->types[constrid] == MAT_TYPE_ZERO ||
        alpha == 0.0) { return retcode; }
    void *data = dsdpSolver->sdpData[blockid]->sdpData[constrid];
    retcode = addMattoIter(data, dsdpSolver->sdpData[blockid]->types[constrid],
                           alpha, dsdpSolver->Scker[blockid], dsdpSolver->symS[blockid]);
    return retcode;
}

extern DSDP_INT addMattodS( HSDSolver *dsdpSolver, DSDP_INT blockid, DSDP_INT constrid, double alpha ) {
    
    DSDP_INT retcode = DSDP_RETCODE_OK;
    if (dsdpSolver->sdpData[blockid]->types[constrid] == MAT_TYPE_ZERO ||
        alpha == 0.0) { return retcode; }
    void *data = dsdpSolver->sdpData[blockid]->sdpData[constrid];
    retcode = addMattoIter(data, dsdpSolver->sdpData[blockid]->types[constrid],
                           alpha, dsdpSolver->dS[blockid], dsdpSolver->symS[blockid]);
    return retcode;
}

/* Free various matrices */
extern DSDP_INT freesdpMat( HSDSolver *dsdpSolver, DSDP_INT blockid, DSDP_INT constrid ) {
    
    DSDP_INT retcode = DSDP_RETCODE_OK;
    void *data = dsdpSolver->sdpData[blockid]->sdpData[constrid];
    switch (dsdpSolver->sdpData[blockid]->types[constrid]) {
        case MAT_TYPE_ZERO  : break;
        case MAT_TYPE_DENSE : denseMatFree(data); break;
        case MAT_TYPE_SPARSE: spsMatFree(data); break;
        case MAT_TYPE_RANKK : rkMatFree(data); break;
        default             : error(etype, "Unknown matrix type. \n"); break;
    }
    return retcode;
}

/* Objective value computer */
extern DSDP_INT getDualObj( HSDSolver *dsdpSolver ) {
    // Compute b' * y
    DSDP_INT retcode = DSDP_RETCODE_OK;
    checkIterProgress(dsdpSolver, ITER_DUAL_OBJ);
    vec *ydata = dsdpSolver->y, *bdata = dsdpSolver->dObj;
    vec_dot(bdata, ydata, &dsdpSolver->dObjVal);
    dsdpSolver->iterProgress[ITER_DUAL_OBJ] = TRUE;
    return retcode;
}

extern DSDP_INT getSDPPrimalObjPhaseB( HSDSolver *dsdpSolver ) {
    // Compute primal objective given primal feasible projection
    DSDP_INT retcode = DSDP_RETCODE_OK;
    if (dsdpSolver->eventMonitor[EVENT_NO_RY] &&
        dsdpSolver->eventMonitor[EVENT_PFEAS_FOUND] &&
        dsdpSolver->eventMonitor[EVENT_IN_PHASE_B]) {
    } else {
        return retcode;
    }
    
    DSDP_INT m = dsdpSolver->m, incx = 1;
    vec *asinv = dsdpSolver->asinv;
    double pObjVal = 0.0, mu = dsdpSolver->mu, *x1 = asinv->x, *d2 = dsdpSolver->d2->x;;
    pObjVal = (double) (dsdpSolver->n + dsdpSolver->m * 2) + dot(&m, x1, &incx, d2, &incx);
    dsdpSolver->pObjVal = pObjVal * mu + dsdpSolver->dObjVal;
    return retcode;
}
