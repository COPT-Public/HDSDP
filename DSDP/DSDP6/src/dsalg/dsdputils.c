#include "dsdputils.h"
#include "dsdpsolver.h"
static char etype[] = "DSDP Utility";

/*
 DSDP utility routines that manage the operations between
 different types of matrices
*/

extern DSDP_INT getSinvASinv( HSDSolver *dsdpSolver, DSDP_INT blockid, DSDP_INT constrid,
                              void *SinvASinv ) {
    
    // Given S and A, the routine computes A, asinv and trace(S, Sinv A Sinv)
    DSDP_INT retcode = DSDP_RETCODE_OK;
    DSDP_INT typeA = dsdpSolver->sdpData[blockid]->types[constrid];
    spsMat *S = dsdpSolver->S[blockid];
    void *A = dsdpSolver->sdpData[blockid]->sdpData[constrid];
    
    double *asinv = NULL;
    double *asinvrysinv = NULL;
    double tracediag = 0.0;
    
    if (constrid < dsdpSolver->m) {
        asinv = &dsdpSolver->asinv->x[constrid];
        asinvrysinv = &dsdpSolver->d4->x[constrid];
    } else {
        asinv = &dsdpSolver->csinv;
        asinvrysinv = &dsdpSolver->csinvrysinv;
    }
    
    if (typeA == MAT_TYPE_RANK1) {
        r1Mat *dataA = (r1Mat *) A;
        r1Mat *dataSinvASinv = (r1Mat *) SinvASinv;
        retcode = spsSinvR1SinvSolve(S, dataA, dataSinvASinv, asinv); checkCode;
    } else if (typeA == MAT_TYPE_SPARSE) {
        spsMat *dataA = (spsMat *) A;
        dsMat *dataSinvASinv = (dsMat *) SinvASinv;
        retcode = spsSinvSpSinvSolve(S, dataA, dataSinvASinv, asinv); checkCode;
    } else if (typeA == MAT_TYPE_DENSE) {
        dsMat *dataA = (dsMat *) A;
        dsMat *dataSinvASinv = (dsMat *) SinvASinv;
        retcode = spsSinvDsSinvSolve(S, dataA, dataSinvASinv, asinv); checkCode;
    } else {
        error(etype, "Invalid matrix type. \n");
    }
    
    switch (typeA) {
        case MAT_TYPE_RANK1:
            retcode = r1MatdiagTrace(SinvASinv, dsdpSolver->Ry, &tracediag);
            break;
        default:
            retcode = denseDiagTrace(SinvASinv, dsdpSolver->Ry, &tracediag);
            break;
    }
    
    checkCode;
    *asinvrysinv += tracediag;
    
    return retcode;
}

extern DSDP_INT getTraceASinvASinv( HSDSolver *dsdpSolver, DSDP_INT blockid, DSDP_INT constrid,
                                    DSDP_INT mattype, DSDP_INT constrid2, void *SinvASinv ) {
    
    // Compute trace between SinvASinv and some A
    // constrid is the position of (A in A * (.)) and constrid2 is the position of (A in (.) * SinvASinv)
    DSDP_INT retcode = DSDP_RETCODE_OK;
    DSDP_INT Atype = dsdpSolver->sdpData[blockid]->types[constrid];
    DSDP_INT m = dsdpSolver->m;
    void *A = dsdpSolver->sdpData[blockid]->sdpData[constrid];
    double trace = 0.0;
    double *M = dsdpSolver->Msdp->array;
    
    assert( constrid <= constrid2 );
    assert( constrid <= m && constrid2 <= m);
    
    if (mattype == MAT_TYPE_ZERO || Atype == MAT_TYPE_ZERO) {
        return retcode;
    }
    
    if (mattype == MAT_TYPE_RANK1) {
        switch (Atype) {
            case MAT_TYPE_DENSE:
                retcode = r1MatdenseTrace(SinvASinv, A, &trace);
                break;
            case MAT_TYPE_SPARSE:
                retcode = r1MatspsTrace(SinvASinv, A, &trace);
                break;
            case MAT_TYPE_RANK1:
                retcode = r1Matr1Trace(SinvASinv, A, &trace);
                break;
            default:
                error(etype, "Invalid matrix type. \n");
                break;
        }
    } else if (mattype == MAT_TYPE_DENSE) {
        switch (Atype) {
            case MAT_TYPE_DENSE:
                retcode = denseDsTrace(SinvASinv, A, &trace);
                break;
            case MAT_TYPE_SPARSE:
                retcode = denseSpsTrace(SinvASinv, A, &trace);
                break;
            case MAT_TYPE_RANK1:
                retcode = r1MatdenseTrace(A, SinvASinv, &trace);
                break;
            default:
                error(etype, "Invalid matrix type. \n");
                break;
        }
    } else {
        error(etype, "Invalid matrix type. \n");
    }
    
    checkCode;
    
    // Update Schur/auxiliary vectors
    if (constrid2 == m) {
        // The first A is C
        if (constrid == m) {
            dsdpSolver->csinvcsinv += trace;
        } else {
            dsdpSolver->u->x[constrid] += trace;
        }
    } else {
        packIdx(M, m, constrid2, constrid) += trace;
    }
    
    return retcode;
}

/* Retrieve S = - Ry + C * tau - ATy across all the blocks */
extern DSDP_INT getPhaseAS( HSDSolver *dsdpSolver, double *y, double tau ) {
    
    DSDP_INT retcode = DSDP_RETCODE_OK;
    spsMat *S = NULL;
    DSDP_INT m = dsdpSolver->m;
    
    for (DSDP_INT i = 0; i < dsdpSolver->nBlock; ++i) {
        S = dsdpSolver->S[i];
        retcode = spsMatReset(S);
        for (DSDP_INT j = 0; j < m; ++j) {
            retcode = addMattoS(dsdpSolver, i, j, - y[j]);
        }
        retcode = addMattoS(dsdpSolver, i, m + 1, tau);
    }
    
    retcode = spsMatAdddiag(S, - dsdpSolver->Ry);
    checkCode;
    return retcode;
}

/* Retrieve dS = Ry + C * dtau - ATdy across all the blocks */
extern DSDP_INT getPhaseAdS( HSDSolver *dsdpSolver, double *dy, double dtau ) {
    
    DSDP_INT retcode = DSDP_RETCODE_OK;
    spsMat *dS = NULL;
    DSDP_INT m = dsdpSolver->m;
    
    for (DSDP_INT i = 0; i < dsdpSolver->nBlock; ++i) {
        dS = dsdpSolver->dS[i];
        retcode = spsMatReset(dS);
        for (DSDP_INT j = 0; j < m; ++j) {
            retcode = addMattodS(dsdpSolver, i, j, - dy[j]);
        }
        retcode = addMattodS(dsdpSolver, i, m + 1, dtau);
    }
    
    retcode = spsMatAdddiag(dS, dsdpSolver->Ry);
    checkCode;
    return retcode;
}

/* DSDP routine for checking positive definite-ness of matrix */
extern DSDP_INT dsdpInCone( HSDSolver *dsdpSolver, DSDP_INT *ispsd ) {
    // Determine whether the current dual variable lies in the cone
    DSDP_INT retcode = DSDP_RETCODE_OK;
    DSDP_INT incone = FALSE;
    
    for (DSDP_INT i = 0; i < dsdpSolver->nBlock; ++i) {
        spsMatIspd(dsdpSolver->S[i], &incone);
        if (!incone) {
            break;
        }
    }
    
    *ispsd = incone;
    return retcode;
}

/* Coefficient norm computer */
extern DSDP_INT getMatnrm( HSDSolver *dsdpSolver, DSDP_INT blockid, DSDP_INT constrid, double *nrm ) {
    
    DSDP_INT retcode = DSDP_RETCODE_OK;
    void *data = dsdpSolver->sdpData[blockid]->sdpData[constrid];
    
    switch (dsdpSolver->sdpData[blockid]->types[constrid]) {
        case MAT_TYPE_ZERO:
            *nrm = 0.0;
            break;
        case MAT_TYPE_DENSE:
            retcode = denseMatFnorm(data, nrm);
            checkCode;
            break;
        case MAT_TYPE_SPARSE:
            retcode = spsMatFnorm(data, nrm);
            checkCode;
            break;
        case MAT_TYPE_RANK1:
            retcode = r1MatFnorm(data, nrm);
            checkCode;
            break;
        default:
            error(etype, "Unknown matrix type. \n");
            break;
    }
    
    return retcode;
}

/* Compute S + alpha * some matrix */
extern DSDP_INT addMattoS( HSDSolver *dsdpSolver, DSDP_INT blockid, DSDP_INT constrid, double alpha ) {
    
    DSDP_INT retcode = DSDP_RETCODE_OK;
    void *data = dsdpSolver->sdpData[blockid]->sdpData[constrid];
    
    switch (dsdpSolver->sdpData[blockid]->types[constrid]) {
        case MAT_TYPE_ZERO:
            break;
        case MAT_TYPE_DENSE:
            retcode = spsMatAddds(dsdpSolver->S[blockid], alpha, data);
            checkCode;
            break;
        case MAT_TYPE_SPARSE:
            retcode = spsMataXpbY(alpha, data, 1.0, dsdpSolver->S[blockid]);
            checkCode;
            break;
        case MAT_TYPE_RANK1:
            retcode = spsMatAddr1(dsdpSolver->S[blockid], alpha, data);
            checkCode;
            break;
        default:
            error(etype, "Unknown matrix type. \n");
            break;
    }
    
    return retcode;
}

extern DSDP_INT addMattodS( HSDSolver *dsdpSolver, DSDP_INT blockid, DSDP_INT constrid, double alpha ) {
    
    DSDP_INT retcode = DSDP_RETCODE_OK;
    void *data = dsdpSolver->sdpData[blockid]->sdpData[constrid];
    
    switch (dsdpSolver->sdpData[blockid]->types[constrid]) {
        case MAT_TYPE_ZERO:
            break;
        case MAT_TYPE_DENSE:
            retcode = spsMatAddds(dsdpSolver->dS[blockid], alpha, data);
            checkCode;
            break;
        case MAT_TYPE_SPARSE:
            retcode = spsMataXpbY(alpha, data, 1.0, dsdpSolver->dS[blockid]);
            checkCode;
            break;
        case MAT_TYPE_RANK1:
            retcode = spsMatAddr1(dsdpSolver->dS[blockid], alpha, data);
            checkCode;
            break;
        default:
            error(etype, "Unknown matrix type. \n");
            break;
    }
    
    return retcode;
}

/* Free various matrices */
extern DSDP_INT freesdpMat( HSDSolver *dsdpSolver, DSDP_INT blockid, DSDP_INT constrid ) {
    
    DSDP_INT retcode = DSDP_RETCODE_OK;
    void *data = dsdpSolver->sdpData[blockid]->sdpData[constrid];
    switch (dsdpSolver->sdpData[blockid]->types[constrid]) {
        case MAT_TYPE_ZERO:
            break;
        case MAT_TYPE_DENSE:
            retcode = denseMatFree(data); checkCode;
            break;
        case MAT_TYPE_SPARSE:
            retcode = spsMatFree(data); checkCode;
            checkCode;
            break;
        case MAT_TYPE_RANK1:
            retcode = r1MatFree(data); checkCode;
            break;
        default:
            error(etype, "Unknown matrix type. \n");
            break;
    }
    
    return retcode;
}

/* Objective value computer */
extern DSDP_INT getDualObj( HSDSolver *dsdpSolver ) {
    
    // Compute b' * y
    DSDP_INT retcode = DSDP_RETCODE_OK;
    retcode = checkIterProgress(dsdpSolver, ITER_DUAL_OBJ);
    checkCode;
    vec *ydata = dsdpSolver->y;
    vec *bdata = dsdpSolver->dObj;
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
    
    double pObjVal  = 0.0;
    DSDP_INT m    = dsdpSolver->m;
    double  mu    = dsdpSolver->mu;
    DSDP_INT incx = 1;
    
    vec *asinv = dsdpSolver->asinv;

    double *x1 = asinv->x;
    double *d2 = dsdpSolver->d2->x;
        
    pObjVal = (double) dsdpSolver->n + dot(&m, x1, &incx, d2, &incx);
    dsdpSolver->pObjVal = pObjVal * mu + dsdpSolver->dObjVal;
    
    return retcode;
}
