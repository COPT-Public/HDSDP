#include "dsdppsol.h"
#include "dsdputils.h"

extern DSDP_INT computePrimalX( HSDSolver *dsdpSolver ) {
    
    // Extract primal solution of the current block
    DSDP_INT retcode = DSDP_RETCODE_OK;
    
    assert( dsdpSolver->solStatus == DSDP_OPTIMAL );
    
    DSDP_INT nblock = dsdpSolver->nBlock, dim = 0;
    vec *dymaker = dsdpSolver->dymaker;
    vec *ymaker = dsdpSolver->ymaker;
    spsMat *Smaker = NULL, *bnmaker = NULL;
    double mumaker = dsdpSolver->mumaker;
    double *Xtmp = NULL;
    
    // getPhaseBS(dsdpSolver, dsdpSolver->y->x);
    // Smaker = C - dsdpgetATy(A, ymaker);
    getPhaseBCheckerS(dsdpSolver, ymaker->x);
    // bnmaker = dsdpgetATy(A, dymaker);
    getPhaseBdS(dsdpSolver, -1.0, dymaker->x, 0.0);
    
    dsMat  *dsaux = NULL;
    rkMat  *rkaux = NULL;
    spsMat *spaux = NULL;
    
    for (DSDP_INT blockid = 0; blockid < nblock; ++blockid) {
    
        dsaux   = dsdpSolver->dsaux[blockid];
        rkaux   = dsdpSolver->rkaux[blockid];
        spaux   = dsdpSolver->spaux[blockid];
        Smaker  = dsdpSolver->Scker[blockid];
        bnmaker = dsdpSolver->dS[blockid];
        dim     = bnmaker->dim;
        Xtmp    = (double *) calloc(dim * dim, sizeof(double));
        spsMatFactorize(Smaker);
        spsMatGetX(Smaker, bnmaker, Xtmp);
        
        double *Stmp = (double *) calloc(dim * dim, sizeof(double));
        spsMatFill(dsdpSolver->S[blockid], Stmp);
        
        DSDP_FREE(Stmp);
        
        // X = mu * D * (speye(n) + D' * bnmaker * D) * D';
        for (DSDP_INT i = 0, idx = 0; i < dim; ++i) {
            memcpy(&dsaux->array[idx], &Xtmp[i * dim + i],
                   sizeof(double) * (dim - i));
            idx += dim - i;
        }
        
        for (DSDP_INT i = 0; i < nsym(dim); ++i) {
            dsaux->array[i] *= mumaker;
        }
        
        DSDP_FREE(Xtmp);
    }
    
    return retcode;
}

extern DSDP_INT computeDIMACS( HSDSolver *dsdpSolver,
                               double *err1, double *err2,
                               double *err3, double *err4,
                               double *err5, double *err6 ) {
    
    /* Compute the DIMACS error after solution is setup */
    DSDP_INT retcode = DSDP_RETCODE_OK;
    
    double bnrm, Cnrm, pObj, dObj, dInf, pInf,
           gap, trace, tmp, minEigX, minEigS, compslack;
    DSDP_INT nblock = dsdpSolver->nBlock, m = dsdpSolver->m;
    vec_norm(dsdpSolver->dObj, &bnrm);
    
    DSDP_INT dim = dsdpSolver->dsaux[0]->dim;
    double res = 0.0;
    double *aux1 = (double *) calloc(dim * dim, sizeof(double));
    double *aux2 = (double *) calloc(dim * dim, sizeof(double));
    spsMatFill(dsdpSolver->S[0], aux1);
    denseMatFill(dsdpSolver->dsaux[0], aux2);
    
    for (DSDP_INT i = 0; i < dim; ++i) {
        for (DSDP_INT j = 0; j < dim; ++j) {
            res += aux1[i * dim + j] * aux2[i * dim + j];
        }
    }
    
    DSDP_FREE(aux1);
    DSDP_FREE(aux2);
    
    // TODO: Replace Cnrm by 1-norm
    Cnrm = 0.0;
    for (DSDP_INT i = 0; i < dsdpSolver->nBlock; ++i) {
        getMatnrm(dsdpSolver, i, dsdpSolver->m, &tmp);
        Cnrm += tmp * tmp;
    }
    Cnrm = sqrt(Cnrm);
    
    vec_dot(dsdpSolver->dObj, dsdpSolver->y, &dObj);
    
    /*  DIMACS Error 1    */
    pInf = 0.0;
    for (DSDP_INT i = 0; i < m; ++i) {
        trace = 0.0;
        for (DSDP_INT j = 0; j < nblock; ++j) {
            rkMatdenseTrace(dsdpSolver->sdpData[j]->sdpData[i],
                            dsdpSolver->dsaux[j], &tmp);
            trace += tmp;
        }
        tmp = dsdpSolver->dObj->x[i] - trace;
        pInf += tmp * tmp;
    }
        
    /*  DIMACS Error 2    */
    minEigS = DSDP_INFINITY;
    for (DSDP_INT i = 0; i < dsdpSolver->nBlock; ++i) {
        spsMatMinEig(dsdpSolver->S[i], &tmp);
        minEigS = MIN(minEigS, tmp);
    }
    
    /*  DIMACS Error 3    */
    dInf = 0.0;
    
    pObj = 0.0;
    for (DSDP_INT i = 0; i < nblock; ++i) {
        rkMatdenseTrace(dsdpSolver->sdpData[i]->sdpData[m],
                        dsdpSolver->dsaux[i], &tmp);
        pObj += tmp;
    }
    
    /* DIMACS Error 4     */
    minEigX = DSDP_INFINITY;
    for (DSDP_INT i = 0; i < dsdpSolver->nBlock; ++i) {
        denseMatMinEig(dsdpSolver->dsaux[i], &tmp);
        minEigX = MIN(tmp, minEigX);
    }

    /*  DIMACS Error 5    */
    gap = pObj - dObj;
    
    /*  DIMACS Error 6    */
    compslack = 0.0;
    for (DSDP_INT i = 0; i < nblock; ++i) {
        denseSpsTrace(dsdpSolver->dsaux[i], dsdpSolver->S[i], &tmp);
        compslack += tmp;
    }
    
    // Collect errors
    *err1 = pInf / (1 + bnrm);
    *err2 = MAX(0.0, -minEigX) / (1 + bnrm);
    *err3 = 0.0;
    *err4 = MAX(0.0, -minEigS) / (1 + Cnrm);
    *err5 = gap / (1 + fabs(pObj) + fabs(dObj));
    *err6 = compslack / (1 + fabs(pObj) + fabs(dObj));

    return retcode;
}
