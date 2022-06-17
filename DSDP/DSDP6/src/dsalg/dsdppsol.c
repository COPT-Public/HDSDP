#include "dsdppsol.h"
#include "dsdputils.h"
#include "dsdplog.h"
#include "schurmat.h"

static char etype[] = "Primal solution extraction";

static double computePres( HSDSolver *dsdpSolver ) {
    // Compute ||A * x - b||
    double tmp; vec *AX = dsdpSolver->d1;
    DSDPConic( COPS_GET_AX )(dsdpSolver, AX);
    vec_axpy(-1.0, dsdpSolver->dObj, AX); vec_norm(AX, &tmp);
    return tmp;
}

extern DSDP_INT computePrimalX( HSDSolver *dsdpSolver ) {
    // Extract primal solution of the current block
    DSDP_INT retcode = DSDP_RETCODE_OK;
    DSDP_INT nblock = dsdpSolver->nBlock, dim = 0, i, idx;
    vec *dymaker = dsdpSolver->dymaker, *ymaker = dsdpSolver->ymaker;
    spsMat *Smaker = NULL, *bnmaker = NULL;
    double mumaker = dsdpSolver->mumaker, *Xtmp = NULL, start = my_clock();
    getPhaseBS(dsdpSolver, dsdpSolver->y); getPhaseBLps(dsdpSolver, dsdpSolver->y);
    getPhaseBCheckerS(dsdpSolver, ymaker); getPhaseBLpCheckers(dsdpSolver, ymaker);
    getPhaseBdS(dsdpSolver, -1.0, dymaker, 0.0); getPhaseBLpds(dsdpSolver, -1.0, dymaker, 0.0);
    dsMat *dsaux = NULL; rkMat *rkaux = NULL;
    
    for (DSDP_INT blockid = 0; blockid < nblock; ++blockid) {
        dsaux = dsdpSolver->dsaux[blockid]; rkaux = dsdpSolver->rkaux[blockid];
        Smaker = dsdpSolver->Scker[blockid]; bnmaker = dsdpSolver->dS[blockid];
        dim = bnmaker->dim; Xtmp = (double *) calloc(dim * dim, sizeof(double));
        spsMatFactorize(Smaker); spsMatGetX(Smaker, bnmaker, Xtmp);
        // X = mu * D * (speye(n) + D' * bnmaker * D) * D';
        for (i = 0, idx = 0; i < dim; ++i) {
            memcpy(&dsaux->array[idx], &Xtmp[i * dim + i], sizeof(double) * (dim - i));
            idx += dim - i;
        }
        for (i = 0; i < nsym(dim); ++i) { dsaux->array[i] *= mumaker; }
        DSDP_FREE(Xtmp);
    }
    // LP Cone
    double *lpbnmaker = dsdpSolver->ds->x, *lps = dsdpSolver->scker->x, tmp = 0.0;
    vec_reset(dsdpSolver->x);
    for (DSDP_INT i = 0; i < dsdpSolver->lpDim; ++i) {
        tmp = 1 / lps[i]; dsdpSolver->x->x[i] = tmp + tmp * tmp * lpbnmaker[i];
    }
    vec_scale(dsdpSolver->x, mumaker);
    DSDPStatUpdate(&dsdpSolver->dsdpStats, STAT_GET_X_TIME, my_clock() - start);
    return retcode;
}

extern DSDP_INT computeDIMACS( HSDSolver *dsdpSolver ) {
    
    /* Compute the DIMACS error after solution is setup */
    DSDP_INT retcode = DSDP_RETCODE_OK;
    
    DSDPStats *stat = &dsdpSolver->dsdpStats;
    double bnrm, Cnrm, pObj, dObj, dInf, pInf,
           gap, tmp, minEigX, minEigS, compslack;
    DSDP_INT nblock = dsdpSolver->nBlock;
    
    DSDPGetStats(stat, STAT_ONE_NORM_C, &Cnrm);
    DSDPGetStats(stat, STAT_ONE_NORM_B, &bnrm);
    
    vec_dot(dsdpSolver->dObj, dsdpSolver->y, &dObj);
    dObj *= dsdpSolver->cScaler;
    
    /*  DIMACS Error 1    */
    pInf = computePres(dsdpSolver);
    
    if (pInf / (1 + bnrm) > 1e-02 && dsdpSolver->mumaker > 0) {
        printf("| Bad primal solution. Trying backup Newton step. \n");
        dsdpSolver->mumaker = dsdpSolver->mumaker2;
        vec_copy(dsdpSolver->ymaker2, dsdpSolver->ymaker);
        vec_copy(dsdpSolver->dymaker2, dsdpSolver->dymaker);
        retcode = computePrimalX(dsdpSolver);
        dsdpSolver->mumaker = -1.0; computeDIMACS(dsdpSolver);
        return retcode;
    }
    
    /*  DIMACS Error 3    */
    dInf = sqrt(dsdpSolver->n) * dsdpSolver->dperturb;
    pObj = DSDPConic( COPS_GET_CX )(dsdpSolver); pObj *= dsdpSolver->cScaler;
    
    /* DIMACS Error 4     */
    minEigX = DSDP_INFINITY;
    for (DSDP_INT i = 0; i < dsdpSolver->nBlock; ++i) {
        denseMatMinEig(dsdpSolver->dsaux[i], &tmp);
        minEigX = MIN(tmp, minEigX);
    }

    /*  DIMACS Error 5    */
    gap = pObj - dObj;
    if ((fabs(gap) / (1 + fabs(pObj) + fabs(dObj)) > 1e-02) && dsdpSolver->mumaker > 0) {
        printf("| Bad gap. Trying backup Newton step. \n");
        dsdpSolver->mumaker = dsdpSolver->mumaker2;
        vec_copy(dsdpSolver->ymaker2, dsdpSolver->ymaker);
        vec_copy(dsdpSolver->dymaker2, dsdpSolver->dymaker);
        retcode = computePrimalX(dsdpSolver);
        dsdpSolver->mumaker = -1.0; computeDIMACS(dsdpSolver);
        return retcode;
    }
    
    /*  DIMACS Error 6    */
    compslack = 0.0;
    for (DSDP_INT i = 0; i < nblock; ++i) {
        denseSpsTrace(dsdpSolver->dsaux[i], dsdpSolver->S[i], &tmp);
        compslack += tmp;
    }
    compslack *= dsdpSolver->cScaler;
    
    /*  DIMACS Error 2    */
    minEigS = DSDP_INFINITY;
    for (DSDP_INT i = 0; i < dsdpSolver->nBlock; ++i) {
        denseMatReset(dsdpSolver->dsaux[i]);
        spsMatFillLower2(dsdpSolver->S[i], dsdpSolver->dsaux[i]);
        denseMatMinEig(dsdpSolver->dsaux[i], &tmp);
        minEigS = MIN(minEigS, tmp);
    }
    
    if (minEigS > 0) {
        dInf -= minEigS * sqrt(dsdpSolver->n); dInf = MAX(dInf, 0.0);
    }
    
    // Collect errors
    DSDPStatUpdate(stat, STAT_DIMACS_ERR1, pInf / (1 + bnrm));
    DSDPStatUpdate(stat, STAT_DIMACS_ERR2, MAX(0.0, -minEigX) / (1 + bnrm));
    DSDPStatUpdate(stat, STAT_DIMACS_ERR3, dInf * dsdpSolver->cScaler / (1 + Cnrm));
    DSDPStatUpdate(stat, STAT_DIMACS_ERR4, MAX(0.0, -minEigS) * dsdpSolver->cScaler / (1 + Cnrm));
    DSDPStatUpdate(stat, STAT_DIMACS_ERR5, gap / (1 + fabs(pObj) + fabs(dObj)));
    DSDPStatUpdate(stat, STAT_DIMACS_ERR6, compslack / (1 + fabs(pObj) + fabs(dObj)));
    
    printf("| Final pObj: %10.5e   dObj: %10.5e \n", pObj, dObj);
    dsdpshowdash();
                   
    return retcode;
}
