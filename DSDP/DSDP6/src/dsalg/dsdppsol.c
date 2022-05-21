#include "dsdppsol.h"
#include "dsdputils.h"
#include "dsdplog.h"

static char etype[] = "Primal solution extraction";

extern DSDP_INT computePrimalX( HSDSolver *dsdpSolver ) {
    
    // Extract primal solution of the current block
    DSDP_INT retcode = DSDP_RETCODE_OK;
        
    DSDP_INT nblock = dsdpSolver->nBlock, dim = 0;
    vec *dymaker = dsdpSolver->dymaker, *ymaker = dsdpSolver->ymaker;
    spsMat *Smaker = NULL, *bnmaker = NULL;
    double mumaker = dsdpSolver->mumaker;
    double *Xtmp = NULL;
    
    double start = my_clock();
    
    // getPhaseBS(dsdpSolver, dsdpSolver->y->x);
    // Smaker = C - dsdpgetATy(A, ymaker);
    getPhaseBCheckerS(dsdpSolver, ymaker);
    // bnmaker = dsdpgetATy(A, dymaker);
    getPhaseBdS(dsdpSolver, -1.0, dymaker, 0.0);
    
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
    
    DSDPStatUpdate(&dsdpSolver->dsdpStats, STAT_GET_X_TIME,
                   my_clock() - start);
    
    return retcode;
}

extern DSDP_INT computeDIMACS( HSDSolver *dsdpSolver ) {
    
    /* Compute the DIMACS error after solution is setup */
    DSDP_INT retcode = DSDP_RETCODE_OK;
    
    DSDPStats *stat = &dsdpSolver->dsdpStats;
    double bnrm, Cnrm, pObj, dObj, dInf, pInf,
           gap, trace, tmp, minEigX, minEigS, compslack;
    DSDP_INT nblock = dsdpSolver->nBlock, m = dsdpSolver->m;
    
    DSDPGetStats(stat, STAT_ONE_NORM_C, &Cnrm);
    DSDPGetStats(stat, STAT_ONE_NORM_B, &bnrm);
    
    vec_dot(dsdpSolver->dObj, dsdpSolver->y, &dObj);
    dObj *= dsdpSolver->cScaler;
    
    /*  DIMACS Error 1    */
    pInf = 0.0;
    for (DSDP_INT i = 0; i < m; ++i) {
        trace = 0.0;
        for (DSDP_INT j = 0; j < nblock; ++j) {
            switch (dsdpSolver->sdpData[j]->types[i]) {
                case MAT_TYPE_ZERO: tmp = 0.0; break;
                case MAT_TYPE_DENSE:
                    denseDsTrace(dsdpSolver->dsaux[j],
                                 dsdpSolver->sdpData[j]->sdpData[i],
                                 &tmp);
                    break;
                case MAT_TYPE_SPARSE:
                    denseSpsTrace(dsdpSolver->dsaux[j],
                                  dsdpSolver->sdpData[j]->sdpData[i],
                                  &tmp);
                    break;
                case MAT_TYPE_RANKK:
                    rkMatdenseTrace(dsdpSolver->sdpData[j]->sdpData[i],
                                    dsdpSolver->dsaux[j], &tmp);
                    break;
                default:
                    error(etype, "Invalid matrix type. \n");
                    break;
            }
            trace += tmp;
        }
        tmp = dsdpSolver->dObj->x[i] - trace;
        pInf += tmp * tmp;
    }
    
    pInf = sqrt(pInf);
    
    if (pInf / (1 + bnrm) > 1e-02 && dsdpSolver->mumaker > 0) {
        printf("| Bad primal solution. Trying backup Newton step. \n");
        dsdpSolver->mumaker = dsdpSolver->mumaker2;
        vec_copy(dsdpSolver->ymaker2, dsdpSolver->ymaker);
        vec_copy(dsdpSolver->dymaker2, dsdpSolver->dymaker);
        retcode = computePrimalX(dsdpSolver);
        dsdpSolver->mumaker = -1.0;
        computeDIMACS(dsdpSolver);
        return retcode;
    }
    
    /*  DIMACS Error 3    */
    dInf = sqrt(dsdpSolver->n) * dsdpSolver->dperturb;
    
    pObj = 0.0; tmp = 0.0;
    for (DSDP_INT i = 0; i < nblock; ++i) {
        switch (dsdpSolver->sdpData[i]->types[m]) {
            case MAT_TYPE_ZERO: tmp = 0.0; break;
            case MAT_TYPE_DENSE:
                denseDsTrace(dsdpSolver->dsaux[i],
                             dsdpSolver->sdpData[i]->sdpData[m],
                             &tmp);
                break;
            case MAT_TYPE_SPARSE:
                denseSpsTrace(dsdpSolver->dsaux[i],
                              dsdpSolver->sdpData[i]->sdpData[m],
                              &tmp);
                break;
            case MAT_TYPE_RANKK:
                rkMatdenseTrace(dsdpSolver->sdpData[i]->sdpData[m],
                                dsdpSolver->dsaux[i],
                                &tmp);
                break;
            default:
                break;
        }
        pObj += tmp;
    }
    
    pObj *= dsdpSolver->cScaler;
    
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
    compslack *= dsdpSolver->cScaler;
    
    /*  DIMACS Error 2    */
    minEigS = DSDP_INFINITY;
    for (DSDP_INT i = 0; i < dsdpSolver->nBlock; ++i) {
        denseMatReset(dsdpSolver->dsaux[i]);
        spsMatFillLower2(dsdpSolver->S[i], dsdpSolver->dsaux[i]);
        denseMatMinEig(dsdpSolver->dsaux[i], &tmp);
        
        // spsMatMinEig(dsdpSolver->S[i], &tmp);
        minEigS = MIN(minEigS, tmp);
    }
    
    // Collect errors
    DSDPStatUpdate(stat, STAT_DIMACS_ERR1, pInf / (1 + bnrm));
    DSDPStatUpdate(stat, STAT_DIMACS_ERR2, MAX(0.0, -minEigX) / (1 + bnrm));
    DSDPStatUpdate(stat, STAT_DIMACS_ERR3, dInf * dsdpSolver->cScaler / (1 + Cnrm));
    DSDPStatUpdate(stat, STAT_DIMACS_ERR4, MAX(0.0, -minEigS) * dsdpSolver->cScaler / (1 + Cnrm));
    DSDPStatUpdate(stat, STAT_DIMACS_ERR5, gap / (1 + fabs(pObj) + fabs(dObj)));
    DSDPStatUpdate(stat, STAT_DIMACS_ERR6, compslack / (1 + fabs(pObj) + fabs(dObj)));
                   
    return retcode;
}
