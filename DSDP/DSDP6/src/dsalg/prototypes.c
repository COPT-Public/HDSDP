#include <stdio.h>
#include "dsdphsd.h"
#include "dsdputils.h"
#include "dsdpsolver.h"

static char etype[] = "Prototype methods";
/* Recording prototype methods that are no longer in use */

extern DSDP_INT dsdpGetAlphaLS( spsMat *S, spsMat *dS, spsMat *Scker,
                                double alphamax, double *alpha, DSDP_INT *sumHash ) {
    // Get the maximum alpha such that S + alpha * dS is PSD by line-search
    DSDP_INT retcode = DSDP_RETCODE_OK;
    
    double step = 1 / alphamax, *src = NULL;
    DSDP_INT ispsd = FALSE;
    spsMat *buffer = NULL;
    
    if (Scker) {
        src = Scker->x; buffer = Scker;
    } else {
        assert( FALSE );
        src = (double *) calloc(S->nnz, sizeof(double));
        memcpy(src, S->x, sizeof(double) * S->nnz); buffer = S;
    }
    
    for (DSDP_INT i = 0; ; ++i) {
        if (step <= 1e-05) {
            *alpha = 0.0; break;
        }
        memcpy(src, S->x, sizeof(double) * S->nnz);
        spsMataXpbY(step, dS, 1.0, buffer, sumHash); spsMatIspd(buffer, &ispsd);
        if (ispsd) { *alpha = step; break; }
        step *= 0.8;
    }

    if (!Scker) { DSDP_FREE(src); }
    return retcode;
}

static DSDP_INT packFactorize( dsMat *S ) {
    
    /* Factorize the dsMat matrix */
    DSDP_INT retcode = DSDP_RETCODE_OK;
    if (S->isFactorized) {
        error(etype, "Matrix is already factorized. \n");
    }
    
    DSDP_INT n = S->dim, info = 0; char uplo = DSDP_MAT_LOW;
    memcpy(S->lfactor, S->array, sizeof(double) * nsym(n));
    
    if (!S->isillCond) {
        packchol(&uplo, &n, S->lfactor, &info);
        if (info > 0) {
            S->isillCond = TRUE;
            packldl(&uplo, &n, S->lfactor, S->ipiv, &info);
        }
    } else {
        packldl(&uplo, &n, S->lfactor, S->ipiv, &info);
    }
    
    if (info < 0) {
        error(etype, "Illegal value detected in packed dense format. \n");
    }
    
    S->isFactorized = TRUE;
    return retcode;
}

static DSDP_INT packSolve( dsMat *S, DSDP_INT nrhs, double *B, double *X ) {
    /* Solve the linear system S * X = B using Lapack packed format */
    DSDP_INT retcode = DSDP_RETCODE_OK;
    char uplo = DSDP_MAT_LOW;
    DSDP_INT n = S->dim, ldb = S->dim, info = 0;
    
    // Copy solution data
    memcpy(X, B, sizeof(double) * nrhs * n);
    if (S->isillCond) {
        ldlsolve(&uplo, &n, &nrhs, S->lfactor, S->ipiv, X, &ldb, &info);
    } else {
        packsolve(&uplo, &n, &nrhs, S->lfactor, X, &ldb, &info);
    }
    if (info < 0) {
        error(etype, "Packed linear system solution failed. \n");
    }
    return retcode;
}

static DSDP_INT isConstant( HSDSolver *dsdpSolver, DSDP_INT *isConstant ) {
    // Check whether C is a constant
    DSDP_INT retcode = DSDP_RETCODE_OK;
    *isConstant = FALSE;
    
    DSDP_INT m = dsdpSolver->m, dim = 0;
    DSDP_INT isCons = TRUE, isRank1 = FALSE;
    r1Mat *C = NULL; double num = 0.0;
    
    for (DSDP_INT i = 0; i < dsdpSolver->nBlock; ++i) {
        isCons = TRUE;
        if (dsdpSolver->sdpData[i]->types[m] == MAT_TYPE_RANKK) {
            rkMatisRank1(dsdpSolver->sdpData[i]->sdpData[m], &isRank1);
            if (isRank1) {
                C = (r1Mat *) dsdpSolver->sdpData[i]->sdpData[m];
                dim = C->dim; num = C->x[0];
                
                for (DSDP_INT j = 1; j < dim; ++j) {
                    if (fabs(num - C->x[j]) > 1e-03) {
                        isCons = FALSE; break;
                    }
                }
            } else {
                isCons = FALSE;
            }
        } else {
            isCons = FALSE;
        }
        
        if (isCons) {
            *isConstant = TRUE; break;
        }
    }
    
    return retcode;
}

/* Check if a dense matrix is rank-one
   Depreciated due to inaccuracy
 */
static DSDP_INT isDenseRank1InAcc( dsMat *dataMat, DSDP_INT *isRank1 ) {
    // Check if a dense matrix is rank-one
    // This version is fast but not accurate due to potential numerical error
    DSDP_INT retcode = DSDP_RETCODE_OK;
    double *A    = dataMat->array;
    DSDP_INT n   = dataMat->dim;
    DSDP_INT r1  = TRUE;
    DSDP_INT col = 0;
    
    double benchCol  = 0.0;
    double scaleCol  = 0.0;
    double benchCol2 = 0.0;
    double scaleCol2 = 0.0;
    double diff      = 0.0;
    
    // Get the first column that contains non-zero elements
    for (DSDP_INT i = 0; i < n; ++i) {
        col = i;
        if (packIdx(A, n, i, i) != 0) {
            break;
        }
    }
    
    assert( col != n - 1 ); // or it is a zero matrix
    
    // Check the scaling coefficient
    for (DSDP_INT i = col + 1; i < n; ++i) {
        scaleCol = packIdx(A, n, i, col);
        for (DSDP_INT j = col; j < n; ++j) {
            benchCol2 = packIdx(A, n, j, col);
            if (i <= j) {
                scaleCol2 = packIdx(A, n, j, i);
            } else {
                scaleCol2 = packIdx(A, n, i, j);
            }
            diff = benchCol * scaleCol2 - benchCol2 * scaleCol;
            if (fabs(diff) > 1e-04 * MAX(1, benchCol)) {
                r1 = FALSE;
                break;
            }
        }
    }
    
    *isRank1 = r1;
    
    return retcode;
}

/* Set up Schur matrix by simple loop
   Depreciated due to low speed
 */

static DSDP_INT setupSDPSchurBlock( HSDSolver *dsdpSolver, DSDP_INT blockid ) {
    /* Set up the schur matrix for a given block. This routine updates
     
     1. Schur matrix Msdp
     2. asinv
     3. u
     4. d4 ( = ASinvRyASinv by the event indicator)
    */
    
    DSDP_INT retcode = DSDP_RETCODE_OK;
    assert( blockid < dsdpSolver->nBlock );
    
    DSDP_INT m       = dsdpSolver->m;
    DSDP_INT mattype = MAT_TYPE_UNKNOWN;
    sdpMat *sdpData = dsdpSolver->sdpData[blockid];
    
    // Temporary storage array for SinvASinv
    rkMat *rkdata = dsdpSolver->rkaux[blockid];
    dsMat *dsdata = dsdpSolver->dsaux[blockid];
    
    DSDP_INT dim = 0;
    
    if (dsdpSolver->eventMonitor[EVENT_IN_PHASE_B]) {
        dim = m;
    } else {
        dim = m + 1;
    }
    
    void *data = NULL;
    DSDP_INT i, j;
    for (i = 0; i < dim; ++i) {
        // Compute SinvASinv
        mattype = sdpData->types[i];
        
        if (mattype == MAT_TYPE_ZERO) {
            continue;
        } else if (mattype == MAT_TYPE_RANKK) {
            retcode = getSinvASinv(dsdpSolver, blockid, i, rkdata);
            data = (void *) rkdata;
            mattype = MAT_TYPE_RANKK;
            checkCode;
        } else {
            retcode = getSinvASinv(dsdpSolver, blockid, i, dsdata);
            data = (void *) dsdata;
            mattype = MAT_TYPE_DENSE;
            checkCode;
        }
        
        for (j = 0; j <= i; ++j) {
            getTraceASinvASinv(dsdpSolver, blockid, j, mattype, i, data);
        }
    }
    
    if (dsdpSolver->eventMonitor[EVENT_IN_PHASE_B] &&
        !dsdpSolver->eventMonitor[EVENT_INVALID_GAP]) {
        
        double maxdiag = 0.0;
        for (DSDP_INT i = 0; i < m; ++i) {
            maxdiag = MAX(packIdx(dsdpSolver->Msdp->array, m, i, i), maxdiag);
        }
        
        dsdpSolver->Mscaler = maxdiag;
        
        for (DSDP_INT i = 0; i < m; ++i) {
            packIdx(dsdpSolver->Msdp->array, m, i, i) += \
            MIN(maxdiag * 1e-06, 1e-08);
        }
    }
    
    // Perturb in Phase A
//    if (dsdpSolver->eventMonitor[EVENT_IN_PHASE_A]) {
//        for (DSDP_INT i = 0; i < m; ++i) {
//            packIdx(dsdpSolver->Msdp->array, m, i, i) += \
//            MIN(maxdiag * 1e-05, 1e-05);
//        }
//    }
    
    assert( retcode == DSDP_RETCODE_OK );
    return retcode;
}

/* Evaluate a single element of the Schur matrix
   Depreciated due to low speed
 */
extern DSDP_INT getTraceASinvASinvSlow( HSDSolver *dsdpSolver, DSDP_INT blockid, DSDP_INT constrid,
                                    DSDP_INT mattype, DSDP_INT constrid2, void *SinvASinv ) {
    
    // Compute trace between SinvASinv and some A
    // constrid is the position of (A in A * (.)) and constrid2 is the position of (A in (.) * SinvASinv)
    DSDP_INT retcode = DSDP_RETCODE_OK;
    DSDP_INT Atype = dsdpSolver->sdpData[blockid]->types[constrid];
    DSDP_INT m = dsdpSolver->m;
    void *A = dsdpSolver->sdpData[blockid]->sdpData[constrid];
    double trace = 0.0;
    double *M = dsdpSolver->Msdp->array;
    
    // assert( constrid <= constrid2 );
    // assert( constrid <= m && constrid2 <= m);
    
    if (mattype == MAT_TYPE_ZERO || Atype == MAT_TYPE_ZERO) {
        return retcode;
    }
    
    retcode = rkMatrkTrace((rkMat *) SinvASinv, (rkMat *) A, &trace);
    
    if (mattype == MAT_TYPE_RANKK) {
        switch (Atype) {
            case MAT_TYPE_DENSE:
                retcode = rkMatdenseTrace((rkMat *) SinvASinv, (dsMat *) A, &trace);
                // retcode = r1MatdenseTrace((r1Mat *) SinvASinv, (dsMat *) A, &trace);
                break;
            case MAT_TYPE_SPARSE:
                retcode = rkMatspsTrace((rkMat *) SinvASinv, (spsMat *) A, &trace);
                // retcode = r1MatspsTrace((r1Mat *) SinvASinv, (spsMat *) A, &trace);
                break;
            case MAT_TYPE_RANKK:
                retcode = rkMatrkTrace((rkMat *) SinvASinv, (rkMat *) A, &trace);
                // retcode = r1Matr1Trace((r1Mat *) SinvASinv, (r1Mat *) A, &trace);
                break;
            default:
                error(etype, "Invalid matrix type. \n");
                break;
        }
    } else if (mattype == MAT_TYPE_DENSE) {
        switch (Atype) {
            case MAT_TYPE_DENSE:
                retcode = denseDsTrace((dsMat *) SinvASinv, (dsMat *) A, &trace);
                break;
            case MAT_TYPE_SPARSE:
                retcode = denseSpsTrace((dsMat *) SinvASinv, (spsMat *) A, &trace);
                break;
            case MAT_TYPE_RANKK:
                retcode = rkMatdenseTrace((rkMat *) A, (dsMat *) SinvASinv, &trace);
                // retcode = r1MatdenseTrace((r1Mat *) A, (dsMat *) SinvASinv, &trace);
                break;
            default:
                error(etype, "Invalid matrix type. \n");
                break;
        }
    } else {
        error(etype, "Invalid matrix type. \n");
    }
    
    checkCode;
    
    // Perturbation of M diagonal element
    if (constrid == constrid2 &&
        dsdpSolver->eventMonitor[EVENT_IN_PHASE_B]) {
        if (dsdpSolver->mu < 1) {
            dsdpSolver->Msdp->isillCond = TRUE;
        }
        
        if (dsdpSolver->mu < 1e-05) {
            trace += 1e-08;
        }
    }
    
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

static DSDP_INT setupLPSchur( HSDSolver *dsdpSolver ) {
    // Set up the schur matrix for LP
    // After calling this routine, Msdp, asinv, u, csinv and d3 will be updated
    DSDP_INT retcode = DSDP_RETCODE_OK;
    
    assert( !dsdpSolver->iterProgress[ITER_SCHUR] );
    if (!dsdpSolver->iterProgress[ITER_SCHUR]) {
        error(etype, "Schur matrix for LP is already setup. \n");
    }
    
    DSDP_INT m = dsdpSolver->m;
    DSDP_INT n = dsdpSolver->lpDim;
    
    vec *c  = dsdpSolver->lpObj, *ry = dsdpSolver->ry, *u  = dsdpSolver->u, *d3 = dsdpSolver->d3;
    
    cs *A = dsdpSolver->lpData->lpdata;
    double *M = dsdpSolver->Msdp->array;
        
    // Setup the LP schur matrix AD^-2AT
    cs *AT = NULL;
    AT = cs_transpose(A, 0);
    
    DSDP_INT *ATp = AT->p;
    DSDP_INT *ATi = AT->i;
    double *ATx = AT->x;
    vec *s = dsdpSolver->s;
    double *sdata = s->x;
    double *cdata = c->x;
    double *rydata = ry->x;
    double sk = 0.0;
    
    double *Arow = NULL;
    Arow = (double *) calloc(n, sizeof(double));
    
    double Mij = 0.0;
    
    // Only compute the lower triangular
    for (DSDP_INT i = 0; i < m; ++i) {
        Mij = 0.0;
        memset(Arow, 0, sizeof(double) * n);
        cs_scatter2(AT, i, Arow);
        for (DSDP_INT j = 0; j <= i; ++j) {
            for (DSDP_INT k = ATp[j]; k < ATp[j + 1]; ++k) {
                sk = sdata[ATi[k]];
                Mij += ATx[k] * Arow[ATi[k]] * sk * sk;
            }
            packIdx(M, m, i, j) += Mij;
        }
    }
    
    DSDP_FREE(Arow);
    cs_spfree(AT);
    
    // Update asinv
    vec *sinv = (vec *) calloc(1, sizeof(vec));
    retcode = vec_init(sinv);
    retcode = vec_inv(sinv, s);
    cs_gaxpy(A, sinv->x, dsdpSolver->asinv->x);
    
    // Update csinv
    DSDP_INT one = 1;
    dsdpSolver->csinv += dot(&n, cdata, &one, sinv->x, &one);
    
    // Update u
    retcode = vec_invsqr(sinv, s);
    
    register DSDP_INT i;
    
    for (i = 0; i < n; ++i) {
        sdata[i] = sdata[i] * cdata[i];
    }
    
    cs_gaxpy(A, sdata, u->x);
    
    // Update csinvrycsinv
    dsdpSolver->csinvrysinv += dot(&n, sdata, &one, ry->x, &one);
    
    // Update d3
    retcode = vec_invsqr(sinv, s);
    
    for (i = 0; i < n; ++i) {
        sdata[i] = sdata[i] * rydata[i];
    }
    
    cs_gaxpy(A, sdata, d3->x);
    
    // Free allocated memory
    retcode = vec_free(sinv);
    DSDP_FREE(sinv);
    
    dsdpSolver->iterProgress[ITER_SCHUR] = TRUE;
    
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

static DSDP_INT getLPsStep( HSDSolver *dsdpSolver, double *sStep ) {
    // Compute the maixmum step size to take at s for LP
    DSDP_INT retcode = DSDP_RETCODE_OK;
    
    vec *s  = dsdpSolver->s, *ds = dsdpSolver->ds;
    
    double step = 100.0;
    
    for (DSDP_INT i = 0; i < s->dim; ++i) {
        step = MIN((ds->x[i] / s->x[i]), step);
    }
    
    if (step < 0.0) {
        *sStep = fabs(0.995 / step);
    }
    
    return retcode;
}

static DSDP_INT takelpsStep( HSDSolver *dsdpSolver ) {
    // Take step in LP s
    double step = dsdpSolver->alpha;
    vec *s  = dsdpSolver->s, *ds = dsdpSolver->ds;
    vec_axpy(step, ds, s);
    return DSDP_RETCODE_OK;
}

extern DSDP_INT dsdpgetPhaseAProxMeasure( HSDSolver *dsdpSolver, double newmu ) {
    
    // Check primal feasibility and proximity measure given a new mu parameter
    DSDP_INT retcode = DSDP_RETCODE_OK;
    
    DSDP_INT ispfeas = FALSE;
    double taudenom = 0.0, dtaudelta = 0.0, tau = dsdpSolver->tau,
           dObj = dsdpSolver->dObjVal, utd2 = 0.0;
    
    // Use b1 as auxiliary
    vec *dydelta = dsdpSolver->b1;
    vec_dot(dsdpSolver->u, dsdpSolver->d12, &taudenom);
    
    taudenom = - taudenom + dsdpSolver->csinvcsinv + 1 / (tau * tau);
    
    if (fabs(taudenom) < 1e-15) {
        dtaudelta = 0.0;
    } else {
        // dtaudelta = - dObj + mu / tau + mu * csinv + tau * u' * d2 - mu * u' * d3;
        vec_dot(dsdpSolver->u, dsdpSolver->d2, &utd2);
        dtaudelta = - dObj + newmu * (1 / tau + dsdpSolver->csinv) + tau * utd2;
        vec_dot(dsdpSolver->u, dsdpSolver->d3, &utd2);
        dtaudelta -= utd2 * newmu;
        
        // dtaudelta = (1 / mu) *  dtaudelta / taudenom;
        dtaudelta = dtaudelta / (newmu * taudenom);
    }
    
    // dydelta = d12 * dtaudelta + d2 * tau / mu - d3;
    vec_zaxpby(dydelta, dtaudelta, dsdpSolver->d12, tau / newmu, dsdpSolver->d2);
    vec_axpy(-1.0, dsdpSolver->d3, dydelta);
    vec_dot(dsdpSolver->u, dydelta, &utd2);
    
    dsdpSolver->Pnrm = (dsdpSolver->csinvcsinv + 1 / (tau * tau)) * (dtaudelta * dtaudelta);
    dsdpSolver->Pnrm -= 2 * utd2 * dtaudelta;
    utd2 = denseMatxTAx(dsdpSolver->Msdp, dsdpSolver->M->schurAux, dydelta->x);
    dsdpSolver->Pnrm = sqrt(dsdpSolver->Pnrm + dsdpSolver->Mscaler * utd2);
    
    retcode = dsdpCheckPhaseAPfeas(dsdpSolver, dtaudelta, dydelta, &ispfeas);
    
    if (ispfeas) {
        dsdpSolver->eventMonitor[EVENT_PFEAS_FOUND] = TRUE;
        utd2 = newmu / (tau * tau) * (tau - dtaudelta);
        dsdpSolver->pObjVal = dObj - utd2;
        dsdpSolver->pObjVal = dsdpSolver->pObjVal / dsdpSolver->tau;
    }
    
    return retcode;
}

static void SDPConic( COPS_CONSTR_EXPR_OLD )
( HSDSolver *dsdpSolver, DSDP_INT type,
  double ycoef, vec *y, double tau, double r ) {
    // r * Ry + tau * C + ycoef * ATy
    spsMat **targets = (type == DUALVAR) ? \
    dsdpSolver->S : dsdpSolver->Scker;
    void (*f)(HSDSolver*, DSDP_INT, DSDP_INT, double) = \
    (type == DUALVAR) ? addMattoS : addMattoChecker;
    double dperturb = dsdpSolver->dperturb;
    if (type == DELTAS) {
        targets = dsdpSolver->dS;
        f = addMattodS;
        dperturb = 0.0;
    }
    
    dperturb += r * dsdpSolver->Ry;
    spsMat *target = NULL; DSDP_INT m = dsdpSolver->m, i, j;
    
    for (i = 0; i < dsdpSolver->nBlock; ++i) {
        target = targets[i]; spsMatReset(target);
        if (ycoef) {
            for (j = 0; j < m; ++j) {
                f(dsdpSolver, i, j, ycoef * y->x[j]);
            }
        }
        if (tau) { f(dsdpSolver, i, m, tau); }
        if (dperturb)   {spsMatAdddiag(target, dperturb, dsdpSolver->symS[i]);}
    }
    
    return;
}

/*
extern DSDP_INT getSinvASinvSlow( HSDSolver *dsdpSolver, DSDP_INT blockid, DSDP_INT constrid,
                              void *SinvASinv ) {
    
    // Given S and A, the routine computes A, asinv and trace(S, Sinv A Sinv)
    DSDP_INT retcode, typeA = dsdpSolver->sdpData[blockid]->types[constrid];
    spsMat *S = dsdpSolver->S[blockid];
    void *A = dsdpSolver->sdpData[blockid]->sdpData[constrid];
    double *asinv = NULL, *asinvrysinv = NULL, tracediag = 0.0;
    
    if (constrid < dsdpSolver->m) {
        asinv = &dsdpSolver->asinv->x[constrid];
        asinvrysinv = &dsdpSolver->d4->x[constrid];
    } else {
        asinv = &dsdpSolver->csinv;
        asinvrysinv = &dsdpSolver->csinvrysinv;
    }
    
    switch (typeA) {
        case MAT_TYPE_RANKK:
            retcode = spsSinvRkSinvSolve(S, A, SinvASinv, &tracediag);
            break;
        case MAT_TYPE_SPARSE:
            retcode = spsSinvSpSinvSolve(S, A, SinvASinv, &tracediag);
            break;
        case MAT_TYPE_DENSE:
            retcode = spsSinvDsSinvSolve(S, A, SinvASinv, &tracediag);
            break;
        default:
            error(etype, "Invalid matrix type. \n");
            break;
    }
    
    *asinv += tracediag;
    
    if (dsdpSolver->eventMonitor[EVENT_IN_PHASE_A]) {
        switch (typeA) {
            case MAT_TYPE_RANKK:
                retcode = rkMatdiagTrace(SinvASinv, dsdpSolver->Ry, &tracediag);
                break;
            default:
                retcode = denseDiagTrace(SinvASinv, dsdpSolver->Ry, &tracediag);
                break;
        }
        *asinvrysinv += tracediag;
    }
    return retcode;
}
 
 */

/*
static DSDP_INT setupSDPSchurBlockA( HSDSolver *dsdpSolver, DSDP_INT blockid ) {
    
    DSDP_INT retcode = DSDP_RETCODE_OK;
    assert( blockid < dsdpSolver->nBlock );
    
    DSDP_INT mattype = MAT_TYPE_UNKNOWN;
    sdpMat *sdpData = dsdpSolver->sdpData[blockid];
    DSDP_INT m = dsdpSolver->m, n = sdpData->dimS, useM1 = TRUE;
    
    // Temporary storage array for SinvASinv
    rkMat *rkaux = dsdpSolver->rkaux[blockid], *rkdata = NULL;
    dsMat *dsaux = dsdpSolver->dsaux[blockid];
        
    double coeff = 0.0, res = 0.0, *M = dsdpSolver->Msdp->array, *Ax = NULL;
    DSDP_INT i, j, r, rank;
    
    for (i = 0; i < m; ++i) {
        // Compute SinvASinv
        mattype = sdpData->types[i];
        
        if (mattype == MAT_TYPE_ZERO) {
            continue;
        } else {
            retcode = getSinvASinv(dsdpSolver, blockid, i, rkaux);
            mattype = MAT_TYPE_RANKK;
            checkCode;
        }
        
        useM1 = TRUE;
        rank = rkaux->rank;
        
        if (rank <= M1Threshold * n) {
            useM1 = FALSE;
        }
        
        if (useM1) {
            // M1 Technique
            denseMatReset(dsaux);
            rkMatdenseUpdate(dsaux, rkaux);
            
            for (j = 0; j <= i; ++j) {
                
                if (sdpData->types[j] == MAT_TYPE_ZERO) {
                    continue;
                }
                
                rkdata = sdpData->sdpData[j];
                switch (rkdata->mattype) {
                    case MAT_TYPE_SPARSE:
                        denseSpsTrace(dsaux, rkdata->origdata, &res);
                        break;
                    case MAT_TYPE_DENSE:
                        denseDsTrace(dsaux, rkdata->origdata, &res);
                        break;
                    case MAT_TYPE_RANKK:
                        r1MatdenseTrace(rkdata->data[0], dsaux, &res);
                        break;
                    default:
                        error(etype, "Invalid matrix type");
                        break;
                }
                packIdx(M, m, i, j) += res;
            }
            
            if (sdpData->types[m] == MAT_TYPE_ZERO) {
                continue;
            }
            
            rkdata = sdpData->sdpData[m];
            switch (rkdata->mattype) {
                case MAT_TYPE_SPARSE:
                    denseSpsTrace(dsaux, rkdata->origdata, &res);
                    break;
                case MAT_TYPE_DENSE:
                    denseDsTrace(dsaux, rkdata->origdata, &res);
                    break;
                case MAT_TYPE_RANKK:
                    r1MatdenseTrace(rkaux->data[0], dsaux, &res);
                    break;
                default:
                    error(etype, "Invalid matrix type");
                    break;
            }
            
            dsdpSolver->u->x[i] += res;
            
        } else {
            // M2 Technique
            for (r = 0; r < rank; ++r) {
                Ax = rkaux->data[r]->x;
                coeff = rkaux->data[r]->sign;
                
                for (j = 0; j <= i; ++j) {
                    
                    if (sdpData->types[j] == MAT_TYPE_ZERO) {
                        continue;
                    }
                    
                    rkdata = sdpData->sdpData[j];
                    switch (rkdata->mattype) {
                        case MAT_TYPE_SPARSE:
                            spsMatxTAx(rkdata->origdata, Ax, &res);
                            res = coeff * res;
                            break;
                        case MAT_TYPE_DENSE:
                            denseMatxTAx(rkdata->origdata, Ax, &res);
                            res = coeff * res;
                            break;
                        case MAT_TYPE_RANKK:
                            r1Matr1Trace(rkdata->data[0], rkaux->data[r], &res);
                            break;
                        default:
                            error(etype, "Invalid matrix type");
                            break;
                    }
                    packIdx(M, m, i, j) += res;
                }
                
                if (sdpData->types[m] == MAT_TYPE_ZERO) {
                    continue;
                }
                
                rkdata = sdpData->sdpData[m];
                switch (rkdata->mattype) {
                    case MAT_TYPE_SPARSE:
                        spsMatxTAx(rkdata->origdata, Ax, &res);
                        res *= coeff;
                        break;
                    case MAT_TYPE_DENSE:
                        denseMatxTAx(rkdata->origdata, Ax, &res);
                        res *= coeff;
                        break;
                    case MAT_TYPE_RANKK:
                        r1Matr1Trace(rkdata->data[0], rkaux->data[r], &res);
                        break;
                    default:
                        error(etype, "Invalid matrix type");
                        break;
                }
                
                dsdpSolver->u->x[i] += res;
            }
        }
    }
        
    // Compute csinvcsinv
    mattype = sdpData->types[m];
    
    if (mattype == MAT_TYPE_ZERO) {
        return retcode;
    }
    
    retcode = getSinvASinv(dsdpSolver, blockid, m, rkaux);
    rank = rkaux->rank;
    for (r = 0; r < rank; ++r) {
        Ax = rkaux->data[r]->x;
        coeff = rkaux->data[r]->sign;
        rkdata = sdpData->sdpData[m];
        switch (rkdata->mattype) {
            case MAT_TYPE_SPARSE:
                spsMatxTAx(rkdata->origdata, Ax, &res);
                res *= coeff;
                break;
            case MAT_TYPE_DENSE:
                denseMatxTAx(rkdata->origdata, Ax, &res);
                res *= coeff;
                break;
            case MAT_TYPE_RANKK:
                r1Matr1Trace(rkdata->data[0], rkaux->data[r], &res);
                break;
            default:
                error(etype, "Invalid matrix type");
                break;
        }
        dsdpSolver->csinvcsinv += res;
    }
    
    assert( retcode == DSDP_RETCODE_OK );
    return retcode;
}

static DSDP_INT setupSDPSchurBlockB( HSDSolver *dsdpSolver, DSDP_INT blockid ) {
    
    DSDP_INT retcode = DSDP_RETCODE_OK;
    assert( blockid < dsdpSolver->nBlock );
    
    DSDP_INT mattype = MAT_TYPE_UNKNOWN;
    sdpMat *sdpData = dsdpSolver->sdpData[blockid];
    DSDP_INT m = dsdpSolver->m, n = sdpData->dimS, useM1 = TRUE;
    
    // Temporary storage array for SinvASinv
    rkMat *rkaux = dsdpSolver->rkaux[blockid], *rkdata = NULL;
    dsMat *dsaux = dsdpSolver->dsaux[blockid];
    
    DSDP_INT shift = 0;
    
    if (dsdpSolver->eventMonitor[EVENT_IN_PHASE_B]) {
        shift = 0;
    } else {
        shift = 1;
    }
    
    void *data = NULL;
    double coeff = 0.0, res = 0.0, *M = dsdpSolver->Msdp->array, *Ax = NULL;
    DSDP_INT i, j, r, rank;
    
    for (i = 0; i < m; ++i) {
        // Compute SinvASinv
        mattype = sdpData->types[i];
        
        if (mattype == MAT_TYPE_ZERO) {
            continue;
        } else if (mattype == MAT_TYPE_RANKK) {
            retcode = getSinvASinv(dsdpSolver, blockid, i, rkaux);
            data = (void *) rkaux;
            mattype = MAT_TYPE_RANKK;
            checkCode;
        }
        
        useM1 = TRUE;
        rank = rkaux->rank;
        
        // Add more choice heuristics
        if ( rank <= M1Threshold * n ) {
            useM1 = FALSE;
        }
        
        if (rank >= M1Threshold * n) {
            // M1 Technique
            denseMatReset(dsaux);
            rkMatdenseUpdate(dsaux, rkaux);
            
            for (j = 0; j <= i; ++j) {
                
                if (sdpData->types[j] == MAT_TYPE_ZERO) {
                    continue;
                }
                
                rkdata = sdpData->sdpData[j];
                
                switch (rkdata->mattype) {
                    case MAT_TYPE_SPARSE:
                        denseSpsTrace(dsaux, rkdata->origdata, &res);
                        break;
                    case MAT_TYPE_DENSE:
                        denseDsTrace(dsaux, rkdata->origdata, &res);
                        break;
                    case MAT_TYPE_RANKK:
                        r1MatdenseTrace(rkdata->data[0], dsaux, &res);
                        break;
                    default:
                        error(etype, "Invalid matrix type");
                        break;
                }
                packIdx(M, m, i, j) += res;
            }
            
        } else {
            // M2 Technique
            for (r = 0; r < rank; ++r) {
                Ax = rkaux->data[r]->x;
                coeff = rkaux->data[r]->sign;
                
                for (j = 0; j <= i; ++j) {
                    
                    if (sdpData->types[j] == MAT_TYPE_ZERO) {
                        continue;
                    }
                    
                    rkdata = sdpData->sdpData[j];
                    switch (rkdata->mattype) {
                        case MAT_TYPE_SPARSE:
                            spsMatxTAx(rkdata->origdata, Ax, &res);
                            res *= coeff;
                            break;
                        case MAT_TYPE_DENSE:
                            denseMatxTAx(rkdata->origdata, Ax, &res);
                            res *= coeff;
                            break;
                        case MAT_TYPE_RANKK:
                            r1Matr1Trace(rkdata->data[0], rkaux->data[r], &res);
                            break;
                        default:
                            error(etype, "Invalid matrix type. \n");
                            break;
                    }
                    packIdx(M, m, i, j) += res;
                }
            }
        }
    }
    
    assert( retcode == DSDP_RETCODE_OK );
    return retcode;
}
*/
