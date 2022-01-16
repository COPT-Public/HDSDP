#include "schurmat.h"

static char etype[] = "Schur matrix setup";
/* Setup the auxiliary arrays as well as the schur matrix
 
 In DSDP, we set up the schur matrix by
 
 
 for k = 1 : nblock
 
    Sk = S[k]
    Ak = A[k]
 
    for i = 1 : m + 1
 
        X = Skinv * Ak[i] * Skinv
        asinv[k][i] = trace(A[i] * Skinv)
 
        for j = 1 : i
            M[i, j] = trace(X * A[j])
 
        u[i] = trace(X * A[m + 1])
     
*/

// A helping routine to build SinvASinv
static DSDP_INT SinvASinv( spsMat *S, DSDP_INT typeA, void *A, double *asinv, void *SinvASinv ) {
    
    // Given S and A, the routine computes A, asinv and trace(S, Sinv A Sinv)
    DSDP_INT retcode = DSDP_RETCODE_OK;
    
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
    
    return retcode;
}

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
    DSDP_INT useRy   = TRUE;
    double trace     = 0.0;
    
    sdpMat *sdpData = dsdpSolver->sdpData[blockid];
    spsMat *S       = dsdpSolver->S[blockid];
    vec    *asinv   = dsdpSolver->asinv;
    double Ry       = dsdpSolver->Ry;
    
    dsMat  *Msdp    = dsdpSolver->Msdp;
    vec    *u       = dsdpSolver->u;
    vec    *d4      = dsdpSolver->d4;
    
    if (dsdpSolver->eventMonitor[EVENT_SDP_NO_RY]) {
        useRy = FALSE;
    }
    
    void   **blockdata = NULL;
    double *Mdata      = Msdp->array;
    double csinv       = dsdpSolver->csinv;
    double csinvcsinv  = dsdpSolver->csinvcsinv;
    double csinvrysinv = dsdpSolver->csinvrysinv;
    
    // Temporary storage array for SinvASinv
    r1Mat *r1data = (r1Mat *) calloc(1, sizeof(r1Mat));
    dsMat *dsdata = (dsMat *) calloc(1, sizeof(dsMat));
    
    for (DSDP_INT i = 0; i < m + 1; ++i) {
        
        // Compute SinvASinv
        mattype = sdpData->types[i];
        if (mattype == MAT_TYPE_ZERO) {
            continue;
        } else if (mattype == MAT_TYPE_RANK1) {
            retcode = SinvASinv(S, mattype, sdpData->sdpData[i], &trace, (void *) r1data);
            checkCodeFree;
            // retcode = r1MatCountNnz(r1data);
        } else {
            retcode = SinvASinv(S, mattype, sdpData->sdpData[i], &trace, (void *) dsdata);
            checkCodeFree;
        }
        
        if (i < m) {
            asinv->x[i] = trace;
        } else {
            csinv += trace;
        }
        
        if (mattype == MAT_TYPE_RANK1) {
            
            if (useRy) {
                
                // Setup d4: dual infeasibility
                retcode = r1MatdiagTrace(r1data, Ry, &trace);
                checkCodeFree;
                
                if (i < m) {
                    d4->x[i] += trace;
                } else {
                    csinvrysinv += trace;
                }
            }
            
            // Set up M
            for (DSDP_INT j = 0; j <= i; ++i) {
                switch (sdpData->types[i]) {
                    case MAT_TYPE_ZERO:
                        break;
                    case MAT_TYPE_RANK1:
                        retcode = r1Matr1Trace(r1data, (r1Mat *) blockdata[j], &trace);
                        checkCodeFree;
                        break;
                    case MAT_TYPE_SPARSE:
                        retcode = r1MatspsTrace(r1data, (spsMat *) blockdata[j], &trace);
                        checkCodeFree;
                        break;
                    case MAT_TYPE_DENSE:
                        retcode = r1MatdenseTrace(r1data, (dsMat *) blockdata[j], &trace);
                        checkCodeFree;
                        break;
                    default:
                        error(etype, "Unknown matrix type. \n");
                        break;
                }
                
                if (i < m) {
                    packIdx(Mdata, m, i, j) += trace;
                } else {
                    // Set up csinv
                    if (j == m) {
                        csinvcsinv += trace;
                    } else {
                        // Setup csinvasinv
                        u->x[j] += trace;
                    }
                }
            }
        } else {
            
            if (useRy) {
                // Set up d4: dual infeasibility
                retcode = denseDiagTrace(dsdata, Ry, &trace);
                checkCodeFree;
                if (i < m) {
                    d4->x[i] += trace;
                } else {
                    csinvrysinv += trace;
                }
            }
            
            // Set up M
            for (DSDP_INT j = 0; j <= i; ++i) {
                switch (sdpData->types[i]) {
                    case MAT_TYPE_ZERO:
                        break;
                    case MAT_TYPE_RANK1:
                        retcode = r1MatdenseTrace((r1Mat *) blockdata[j], dsdata, &trace);
                        checkCodeFree;
                        break;
                    case MAT_TYPE_SPARSE:
                        retcode = denseSpsTrace(dsdata, (spsMat *) blockdata[j], &trace);
                        checkCodeFree;
                        break;
                    case MAT_TYPE_DENSE:
                        retcode = denseDsTrace(dsdata, (dsMat *) blockdata[j], &trace);
                        checkCodeFree;
                        break;
                    default:
                        error(etype, "Unknown matrix type. \n");
                        break;
                }
                if (i < m) {
                    packIdx(Mdata, m, i, j) += trace;
                } else {
                    if (j == m) {
                        csinvcsinv += trace;
                    } else {
                        u->x[j] += trace;
                    }
                }
            }
        }
    }
    
    dsdpSolver->csinv = csinv;
    dsdpSolver->csinvcsinv = csinvcsinv;
    dsdpSolver->csinvrysinv = csinvrysinv;
    
clean_up:
    
    retcode = r1MatFree(r1data);
    retcode = denseMatFree(dsdata);
    return retcode;
}

static DSDP_INT setupSDPSchur( HSDSolver *dsdpSolver ) {
    // Set up the schur matrix for SDP
    // After calling this routine, Msdp, asinv, u for SDP and csinv will be filled
    DSDP_INT retcode = DSDP_RETCODE_OK;
    
    assert( !dsdpSolver->iterProgress[ITER_SDP_SCHUR] );
    if (!dsdpSolver->iterProgress[ITER_SDP_SCHUR]) {
        error(etype, "Schur matrix is already setup. \n");
    }
    
    DSDP_INT m = dsdpSolver->m;
    DSDP_INT nblock = dsdpSolver->nBlock;
    vec *u      = dsdpSolver->u;
    vec *d4     = dsdpSolver->d4;
    dsMat *Msdp = dsdpSolver->Msdp;
    
    // Clear the Schur matrix and other arrays
    memset(Msdp->array, 0, sizeof(double) * nsym(m));
    
    dsdpSolver->csinv = 0.0;
    dsdpSolver->csinvrysinv = 0.0;
    dsdpSolver->csinvcsinv  = 0.0;
    retcode = vec_reset(u);
    retcode = vec_reset(d4);
    
    // Start setting up the Schur matrix
    for (DSDP_INT k = 0; k < nblock; ++k) {
        retcode = setupSDPSchurBlock(dsdpSolver, k); checkCode;
    }
    
    dsdpSolver->iterProgress[ITER_SDP_SCHUR] = TRUE;
    return retcode;
}

static DSDP_INT setupLPSchur( HSDSolver *dsdpSolver ) {
    // Set up the schur matrix for LP
    // After calling this routine, Msdp, asinv, u, csinv and d3 will be updated
    DSDP_INT retcode = DSDP_RETCODE_OK;
    
    assert( !dsdpSolver->iterProgress[ITER_LP_SCHUR] );
    if (!dsdpSolver->iterProgress[ITER_LP_SCHUR]) {
        error(etype, "Schur matrix for LP is already setup. \n");
    }
    
    DSDP_INT m = dsdpSolver->m;
    DSDP_INT n = dsdpSolver->lpDim;
    
    vec *c  = dsdpSolver->lpObj;
    vec *ry = dsdpSolver->ry;
    vec *u  = dsdpSolver->u;
    vec *d3 = dsdpSolver->d3;
    
    cs *A = dsdpSolver->lpData->lpdata;
    double *M = dsdpSolver->Msdp->array;
    
    
    // TODO: Replace the SuiteSparse routines using MKL Sparse routine library
    
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
    for (DSDP_INT i = 0; i < n - n % 4; i+=4) {
        sdata[i    ] = sdata[i    ] * cdata[i    ];
        sdata[i + 1] = sdata[i + 1] * cdata[i + 1];
        sdata[i + 2] = sdata[i + 2] * cdata[i + 2];
        sdata[i + 3] = sdata[i + 3] * cdata[i + 3];
    }
    
    for (DSDP_INT i = n - n % 4; i < n; ++i) {
        sdata[i] = sdata[i] * cdata[i];
    }
    
    cs_gaxpy(A, sdata, u->x);
    
    // Update csinvrycsinv
    dsdpSolver->csinvrysinv += dot(&n, sdata, &one, ry->x, &one);
    
    // Update d3
    retcode = vec_invsqr(sinv, s);
    
    if (n > 64) {
        for (DSDP_INT i = 0; i < n - n % 4; i+=4) {
            sdata[i    ] = sdata[i    ] * rydata[i    ];
            sdata[i + 1] = sdata[i + 1] * rydata[i + 1];
            sdata[i + 2] = sdata[i + 2] * rydata[i + 2];
            sdata[i + 3] = sdata[i + 3] * rydata[i + 3];
        }
        
        for (DSDP_INT i = n - n % 4; i < n; ++i) {
            sdata[i] = sdata[i] * rydata[i];
        }
    } else {
        for (DSDP_INT i = 0; i < n; ++i) {
            sdata[i] = sdata[i] * rydata[i];
        }
    }
    
    cs_gaxpy(A, sdata, d3->x);
    
    // Free allocated memory
    retcode = vec_free(sinv);
    DSDP_FREE(sinv);
    
    dsdpSolver->iterProgress[ITER_LP_SCHUR] = TRUE;
    
    return retcode;
}

static DSDP_INT setupRM( HSDSolver *dsdpSolver, vec *RM ) {
    // Set up the auxiliary vector RM
    // RM = b * tau - mu * asinv + mu * d3
    DSDP_INT retcode = DSDP_RETCODE_OK;
    
    DSDP_INT m = RM->dim;
    vec *asinv = dsdpSolver->asinv;
    vec *d3    = dsdpSolver->d3;
    double mu  = dsdpSolver->mu;
    double tau = dsdpSolver->tau;
    
    assert( dsdpSolver->m == m );
    
    retcode = vec_reset(RM);
    retcode = vec_copy(d3, RM);
    retcode = vec_axpby(tau, dsdpSolver->dObj, mu, RM);
    retcode = vec_axpby(- mu, asinv, 1.0, RM);
    
    return retcode;
}

static DSDP_INT setuprM( HSDSolver *dsdpSolver, double *rM ) {
    // Set up the auxiliary vector rM
    // rM = -dObj + mu / tau + mu * csinv - mu * csinvrysinv
    DSDP_INT retcode = DSDP_RETCODE_OK;
    
    double mu = dsdpSolver->mu;
    *rM = - dsdpSolver->dObjVal + mu * (1 / dsdpSolver->tau
                                        + dsdpSolver->csinv
                                        - dsdpSolver->csinvrysinv);
    
    return retcode;
}

static DSDP_INT setupb( HSDSolver *dsdpSolver, vec *b1, vec *b2 ) {
    // Set up the auxiliary vector b1 = b + mu * u and b2 = b - mu * u
    DSDP_INT retcode = DSDP_RETCODE_OK;
    
    vec *b = dsdpSolver->lpObj;
    vec *u = dsdpSolver->u;
    double mu = dsdpSolver->mu;
    
    retcode = vec_copy(b, b1);
    retcode = vec_copy(b, b2);
    retcode = vec_axpy(  mu, u, b1);
    retcode = vec_axpy(- mu, u, b1);
    
    return retcode;
}

extern DSDP_INT setupFactorize( HSDSolver *dsdpSolver ) {
    // Factorize all the dual solutions
    DSDP_INT retcode = DSDP_RETCODE_OK;
    
    retcode = checkIterProgress(dsdpSolver, ITER_DUAL_FACTORIZE);
    assert( !dsdpSolver->iterProgress[ITER_DUAL_FACTORIZE] );
    
    if (dsdpSolver->iterProgress[ITER_DUAL_FACTORIZE]) {
        error(etype, "Dual variables have been factorized. \n");
    }
    
    DSDP_INT nblock = dsdpSolver->nBlock;
    
    for (DSDP_INT i = 0; i < nblock; ++i) {
        retcode = spsMatFactorize(dsdpSolver->S[i]);
        checkCode
    }
    
    return retcode;
}

extern DSDP_INT setupSchur( HSDSolver *dsdpSolver ) {
    // Setup the schur matrix Msdp and some of the temporary arrays
    DSDP_INT retcode = DSDP_RETCODE_OK;
    
    retcode = setupSDPSchur(dsdpSolver);
    retcode = setupLPSchur(dsdpSolver);
    
    return retcode;
}

extern DSDP_INT setupAux( HSDSolver *dsdpSolver,
                          vec       *RM,
                          double    *rM,
                          vec       *b1,
                          vec       *b2 ) {
    
    // Set up auxiliary arrays rm and RM
    DSDP_INT retcode = DSDP_RETCODE_OK;
    
    retcode = checkIterProgress(dsdpSolver, ITER_AUX_ARRAY);
    assert( !dsdpSolver->iterProgress[ITER_AUX_ARRAY] );
    
    if (dsdpSolver->iterProgress[ITER_AUX_ARRAY]) {
        error(etype, "Auxiliary arrays have been set up. \n");
    }
    
    retcode = setupRM(dsdpSolver, RM); checkCode;
    retcode = setuprM(dsdpSolver, rM); checkCode;
    retcode = setupb(dsdpSolver, b1, b2); checkCode;
    
    dsdpSolver->iterProgress[ITER_AUX_ARRAY] = TRUE;
    
    return retcode;
}

extern DSDP_INT schurMatSolve( HSDSolver *dsdpSolver, vec *b1, vec *RM ) {
    // Solve the internal system for getting the directions
    // After this routine, d1 and d2 will be filled
    DSDP_INT retcode = DSDP_RETCODE_OK;
    
    retcode = checkIterProgress(dsdpSolver, ITER_SCHUR_SOLVE);
    
    assert( !dsdpSolver->iterProgress[ITER_SCHUR_SOLVE] );
    
    if (dsdpSolver->iterProgress[ITER_SCHUR_SOLVE]) {
        error(etype, "Schur matrix has been factorized. \n");
    }
    
    dsMat *M = dsdpSolver->Msdp;
    assert( !M->isFactorized );
    
    DSDP_INT m = M->dim;
    assert( b1->dim == m );
    assert( RM->dim == m );
    
    // Factorize then solve
    retcode = denseMatFactorize(M); checkCode;
    double *b1RM  = (double *) calloc(2 * m, sizeof(double));
    
    memcpy(b1RM, b1->x, sizeof(double) * m);
    memcpy(&b1RM[m], RM->x, sizeof(double) * m);
    
    retcode = denseArrSolveInp(M, 2, b1RM);
    
    // Collect solution
    memcpy(dsdpSolver->d1->x, b1RM, sizeof(double) * m);
    memcpy(dsdpSolver->d2->x, &b1RM[m], sizeof(double) * m);
    
    DSDP_FREE(b1RM);
    
    dsdpSolver->iterProgress[ITER_SCHUR_SOLVE] = TRUE;
    
    return retcode;
}
