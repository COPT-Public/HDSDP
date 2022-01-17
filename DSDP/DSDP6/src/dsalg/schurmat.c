#include "schurmat.h"
#include "dsdputils.h"

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
    r1Mat *r1data = (r1Mat *) calloc(1, sizeof(r1Mat));
    dsMat *dsdata = (dsMat *) calloc(1, sizeof(dsMat));
    void *data = NULL;
    
    for (DSDP_INT i = 0; i < m + 1; ++i) {
        // Compute SinvASinv
        mattype = sdpData->types[i];
        if (mattype == MAT_TYPE_ZERO) {
            continue;
        } else if (mattype == MAT_TYPE_RANK1) {
            retcode = getSinvASinv(dsdpSolver, blockid, i, r1data);
            data = (void *) r1data;
            checkCodeFree;
            // retcode = r1MatCountNnz(r1data);
        } else {
            retcode = getSinvASinv(dsdpSolver, blockid, i, dsdata);
            data = (void *) dsdata;
            checkCodeFree;
        }
        
        for (DSDP_INT j = 0; j <= i; ++j) {
            getTraceASinvASinv(dsdpSolver, blockid, j, mattype, i, data);
        }
    }
    
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

static DSDP_INT setupPhaseAdvecs( HSDSolver *dsdpSolver ) {
    // Set up the auxiliary vector
    /*
     d2_12_3_4 = M \ [b, u, asinv, asinvrysinv];
    */
    DSDP_INT retcode = DSDP_RETCODE_OK;

    retcode = vec_copy(dsdpSolver->dObj, dsdpSolver->d2);
    retcode = vec_copy(dsdpSolver->u, dsdpSolver->d12);
    retcode = vec_copy(dsdpSolver->asinv, dsdpSolver->d3);
    retcode = vec_copy(dsdpSolver->asinv, dsdpSolver->d4);
    
    return retcode;
}

static DSDP_INT setupFactorize( HSDSolver *dsdpSolver ) {
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

extern DSDP_INT schurPhaseAMatSolve( HSDSolver *dsdpSolver ) {
    // Solve the internal system for getting the directions
    // After this routine, d1 and d2 will be filled
    DSDP_INT retcode = DSDP_RETCODE_OK;
    retcode = checkIterProgress(dsdpSolver, ITER_SCHUR_SOLVE);
    
    assert( !dsdpSolver->iterProgress[ITER_SCHUR_SOLVE] );
    
    if (dsdpSolver->iterProgress[ITER_SCHUR_SOLVE]) {
        error(etype, "Schur matrix has been factorized. \n");
    }
    
    retcode = setupPhaseAdvecs(dsdpSolver); checkCode;
    
    dsMat *M = dsdpSolver->Msdp;
    assert( !M->isFactorized );
    
    denseArrSolveInp(M, 1, dsdpSolver->d2->x);
    denseArrSolveInp(M, 1, dsdpSolver->d12->x);
    denseArrSolveInp(M, 1, dsdpSolver->d3->x);
    denseArrSolveInp(M, 1, dsdpSolver->d4->x);
    
    
    dsdpSolver->iterProgress[ITER_SCHUR_SOLVE] = TRUE;
    
    return retcode;
}

extern DSDP_INT setupPhaseASchur( HSDSolver *dsdpSolver ) {
    // Setup the schur matrix Msdp and some of the temporary arrays
    DSDP_INT retcode = DSDP_RETCODE_OK;
    
    retcode = setupFactorize(dsdpSolver);
    retcode = setupSDPSchur(dsdpSolver);
    retcode = setupLPSchur(dsdpSolver);
    retcode = schurPhaseAMatSolve(dsdpSolver);
    
    return retcode;
}

