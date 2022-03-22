#include "schurmat.h"
#include "dsdputils.h"

#define M1Threshold 0.7
static char etype[] = "Schur matrix setup";

static DSDP_INT setupSDPSchur( HSDSolver *dsdpSolver ) {
    // Set up the schur matrix for SDP
    // After calling this routine, Msdp, asinv, u for SDP and csinv will be filled
    DSDP_INT retcode = DSDP_RETCODE_OK;
    
    assert( !dsdpSolver->iterProgress[ITER_SCHUR] );
    if (dsdpSolver->iterProgress[ITER_SCHUR]) {
        error(etype, "Schur matrix is already setup. \n");
        return retcode;
    }
    
    dsdpSolver->schurmu = dsdpSolver->mu;
    retcode = DSDPSchurSetup(dsdpSolver->M);
    dsdpSolver->iterProgress[ITER_SCHUR] = TRUE;
    return retcode;
}

static DSDP_INT setupBoundYSchur( HSDSolver *dsdpSolver ) {
    // Set up the Schur matrix for the bound
    DSDP_INT retcode = DSDP_RETCODE_OK;
    
    double *M = dsdpSolver->Msdp->array, *asinv = dsdpSolver->asinv->x, s = 0.0;
    double *sl = dsdpSolver->sl->x, *su = dsdpSolver->su->x;
    DSDP_INT m = dsdpSolver->m;
    
    for (DSDP_INT i = 0, idx = 0; i < m; ++i) {
        s = su[i]; M[idx] += 1.0 / (s * s); asinv[i] += 1.0 / s;
        s = sl[i]; M[idx] += 1.0 / (s * s); asinv[i] -= 1.0 / s;
        idx += m - i;
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
    
    vec *c  = dsdpSolver->lpObj;
    vec *ry = dsdpSolver->ry;
    vec *u  = dsdpSolver->u;
    vec *d3 = dsdpSolver->d3;
    
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

static DSDP_INT schurMatPerturb( HSDSolver *dsdpSolver ) {
    
    // Perturb the Schur matrix
    DSDP_INT retcode = DSDP_RETCODE_OK;
    DSDP_INT m = dsdpSolver->m, i;
    double perturb = 0.0, maxdiag = 0.0;
    denseMatGetdiag(dsdpSolver->Msdp, dsdpSolver->Mdiag);
    
    if (!dsdpSolver->eventMonitor[EVENT_INVALID_GAP]) {
        
        if (dsdpSolver->mu < 1e-05) {
             perturb += 1e-13;
        }
        
        for (DSDP_INT i = 0; i < m; ++i) {
            maxdiag = MAX(dsdpSolver->Mdiag->x[i], maxdiag);
        }
        
        dsdpSolver->Mscaler = maxdiag;
        
        if (dsdpSolver->m < 100) {
            perturb += MIN(maxdiag * 1e-08, 1e-12);
        } else if (dsdpSolver->m < 1000) {
            perturb += MIN(maxdiag * 1e-08, 1e-11);
        } else {
            perturb += MIN(maxdiag * 1e-08, 1e-09);
        }
        
        double invalid;
        DSDPGetStats(&dsdpSolver->dsdpStats, STAT_GAP_BROKEN, &invalid);
        
        if (dsdpSolver->eventMonitor[EVENT_IN_PHASE_B] && !invalid) {
            for (i = 0; i < m; ++i) {
                packIdx(dsdpSolver->Msdp->array, m, i, i) += perturb;
            }
        }
    }
    
    return retcode;
}

static DSDP_INT schurMatscale( HSDSolver *dsdpSolver ) {
    
    // Scale Schur matrix if necessary
    
    dsdpSolver->Mscaler = 1.0;
    return DSDP_RETCODE_OK;
    dsMat *M = dsdpSolver->Msdp;
    double scaler = sqrt(dsdpSolver->Mscaler);
    
    if (dsdpSolver->eventMonitor[EVENT_IN_PHASE_A]) {
        if (scaler > 1e+06 || scaler <= 1e-04) {
            denseMatRscale(M, scaler);
            vec_rscale(dsdpSolver->d2, scaler);
            vec_rscale(dsdpSolver->d12, scaler);
            vec_rscale(dsdpSolver->d3, scaler);
            vec_rscale(dsdpSolver->d4, scaler);
        } else {
            dsdpSolver->Mscaler = 1.0;
        }
    } else {
        if (scaler > 1e+06 || scaler <= 1e-04) {
            denseMatRscale(M, scaler);
            vec_rscale(dsdpSolver->d1, scaler);
            vec_rscale(dsdpSolver->d2, scaler);
        } else {
            dsdpSolver->Mscaler = 1.0;
        }
    }
    
    return DSDP_RETCODE_OK;
}

static DSDP_INT schurCGSetup( HSDSolver *dsdpSolver ) {
    // Set CG tolerance and maxiteration based on simple heuristic
    DSDP_INT retcode = DSDP_RETCODE_OK;
    
    double tol = 0.0; DSDP_INT cgiter = 20;
    CGSolver *cgsolver = dsdpSolver->cgSolver;
    
    // Set parameters
    if (dsdpSolver->mu > 1.0) {
        tol = 1e-04;
    } else if (dsdpSolver->mu > 1e-02) {
        tol = 1e-05;
    } else if (dsdpSolver->mu > 1e-05){
        tol = 1e-06;
    } else if (dsdpSolver->mu > 1e-07) {
        tol = 1e-07;
    } else {
        tol = 1e-08;
    }
    
    
    cgsolver->status = CG_STATUS_INDEFINITE;
    // dsdpCGSetPreReuse(cgsolver, 0);
    
    if (dsdpSolver->m > 20000) {
        cgiter = 500; tol *= 10.0;
    } else if (dsdpSolver->m > 10000) {
        cgiter = 400; tol *= 1.2;
    } else if (dsdpSolver->m > 5000) {
        cgiter = 250; tol *= 1.1;
    } else {
        cgiter = MAX(dsdpSolver->m / 50, 30);
    }
    
    dsdpCGSetTol(cgsolver, tol);
    dsdpCGSetMaxIter(cgsolver, cgiter);
    dsdpCGprepareP(cgsolver);
    
    return retcode;
}

static DSDP_INT setupPhaseAdvecs( HSDSolver *dsdpSolver ) {
    // Set up the auxiliary vectors
    /*
     d2_12_3_4 = M \ [b, u, asinv, asinvrysinv];
    */
    DSDP_INT retcode = DSDP_RETCODE_OK;
    retcode = vec_copy(dsdpSolver->dObj,  dsdpSolver->d2);
    retcode = vec_copy(dsdpSolver->u,     dsdpSolver->d12);
    retcode = vec_copy(dsdpSolver->asinv, dsdpSolver->d3);
    
    return retcode;
}

static DSDP_INT setupPhaseBdvecs( HSDSolver *dsdpSolver ) {
    // Set up the auxiliary vectors
    /*
     dy1dy2 = Mhat \ [b, asinv];
     dy1 is d1
     dy2 is d2
     */
    
    DSDP_INT retcode = DSDP_RETCODE_OK;
    
    retcode = vec_copy(dsdpSolver->dObj, dsdpSolver->d1);
    retcode = vec_copy(dsdpSolver->asinv, dsdpSolver->d2);
    
    return retcode;
}

static DSDP_INT cgSolveCheck( CGSolver *cgSolver, vec *b ) {
    
    DSDP_INT retcode = DSDP_RETCODE_OK;
    DSDP_INT status;
    
    double *tmp = (double *) calloc(b->dim, sizeof(double));
    memcpy(tmp, b->x, sizeof(double) * b->dim);
    dsdpCGSolve(cgSolver, b, NULL);
    dsdpCGGetStatus(cgSolver, &status);
    
    if (status != CG_STATUS_SOLVED && status != CG_STATUS_INDEFINITE) {
        
        if (cgSolver->pType == CG_PRECOND_DIAG) {
            dsdpCGSetPType(cgSolver, CG_PRECOND_CHOL);
        } else {
            cgSolver->nfailed += 1;
        }
        
        cgSolver->nused = 1024;
        if (cgSolver->nfailed >= 5) { dsdpCGSetPreReuse(cgSolver, MIN(cgSolver->reuse, 5));}
        if (cgSolver->nfailed >= 10) { dsdpCGSetPreReuse(cgSolver, MIN(cgSolver->reuse, 3));}
        if (cgSolver->nfailed >= 15) { dsdpCGSetPreReuse(cgSolver, MIN(cgSolver->reuse, 1));}
        if (cgSolver->nfailed >= 50) { cgSolver->M->isillCond = TRUE;}
        
        // Restart from the last unfinished solution
        dsdpCGprepareP(cgSolver);
        memcpy(b->x, tmp, sizeof(double) * b->dim);
        dsdpCGSolve(cgSolver, b, cgSolver->x);
        dsdpCGGetStatus(cgSolver, &status);
        
        if (status != CG_STATUS_SOLVED && status != CG_STATUS_INDEFINITE) {
            error(etype, "Invalid CG iteration. \n");
        }
    }
    
    DSDP_FREE(tmp);
    return retcode;
}

extern DSDP_INT setupFactorize( HSDSolver *dsdpSolver ) {
    // Factorize the dual solution and compute slack
    DSDP_INT retcode = DSDP_RETCODE_OK;
    retcode = checkIterProgress(dsdpSolver, ITER_DUAL_FACTORIZE);
    assert( !dsdpSolver->iterProgress[ITER_DUAL_FACTORIZE] );
    if (dsdpSolver->iterProgress[ITER_DUAL_FACTORIZE]) {
        error(etype, "Dual variables have been factorized. \n");
    }
    
    DSDP_INT nblock = dsdpSolver->nBlock;
    double iterA = 0.0;
    
    DSDPGetStats(&dsdpSolver->dsdpStats, STAT_PHASE_A_ITER, &iterA);
    
    if (iterA == 0.0) {
        for (DSDP_INT i = 0; i < nblock; ++i) {
            retcode = spsMatSymbolic(dsdpSolver->S[i]);
            retcode = spsMatFactorize(dsdpSolver->S[i]);
            retcode = spsMatSymbolic(dsdpSolver->Scker[i]);
            checkCode;
        }
    } else {
        for (DSDP_INT i = 0; i < nblock; ++i) {
            retcode = spsMatFactorize(dsdpSolver->S[i]);
            checkCode;
        }
    }
    
    // Get dual slack
    double bound = 0.0; DSDPGetDblParam(dsdpSolver, DBL_PARAM_PRLX_PENTALTY, &bound);
    vec_lslack(dsdpSolver->y, dsdpSolver->sl, -bound);
    vec_uslack(dsdpSolver->y, dsdpSolver->su, bound);
    
    dsdpSolver->iterProgress[ITER_DUAL_FACTORIZE] = TRUE;
    return retcode;
}

extern DSDP_INT schurPhaseAMatSolve( HSDSolver *dsdpSolver ) {
    // Solve the internal system to get the directions
    // After this routine, d1, d12, d3 and d4 will be filled
    DSDP_INT retcode = DSDP_RETCODE_OK;
    retcode = checkIterProgress(dsdpSolver, ITER_SCHUR_SOLVE);
    assert( !dsdpSolver->iterProgress[ITER_SCHUR_SOLVE] );
    if (dsdpSolver->iterProgress[ITER_SCHUR_SOLVE]) {
        error(etype, "Schur matrix has been factorized. \n");
    }
    // Currently no CG warmstart
    retcode = cgSolveCheck(dsdpSolver->cgSolver, dsdpSolver->d2);
    retcode = cgSolveCheck(dsdpSolver->cgSolver, dsdpSolver->d12);
    retcode = cgSolveCheck(dsdpSolver->cgSolver, dsdpSolver->d3);
    retcode = cgSolveCheck(dsdpSolver->cgSolver, dsdpSolver->d4);
    
    dsdpSolver->iterProgress[ITER_SCHUR_SOLVE] = TRUE;
    
    return retcode;
}

extern DSDP_INT schurPhaseBMatSolve( HSDSolver *dsdpSolver ) {
    // Solve the internal system to get the directions
    // After this routine, d1 and d2 will be filled
    DSDP_INT retcode = DSDP_RETCODE_OK;
    retcode = checkIterProgress(dsdpSolver, ITER_SCHUR_SOLVE);
    assert( !dsdpSolver->iterProgress[ITER_SCHUR_SOLVE] );
    if (dsdpSolver->iterProgress[ITER_SCHUR_SOLVE]) {
        error(etype, "Schur matrix has been factorized. \n");
    }
    cgSolveCheck(dsdpSolver->cgSolver, dsdpSolver->d1);
    cgSolveCheck(dsdpSolver->cgSolver, dsdpSolver->d2);
    dsdpSolver->iterProgress[ITER_SCHUR_SOLVE] = TRUE;
    
    return retcode;
}

extern DSDP_INT setupSchur( HSDSolver *dsdpSolver ) {
    // Setup the schur matrix Msdp and some of the temporary arrays
    DSDP_INT retcode = DSDP_RETCODE_OK;
    
    retcode = setupSDPSchur(dsdpSolver);
    
    if (dsdpSolver->eventMonitor[EVENT_IN_PHASE_A]) {
        // TODO: Add primal relaxation to this phase ?
    } else {
        retcode = setupBoundYSchur(dsdpSolver);
    }
    
    checkCode;
    dsdpSolver->iterProgress[ITER_SCHUR] = TRUE;
    // retcode = setupLPSchur(dsdpSolver);
    schurMatPerturb(dsdpSolver);
    if (dsdpSolver->eventMonitor[EVENT_IN_PHASE_A]) {
        retcode = setupPhaseAdvecs(dsdpSolver); checkCode;
    } else {
        retcode = setupPhaseBdvecs(dsdpSolver); checkCode;
    }
    schurMatscale(dsdpSolver);
    schurCGSetup(dsdpSolver);
    if (dsdpSolver->eventMonitor[EVENT_IN_PHASE_A]) {
        retcode = schurPhaseAMatSolve(dsdpSolver);
    } else {
        retcode = schurPhaseBMatSolve(dsdpSolver);
    }
    
    return retcode;
}
