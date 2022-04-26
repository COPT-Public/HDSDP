#include "schurmatops.h"
#include "dsdputils.h"
#include "schurmat.h"

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
    
    if (dsdpSolver->ybound == DSDP_INFINITY) {
        return retcode;
    }
    
    double **Mdiag = dsdpSolver->Msdp->diag, *asinv = dsdpSolver->asinv->x, s = 0.0;
    double *sl = dsdpSolver->sl->x, *su = dsdpSolver->su->x, bound = dsdpSolver->ybound;
    DSDP_INT m = dsdpSolver->m;
    
    if (dsdpSolver->eventMonitor[EVENT_IN_PHASE_A]) {
        double *u = dsdpSolver->u->x;
        double sinvsqr, csinvsuml = 0.0, csinvsumu = 0.0, cscssuml = 0.0, cscssumu = 0.0;
        
        for (DSDP_INT i = 0; i < m; ++i) {
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
        
        dsdpSolver->csinv      += bound * csinvsumu - bound * csinvsuml;
        dsdpSolver->csinvcsinv += bound * bound * (cscssumu + cscssuml);
        
    } else {
        for (DSDP_INT i = 0; i < m; ++i) {
            s = su[i]; *Mdiag[i] += 1.0 / (s * s); asinv[i] += 1.0 / s;
            s = sl[i]; *Mdiag[i] += 1.0 / (s * s); asinv[i] -= 1.0 / s;
        }
    }
    
    return retcode;
}

static DSDP_INT schurMatPerturb( HSDSolver *dsdpSolver ) {
    
    // Perturb the Schur matrix
    DSDP_INT retcode = DSDP_RETCODE_OK;
    double perturb = 0.0, maxdiag = 0.0, iterB;
    
    schurMatGetdiag(dsdpSolver->Msdp, dsdpSolver->Mdiag);
    DSDPGetStats(&dsdpSolver->dsdpStats, STAT_PHASE_B_ITER, &iterB);
    
    if (!dsdpSolver->eventMonitor[EVENT_INVALID_GAP]) {
        if (dsdpSolver->mu < 1e-05) { perturb += 1e-13; }
        maxdiag = vec_infnorm(dsdpSolver->Mdiag);
        dsdpSolver->Mscaler = maxdiag;
        
        if (dsdpSolver->m < 100) {
            perturb += MIN(maxdiag * 1e-08, 1e-14);
        } else if (dsdpSolver->m < 1000) {
            perturb += MIN(maxdiag * 1e-08, 1e-13);
        } else {
            perturb += MIN(maxdiag * 1e-08, 1e-12);
        }
        
        double invalid;
        DSDPGetStats(&dsdpSolver->dsdpStats, STAT_GAP_BROKEN, &invalid);
        if (dsdpSolver->eventMonitor[EVENT_IN_PHASE_B] && !invalid) {
            schurMatAdddiag(dsdpSolver->Msdp, perturb);
        }
        schurMatGetdiag(dsdpSolver->Msdp, dsdpSolver->Mdiag);
    }
    
    return retcode;
}

static DSDP_INT schurMatscale( HSDSolver *dsdpSolver ) {
    
    // Scale Schur matrix if necessary
    dsdpSolver->Mscaler = 1.0;
    if ((TRUE)) { return DSDP_RETCODE_OK; }
    double scaler = dsdpSolver->Mscaler;
    if (dsdpSolver->eventMonitor[EVENT_IN_PHASE_A]) {
        if (scaler > 1e+06 || scaler <= 1e-04) {
            vec_scale(dsdpSolver->d2, scaler);
            vec_scale(dsdpSolver->d12, scaler);
            vec_scale(dsdpSolver->d3, scaler);
            vec_scale(dsdpSolver->d4, scaler);
        } else {
            dsdpSolver->Mscaler = 1.0;
        }
    } else {
        if (scaler > 1e+06 || scaler <= 1e-04) {
            vec_scale(dsdpSolver->d1, scaler);
            vec_scale(dsdpSolver->d2, scaler);
        } else {
            dsdpSolver->Mscaler = 1.0;
        }
    }
    
    return DSDP_RETCODE_OK;
}

static DSDP_INT schurCGSetup( HSDSolver *dsdpSolver ) {
    // Set CG tolerance and maxiteration based on simple heuristic
    DSDP_INT retcode = DSDP_RETCODE_OK;
    
    double tol = 0.0; DSDP_INT cgiter;
    CGSolver *cgsolver = dsdpSolver->cgSolver;
    
    // Set parameters
    if (dsdpSolver->mu > 1.0) {
        tol = 1e-05;
    } else if (dsdpSolver->mu > 1e-02) {
        tol = 1e-06;
    } else if (dsdpSolver->mu > 1e-05){
        tol = 5e-07;
    } else if (dsdpSolver->mu > 1e-07) {
        tol = 1e-07;
    } else {
        tol = 1e-07;
    }
    
#ifdef compareMode
    cgsolver->status = CG_STATUS_INDEFINITE;
#endif
    
    if (dsdpSolver->m > 20000) {
        cgiter = 500; tol *= 10.0;
    } else if (dsdpSolver->m > 15000) {
        cgiter = 400; tol *= 10.0;
    } else if (dsdpSolver->m > 5000) {
        cgiter = 250; tol *= 2.0;
    } else {
        cgiter = MAX(dsdpSolver->m / 50, 30);
    }
    
    dsdpCGSetTol(cgsolver, tol);
    dsdpCGSetMaxIter(cgsolver, cgiter);
    dsdpCGprepareP(cgsolver);
    
    return retcode;
}

static DSDP_INT setupPhaseAdvecs( HSDSolver *dsdpSolver ) {
    vec_copy(dsdpSolver->dObj,  dsdpSolver->d2);
    vec_copy(dsdpSolver->u,     dsdpSolver->d12);
    vec_copy(dsdpSolver->asinv, dsdpSolver->d3);
    vec_copy(dsdpSolver->asinvrysinv, dsdpSolver->d4);
    return DSDP_RETCODE_OK;
}

static DSDP_INT setupPhaseBdvecs( HSDSolver *dsdpSolver ) {
    vec_copy(dsdpSolver->dObj,  dsdpSolver->d1);
    vec_copy(dsdpSolver->asinv, dsdpSolver->d2);
    return DSDP_RETCODE_OK;
}

static DSDP_INT cgSolveCheck( CGSolver *cgSolver, vec *b ) {
    
    DSDP_INT retcode = DSDP_RETCODE_OK, status;
    
    dsdpCGStoreRHS(cgSolver, b);
    dsdpCGSolve(cgSolver, b, NULL);
    dsdpCGGetStatus(cgSolver, &status);
    
    if (status != CG_STATUS_SOLVED && status != CG_STATUS_INDEFINITE) {
        if (cgSolver->pType == CG_PRECOND_DIAG) {
            dsdpCGSetPType(cgSolver, CG_PRECOND_CHOL);
        } else {
            cgSolver->nfailed += 1;
        }
        cgSolver->nused = 1024;
        if (cgSolver->nfailed >= 5) {
            dsdpCGSetPreReuse(cgSolver, MIN(cgSolver->reuse, 5));
            if (cgSolver->nfailed >= 10) { dsdpCGSetPreReuse(cgSolver, MIN(cgSolver->reuse, 3));}
            if (cgSolver->nfailed >= 15) { dsdpCGSetPreReuse(cgSolver, MIN(cgSolver->reuse, 2));}
            if (cgSolver->nfailed >= 50) { cgSolver->M->isillCond = TRUE;}
        }
        dsdpCGprepareP(cgSolver);
        dsdpCGRestoreRHS(cgSolver, b);
        dsdpCGSolve(cgSolver, b, cgSolver->x);
        // dsdpCGGetStatus(cgSolver, &status);
    }
    
    return retcode;
}

static DSDP_INT setupSDPCones( HSDSolver *dsdpSolver ) {
    
    DSDP_INT retcode = DSDP_RETCODE_OK;
    DSDP_INT nblock = dsdpSolver->nBlock; double iterA = 0.0;
    DSDPGetStats(&dsdpSolver->dsdpStats, STAT_PHASE_A_ITER, &iterA);
    
    if (iterA == 0.0) {
        for (DSDP_INT i = 0; i < nblock; ++i) {
            retcode = spsMatSymbolic(dsdpSolver->S[i]);
            retcode = spsMatFactorize(dsdpSolver->S[i]);
            retcode = spsMatSymbolic(dsdpSolver->Scker[i]);
            checkCode;
        }
    } else {
//        for (DSDP_INT i = 0; i < nblock; ++i) {
//            retcode = spsMatFactorize(dsdpSolver->S[i]);
//            checkCode;
//        }
    }
    return retcode;
}

static void setupBoundYCones( HSDSolver *dsdpSolver ) {
    // Get dual slack
    double bound = dsdpSolver->ybound;
    if (bound != DSDP_INFINITY) {
        vec_lslack(dsdpSolver->y, dsdpSolver->sl, -bound * dsdpSolver->tau);
        vec_uslack(dsdpSolver->y, dsdpSolver->su,  bound * dsdpSolver->tau);
    }
}


extern DSDP_INT setupFactorize( HSDSolver *dsdpSolver ) {
    // Factorize the dual solution and compute slack
    DSDP_INT retcode = DSDP_RETCODE_OK;
    retcode = checkIterProgress(dsdpSolver, ITER_DUAL_FACTORIZE);
    assert( !dsdpSolver->iterProgress[ITER_DUAL_FACTORIZE] );
    if (dsdpSolver->iterProgress[ITER_DUAL_FACTORIZE]) {
        error(etype, "Dual variables have been factorized. \n");
    }
    
    retcode = setupSDPCones(dsdpSolver); checkCode;
    setupBoundYCones(dsdpSolver);
    
    dsdpSolver->iterProgress[ITER_DUAL_FACTORIZE] = TRUE;
    return retcode;
}

extern DSDP_INT schurPhaseAMatSolve( HSDSolver *dsdpSolver ) {
    // Solve the internal system to get the directions
    // After this routine, d1, d12, d3 and d4 will be filled
    DSDP_INT retcode = DSDP_RETCODE_OK;
    retcode = checkIterProgress(dsdpSolver, ITER_SCHUR_SOLVE);
    retcode = cgSolveCheck(dsdpSolver->cgSolver, dsdpSolver->d2);
    if (dsdpSolver->eventMonitor[EVENT_HSD_UPDATE]) {
        retcode = cgSolveCheck(dsdpSolver->cgSolver, dsdpSolver->d12);
    }
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
    cgSolveCheck(dsdpSolver->cgSolver, dsdpSolver->d1);
    cgSolveCheck(dsdpSolver->cgSolver, dsdpSolver->d2);
    dsdpSolver->iterProgress[ITER_SCHUR_SOLVE] = TRUE;
    return retcode;
}

extern void schurMatPrint( schurMat *sMat ) {
    if (sMat->stype == SCHUR_TYPE_DENSE) {
        for (DSDP_INT i = 0, j; i < sMat->m; ++i) {
            for (j = 0; j < sMat->m; ++j) {
                printf("%6.3e, ", fullIdx(sMat->denseM->array, sMat->m, i, j));
            }
            printf("\n");
        }
    }
}

/*
 
 parray 10 dsdpSolver->Msdp->denseM->array
 
 
 */
extern DSDP_INT setupSchur( HSDSolver *dsdpSolver ) {
    // Setup the schur matrix Msdp and some of the temporary arrays
    DSDP_INT retcode = DSDP_RETCODE_OK;
    
    setupSDPSchur(dsdpSolver);
#ifdef compareMode
        assert( TRUE );
#endif
    vec_copy(dsdpSolver->asinv, dsdpSolver->d12);
    setupBoundYSchur(dsdpSolver);
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
    
#ifdef compareMode
        assert( TRUE );
#endif
    return retcode;
}
