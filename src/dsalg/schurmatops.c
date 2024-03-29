#include "schurmatops.h"
#include "dsdplapack.h"
#include "dsdplog.h"
#include "adpcg.h"
#include "dsdputils.h"
#include "schurmat.h"
#include "vec.h"

#define M1Threshold 0.7

static DSDP_INT schurMatPerturb( HSDSolver *dsdpSolver ) {
    
    // Perturb the Schur matrix
    DSDP_INT retcode = DSDP_RETCODE_OK;
    double perturb = 0.0, maxdiag = 0.0, iterB;
    
    schurMatGetdiag(dsdpSolver->Msdp, dsdpSolver->Mdiag);
    DSDPGetStats(&dsdpSolver->dsdpStats, STAT_PHASE_B_ITER, &iterB);
    
    if (!dsdpSolver->eventMonitor[EVENT_INVALID_GAP] &&
        dsdpSolver->eventMonitor[EVENT_IN_PHASE_A]) {
        
        if (dsdpSolver->mu < 1e-05) { perturb += 1e-12; }
        maxdiag = vec_infnorm(dsdpSolver->Mdiag);
        dsdpSolver->Mscaler = maxdiag;
        
        if (dsdpSolver->m < 100) {
            perturb += MIN(maxdiag * 1e-08, 1e-14);
        } else if (dsdpSolver->m < 1000) {
            perturb += MIN(maxdiag * 1e-08, 1e-13);
        } else {
            perturb += MIN(maxdiag * 1e-08, 1e-12);
        }
        
//        if (dsdpSolver->Msdp->isillCond ||
//            dsdpSolver->cgSolver->status == CG_STATUS_INDEFINITE) {
//            perturb = 0.0;
//        }
        
        // if (dsdpSolver->mu < 1e-05) { perturb += maxdiag * 1e-05; }
        double invalid;
        DSDPGetStats(&dsdpSolver->dsdpStats, STAT_GAP_BROKEN, &invalid);
        if (dsdpSolver->eventMonitor[EVENT_IN_PHASE_B] && !invalid) {
            schurMatAdddiag(dsdpSolver->Msdp, perturb);
        }
        schurMatGetdiag(dsdpSolver->Msdp, dsdpSolver->Mdiag);
    }
    
    return retcode;
}

static DSDP_INT schurCGSetup( HSDSolver *dsdpSolver ) {
    // Set CG tolerance and maxiteration based on simple heuristic
    DSDP_INT retcode = DSDP_RETCODE_OK;
    
    double tol = 0.0; DSDP_INT cgiter;
    adpcg *cg = dsdpSolver->cg;
    
    retcode = cg_start(cg);
    
    // Set parameters
    if (dsdpSolver->mu > 1000) {
        tol = 1e-05;
    } else if (dsdpSolver->mu > 1.0) {
        tol = 1e-07;
    } else if (dsdpSolver->mu > 1e-02) {
        tol = 1e-08;
    } else if (dsdpSolver->mu > 1e-05){
        tol = 1e-09;
    } else if (dsdpSolver->mu > 1e-07) {
        tol = 1e-10;
    } else {
        tol = 1e-12;
    }
    
    if (dsdpSolver->m > 20000) {
        cgiter = 500; tol *= 10.0;
    } else if (dsdpSolver->m > 15000) {
        cgiter = 450; tol *= 10.0;
    } else if (dsdpSolver->m > 5000) {
        cgiter = 120; tol *= 2.0;
    } else {
        cgiter = MAX(dsdpSolver->m / 50, 30);
    }
    
    cg_setparam(cg, tol, -1, cgiter, -1);
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

extern DSDP_INT setupFactorize( HSDSolver *dsdpSolver ) {
    DSDP_INT retcode = DSDP_RETCODE_OK;
    checkIterProgress(dsdpSolver, ITER_DUAL_FACTORIZE);
    double iterA; DSDPGetStats(&dsdpSolver->dsdpStats, STAT_PHASE_A_ITER, &iterA);
    if (iterA == 0.0) {
        DSDPConic( COPS_GET_SLACK )(dsdpSolver, DUALVAR);
        DSDPConic( COPS_SYMFAC )(dsdpSolver);
        if (!DSDPConic( COPS_CHECK_INCONE ) (dsdpSolver, DUALVAR)) {
            printf("| Invalid starting point. \n");
            retcode = DSDP_RETCODE_FAILED;
        }
    }
    dsdpSolver->iterProgress[ITER_DUAL_FACTORIZE] = TRUE;
    return retcode;
}

extern DSDP_INT schurPhaseAMatSolve( HSDSolver *dsdpSolver ) {
    // Solve the internal system to get the directions
    // After this routine, d1, d12, d3 and d4 will be filled
    DSDP_INT retcode = DSDP_RETCODE_OK;
    checkIterProgress(dsdpSolver, ITER_SCHUR_SOLVE);
    retcode = cg_solve(dsdpSolver->cg, dsdpSolver->d2, NULL);
    if ( retcode != DSDP_RETCODE_OK ) {
        return retcode;
    }
    
    retcode = cg_solve(dsdpSolver->cg, dsdpSolver->d3, NULL);
    if ( retcode != DSDP_RETCODE_OK ) {
        return retcode;
    }
    
    retcode = cg_solve(dsdpSolver->cg, dsdpSolver->d4, NULL);
    if ( retcode != DSDP_RETCODE_OK ) {
        return retcode;
    }
    dsdpSolver->iterProgress[ITER_SCHUR_SOLVE] = TRUE;
    return retcode;
}

extern DSDP_INT schurPhaseBMatSolve( HSDSolver *dsdpSolver ) {
    // Solve the internal system to get the directions
    // After this routine, d1 and d2 will be filled
    DSDP_INT retcode = DSDP_RETCODE_OK;
    checkIterProgress(dsdpSolver, ITER_SCHUR_SOLVE);
    cg_solve(dsdpSolver->cg, dsdpSolver->d1, NULL);
    cg_solve(dsdpSolver->cg, dsdpSolver->d2, NULL);
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

extern DSDP_INT setupSchur( HSDSolver *dsdpSolver ) {
    // Setup the schur matrix Msdp and some of the temporary arrays
    DSDP_INT retcode = DSDP_RETCODE_OK;
    DSDPConic( COPS_GET_SCHUR ) (dsdpSolver);
    dsdpSolver->iterProgress[ITER_SCHUR] = TRUE;
    schurMatPerturb(dsdpSolver); dsdpSolver->Mscaler = 1.0;
    schurCGSetup(dsdpSolver);
    if (dsdpSolver->eventMonitor[EVENT_IN_PHASE_A]) {
        setupPhaseAdvecs(dsdpSolver); schurPhaseAMatSolve(dsdpSolver);
    } else {
        setupPhaseBdvecs(dsdpSolver); schurPhaseBMatSolve(dsdpSolver);
    }
    return retcode;
}
