#include "dsdpdinfeas.h"
#include "dsdputils.h"
#include "dsdpinitializer.h"
#include "schurmatops.h"
#include "residualsetup.h"
#include "stepdirection.h"
#include "stepheur.h"
#include "dsdpproxmeasure.h"
#include "dsdpcorrector.h"
#include "dsdplog.h"

// Implement the phase A of the dual scaling algorithm: Dual infeasibility elimination
static char etype[] = "DSDP Dual infeasibility elimination";

/* Driver routine for dual infeasibility elimination
 
 This routine tries to elimininate dual infeasibility and returns a
 solution to start the second phase

 */
static const double nomu[1]    = {1.0};
static const double consvmu[2] = {0.3, 1.0};
static const double aggmu[4]   = {0.1, 0.4, 0.7, 1.0};

extern void funcA( dsMat *A, DSDP_INT i, DSDP_INT j) {
    assert( i >= j );
    printf("%30.20e \n", packIdx(A->array, A->dim, i, j));
}

extern DSDP_INT DSDPDInfeasEliminator( HSDSolver *dsdpSolver ) {
    
    DSDP_INT retcode = DSDP_RETCODE_OK;
    DSDP_INT goOn = TRUE;
    DSDPStats *stat = &dsdpSolver->dsdpStats;
    dsdpSolver->pObjVal = DSDP_INFINITY;
    
    // Initialize
    double muprimal, tol, sigma;
    DSDP_INT attempt, agiter;
    
    retcode = DSDPGetDblParam(dsdpSolver, DBL_PARAM_INIT_MU,    &muprimal);
    retcode = DSDPGetDblParam(dsdpSolver, DBL_PARAM_ABS_OPTTOL, &tol     );
    retcode = DSDPGetDblParam(dsdpSolver, DBL_PARAM_ASIGMA,     &sigma   );
    retcode = DSDPGetIntParam(dsdpSolver, INT_PARAM_AATTEMPT,   &attempt );
    retcode = DSDPGetIntParam(dsdpSolver, INT_PARAM_AMAXITER,   &agiter  );
    agiter  = MIN(agiter, 30);
    
    double trymu = 0.0;
    double time = 0.0;
    
    DSDP_INT ntry   = 0;
    double start    = my_clock();
    
    const double *newmu = NULL;
    
    if (attempt == AATEMT_CONSERVATIVE) {
        ntry = 1; newmu = nomu;
    } else if (attempt == AATEMPT_MILD) {
        ntry = 2; newmu = consvmu;
    } else {
        ntry = 4; newmu = aggmu;
    }
        
    dsdpSolver->eventMonitor[EVENT_IN_PHASE_A] = TRUE;
    dsdpSolver->eventMonitor[EVENT_IN_PHASE_B] = FALSE;
    
    /* Start Phase A algorithm */
    dsdpshowdash();
    retcode = dsdpInitializeA(dsdpSolver); checkCode;
    
    /* Print algorithm header */
    dsdpprintPhaseAheader();
    
    DSDP_INT i;
    for (i = 0; ; ++i) {
        
        // Start iteration
        DSDPStatUpdate(stat, STAT_PHASE_A_ITER, (double) i);
        dsdpSolver->iterProgress[ITER_LOGGING] = FALSE;
        dsdpSolver->iterProgress[ITER_DUAL_OBJ] = FALSE;
        
        // Check NaN
        dsdpCheckNan(dsdpSolver);
        // Check algorithm convergence
        DSDPCheckPhaseAConvergence(dsdpSolver, &goOn);
        // Compute dual objective
        getDualObj(dsdpSolver);
        // Logging
        DSDPPhaseALogging(dsdpSolver);
        // Reset monitor
        DSDPResetPhaseAMonitor(dsdpSolver);
        
        if (i > 90) {
            dsdpSolver->eventMonitor[EVENT_HSD_UPDATE] = TRUE;
        }
        
        if (goOn) {
            if (i >= 1) {
                dsdpSolver->ybound = (dsdpSolver->ybound == DSDP_INFINITY) ? 1e+05 * dsdpSolver->tau : dsdpSolver->ybound;
                DSDPSetIntParam(dsdpSolver, INT_PARAM_ACORRECTOR, 2);
                retcode = dInfeasCorrectorStep(dsdpSolver, TRUE); checkCode;
            }
            break;
        }
        
        // Factorize dual matrices
        retcode = setupFactorize(dsdpSolver); checkCode;

#ifdef compareMode
        assert( TRUE );
#endif
        
        // Set up Schur matrix and solve the system
        retcode = setupSchur(dsdpSolver);
        
        // Get proximity and check primal feasibility
        trymu = muprimal;
        retcode = dsdpgetPhaseAProxMeasure(dsdpSolver, trymu); checkCode;
        if (dsdpSolver->eventMonitor[EVENT_PFEAS_FOUND]) {
            muprimal = trymu;
        }
        trymu = (dsdpSolver->pObjVal - dsdpSolver->dObjVal -
                 dsdpSolver->Ry * 1e+10) / (5.0 * dsdpSolver->nall);
        
        
        if (dsdpSolver->Pnrm < 5) {
            dsdpSolver->mu = MAX(dsdpSolver->mu * 0.01, trymu);
        } else if (dsdpSolver->mu < 10) {
            dsdpSolver->mu = MAX(dsdpSolver->mu * 0.1, trymu);
        } else {
            dsdpSolver->mu = MAX(dsdpSolver->mu * 0.95, trymu);
        }
        
        
        dsdpSolver->iterProgress[ITER_PROX_POBJ] = TRUE;
        // Adaptive dual infeasibility
#ifdef compareMode
        assert( TRUE );
#endif
        retcode = computeAdaptivedRate(dsdpSolver); checkCode;
        // Setup directions
        retcode = getStepDirs(dsdpSolver); checkCode;
        dsdpSolver->iterProgress[ITER_STEP_DIRECTION] = TRUE;
        // Compute maximum available stepsize
        retcode = getMaxStep(dsdpSolver); checkCode;
        dsdpSolver->iterProgress[ITER_COMPUTE_STEP] = TRUE;
        // Compute residual
        retcode = setupRes(dsdpSolver); checkCode;
        dsdpSolver->iterProgress[ITER_RESIDUAL] = TRUE;
        // Take step
        retcode = takeStep(dsdpSolver); checkCode;
        dsdpSolver->iterProgress[ITER_TAKE_STEP] = TRUE;
#ifdef compareMode
        assert( TRUE );
#endif
        // Corrector
        retcode = dInfeasCorrectorStep(dsdpSolver, FALSE); checkCode;

        // Decrease mu with sufficient proximity
        checkIterProgress(dsdpSolver, ITER_DECREASE_MU);
        if (dsdpSolver->Pnrm < 0.1 &&
            !dsdpSolver->eventMonitor[EVENT_MU_QUALIFIES]) {
            // dsdpSolver->mu *= 0.1;
            // muprimal *= 0.1;
        }
        
        dsdpSolver->iterProgress[ITER_DECREASE_MU] = TRUE;
        checkIterProgress(dsdpSolver, ITER_NEXT_ITERATION);
        time = my_clock() - start;
    }
    
    dsdpSolver->mu = muprimal;
    
    DSDPStatUpdate(stat, STAT_PHASE_A_TIME, (double) time);
    DSDPStatUpdate(stat, STAT_PHASE_A_ITER, (double) i + 1);
    printPhaseASummary(dsdpSolver);
    
    return retcode;
}
