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
extern DSDP_INT DSDPDInfeasEliminator( HSDSolver *dsdpSolver ) {
    
    DSDP_INT retcode = DSDP_RETCODE_OK;
    DSDP_INT goOn = TRUE;
    DSDPStats *stat = &dsdpSolver->dsdpStats;
    dsdpSolver->pObjVal = DSDP_INFINITY;
    
    // Initialize
    double muprimal, tol, sigma, initpObj;
    double trymu = 0.0, time = 0.0, start = my_clock();
    DSDP_INT attempt;
    
    DSDPGetDblParam(dsdpSolver, DBL_PARAM_INIT_MU,    &muprimal);
    DSDPGetDblParam(dsdpSolver, DBL_PARAM_ABS_OPTTOL, &tol     );
    DSDPGetDblParam(dsdpSolver, DBL_PARAM_ASIGMA,     &sigma   );
    DSDPGetIntParam(dsdpSolver, INT_PARAM_AATTEMPT,   &attempt );
    dsdpSolver->eventMonitor[EVENT_IN_PHASE_A] = TRUE;
    dsdpSolver->eventMonitor[EVENT_IN_PHASE_B] = FALSE;
    
    /* Start Phase A algorithm */
    dsdpshowdash();
    dsdpInitializeA(dsdpSolver);
    initpObj = dsdpSolver->pObjVal;
    
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
        DSDPConic( COPS_GET_DOBJ )(dsdpSolver);
        // Logging
        DSDPPhaseALogging(dsdpSolver);
        // Reset monitor
        DSDPResetPhaseAMonitor(dsdpSolver);
        
        if (i > 100) {
            dsdpSolver->eventMonitor[EVENT_HSD_UPDATE] = TRUE;
        }
        
        if (goOn) {
//            if (i >= 1) {
//                DSDPSetIntParam(dsdpSolver, INT_PARAM_ACORRECTOR, 2);
//                dInfeasCorrectorStep(dsdpSolver, TRUE);
//            }
            break;
        }
        
        // Factorize dual matrices
        setupFactorize(dsdpSolver);
        // Set up Schur matrix and solve the system
        setupSchur(dsdpSolver);
        // Get proximity and check primal feasibility
        trymu = muprimal;
        dsdpgetPhaseAProxMeasure(dsdpSolver, trymu);
        if (dsdpSolver->eventMonitor[EVENT_PFEAS_FOUND]) {
            muprimal = trymu;
        }
        trymu = (dsdpSolver->pObjVal - dsdpSolver->dObjVal -
                 dsdpSolver->Ry * 1e+08) / (5.0 * dsdpSolver->nall);
        
        if (dsdpSolver->Pnrm < 5) {
            dsdpSolver->mu = MAX(dsdpSolver->mu * 0.01, trymu * 0.1);
        } else if (dsdpSolver->mu < 10) {
            dsdpSolver->mu = MAX(dsdpSolver->mu * 0.1, trymu * 0.8);
        } else {
            dsdpSolver->mu = MAX(dsdpSolver->mu * 0.95, trymu);
        }
        
        dsdpSolver->iterProgress[ITER_PROX_POBJ] = TRUE;
        // Adaptive dual infeasibility
        computeAdaptivedRate(dsdpSolver);
        // Setup directions
        getStepDirs(dsdpSolver);
        // Compute maximum available stepsize
        getMaxStep(dsdpSolver);
        // Compute residual
        setupRes(dsdpSolver);
        if (i == 3 && fabs(dsdpSolver->pObjVal - initpObj) < 1e-10) {
            dsdpSolver->Ry = - MIN(fabs(dsdpSolver->Ry) * 1e+10, 1e+15);
        }
        // Take step
        takeStep(dsdpSolver);
        dsdpSolver->iterProgress[ITER_TAKE_STEP] = TRUE;
        // Corrector
        dInfeasCorrectorStep(dsdpSolver, FALSE);
        // Decrease mu with sufficient proximity
        checkIterProgress(dsdpSolver, ITER_DECREASE_MU);
        if (dsdpSolver->Pnrm < 0.1 &&
            !dsdpSolver->eventMonitor[EVENT_MU_QUALIFIES]) {
            dsdpSolver->mu *= 0.1;
        }
        
        dsdpSolver->iterProgress[ITER_DECREASE_MU] = TRUE;
        checkIterProgress(dsdpSolver, ITER_NEXT_ITERATION);
        time = my_clock() - start;
    }
    
    DSDPStatUpdate(stat, STAT_PHASE_A_TIME, (double) time);
    DSDPStatUpdate(stat, STAT_PHASE_A_ITER, (double) i + 1);
    printPhaseASummary(dsdpSolver);
    
    return retcode;
}
