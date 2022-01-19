#include "dsdpdinfeas.h"
#include "dsdputils.h"
#include "dsdpinitializer.h"
#include "schurmat.h"
#include "residualsetup.h"
#include "stepdirection.h"
#include "stepheur.h"
#include "dsdpproxmeasure.h"
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

extern DSDP_INT DSDPDInfeasEliminator( HSDSolver *dsdpSolver ) {
    
    DSDP_INT retcode = DSDP_RETCODE_OK;
    DSDP_INT goOn = TRUE;
    
    // Initialize
    double muprimal = dsdpSolver->param->initMu;
    double tol      = dsdpSolver->param->absOptTol;
    double sigma    = dsdpSolver->param->Asigma;
    double attempt  = dsdpSolver->param->Aattempt;
    double trymu    = 0.0;
    double time     = 0.0;
    DSDP_INT agiter = MIN(dsdpSolver->param->AmaxIter, 30);
    DSDP_INT ntry   = 0;
    clock_t start   = clock();
    
    const double *newmu = NULL;
    
    if (attempt == DSDP_ATTEMPT_NO) {
        ntry = 1;
        newmu = nomu;
    } else if (attempt == DSDP_ATTEMPT_CONSV) {
        ntry = 2;
        newmu = consvmu;
    } else {
        ntry = 4;
        newmu = aggmu;
    }
        
    /* Start Phase A algorithm */
    dsdpshowdash();
    retcode = dsdpInitialize(dsdpSolver); checkCode;
    dsdpshowdash();
    
    /* Print algorithm header */
    dsdpprintPhaseAheader();
    
    for (DSDP_INT i = 0; ; ++i) {
        
        dsdpSolver->iterA = i;
        
        // Check NaN
        dsdpCheckNan(dsdpSolver);
        // Check algorithm convergence
        DSDPCheckPhaseAConvergence(dsdpSolver, &goOn);
        // Logging
        DSDPPhaseALogging(dsdpSolver);
        
        if (!goOn) {
            break;
        }
        
        // Reset monitor
        DSDPResetPhaseAMonitor(dsdpSolver);
        // Compute dual objective
        retcode = getDualObj(dsdpSolver); checkCode;
        // Factorize dual matrices
        retcode = setupFactorize(dsdpSolver); checkCode;
        // Set up schur matrix and solve the system
        retcode = setupPhaseASchur(dsdpSolver);
        // Get proximity and check primal feasibility
        for (DSDP_INT j = 0; j < ntry; ++j) {
            trymu = newmu[j] * muprimal;
            retcode = dsdpgetPhaseAProxMeasure(dsdpSolver, trymu); checkCode;
            if (dsdpSolver->eventMonitor[EVENT_PFEAS_FOUND]) {
                if (trymu > tol * tol) {
                    muprimal = trymu;
                }
                dsdpSolver->mu = MAX(muprimal, dsdpSolver->mu * sigma);
                break;
            }
        }
        dsdpSolver->iterProgress[ITER_PROX_POBJ] = TRUE;
        
        // Setup directions
        retcode = getStepDirs(dsdpSolver); checkCode;
        // Compute maximum available stepsize
        retcode = getMaxStep(dsdpSolver); checkCode;
        // Take step
        retcode = takeStep(dsdpSolver); checkCode;
        // Compute residual
        retcode = setupRes(dsdpSolver); checkCode;
        
        // Decrease mu with sufficient proximity
        if (dsdpSolver->Pnrm < 0.1) {
            dsdpSolver->mu *= 0.1;
            muprimal *= 0.1;
        }
        
        dsdpSolver->iterProgress[ITER_DECREASE_MU] = TRUE;
        
        // Be more aggressive if
        if (i == agiter) {
            sigma = 0.1;
            dsdpSolver->param->Aalpha = 0.2;
        }
    }
    
    time = (double) (clock() - start) / CLOCKS_PER_SEC;
    printPhaseASummary(dsdpSolver, time);
    
    return retcode;
}
