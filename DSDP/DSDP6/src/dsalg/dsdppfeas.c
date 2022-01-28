#include "dsdppfeas.h"
#include "dsdputils.h"
#include "dsdpinitializer.h"
#include "schurmat.h"
#include "residualsetup.h"
#include "stepdirection.h"
#include "stepheur.h"
#include "dsdpproxmeasure.h"
#include "dsdplog.h"
#include "dsdppfeascheck.h"

static char etype[] = "DSDP Primal Feasibility Certificate";

/* Implement DSDP Primal heuristic: primal feasibility phase */
extern DSDP_INT DSDPPFeasPhase( HSDSolver *dsdpSolver ) {
    
    /* Driver routine for primal feasibility/optimality certificate
     This routine applies the original DSDP dual scaling algorithm to solve the
     problem. This phase certificates primal infeasibility / primal-dual optimality.
     */
    
    DSDP_INT retcode = DSDP_RETCODE_OK;

    // Switch to phase B
    dsdpSolver->eventMonitor[EVENT_IN_PHASE_A] = FALSE;
    dsdpSolver->eventMonitor[EVENT_IN_PHASE_B] = TRUE;
    dsdpSolver->iterProgress[ITER_RESIDUAL] = TRUE;
    
    retcode = dsdpInitializeB(dsdpSolver);
    dsdpprintPhaseBheader();
    
    // Start dual scaling
    DSDP_INT stop = FALSE;
    double start = my_clock();
    
    double newmu = 0.0;
    double muub  = 0.0;
    double mulb  = 0.0;
    double time  = 0.0;
    
    for (DSDP_INT i = 0; ; ++i) {

        // Start iteration
        dsdpSolver->iterB = i;
        dsdpSolver->iterProgress[ITER_LOGGING] = FALSE;

        // Check NAN
        dsdpCheckNan(dsdpSolver);
        // Check Algorithm convergence
        DSDPCheckPhaseBConvergence(dsdpSolver, &stop);
        // Compute dual objective value
        if (i > 0) {
            dsdpSolver->iterProgress[ITER_DUAL_OBJ] = FALSE;
            retcode = getDualObj(dsdpSolver); checkCode;
        }
        
        // Logging
        DSDPPhaseBLogging(dsdpSolver);
        
        // Reset monitor
        DSDPResetPhaseBMonitor(dsdpSolver);
        
        if (stop) {
            break;
        }
        
        // Factorize dual matrices
        retcode = setupFactorize(dsdpSolver); checkCode;
        // Set up Schur matrix and solve the system
        retcode = setupSchur(dsdpSolver); checkCode;
        // Get proximity and check primal feasibility
        retcode = dsdpgetPhaseBProxMeasure(dsdpSolver, &muub, &mulb); checkCode;
        // Select new mu
        retcode = selectMu(dsdpSolver, &newmu); checkCode;
        newmu = MIN(newmu, muub);
        newmu = MAX(newmu, mulb);
        dsdpSolver->mu = newmu;
        
        if (dsdpSolver->Pnrm < 0.1) {
            dsdpSolver->mu *= 0.1;
        }
        
        dsdpSolver->iterProgress[ITER_PROX_POBJ] = TRUE;
        // Get step direction
        retcode = getStepDirs(dsdpSolver); checkCode;
        // Verify primal infeasibility
        retcode = dsdpCheckPrimalInfeas(dsdpSolver); checkCode;
        // Potential reduction
        retcode = dualPotentialReduction(dsdpSolver); checkCode;
        // Corrector step
        dsdpSolver->iterProgress[ITER_CORRECTOR] = TRUE;
        dsdpSolver->iterProgress[ITER_DECREASE_MU] = TRUE;
        
        checkIterProgress(dsdpSolver, ITER_NEXT_ITERATION);

        time =  my_clock() - start;
    }
    
    printPhaseBSummary(dsdpSolver, time);
    
    return retcode;
}