#include "dsdppfeas.h"
#include "dsdputils.h"
#include "dsdpinitializer.h"
#include "schurmat.h"
#include "residualsetup.h"
#include "stepdirection.h"
#include "stepheur.h"
#include "dsdpproxmeasure.h"
#include "dsdplog.h"

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
    
    retcode = dsdpInitializeB(dsdpSolver);
    dsdpprintPhaseBheader();
    
    // Start dual scaling
    double tol = dsdpSolver->param->absOptTol;
    DSDP_INT maxiter = dsdpSolver->param->BmaxIter;
    DSDP_INT stop = FALSE;
    clock_t start = clock();
    
    double muub = 0.0;
    double mulb = 0.0;
    
    for (DSDP_INT i = 0; ; ++i) {
        
        dsdpSolver->iterB = i;

        // Check NAN
        dsdpCheckNan(dsdpSolver);
        // Check Algorithm convergence
        DSDPCheckPhaseBConvergence(dsdpSolver, &stop);
        // Compute dual objective value
        
        if (i > 0) {
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
        retcode = setupSchur(dsdpSolver);
        // Get proximity and check primal feasibility
        retcode = dsdpgetPhaseBProxMeasure(dsdpSolver, &muub, &mulb);
        // Select new mu
        
        
        
    }
    
    return retcode;
}
