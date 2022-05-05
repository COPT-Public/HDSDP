#include "dsdppfeas.h"
#include "dsdputils.h"
#include "dsdpinitializer.h"
#include "schurmatops.h"
#include "residualsetup.h"
#include "stepdirection.h"
#include "stepheur.h"
#include "dsdpproxmeasure.h"
#include "dsdplog.h"
#include "dsdppfeascheck.h"
#include "dsdpcorrector.h"

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
    DSDP_INT stop = FALSE, usegold = FALSE;
    double start = my_clock();
    double newmu = 0.0, muub = 0.0, mulb = 0.0, time = 0.0;
    double rhon, tol, sigma, approxpObj;

    DSDPStats *stat = &dsdpSolver->dsdpStats;
    
    DSDPGetIntParam(dsdpSolver, INT_PARAM_GOLDSEARCH, &usegold);
    DSDPGetDblParam(dsdpSolver, DBL_PARAM_RHON, &rhon);
    DSDPGetDblParam(dsdpSolver, DBL_PARAM_ABS_OPTTOL, &tol); tol *= 0.1;
    DSDPGetDblParam(dsdpSolver, DBL_PARAM_BSIGMA, &sigma);
    
    DSDP_INT i;
    for (i = 0; ; ++i) {
        
        if (i > 300) { dsdpSolver->dperturb = 1e-07; }
        if (i > 400) { dsdpSolver->dperturb = 1e-06; }
        
        // Start iteration
        DSDPStatUpdate(&dsdpSolver->dsdpStats, STAT_PHASE_B_ITER, (double) i);
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
        
        // Decrease primal objective via golden search
        if (usegold && dsdpSolver->mu < 1e+07) {
            retcode = searchpObj(dsdpSolver, &approxpObj);
            approxpObj = MAX(dsdpSolver->pObjVal * 0.8, approxpObj);
            dsdpSolver->pObjVal = MIN(dsdpSolver->pObjVal, approxpObj);
        }

        // Get proximity and check primal feasibility
        retcode = dsdpgetPhaseBProxMeasure(dsdpSolver, &muub, &mulb); checkCode;
        
        // Select new mu
        if (dsdpSolver->mu > 1e-10 && i <= 480) {
            retcode = selectMu(dsdpSolver, &newmu); checkCode;
            newmu = MIN(newmu, muub); newmu = MAX(newmu, mulb);
            dsdpSolver->mu = newmu;
        } else {
            dsdpSolver->mu = MAX(dsdpSolver->mu * sigma, 1e-12);
        }

        DSDPSetDblParam(dsdpSolver, DBL_PARAM_RHO,
                        (dsdpSolver->pObjVal - dsdpSolver->dObjVal) / dsdpSolver->mu);
        
        dsdpSolver->iterProgress[ITER_PROX_POBJ] = TRUE;
        // Get step direction
        retcode = getStepDirs(dsdpSolver); checkCode;
        // Verify primal infeasibility
        retcode = dsdpCheckPrimalInfeas(dsdpSolver); checkCode;
        // Potential reduction
        retcode = dualPotentialReduction(dsdpSolver); checkCode;
        // Corrector step
        retcode = dualCorrectorStep(dsdpSolver); checkCode;
        // dsdpSolver->iterProgress[ITER_CORRECTOR] = TRUE;
        dsdpSolver->iterProgress[ITER_DECREASE_MU] = TRUE;
        checkIterProgress(dsdpSolver, ITER_NEXT_ITERATION);
        time = my_clock() - start;
    }
    
    DSDPStatUpdate(stat, STAT_PHASE_B_TIME, time);
    DSDPStatUpdate(stat, STAT_PHASE_B_ITER, i + 1);
    printPhaseBSummary(dsdpSolver);
    
    return retcode;
}
