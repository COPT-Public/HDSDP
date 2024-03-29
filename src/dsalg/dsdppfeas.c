#include "dsdppfeas.h"
#include "dsdplapack.h"
#include "dsdputils.h"
#include "dsdpinitializer.h"
#include "schurmatops.h"
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
    
    dsdpInitializeB(dsdpSolver);
    printPhaseBheader();
    
    // Start dual scaling
    DSDP_INT stop = FALSE, usegold = FALSE, nopfeasIter = 0;
    double start = my_clock();
    double newmu = 0.0, muub = 0.0, mulb = 0.0, time = 0.0;
    double rhon, tol, sigma = 0.7, initpObj = dsdpSolver->pObjVal;

    DSDPStats *stat = &dsdpSolver->dsdpStats;
    DSDPGetIntParam(dsdpSolver, INT_PARAM_GOLDSEARCH, &usegold);
    DSDPGetDblParam(dsdpSolver, DBL_PARAM_RHON, &rhon);
    DSDPGetDblParam(dsdpSolver, DBL_PARAM_ABS_OPTTOL, &tol); tol *= 0.1;
    
    DSDP_INT i;
    for (i = 0; ; ++i) {
        
        if (i == 350) { dsdpSolver->dperturb += 1e-06; }
        if (i == 400) { dsdpSolver->dperturb += 1e-05; }
        if (i == 450) { dsdpSolver->dperturb += 1e-04; }
        if (i >= 480) { dsdpSolver->dperturb += 5e-05; }
        
        // Start iteration
        DSDPStatUpdate(&dsdpSolver->dsdpStats, STAT_PHASE_B_ITER, (double) i);
        dsdpSolver->iterProgress[ITER_LOGGING] = FALSE;
        // Check NAN
        checkNan(dsdpSolver);
        // Check Algorithm convergence
        checkPhaseBConvergence(dsdpSolver, &stop);
        // Compute dual objective value
        dsdpSolver->iterProgress[ITER_DUAL_OBJ] = FALSE;
        DSDPConic( COPS_GET_DOBJ )(dsdpSolver);
        // Logging
        phaseBLogging(dsdpSolver);
        // Reset monitor
        resetPhaseBMonitor(dsdpSolver);
        if (stop) break;
        
        // Factorize dual matrices
        dsdpSolver->iterProgress[ITER_DUAL_FACTORIZE] = TRUE;
        // setupFactorize(dsdpSolver);
        // Set up Schur matrix and solve the system
        
        getBslack(dsdpSolver, dsdpSolver->y, DUALVAR);
        setupSchur(dsdpSolver);
        
        dsdpgetPhaseBProxMeasure(dsdpSolver, &muub, &mulb); checkCode;
        if (dsdpSolver->eventMonitor[EVENT_PFEAS_FOUND]) {
            nopfeasIter = 0;
        } else {
            nopfeasIter += 1;
        }
        
        // Select new mu
        if ((dsdpSolver->mu > 1e-12 && i <= 480 && nopfeasIter < 10) || dsdpSolver->pObjVal == initpObj) {
            selectMu(dsdpSolver, &newmu); checkCode;
            newmu = MIN(newmu, muub); newmu = MAX(newmu, mulb);
            dsdpSolver->mu = newmu;
            if (nopfeasIter >= 8) {
                dsdpSolver->mu *= 100;
            }
        } else {
            if (dsdpSolver->mu < 1e-13 || nopfeasIter >= 15) {
                printf("| Hard to find a new primal solution. Give up. \n"); break;
            }
            dsdpSolver->mu = MAX(dsdpSolver->mu * sigma, 1e-13);
        }
        DSDPSetDblParam(dsdpSolver, DBL_PARAM_RHO,
                        (dsdpSolver->pObjVal - dsdpSolver->dObjVal) / dsdpSolver->mu);
        
        dsdpSolver->iterProgress[ITER_PROX_POBJ] = TRUE;
        // Get step direction
        getStepDirs(dsdpSolver);
        // Verify primal infeasibility
        dsdpCheckPrimalInfeas(dsdpSolver);
        // Potential reduction
        dualPotentialReduction(dsdpSolver);
        // Corrector step
        dualCorrectorStep(dsdpSolver);
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
