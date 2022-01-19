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

static void dsdpCheckNan( HSDSolver *dsdpSolver ) {
    
    // Check nan in some iterations
    if (isnan(dsdpSolver->y->x[0])) {
        dsdpSolver->eventMonitor[EVENT_NAN_IN_ITER] = TRUE;
        return;
    }
    
    if (isnan(dsdpSolver->tau)   || isnan(dsdpSolver->kappa) ||
        isnan(dsdpSolver->alpha) || isnan(dsdpSolver->Ry)    ||
        isnan(dsdpSolver->dObjVal|| isnan(dsdpSolver->pObjVal))) {
        dsdpSolver->eventMonitor[EVENT_NAN_IN_ITER] = TRUE;
    }
}

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
        
        // Compute dual objective
        retcode = getDualObj(dsdpSolver); checkCode;
        
        // Set up schur matrix
        retcode = setupPhaseASchur(dsdpSolver);
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
        
        retcode = getStepDirs(dsdpSolver); checkCode;
        retcode = getMaxStep(dsdpSolver); checkCode;
        retcode = takeStep(dsdpSolver); checkCode;
        retcode = setupRes(dsdpSolver); checkCode;
        
        if (dsdpSolver->Pnrm < 0.1) {
            dsdpSolver->mu *= 0.1;
            muprimal *= 0.1;
        }
        
        if (i >= 30) {
            sigma = 0.1;
            dsdpSolver->param->Aalpha = 0.2;
        }
    }
    
    time = (double) (clock() - start) / CLOCKS_PER_SEC;
    printPhaseASummary(dsdpSolver, time);
    
    return retcode;
}
