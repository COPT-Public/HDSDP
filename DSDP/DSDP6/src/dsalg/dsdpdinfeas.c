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
extern DSDP_INT DSDPCheckPhaseAConvergence( HSDSolver *dsdpSolver, DSDP_INT *isOK ) {
    /* Check convergence of DSDP */
    DSDP_INT retcode = DSDP_RETCODE_OK;
    
    *isOK = FALSE;
    
    if (dsdpSolver->Ry < dsdpSolver->param->absOptTol) {
        dsdpSolver->solStatus = DSDP_PD_FEASIBLE;
        *isOK = TRUE;
    }
    
    if (dsdpSolver->tau < 0.0001 * dsdpSolver->kappa &&
        dsdpSolver->mu < dsdpSolver->param->absOptTol) {
        dsdpSolver->solStatus = DSDP_PUNKNOWN_DFEAS;
        *isOK = TRUE;
    }

    if (dsdpSolver->iterA >= dsdpSolver->param->AmaxIter) {
        dsdpSolver->solStatus = DSDP_MAXITER;
        *isOK = TRUE;
    }
    
    if (dsdpSolver->alpha < 1e-05) {
        dsdpSolver->eventMonitor[EVENT_STEP_TOO_SMALL] = TRUE;
        dsdpSolver->solStatus = DSDP_UNKNOWN;
        *isOK = TRUE;
    }
    
    return retcode;
}

extern DSDP_INT DSDPDInfeasEliminator( HSDSolver *dsdpSolver ) {
    
    DSDP_INT retcode = DSDP_RETCODE_OK;
    DSDP_INT goOn = TRUE;
    double muprimal = dsdpSolver->param->initMu;
    double tol = dsdpSolver->param->absOptTol;
    double trymu = 0.0;
    double sigma = dsdpSolver->param->Asigma;
    double newmu[4] = {0.1, 0.3, 0.7, 1.0};
    
    /* Start Phase A algorithm */
    retcode = dsdpInitialize(dsdpSolver); checkCode;
    
    for (DSDP_INT i = 0;; ++i) {
        dsdpSolver->iterA = i;
        DSDPCheckPhaseAConvergence(dsdpSolver, &goOn);
        if (!goOn) {
            break;
        }
        
        retcode = setupPhaseASchur(dsdpSolver);
        
        for (DSDP_INT j = 0; j < 4; ++j) {
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
    return retcode;
}
