#include <stdio.h>
#include "dsdplog.h"
#include "hsd.h"
#define dsdplog printf

#ifndef isnan
#define isnan(x) ((x) != (x))
#endif

#ifndef EVENT
#define E(x) dsdpSolver->eventMonitor(x);
#endif

static void dsdpstatus( DSDP_INT code, char *word ) {
    // Convert status code
    switch (code) {
        case DSDP_UNKNOWN:
            strcpy(word, "DSDP_UKNOWN");
            break;
        case DSDP_OPTIMAL:
            strcpy(word, "DSDP_OPTIMAL");
            break;
        case DSDP_MAXITER:
            strcpy(word, "DSDP_MAXITER");
            break;
        case DSDP_INTERNAL_ERROR:
            strcpy(word, "DSDP_INTERNAL_ERROR");
            break;
        case DSDP_PD_FEASIBLE:
            strcpy(word, "DSDP_PRIMAL_DUAL_FEASIBLE");
            break;
        case DSDP_PFEAS_DINFEAS:
            strcpy(word, "DSDP_PFEASIBLE_DINFEASIBLE");
            break;
        case DSDP_PUNKNOWN_DFEAS:
            strcpy(word, "DSDP_PUNKNOWN_DFEASIBLE");
            break;
        case DSDP_PUNBOUND_DINFEAS:
            strcpy(word, "DSDP_PUNBOUND_DINFEASIBLE");
            break;
        default:
            break;
    }
}

extern void dsdpshowdash(void) {
    printf("_______________________________________"
           "_______________________________________"
           "______________________\n");
}

extern void dsdpCheckNan( HSDSolver *dsdpSolver ) {
    
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

extern void dsdpprintPhaseAheader(void) {
    
    dsdpshowdash();
    /*      ________________________________________________________________________________*/
    /*      | niter |   pObj   |   dObj   |  dInf  |  k/t  |  mu  |  alpha  |  pNrm  |  E  |*/
    /*      --------------------------------------------------------------------------------*/
    printf("| %4s | %12s | %12s | %8s | %8s | %8s | %8s | %8s | %3s |\n",
            "iter", "pObj", "dObj", "dInf", "k/t", "mu", "alpha", "pNrm", "E");
    printf("_______________________________________"
           "_______________________________________"
           "______________________\n");
}

extern void DSDPResetPhaseAMonitor( HSDSolver *dsdpSolver ) {
    
    // Reset event and iteration monitor at the end of each iteration
    DSDP_INT logged = dsdpSolver->iterProgress[ITER_LOGGING];
    memset(dsdpSolver->iterProgress, 0, sizeof(DSDP_INT) * IterStep);
    dsdpSolver->iterProgress[ITER_INITIALIZE] = TRUE;
    dsdpSolver->iterProgress[ITER_LOGGING] = logged;
    
    DSDP_INT *monitor = dsdpSolver->eventMonitor;
    assert( monitor[EVENT_IN_PHASE_A] );
    monitor[EVENT_PFEAS_FOUND] = FALSE;
    monitor[EVENT_LARGE_NORM] = FALSE;
    monitor[EVENT_SMALL_STEP] = FALSE;
}

extern DSDP_INT DSDPPhaseALogging( HSDSolver *dsdpSolver ) {
    
    /*  _______________________________________________________________________________*/
    /*  | niter |   pObj   |   dObj   |  dInf  |  k/t  |  mu  |  alpha  |  pNrm  |  E  |*/
    /*  -------------------------------------------------------------------------------*/
    
    // Implement the logging procedure of DSDP phase A
    DSDP_INT retcode = DSDP_RETCODE_OK;
    retcode = checkIterProgress(dsdpSolver, ITER_LOGGING);
    
    char event[3] = " ";
    double nRy = fabs(dsdpSolver->Ry) * sqrt(dsdpSolver->n);
    double kovert = dsdpSolver->kappa / dsdpSolver->tau;
    DSDP_INT *moniter = dsdpSolver->eventMonitor;
    
    assert( moniter[EVENT_IN_PHASE_A] );
    
    if (moniter[EVENT_PFEAS_FOUND]) {
        strcpy(&event[0], "H");
    }
    
    if (moniter[EVENT_NO_RY]) {
        strcpy(&event[0], "*");
        goto print_log;
    }
    
    if (moniter[EVENT_MU_QUALIFIES] &&
        moniter[EVENT_KT_QUALIFIES]) {
        strcpy(&event[0], "I");
        goto print_log;
    }
    
    if (dsdpSolver->smallIter == 3 ||
        moniter[EVENT_NAN_IN_ITER] ||
        moniter[EVENT_BAD_SCHUR]) {
        strcpy(&event[0], "F");
        goto print_log;
    }
    
    if (moniter[EVENT_LARGE_NORM]) {
        strcpy(&event[1], "N");
    }
    
    if (moniter[EVENT_MAX_ITERATION]) {
        strcpy(&event[1], "M");
    }
    
print_log:
    printf("| %4d | %12.3e | %12.3e | %8.2e | %8.2e | %8.2e | %8.2e | %8.2e | %3s |\n",
            dsdpSolver->iterA + 1, dsdpSolver->pObjVal, dsdpSolver->dObjVal,
            nRy, kovert, dsdpSolver->mu, dsdpSolver->alpha, dsdpSolver->Pnrm,
            event);
    
    dsdpSolver->iterProgress[ITER_LOGGING] = TRUE;
    return retcode;
}

extern DSDP_INT DSDPCheckPhaseAConvergence( HSDSolver *dsdpSolver, DSDP_INT *isOK ) {
    /* Check convergence of DSDP */
    DSDP_INT retcode = DSDP_RETCODE_OK;
    DSDP_INT *monitor = dsdpSolver->eventMonitor;
    *isOK = FALSE;
    
    // Dual infeasibility eliminated
    if (fabs(dsdpSolver->Ry) * sqrt(dsdpSolver->n) < dsdpSolver->param->absOptTol * 0.1) {
        monitor[EVENT_NO_RY] = TRUE;
        
        if (dsdpSolver->pObjVal != dsdpSolver->param->initpObj) {
            dsdpSolver->solStatus = DSDP_PD_FEASIBLE;
        } else {
            dsdpSolver->solStatus = DSDP_PUNKNOWN_DFEAS;
        }
        
        *isOK = TRUE;
        return retcode;
    }
    
    // Dual infeasibility certificate found
    if (dsdpSolver->tau < 0.001 * dsdpSolver->kappa) {
        monitor[EVENT_KT_QUALIFIES] = TRUE;
    }
    
    if (dsdpSolver->mu < dsdpSolver->param->absOptTol) {
        monitor[EVENT_MU_QUALIFIES] = TRUE;
    }
    
    if (monitor[EVENT_KT_QUALIFIES] &&
        monitor[EVENT_MU_QUALIFIES]) {
        if (dsdpSolver->pObjVal != dsdpSolver->param->initpObj) {
            dsdpSolver->solStatus = DSDP_PFEAS_DINFEAS;
        } else {
            dsdpSolver->solStatus = DSDP_PUNBOUND_DINFEAS;
        }
        *isOK = TRUE;
        return retcode;
    }

    // Maximum iteration
    if (dsdpSolver->iterA >= dsdpSolver->param->AmaxIter) {
        monitor[EVENT_MAX_ITERATION] = TRUE;
        dsdpSolver->solStatus = DSDP_MAXITER;
        *isOK = TRUE;
        return retcode;
    }
    
    // Small step
    if (dsdpSolver->iterA > 0 && dsdpSolver->alpha < 1e-04 &&
        !monitor[EVENT_PFEAS_FOUND]) {
        monitor[EVENT_SMALL_STEP] = TRUE;
        dsdpSolver->smallIter += 1;
    }
    
    // NAN in iteration
    if (monitor[EVENT_NAN_IN_ITER]) {
        dsdpSolver->solStatus = DSDP_INTERNAL_ERROR;
        *isOK = TRUE;
        return retcode;
    }
    
    if (monitor[EVENT_BAD_SCHUR]) {
        dsdpSolver->solStatus = DSDP_INTERNAL_ERROR;
        *isOK = TRUE;
        return retcode;
    }

    return retcode;
}

extern DSDP_INT printPhaseASummary( HSDSolver *dsdpSolver, double time ) {
    // Summarize Phase A iteration
    DSDP_INT retcode = DSDP_RETCODE_OK;
    char sAlog[100] = "DSDP Phase A ends with status: ";
    dsdpshowdash();
    dsdpstatus(dsdpSolver->solStatus, &sAlog[31]);
    printf("%s \n", sAlog);
    printf("Elapsed Time: %5f seconds \n", time);
    dsdpshowdash();
    return retcode;
}
