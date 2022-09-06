#include <sys/time.h>
#include "dsdplapack.h"
#include "dsdplog.h"

#ifndef isnan
#define isnan(x) ((x) != (x))
#endif

static char etype[] = "DSDP logging";

static void convertStatusCode( DSDP_INT code, char *word ) {
    // Convert status code
    switch (code) {
        case DSDP_UNKNOWN         : strcpy(word, "DSDP_UKNOWN"); break;
        case DSDP_OPTIMAL         : strcpy(word, "DSDP_OPTIMAL"); break;
        case DSDP_MAXITER         : strcpy(word, "DSDP_MAXITER"); break;
        case DSDP_INTERNAL_ERROR  : strcpy(word, "DSDP_INTERNAL_ERROR"); break;
        case DSDP_PD_FEASIBLE     : strcpy(word, "DSDP_PRIMAL_DUAL_FEASIBLE"); break;
        case DSDP_PFEAS_DINFEAS   : strcpy(word, "DSDP_PFEASIBLE_DINFEASIBLE"); break;
        case DSDP_PUNKNOWN_DFEAS  : strcpy(word, "DSDP_PUNKNOWN_DFEASIBLE"); break;
        case DSDP_PUNKNOWN_DINFEAS: strcpy(word, "DSDP_PUNKNOWN_DINFEASIBLE"); break;
        case DSDP_PINFEAS_DFEAS   : strcpy(word, "DSDP_PINFEASIBLE_DFEASIBLE"); break;
        case DSDP_TIMELIMIT       : strcpy(word, "DSDP_TIMELIMIT"); break;
        default: break;
    }
}

extern double my_clock( void ) {
    
    struct timeval t;
    gettimeofday(&t, NULL);
    return (1.0e-6 * t.tv_usec + t.tv_sec);
}

extern void showBeautifulDashlines( void ) {
    
    printf("---------------------------------------"
           "---------------------------------------"
           "--------------------\n");
}

extern DSDP_INT checkIterProgress( HSDSolver *dsdpSolver, DSDP_INT iter ) {
    
    DSDP_INT retcode = DSDP_RETCODE_OK, npassed = 0, i;
    for (i = 0; i < iter; ++i) { npassed += dsdpSolver->iterProgress[i]; }
    if (npassed < iter) { error(etype, "Pre-requisite steps are not completed. \n"); }
    return retcode;
}

extern void checkNan( HSDSolver *dsdpSolver ) {
    // Check nan in some iterations
    if (isnan(dsdpSolver->tau)  || isnan(dsdpSolver->kappa)   || isnan(dsdpSolver->alpha)   ||
        isnan(dsdpSolver->Pnrm) || isnan(dsdpSolver->dObjVal) || isnan(dsdpSolver->pObjVal) ||
        isnan(dsdpSolver->y->x[0])) {
        dsdpSolver->eventMonitor[EVENT_NAN_IN_ITER] = TRUE;
    }
}

extern void printPhaseAheader( void ) {
    
    printf("| DSDP Phase A starts. Eliminating dual infeasibility %45s\n", "");
    showBeautifulDashlines();
    printf("\nPhase A Log: 'P': Primal solution found. '*': Phase A ends. 'F': Error happens. 'M': Max iteration\n");
    showBeautifulDashlines();
    printf("| %4s | %12s | %12s | %9s | %8s | %8s | %6s | %8s | %3s |\n",
            "Iter", "pObj", "dObj", "dInf", "k/t", "mu", "step", "Pnorm", "E");
    showBeautifulDashlines();
}

extern void printPhaseBheader( void ) {
    
    showBeautifulDashlines();
    printf("\nPhase B Log: 'P': Primal solution found. '*': Optimal. 'F': Error happens. 'M': Max iteration\n");
    showBeautifulDashlines();
    printf("| %4s | %17s | %17s | %10s | %8s | %6s | %8s | %3s |\n",
            "Iter", "pObj", "dObj", "pInf", "mu", "step", "Pnorm", "E");
    showBeautifulDashlines();
}

extern void resetPhaseAMonitor( HSDSolver *dsdpSolver ) {
    // Reset event and iteration monitor in Phase A
    DSDP_INT islogged = dsdpSolver->iterProgress[ITER_LOGGING];
    DSDP_INT isdobj = dsdpSolver->iterProgress[ITER_DUAL_OBJ];
    memset(dsdpSolver->iterProgress, 0, sizeof(DSDP_INT) * IterStep);
    dsdpSolver->iterProgress[ITER_INITIALIZE  ] = TRUE;
    dsdpSolver->iterProgress[ITER_LOGGING     ] = islogged;
    dsdpSolver->iterProgress[ITER_DUAL_OBJ    ] = isdobj;
    dsdpSolver->eventMonitor[EVENT_PFEAS_FOUND] = FALSE;
    dsdpSolver->eventMonitor[EVENT_LARGE_NORM ] = FALSE;
    dsdpSolver->eventMonitor[EVENT_SMALL_STEP ] = FALSE;
}

extern void resetPhaseBMonitor( HSDSolver *dsdpSolver ) {
    // Reset event and iteration monitor in Phase B
    DSDP_INT islogged = dsdpSolver->iterProgress[ITER_LOGGING];
    DSDP_INT isdobj = dsdpSolver->iterProgress[ITER_DUAL_OBJ];
    memset(dsdpSolver->iterProgress, 0, sizeof(DSDP_INT) * IterStep);
    dsdpSolver->iterProgress[ITER_INITIALIZE  ] = TRUE;
    dsdpSolver->iterProgress[ITER_RESIDUAL    ] = TRUE;
    dsdpSolver->iterProgress[ITER_LOGGING     ] = islogged;
    dsdpSolver->iterProgress[ITER_DUAL_OBJ    ] = isdobj;
    dsdpSolver->eventMonitor[EVENT_PFEAS_FOUND] = FALSE;
    dsdpSolver->eventMonitor[EVENT_LARGE_NORM ] = FALSE;
    dsdpSolver->eventMonitor[EVENT_SMALL_STEP ] = FALSE;
}

extern void phaseALogging( HSDSolver *dsdpSolver ) {
    // Implement the logging procedure of DSDP phase A
    char event[3] = " ";
    double nRy = fabs(dsdpSolver->Ry) * sqrt(dsdpSolver->n + dsdpSolver->lpDim);
    double tau = dsdpSolver->tau, kovert = dsdpSolver->kappa / tau, statval, pobj;
    checkIterProgress(dsdpSolver, ITER_LOGGING);
    DSDP_INT *moniter = dsdpSolver->eventMonitor;
    
    pobj = (dsdpSolver->pObjVal == DSDP_INFINITY) ? \
                DSDP_INFINITY : dsdpSolver->pObjVal * dsdpSolver->cScaler;
    // Different events
    if (moniter[EVENT_PFEAS_FOUND]) { strcpy(&event[0], "P"); }
    if (moniter[EVENT_NO_RY]) { strcpy(&event[0], "*"); goto print_log; }
    if (moniter[EVENT_MU_QUALIFIES] && moniter[EVENT_KT_QUALIFIES]) {
        strcpy(&event[0], "I"); goto print_log;
    }
    DSDPGetStats(&dsdpSolver->dsdpStats, STAT_NUM_SMALL_ITER, &statval);
    if ((DSDP_INT) statval >= 3 || moniter[EVENT_NAN_IN_ITER] || moniter[EVENT_BAD_SCHUR]) {
        strcpy(&event[0], "F"); goto print_log;
    }
    if (moniter[EVENT_LARGE_NORM]) { strcpy(&event[1], "N"); }
    if (moniter[EVENT_MAX_ITERATION]) { strcpy(&event[1], "M"); }
    
print_log:
    DSDPGetStats(&dsdpSolver->dsdpStats, STAT_PHASE_A_ITER, &statval);
    if ((DSDP_INT) statval % 1 == 0) {
        /*
                Phase A Log: 'P': Primal solution found. '*': Phase A ends. 'F': Error happens. 'M': Max iteration
                --------------------------------------------------------------------------------------------------
                | Iter |         pObj |         dObj |      dInf |      k/t |       mu |   step |    Pnorm |   E |
                --------------------------------------------------------------------------------------------------
        */
        printf("| %4d | %12.5e | %12.5e | %9.2e | %8.2e | %8.2e | %6.2f | %8.1e | %3s |\n",
                (DSDP_INT) statval + 1, pobj, dsdpSolver->dObjVal / tau * dsdpSolver->cScaler,
                nRy / tau, kovert, dsdpSolver->mu, dsdpSolver->alpha, dsdpSolver->Pnrm, event);
    }
    dsdpSolver->iterProgress[ITER_LOGGING] = TRUE;
}

extern void phaseBLogging( HSDSolver *dsdpSolver ) {
    // Implement the logging procedure of DSDP phase A
    checkIterProgress(dsdpSolver, ITER_LOGGING);
    char event[3] = " "; DSDP_INT *moniter = dsdpSolver->eventMonitor; double statval;
    if (moniter[EVENT_PFEAS_FOUND]) { strcpy(&event[0], "P"); }
    if (dsdpSolver->solStatus == DSDP_OPTIMAL) { strcpy(&event[0], "*"); goto print_log; }
    if (dsdpSolver->solStatus == DSDP_PINFEAS_DFEAS) { strcpy(&event[0], "I"); goto print_log; }
    DSDPGetStats(&dsdpSolver->dsdpStats, STAT_NUM_SMALL_ITER, &statval);
    if (statval >= 2 || moniter[EVENT_NAN_IN_ITER] || moniter[EVENT_BAD_SCHUR] ||
        dsdpSolver->solStatus == DSDP_INTERNAL_ERROR) {
        strcpy(&event[0], "F"); goto print_log;
    }
    if (moniter[EVENT_LARGE_NORM]) { strcpy(&event[1], "N"); }
    if (moniter[EVENT_MAX_ITERATION]) { strcpy(&event[1], "M"); }
    if (moniter[EVENT_INVALID_GAP]) { strcpy(&event[1], "!"); }
    
print_log:
    DSDPGetStats(&dsdpSolver->dsdpStats, STAT_PHASE_B_ITER, &statval);
    if ((DSDP_INT) statval % 1 == 0) {
        /*
                Phase B Log: 'P': Primal solution found. '*': Optimal. 'F': Error happens. 'M': Max iteration
                --------------------------------------------------------------------------------------------------
                | Iter |              pObj |              dObj |       pInf |       mu |   step |    Pnorm |   E |
                --------------------------------------------------------------------------------------------------
         */
        printf("| %4d | %17.10e | %17.10e | %10.3e | %8.2e | %6.2f | %8.1e | %3s |\n",
                (DSDP_INT) statval + 1, dsdpSolver->pObjVal * dsdpSolver->cScaler,
                dsdpSolver->dObjVal * dsdpSolver->cScaler, dsdpSolver->pinfeas, dsdpSolver->mu,
                dsdpSolver->alpha, dsdpSolver->Pnrm, event);
    }
    dsdpSolver->iterProgress[ITER_LOGGING] = TRUE; return;
}

extern void checkPhaseAConvergence( HSDSolver *dsdpSolver, DSDP_INT *isOK ) {
    // Check convergence of DSDP
    DSDP_INT *monitor = dsdpSolver->eventMonitor;
    
    DSDP_INT AmaxIter; *isOK = FALSE;
    double absOptTol, absFeasTol, statval1, statval2, initpObj, tmax;
    DSDPStats *stat = &dsdpSolver->dsdpStats;
    
    DSDPGetDblParam(dsdpSolver, DBL_PARAM_ABS_OPTTOL, &absOptTol);
    DSDPGetDblParam(dsdpSolver, DBL_PARAM_ABS_FEASTOL, &absFeasTol);
    DSDPGetDblParam(dsdpSolver, DBL_PARAM_INIT_POBJ, &initpObj);
    DSDPGetDblParam(dsdpSolver, DBL_PARAM_TIMELIMIT, &tmax);
    DSDPGetIntParam(dsdpSolver, INT_PARAM_AMAXITER, &AmaxIter);
    
    // Dual infeasibility eliminated
    if (fabs(dsdpSolver->Ry) < absFeasTol * dsdpSolver->tau) {
        monitor[EVENT_NO_RY] = TRUE;
        if (dsdpSolver->pObjVal != initpObj) {
            dsdpSolver->solStatus = DSDP_PD_FEASIBLE;
        } else {
            dsdpSolver->solStatus = DSDP_PUNKNOWN_DFEAS;
        }
        *isOK = TRUE; return;
    }
    
    // Dual infeasibility certificate found
    if (dsdpSolver->tau < 0.001 * dsdpSolver->kappa) { monitor[EVENT_KT_QUALIFIES] = TRUE; }
    if (dsdpSolver->mu < absOptTol) { monitor[EVENT_MU_QUALIFIES] = TRUE; }
    if (monitor[EVENT_KT_QUALIFIES] && monitor[EVENT_MU_QUALIFIES]) {
        if (dsdpSolver->pObjVal != DSDP_INFINITY) {
            dsdpSolver->solStatus = DSDP_PFEAS_DINFEAS;
        } else {
            dsdpSolver->solStatus = DSDP_PUNKNOWN_DINFEAS;
        }
        *isOK = TRUE; return;
    }
    if (fabs(dsdpSolver->Ry) * sqrt(dsdpSolver->n) < absOptTol && dsdpSolver->mu < 1e-05) {
        dsdpSolver->solStatus = DSDP_OPTIMAL;
        *isOK = TRUE; return;
    }
    // Maximum iteration
    DSDPGetStats(stat, STAT_PHASE_A_ITER, &statval1);
    if (statval1 >= AmaxIter) {
        monitor[EVENT_MAX_ITERATION] = TRUE;
        dsdpSolver->solStatus = DSDP_MAXITER;
        *isOK = TRUE; return;
    }
    // Small step
    DSDPGetStats(stat, STAT_PHASE_A_ITER, &statval1);
    DSDPGetStats(stat, STAT_NUM_SMALL_ITER, &statval2);
    if (statval1 > 0 && dsdpSolver->alpha < 3e-02 && fabs(dsdpSolver->Ry) < 1e-06 ) {
        monitor[EVENT_SMALL_STEP] = TRUE;
        statval2 += 1.0; DSDPStatUpdate(stat, STAT_NUM_SMALL_ITER, statval2);
        if (statval2 >= 2) { *isOK = TRUE; }
    }
    // NAN in iteration
    if (monitor[EVENT_NAN_IN_ITER]) {
        dsdpSolver->solStatus = DSDP_INTERNAL_ERROR;
        *isOK = TRUE; return;
    }
    if (monitor[EVENT_BAD_SCHUR]) {
        dsdpSolver->solStatus = DSDP_INTERNAL_ERROR;
        *isOK = TRUE; return;
    }
    if (my_clock() - dsdpSolver->startTime > tmax) {
        dsdpSolver->solStatus = DSDP_TIMELIMIT;
        *isOK = TRUE; return;
    }
}

extern void printPhaseASummary( HSDSolver *dsdpSolver ) {
    // Summarize Phase A iteration
    char sAlog[80] = "| DSDP Phase A ends with status: "; double time;
    showBeautifulDashlines(); printf("Phase A Log Ends. \n");
    convertStatusCode(dsdpSolver->solStatus, &sAlog[33]);
    printf("\n"); showBeautifulDashlines();
    DSDPGetStats(&dsdpSolver->dsdpStats, STAT_PHASE_A_TIME, &time);
    printf("%-80s %18s\n", sAlog, "");
    printf("| Elapsed Time: %3.3f seconds %65s \n", time, "");
    showBeautifulDashlines();
}

extern void printPhaseBSummary( HSDSolver *dsdpSolver ) {
    // Summarize Phase B iteration
    char sBlog[80] = "| DSDP Phase B ends with status: "; double time;
    showBeautifulDashlines(); printf("Phase B Log Ends. \n");
    convertStatusCode(dsdpSolver->solStatus, &sBlog[33]);
    printf("\n"); showBeautifulDashlines();
    DSDPGetStats(&dsdpSolver->dsdpStats, STAT_PHASE_B_TIME, &time);
    printf("%-80s %18s\n", sBlog, "");
    printf("| Elapsed Time: %3.3f seconds %65s \n", time, "");
    showBeautifulDashlines();
}

extern void checkPhaseBConvergence( HSDSolver *dsdpSolver, DSDP_INT *isOK ) {
    // Check convergence of DSDP
    DSDP_INT *monitor = dsdpSolver->eventMonitor, BmaxIter;
    DSDPStats *stat = &dsdpSolver->dsdpStats;
    double tmax, statval1, statval2, gap, thresh;
    *isOK = FALSE;
    
    DSDPGetIntParam(dsdpSolver, INT_PARAM_BMAXITER, &BmaxIter);
    DSDPGetDblParam(dsdpSolver, DBL_PARAM_TIMELIMIT, &tmax);
    
    // Primal infeasibility detected
    if (monitor[EVENT_PINFEAS_DETECTED]) {
        dsdpSolver->solStatus = DSDP_PINFEAS_DFEAS; *isOK = TRUE; return;
    }
    DSDPGetStats(stat, STAT_GAP_BROKEN, &statval1);
    if (dsdpSolver->dObjVal > dsdpSolver->pObjVal) {
        if (statval1 == 0.0) {
            statval1 += 1.0; DSDPStatUpdate(stat, STAT_GAP_BROKEN, statval1);
            dsdpSolver->pObjVal = dsdpSolver->dObjVal + 1;
            monitor[EVENT_INVALID_GAP] = TRUE;
        } else {
            dsdpSolver->solStatus = DSDP_INTERNAL_ERROR;
            *isOK = TRUE; return;
        }
    }
    
    // Gap is sufficiently small
    gap = fabs(dsdpSolver->pObjVal - dsdpSolver->dObjVal) * dsdpSolver->cScaler;
    gap = gap / (fabs(dsdpSolver->pObjVal) * dsdpSolver->cScaler + \
                 fabs(dsdpSolver->dObjVal) * dsdpSolver->cScaler + 1);
    
    if (dsdpSolver->dObjVal > 1e+08) {
        dsdpSolver->eventMonitor[EVENT_LARGE_DOBJ] = TRUE;
    }
    
    DSDPGetDblParam(dsdpSolver, DBL_PARAM_ABS_OPTTOL, &thresh);
    
    if (gap < thresh) {
        monitor[EVENT_MU_QUALIFIES] = TRUE;
        dsdpSolver->solStatus = DSDP_OPTIMAL; *isOK = TRUE; return;
    }
    // Maximum iteration
    DSDPGetStats(stat, STAT_PHASE_B_ITER, &statval1);
    if (statval1 >= BmaxIter) {
        monitor[EVENT_MAX_ITERATION] = TRUE;
        dsdpSolver->solStatus = DSDP_MAXITER; *isOK = TRUE;
    }
    // Small step
    DSDPGetStats(stat, STAT_NUM_SMALL_ITER, &statval2);
    if (statval1 > 0 && dsdpSolver->alpha < 3e-02 && gap < 1e-02) {
        statval2 += 1; DSDPStatUpdate(stat, STAT_NUM_SMALL_ITER, statval2);
        monitor[EVENT_SMALL_STEP] = TRUE;
    }
    if (statval2 >= 2) {
        dsdpSolver->solStatus = DSDP_INTERNAL_ERROR;
        *isOK = TRUE; return;
    }
    if (monitor[EVENT_NAN_IN_ITER]) {
        dsdpSolver->solStatus = DSDP_INTERNAL_ERROR;
        *isOK = TRUE; return;
    }
    if (monitor[EVENT_BAD_SCHUR]) {
        dsdpSolver->solStatus = DSDP_INTERNAL_ERROR;
        *isOK = TRUE; return;
    }
    if (my_clock() - dsdpSolver->startTime > tmax) {
        dsdpSolver->solStatus = DSDP_TIMELIMIT;
        *isOK = TRUE; return;
    }
}

extern void printPhaseABConvert( HSDSolver *dsdpSolver, DSDP_INT *goPb ) {
    // Print the conversion logging between phase A and B
    DSDP_INT status = dsdpSolver->solStatus; *goPb = TRUE;
    if (status == DSDP_OPTIMAL) { printf("| DSDP Phase A solves the problem %35s\n", ""); *goPb = FALSE; }
    switch (status) {
        case DSDP_PD_FEASIBLE     : printf("| DSDP Phase A certificates primal-dual feasibility %47s\n", ""); break;
        case DSDP_PUNKNOWN_DFEAS  : printf("| DSDP Phase A certificates dual feasibility %54s\n", ""); break;
        case DSDP_PFEAS_DINFEAS   : printf("| DSDP Phase A certificates dual infeasibility and primal feasibility \n");
            *goPb = FALSE; break;
        case DSDP_PUNKNOWN_DINFEAS: printf("| DSDP Phase A certificates dual infeasibility \n");
            *goPb = FALSE; break;
        case DSDP_INTERNAL_ERROR  : printf("| DSDP Phase A encounters internal error \n");
            *goPb = FALSE; break;
        case DSDP_MAXITER         : printf("| DSDP Phase A reaches maximum iteration \n");
            *goPb = FALSE; break;
        case DSDP_TIMELIMIT       : printf("| DSDP Phase A reaches maximum time limie \n");
            *goPb = FALSE; break;
        default: printf("Invalid status code at the end of Phase A. \n"); break;
    }
}
