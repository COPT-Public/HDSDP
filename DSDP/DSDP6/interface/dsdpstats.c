#include "dsdpstats.h"
#include "dsdplog.h"

// Implement the solver statistics interface for DSDP
static char etype[] = "Solver statistics";

static void printStat( DSDPStats *stats, DSDP_INT sName ) {
    // Print the value of some statistic
    printf("| ");
    switch (sName) {
        case STAT_NUM_DENSE_MAT : printf("Total number of dense matrices: "); break;
        case STAT_NUM_SPARSE_MAT: printf("Total number of sparse matrices: "); break;
        case STAT_NUM_RONE_MAT  : printf("Total number of rankone matrices: "); break;
        case STAT_NUM_ZERO_MAT  : printf("Total number of zero matrices: "); break;
        case STAT_ONE_NORM_C    : printf("|C|: "); break;
        case STAT_ONE_NORM_B    : printf("|b|: "); break;
        case STAT_INF_NORM_Y    : printf("|y|: "); break;
        case STAT_TRACE_S       : printf("tr(S): "); break;
        case STAT_TRACE_X       : printf("tr(X): "); break;
        case STAT_PRESOLVE_TIME : printf("Presolve time: "); break;
        case STAT_PHASE_A_TIME  : printf("Phase A time: "); break;
        case STAT_PHASE_B_TIME  : printf("Phase B time: "); break;
        case STAT_GET_X_TIME    : printf("Time getting X: "); break;
        case STAT_POSTSOLVE_TIME: printf("Postsolve time: "); break;
        case STAT_READ_TIME     : printf("Time reading data: "); break;
        case STAT_SCAL_TIME     : printf("Time scaling data: "); break;
        case STAT_RONE_TIME     : printf("Time detecting rank-1 structure: "); break;
        case STAT_EIG_TIME      : printf("Time in eigen-decomposition: "); break;
        case STAT_MATSTAT_TIME  : printf("Time getting matrix statistics: "); break;
        case STAT_SYMBOLIC_TIME : printf("Time doing symbolic ordering: "); break;
        case STAT_SCHURORD_TIME : printf("Time doing schur ordering: "); break;
        case STAT_SPECIAL_DETECT: printf("Time detecting special structure: "); break;
        case STAT_ONE_NORM_A    : printf("|A|: "); break;
        case STAT_PFEAS_PROBLEM : printf("Primal feasibility problem: "); break;
        case STAT_DFEAS_PROBLEM : printf("Dual feasibility problem: "); break;
        case STAT_LARGEST_BLOCK : printf("Largest SDP block: "); break;
        case STAT_NNZ_OBJ       : printf("Nnzs in C: "); break;
        case STAT_NNZ_SCHUR     : printf("Nnzs in M: "); break;
        case STAT_PHASE_A_ITER  : printf("Phase A iter: "); break;
        case STAT_PHASE_B_ITER  : printf("Phase B time: "); break;
        case STAT_NUM_SMALL_ITER: printf("nIter of small step: "); break;
        case STAT_NO_PINTERIOR  : printf("No primal interior: "); break;
        case STAT_NO_DINTERIOR  : printf("No dual interior: "); break;
        case STAT_IMP_BOUNDX    : printf("Implicit bound on tr(X) "); break;
        case STAT_IMP_UBOUNDY   : printf("Implicit upper bound on y: "); break;
        case STAT_IMP_LBOUNDY   : printf("Implicit lower bound on y: "); break;
        case STAT_GAP_BROKEN    : printf("Gap is broken: "); break;
        case STAT_DIMACS_ERR1   : printf("DIMACS Err 1: "); break;
        case STAT_DIMACS_ERR2   : printf("DIMACS Err 2: "); break;
        case STAT_DIMACS_ERR3   : printf("DIMACS Err 3: "); break;
        case STAT_DIMACS_ERR4   : printf("DIMACS Err 4: "); break;
        case STAT_DIMACS_ERR5   : printf("DIMACS Err 5: "); break;
        case STAT_DIMACS_ERR6   : printf("DIMACS Err 6: "); break;
        default: printf("Invalid Statistic Code %d \n", sName); break;
    }
    return;
}

static void printSolFlag(double err) {
    
    if (err < 1e-04) {
        printf("| Passed Mittelmann's Benchmark Test (*) \n");
    } else if (err < 1e-02) {
        printf("| Passed Mittelmann's Benchmark Test (a) \n");
    } else {
        printf("| Failed to pass Mittelmann's Benchmark Test \n");
    }
}

extern void DSDPStatInit( DSDPStats *stats ) {
    
    memset(stats->stats, 0, sizeof(double) * NUM_STATISTICS);
    DSDPStatUpdate(stats, STAT_DIMACS_ERR1, DSDP_INFINITY);
    DSDPStatUpdate(stats, STAT_DIMACS_ERR2, DSDP_INFINITY);
    DSDPStatUpdate(stats, STAT_DIMACS_ERR3, DSDP_INFINITY);
    DSDPStatUpdate(stats, STAT_DIMACS_ERR4, DSDP_INFINITY);
    DSDPStatUpdate(stats, STAT_DIMACS_ERR5, DSDP_INFINITY);
    DSDPStatUpdate(stats, STAT_DIMACS_ERR6, DSDP_INFINITY);
}

extern void DSDPStatUpdate( DSDPStats *stat, DSDP_INT sName, double val ) {
    
    if (sName > NUM_STATISTICS || sName < 0) {
        printf("| Invalid Statistic Code: %d \n", sName); return;
    }
    stat->stats[sName] = val;
}

extern void DSDPGetStats( DSDPStats *stat, DSDP_INT sName, double *val ) {
    
    if (sName > NUM_STATISTICS || sName < 0) {
        printf("| Invalid Statistic Code: %d \n", sName); return;
    }
    *val = stat->stats[sName];
}

/* DSDP Summary statistics printer */
extern void DSDPDataStatPrint( DSDPStats *stat ) {
    // Print DSDP Data summary statistics. Invoked after pre-solving
    double sps, ds, r1, zero, a, b, c;
    DSDPGetStats(stat, STAT_NUM_SPARSE_MAT, &sps); DSDPGetStats(stat, STAT_NUM_DENSE_MAT, &ds);
    DSDPGetStats(stat, STAT_NUM_RONE_MAT, &r1); DSDPGetStats(stat, STAT_NUM_ZERO_MAT, &zero);
    DSDPGetStats(stat, STAT_ONE_NORM_A, &a); DSDPGetStats(stat, STAT_ONE_NORM_B, &b);
    DSDPGetStats(stat, STAT_ONE_NORM_C, &c);
    showBeautifulDashlines();
    printf("| Matrix statistics [Including C]: \n");
    showBeautifulDashlines();
    printf("| %10s | %10s | %10s | %10s | %8s     | %8s     | %8s     \n",
           "Zero", "Sparse", "Dense", "Rank-1", "|A|", "|b|", "|C|");
    printf("----------------------------------------------------|"
           "---------------------------------------------\n");
    printf("| %10d | %10d | %10d | %10d |  %10.5e |  %10.5e |  %10.5e \n",
           (DSDP_INT) zero, (DSDP_INT) sps, (DSDP_INT) ds, (DSDP_INT) r1, a, b, c);
    return;
}

extern void DSDPDIMACErrorPrint( DSDPStats *stat ) {
    // Print DSDP DIMACS error. Invoked after extracting primal solution (and DIMACS computed)
    double err1, err2, err3, err4, err5, err6, err = 0.0;
    DSDPGetStats(stat, STAT_DIMACS_ERR1, &err1); DSDPGetStats(stat, STAT_DIMACS_ERR2, &err2);
    DSDPGetStats(stat, STAT_DIMACS_ERR3, &err3); DSDPGetStats(stat, STAT_DIMACS_ERR4, &err4);
    DSDPGetStats(stat, STAT_DIMACS_ERR5, &err5); DSDPGetStats(stat, STAT_DIMACS_ERR6, &err6);
    err = MAX(err, fabs(err1)); err = MAX(err, fabs(err2)); err = MAX(err, fabs(err3));
    err = MAX(err, fabs(err4)); err = MAX(err, fabs(err5)); err = MAX(err, fabs(err6));
    
    printf("| DIMACS Error:\n");
    showBeautifulDashlines();
    printf("| Err1: %-6.2e | Err2: %-6.1e | Err3: %-6.1e | Err4: %-6.1e | Err5: %-8.1e | Err6: %-8.1e \n",
           err1, err2, err3, err4, err5, err6);
    showBeautifulDashlines(); printSolFlag(err); showBeautifulDashlines();
}

extern void DSDPBProfilerPrint( DSDPStats *stat ) {
    // Print brief time profiling of HDSDP
    double tpresolve, tphaseA, tphaseB, tgetX, tpostsolve, iterA, iterB;
    DSDPGetStats(stat, STAT_PRESOLVE_TIME, &tpresolve); DSDPGetStats(stat, STAT_PHASE_A_TIME, &tphaseA);
    DSDPGetStats(stat, STAT_PHASE_A_ITER, &iterA); DSDPGetStats(stat, STAT_PHASE_B_TIME, &tphaseB);
    DSDPGetStats(stat, STAT_PHASE_B_ITER, &iterB); DSDPGetStats(stat, STAT_GET_X_TIME, &tgetX);
    DSDPGetStats(stat, STAT_POSTSOLVE_TIME, &tpostsolve);
    printf("| DSDP Time Summary: \n");
    showBeautifulDashlines();
    printf("| %15s | %10s | \n", "Event", "Time(s)");
    showBeautifulDashlines();
    printf("| %15s | %10.3f | \n", "Presolve", tpresolve);
    printf("| %15s | %10.3f | (%d) \n", "Phase A (iter)", tphaseA, (DSDP_INT) iterA);
    printf("| %15s | %10.3f | (%d) \n", "Phase B (iter)", tphaseB, (DSDP_INT) iterB);
    printf("| %15s | %10.3f | \n", "Get X", tgetX);
    printf("| %15s | %10.3f | \n", "Postsolve", tpostsolve);
    printf("| %15s | %10.3f | (%d) \n", "All",
           tpresolve + tphaseA + tphaseB + tgetX + tpostsolve, (DSDP_INT) (iterA + iterB));
}
