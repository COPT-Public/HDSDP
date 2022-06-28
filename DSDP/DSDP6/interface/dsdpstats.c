#include <string.h>
#include "dsdpstats.h"
#include "hsd.h"
#include "dsdplog.h"

// Implement the solver statistics interface for DSDP
static char etype[] = "Solver statistics";

static void printStat( DSDPStats *stats, DSDP_INT sName ) {
    // Print the value of some statistic
    printf("| ");
    
    if (sName > NUM_STATISTICS || sName < 0) {
        printf("Invalid Statistic Code %d \n", sName);
        return;
    }
    
    switch (sName) {
        case STAT_NUM_DENSE_MAT:
            printf("Total number of dense matrices: "); break;
        case STAT_NUM_SPARSE_MAT:
            printf("Total number of sparse matrices: "); break;
        case STAT_NUM_RANKONE_MAT:
            printf("Total number of rankone matrices: "); break;
        case STAT_NUM_ZERO_MAT:
            printf("Total number of zero matrices: "); break;
        case STAT_ONE_NORM_C:
            printf("|C|: "); break;
        case STAT_ONE_NORM_B:
            printf("|b|: "); break;
        case STAT_INF_NORM_Y:
            printf("|y|: "); break;
        case STAT_PRESOLVE_TIME:
            printf("Presolve time: "); break;
        case STAT_PHASE_A_TIME:
            printf("Phase A time: "); break;
        case STAT_PHASE_B_TIME:
            printf("Phase B time: "); break;
        case STAT_GET_X_TIME:
            printf("Time getting X: "); break;
        case STAT_PHASE_A_ITER:
            printf("Phase A iter: "); break;
        case STAT_PHASE_B_ITER:
            printf("Phase B time: "); break;
        case STAT_NUM_SMALL_ITER:
            printf("nIter of small step: "); break;
        case STAT_NO_PINTERIOR:
            printf("No primal interior: "); break;
        case STAT_NO_DINTERIOR:
            printf("No dual interior: "); break;
        case STAT_IMP_BOUNDX:
            printf("Implicit bound on tr(X) "); break;
        case STAT_IMP_UBOUNDY:
            printf("Implicit upper bound on y: "); break;
        case STAT_IMP_LBOUNDY:
            printf("Implicit lower bound on y: "); break;
        case STAT_GAP_BROKEN:
            printf("Gap is broken: "); break;
        case STAT_DIMACS_ERR1:
            printf("DIMACS Err 1: "); break;
        case STAT_DIMACS_ERR2:
            printf("DIMACS Err 2: "); break;
        case STAT_DIMACS_ERR3:
            printf("DIMACS Err 3: "); break;
        case STAT_DIMACS_ERR4:
            printf("DIMACS Err 4: "); break;
        case STAT_DIMACS_ERR5:
            printf("DIMACS Err 5: "); break;
        case STAT_DIMACS_ERR6:
            printf("DIMACS Err 6: "); break;
        default:
            printf("Invalid Statistic Code %d \n", sName);
            break;
    }
    
    return;
}

static DSDP_INT checkStats( DSDP_INT sName, double val ) {
    // Check if a value is valid
    DSDP_INT retcode = DSDP_RETCODE_OK;
    // Currently do nothing
    return retcode;
}

extern DSDP_INT DSDPStatInit( DSDPStats *stats ) {
    
    DSDP_INT retcode = DSDP_RETCODE_OK;
    
    memset(stats->stats, 0, sizeof(double) * NUM_STATISTICS);
    DSDPStatUpdate(stats, STAT_DIMACS_ERR1, DSDP_INFINITY);
    DSDPStatUpdate(stats, STAT_DIMACS_ERR2, DSDP_INFINITY);
    DSDPStatUpdate(stats, STAT_DIMACS_ERR3, DSDP_INFINITY);
    DSDPStatUpdate(stats, STAT_DIMACS_ERR4, DSDP_INFINITY);
    DSDPStatUpdate(stats, STAT_DIMACS_ERR5, DSDP_INFINITY);
    DSDPStatUpdate(stats, STAT_DIMACS_ERR6, DSDP_INFINITY);
    
    return retcode;
}

extern DSDP_INT DSDPStatUpdate( DSDPStats *stat, DSDP_INT sName, double val ) {
    
    DSDP_INT retcode = DSDP_RETCODE_OK;
    if (sName > NUM_STATISTICS || sName < 0) {
        printf("| Invalid Statistic Code: %d \n", sName);
        return DSDP_RETCODE_FAILED;
    }
    
    retcode = checkStats(sName, val); checkCode;
    
#ifdef SHOWALL
    printf("| Before \n");
    printStat(stat, sName);
#endif
    stat->stats[sName] = val;
    
#ifdef SHOWALL
    printf("| After \n")
    printStat(stat, sName);
#endif
    
    return retcode;
}

extern DSDP_INT DSDPGetStats( DSDPStats *stat, DSDP_INT sName, double *val ) {
    
    DSDP_INT retcode = DSDP_RETCODE_OK;
    if (sName > NUM_STATISTICS || sName < 0) {
        printf("| Invalid Statistic Code: %d \n", sName);
        return DSDP_RETCODE_FAILED;
    }
    *val = stat->stats[sName];
    
    return retcode;
}

/* DSDP Summary statistics printer */
extern void DSDPDataStatPrint( DSDPStats *stat ) {
    // Print DSDP Data summary statistics. Invoked after pre-solving
    
    double sps, ds, r1, zero, a, b, c;
    
    DSDPGetStats(stat, STAT_NUM_SPARSE_MAT, &sps);
    DSDPGetStats(stat, STAT_NUM_DENSE_MAT, &ds);
    DSDPGetStats(stat, STAT_NUM_RANKONE_MAT, &r1);
    DSDPGetStats(stat, STAT_NUM_ZERO_MAT, &zero);
    DSDPGetStats(stat, STAT_ONE_NORM_A, &a);
    DSDPGetStats(stat, STAT_ONE_NORM_B, &b);
    DSDPGetStats(stat, STAT_ONE_NORM_C, &c);
    
    dsdpshowdash();
    printf("| Matrix statistics [Including C]: \n");
    dsdpshowdash();
    printf("| %10s | %10s | %10s | %10s | %8s     | %8s     | %8s     \n",
           "Zero", "Sparse", "Dense", "Rank-1", "|A|", "|b|", "|C|");
    printf("----------------------------------------------------|"
           "---------------------------------------------\n");
    printf("| %10d | %10d | %10d | %10d |  %10.5e |  %10.5e |  %10.5e \n",
           (DSDP_INT) zero, (DSDP_INT) sps, (DSDP_INT) ds, (DSDP_INT) r1, a, b, c);
//    dsdpshowdash();
//    printf("| |A|: %10.6e  |b|: %10.6e  |C|: %10.6e "
//           "                                    \n", a, b, c);
    
    return;
}

extern void DSDPDIMACErrorPrint( DSDPStats *stat ) {
    // Print DSDP DIMACS error. Invoked after extracting primal solution (and DIMACS computed)
    double err1, err2, err3, err4, err5, err6, err = 0.0;
    
    DSDPGetStats(stat, STAT_DIMACS_ERR1, &err1);
    DSDPGetStats(stat, STAT_DIMACS_ERR2, &err2);
    DSDPGetStats(stat, STAT_DIMACS_ERR3, &err3);
    DSDPGetStats(stat, STAT_DIMACS_ERR4, &err4);
    DSDPGetStats(stat, STAT_DIMACS_ERR5, &err5);
    DSDPGetStats(stat, STAT_DIMACS_ERR6, &err6);
    
    err = MAX(err, fabs(err1));
    err = MAX(err, fabs(err2));
    err = MAX(err, fabs(err3));
    err = MAX(err, fabs(err4));
    err = MAX(err, fabs(err5));
    err = MAX(err, fabs(err6));
    
    printf("| DIMACS Error:                                  "
           "                                            \n");
    dsdpshowdash();
    printf("| Err1: %-6.2e | Err2: %-6.1e | Err3: %-6.1e | Err4: %-6.1e | Err5: %-8.1e | Err6: %-8.1e \n",
           err1, err2, err3, err4, err5, err6);
    dsdpshowdash();
    if (err < 1e-04) {
        printf("| Passed Mittelmann's Benchmark Test (*) \n");
    } else if (err < 1e-02) {
        printf("| Passed Mittelmann's Benchmark Test (a) \n");
    } else {
        printf("| Failed to pass Mittelmann's Benchmark Test \n");
    }
    
    dsdpshowdash();
    return;
}

extern void DSDPBProfilerPrint( DSDPStats *stat ) {
    // Print brief time profiling of HDSDP
    
    double tpresolve, tphaseA, tphaseB, tgetX, tpostsolve;
    double iterA, iterB;
    
    DSDPGetStats(stat, STAT_PRESOLVE_TIME, &tpresolve);
    DSDPGetStats(stat, STAT_PHASE_A_TIME, &tphaseA);
    DSDPGetStats(stat, STAT_PHASE_A_ITER, &iterA);
    DSDPGetStats(stat, STAT_PHASE_B_TIME, &tphaseB);
    DSDPGetStats(stat, STAT_PHASE_B_ITER, &iterB);
    DSDPGetStats(stat, STAT_GET_X_TIME, &tgetX);
    DSDPGetStats(stat, STAT_POSTSOLVE_TIME, &tpostsolve);
    
    printf("| DSDP Time Summary: \n");
    dsdpshowdash();
    printf("| %15s | %10s | \n", "Event", "Time(s)");
    dsdpshowdash();
    printf("| %15s | %10.3f | \n", "Presolve", tpresolve);
    printf("| %15s | %10.3f | (%d) \n", "Phase A (iter)", tphaseA, (DSDP_INT) iterA);
    printf("| %15s | %10.3f | (%d) \n", "Phase B (iter)", tphaseB, (DSDP_INT) iterB);
    printf("| %15s | %10.3f | \n", "Get X", tgetX);
    printf("| %15s | %10.3f | \n", "Postsolve", tpostsolve);
    printf("| %15s | %10.3f | (%d) \n", "All",
           tpresolve + tphaseA + tphaseB + tgetX + tpostsolve, (DSDP_INT) (iterA + iterB));
    
    return;
}
