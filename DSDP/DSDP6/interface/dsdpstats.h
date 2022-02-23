#ifndef dsdpstats_h
#define dsdpstats_h

/* Implement the statistics interface of DSDP */
#include "dsdphsd.h"

// Data statistics
#define STAT_NUM_DENSE_MAT   0
#define STAT_NUM_SPARSE_MAT  1
#define STAT_NUM_RANKONE_MAT 2
#define STAT_NUM_ZERO_MAT    3

#define STAT_ONE_NORM_C      4
#define STAT_FNORM_C         5
#define STAT_ONE_NORM_B      6
#define STAT_INF_NORM_Y      7
#define STAT_TRACE_S         8
#define STAT_TRACE_X         9

// Time statistics
#define STAT_PRESOLVE_TIME   10
#define STAT_PHASE_A_TIME    11
#define STAT_PHASE_B_TIME    12
#define STAT_GET_X_TIME      13
#define STAT_POSTSOLVE_TIME  14
#define STAT_READ_TIME       15

#define STAT_SCAL_TIME       16
#define STAT_RONE_TIME       17
#define STAT_EIG_TIME        18
#define STAT_MATSTAT_TIME    19
#define STAT_SYMBOLIC_TIME   20

#define STAT_ONE_NORM_A      21

// 20 ~ 50 left for event profiling
#define STAT_PHASE_A_ITER    51
#define STAT_PHASE_B_ITER    52

// Solver heuristic statistics
#define STAT_NUM_SMALL_ITER  53
#define STAT_C_ISCONSTANT    54
#define STAT_B_ISNULL        55
#define STAT_SCHUR_DIAG      56
#define STAT_GAP_BROKEN      57

// DIMACS statistics
#define STAT_DIMACS_ERR1     58
#define STAT_DIMACS_ERR2     59
#define STAT_DIMACS_ERR3     60
#define STAT_DIMACS_ERR4     61
#define STAT_DIMACS_ERR5     62
#define STAT_DIMACS_ERR6     63

#define NUM_STATISTICS       64

typedef struct {
    
    double stats[NUM_STATISTICS];
    
} DSDPStats;

extern DSDP_INT DSDPStatInit    ( DSDPStats *stats );
extern DSDP_INT DSDPStatUpdate  ( DSDPStats *stat, DSDP_INT sName, double  val );
extern DSDP_INT DSDPGetStats    ( DSDPStats *stat, DSDP_INT sName, double *val );
extern void DSDPDataStatPrint   ( DSDPStats *stat );
extern void DSDPDIMACErrorPrint ( DSDPStats *stat );
extern void DSDPBProfilerPrint  ( DSDPStats *stat );


#endif /* dsdpstats_h */
