#ifndef dsdpstats_h
#define dsdpstats_h

/* Implement the statistics interface of DSDP */
#include "dsdphsd.h"

// Data statistics
#define STAT_NUM_DENSE_MAT   0
#define STAT_NUM_SPARSE_MAT  1
#define STAT_NUM_RONE_MAT    2
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
#define STAT_SCHURORD_TIME   21
#define STAT_SPECIAL_DETECT  22

#define STAT_ONE_NORM_A      23

// Problem type
#define STAT_PFEAS_PROBLEM   24
#define STAT_DFEAS_PROBLEM   25

#define STAT_LARGEST_BLOCK   26
#define STAT_NNZ_OBJ         27
#define STAT_NNZ_SCHUR       28

// 20 ~ 50 left for event profiling
#define STAT_PHASE_A_ITER    51
#define STAT_PHASE_B_ITER    52

// Solver heuristic statistics
#define STAT_NUM_SMALL_ITER  53
#define STAT_NO_PINTERIOR    54
#define STAT_NO_DINTERIOR    55
#define STAT_IMP_UBOUNDY     56
#define STAT_IMP_LBOUNDY     57
#define STAT_IMP_BOUNDX      58
#define STAT_GAP_BROKEN      59

// DIMACS statistics
#define STAT_DIMACS_ERR1     60
#define STAT_DIMACS_ERR2     61
#define STAT_DIMACS_ERR3     62
#define STAT_DIMACS_ERR4     63
#define STAT_DIMACS_ERR5     64
#define STAT_DIMACS_ERR6     65

#define NUM_STATISTICS       66

typedef struct {
    double stats[NUM_STATISTICS];
} DSDPStats;


#ifdef __cplusplus
extern "C" {
#endif

extern void DSDPStatInit    ( DSDPStats *stats );
extern void DSDPStatUpdate  ( DSDPStats *stat, DSDP_INT sName, double  val );
extern void DSDPGetStats    ( DSDPStats *stat, DSDP_INT sName, double *val );
extern void DSDPDataStatPrint   ( DSDPStats *stat );
extern void DSDPDIMACErrorPrint ( DSDPStats *stat );
extern void DSDPBProfilerPrint  ( DSDPStats *stat );

#ifdef __cplusplus
}
#endif

#endif /* dsdpstats_h */
