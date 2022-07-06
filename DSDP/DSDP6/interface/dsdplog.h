#ifndef dsdplog_h
#define dsdplog_h
#include "dsdphsd.h"
#include "dsdpsolver.h"

/* Implement problem statistic collection and logging */
// Iteration monitor
#define ITER_INITIALIZE         0  // Initialize y, S, tau, kappa and Ry
#define ITER_DUAL_OBJ           1  // Get b' * y
#define ITER_LOGGING            2  // Print iteration information
#define ITER_DUAL_FACTORIZE     3  // Factorize dual matrices {S}
#define ITER_SCHUR              4  // Set up the schur matrix M and part of auxiliary arrays
#define ITER_SCHUR_SOLVE        5  // Solve the Schur system
#define ITER_PROX_POBJ          6  // Get proximity and verify primal feasibility
#define ITER_STEP_DIRECTION     7  // Recover step directions
#define ITER_COMPUTE_STEP       8  // Compute algorithm stepsize
#define ITER_RESIDUAL           9  // Set up residuals
#define ITER_TAKE_STEP         10  // Take step
#define ITER_CORRECTOR         11  // Take corrector step
#define ITER_DECREASE_MU       12  // Decrease parameter mu
#define ITER_NEXT_ITERATION    13  // Get into next iteration

// Special events
#define EVENT_NO_RY             0  // Dual infeasibility is eliminated
#define EVENT_PFEAS_FOUND       1  // Projection finds a primal feasible solution
#define EVENT_MU_QUALIFIES      2  // Duality gap is qualified by tolerance
#define EVENT_KT_QUALIFIES      3  // Kappa / tau is large enough
#define EVENT_LARGE_NORM        4  // Norm of y is too large
#define EVENT_BAD_SCHUR         5  // Bad schur matrix
#define EVENT_NAN_IN_ITER       6  // Nan detected in iterations
#define EVENT_SMALL_STEP        7  // Step taken is too small
#define EVENT_MAX_ITERATION     8  // Maximum iteration reached
#define EVENT_LARGE_DOBJ        9  // Large dual objective
#define EVENT_PINFEAS_DETECTED 10  // Primal infeasibility is detected
#define EVENT_INVALID_GAP      11  // Dual objective exceeds primal bound
#define EVENT_HSD_UPDATE       12  // Homogeneous self-dual embedding switch
#define EVENT_IN_PHASE_A       13  // We are in phase A
#define EVENT_IN_PHASE_B       14  // We are in phase B

#ifdef __cplusplus
extern "C" {
#endif

/* Timer */
extern double my_clock(void);
/* Utils */
extern DSDP_INT checkIterProgress( HSDSolver *dsdpSolver, DSDP_INT iter );
extern void showBeautifulDashlines          ( void );
extern void printPhaseAheader               ( void );
extern void printPhaseBheader               ( void );
extern void checkNan                        ( HSDSolver *dsdpSolver );

/* Phase A operations */
extern void resetPhaseAMonitor              ( HSDSolver *dsdpSolver );
extern void checkPhaseAConvergence          ( HSDSolver *dsdpSolver, DSDP_INT *isOK );
extern void phaseALogging                   ( HSDSolver *dsdpSolver );
extern void printPhaseASummary              ( HSDSolver *dsdpSolver );
extern void printPhaseABConvert             ( HSDSolver *dsdpSolver, DSDP_INT *goPb );

/* Phase B operations */
extern void resetPhaseBMonitor              ( HSDSolver *dsdpSolver );
extern void checkPhaseBConvergence          ( HSDSolver *dsdpSolver, DSDP_INT *isOK );
extern void phaseBLogging                   ( HSDSolver *dsdpSolver );
extern void printPhaseBSummary              ( HSDSolver *dsdpSolver );

#ifdef __cplusplus
}
#endif

#endif /* dsdplog_h */
