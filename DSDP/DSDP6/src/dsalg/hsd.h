#ifndef hsd_h
#define hsd_h

#include "dsdphsd.h"
#include "dsdpsolver.h"

// Iteration monitor
#define ITER_INITIALIZE         0  // Initialize y, S, tau, kappa and Ry
#define ITER_LOGGING            1  // Print | niter | pObj | dObj | dInf | k/t | mu | alpha | pNrm | E |
#define ITER_DUAL_OBJ           2  // Get b' * y
#define ITER_DUAL_FACTORIZE     3  // Factorize dual matrices {S}
#define ITER_SCHUR              4  // Set up the schur matrix M and part of auxiliary arrays
#define ITER_SCHUR_SOLVE        5  // Solve the Schur system
#define ITER_PROX_POBJ          6  // Get proximity and verify primal feasibility
#define ITER_STEP_DIRECTION     7  // Recover step directions
#define ITER_COMPUTE_STEP       8  // Compute algorithm stepsize
#define ITER_TAKE_STEP          9  // Take step
#define ITER_CORRECTOR         10  // Take corrector step
#define ITER_RESIDUAL          11  // Set up residuals
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
#define EVENT_IN_PHASE_A        9  // We are in phase A
#define EVENT_IN_PHASE_B       10  // We are in phase B

#ifdef __cplusplus
extern "C" {
#endif

extern DSDP_INT checkIterProgress( HSDSolver *dsdpSolver, DSDP_INT iter );

#ifdef __cplusplus
}
#endif

#endif /* hsd_h */
