#ifndef hsd_h
#define hsd_h

#include "dsdphsd.h"
#include "dsdpsolver.h"

// Iteration monitor
#define ITER_LOGGING           0
#define ITER_LP_DUAL_RESIDUAL  1
#define ITER_SDP_DUAL_RESIDUAL 2
#define ITER_DUAL_OBJ          3
#define ITER_AUX_ARRAY         4
#define ITER_LP_SCHUR          5
#define ITER_SDP_SCHUR         6
#define ITER_STEP_DIRECTION    7
#define ITER_RECOVER_LP_DIR    8
#define ITER_RECOVER_SDP_DIR   9
#define ITER_COMPUTE_LP_STEP   10
#define ITER_COMPUTE_SDP_STEP  11
#define ITER_TAKE_STEP         12
#define ITER_PRIMAL_PROJ       13
#define ITER_COMPUTE_POBJ      14
#define ITER_CORRECTOR         15
#define ITER_DECREASE_MU       16
#define ITER_NEXT_ITERATION    17

// Special events
#define EVENT_INIT_ITER    (-1) // Initial iteration
#define EVENT_LP_NO_RY       0  // LP dual infeasibility is eliminated
#define EVENT_SDP_NO_RY      1  // SDP dual infeasibility is eliminated
#define EVENT_NO_RKAPPA      2  // Complementary infeasibility is eliminated
#define EVENT_PFEAS_FOUND    3  // Projection finds a primal feasible solution
#define EVENT_MU_QUALIFIES   4  // Duality gap is qualified by tolerance
#define EVENT_BAD_SCHUR      5  // Bad schur matrix
#define EVENT_NAN_IN_ITER    6  // Nan detected in iterations
#define EVENT_STEP_TOO_SMALL 7  // Step taken is too small
#define EVENT_ITERTION_END   8  // Iteration stops

extern DSDP_INT checkIterProgress( HSDSolver *dsdpSolver, DSDP_INT iter );

#endif /* hsd_h */
