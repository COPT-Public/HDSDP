#ifndef hsd_h
#define hsd_h

#include "dsdphsd.h"
#include "dsdpsolver.h"

// Solver Phase
#define PHASEA_DUAL_INFEAS     0
#define PHASEB_PRIMAL_HEUR     1

// Iteration monitor
#define ITER_LOGGING           0
#define ITER_LP_RESIDUAL       1
#define ITER_SDP_RESIDUAL      2
#define ITER_DUAL_OBJ          3
#define ITER_DUAL_FACTORIZE    4
#define ITER_LP_SCHUR          5
#define ITER_SDP_SCHUR         6
#define ITER_AUX_ARRAY         7
#define ITER_SCHUR_SOLVE       8
#define ITER_STEP_DIRECTION    9
#define ITER_RECOVER_LP_DIR    10
#define ITER_RECOVER_SDP_DIR   11
#define ITER_COMPUTE_STEP      12
#define ITER_TAKE_STEP         13
#define ITER_PRIMAL_PROJ       14
#define ITER_COMPUTE_POBJ      15
#define ITER_CORRECTOR         16
#define ITER_DECREASE_MU       17
#define ITER_NEXT_ITERATION    18

// Special events
#define EVENT_NO_RY          0  // Dual infeasibility is eliminated
#define EVENT_DINFES         1
#define EVENT_PFEAS_FOUND    2  // Projection finds a primal feasible solution
#define EVENT_MU_QUALIFIES   3  // Duality gap is qualified by tolerance
#define EVENT_BAD_SCHUR      4  // Bad schur matrix
#define EVENT_NAN_IN_ITER    5  // Nan detected in iterations
#define EVENT_STEP_TOO_SMALL 6  // Step taken is too small
#define EVENT_ITERTION_END   7  // Iteration stops

#ifdef __cplusplus
extern "C" {
#endif

extern DSDP_INT checkIterProgress( HSDSolver *dsdpSolver, DSDP_INT iter );


#ifdef __cplusplus
}
#endif

#endif /* hsd_h */
