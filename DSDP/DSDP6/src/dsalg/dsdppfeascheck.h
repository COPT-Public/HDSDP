#ifndef dsdppfeascheck_h
#define dsdppfeascheck_h
#include "dsdpdata.h"
#include "structs.h"
#include "dsdpsolver.h"
#include "hsd.h"
/* Implement the primal feasibility check */

#ifdef __cplusplus
extern "C" {
#endif

extern DSDP_INT dsdpCheckPhaseAPfeas   ( HSDSolver *dsdpSolver, double dtaudelta,
                                         vec *dydelta, DSDP_INT *isPfeas );
extern DSDP_INT dsdpCheckBackwardNewton( HSDSolver *dsdpSolver, DSDP_INT *ispfeas );

#ifdef __cplusplus
}
#endif
#endif /* dsdppfeascheck_h */
