#ifndef dsdppfeascheck_h
#define dsdppfeascheck_h
#include "dsdpdata.h"
#include "structs.h"
#include "dsdpsolver.h"
#include "dsdplog.h"
/* Implement the primal feasibility check */

#ifdef __cplusplus
extern "C" {
#endif

extern void dsdpCheckPhaseAPfeas    ( HSDSolver *dsdpSolver, double dtaudelta,
                                          vec *dydelta, DSDP_INT *isPfeas );
extern void dsdpCheckBackwardNewton ( HSDSolver *dsdpSolver, DSDP_INT *ispfeas );
extern void dsdpCheckPrimalInfeas   ( HSDSolver *dsdpSolver );

#ifdef __cplusplus
}
#endif
#endif /* dsdppfeascheck_h */
