#ifndef obj_h
#define obj_h

/* Objective value computation in DSDP */

#include "dsdphsd.h"
#include "hsd.h"

#ifdef __cplusplus
extern "C" {
#endif

extern DSDP_INT getDualObj     ( HSDSolver *dsdpSolver );
extern DSDP_INT getSDPPrimalObj( HSDSolver *dsdpSolver );

#ifdef __cplusplus
}
#endif

#endif /* obj_h */
