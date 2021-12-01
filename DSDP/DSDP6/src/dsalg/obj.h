#ifndef obj_h
#define obj_h

/* Objective value computation in DSDP */

#include "dsdphsd.h"
#include "hsd.h"

extern DSDP_INT getDualObj     ( HSDSolver *dsdpSolver );
extern DSDP_INT getSDPPrimalObj( HSDSolver *dsdpSolver );

#endif /* obj_h */
