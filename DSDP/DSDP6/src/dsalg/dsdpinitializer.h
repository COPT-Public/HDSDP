#ifndef dsdpinitializer_h
#define dsdpinitializer_h
/* Implement the initialization strategy for DSDP
 Currently we initialize DSDP with y = 0 and S =
*/

#include "dsdphsd.h"
#include "dsdpdata.h"
#include "structs.h"
#include "dsdpsolver.h"
#include "dsdplog.h"

#ifdef __cplusplus
extern "C" {
#endif

extern DSDP_INT dsdpInitializeA( HSDSolver *dsdpSolver );
extern DSDP_INT dsdpInitializeB( HSDSolver *dsdpSolver );

#ifdef __cplusplus
}
#endif

#endif /* dsdpinitializer_h */
