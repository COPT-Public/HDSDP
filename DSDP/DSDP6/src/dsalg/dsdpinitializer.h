#ifndef dsdpinitializer_h
#define dsdpinitializer_h
/* Implement the initialization strategy for DSDP
 Currently we initialize DSDP with y = 0 and S =
*/

#include "dsdphsd.h"
#include "dsdpdata.h"
#include "structs.h"
#include "dsdpsolver.h"
#include "hsd.h"

#ifdef __cplusplus
extern "C" {
#endif

extern DSDP_INT dsdpInitialize( HSDSolver *dsdpSolver );

#ifdef __cplusplus
}
#endif

#endif /* dsdpinitializer_h */
