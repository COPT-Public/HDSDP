#ifndef dsdpinitializer_h
#define dsdpinitializer_h

/* Implement the initialization strategy for DSDP */
#include "dsdpsolver.h"

#ifdef __cplusplus
extern "C" {
#endif

extern DSDP_INT dsdpInitializeA( HSDSolver *dsdpSolver );
extern DSDP_INT dsdpInitializeB( HSDSolver *dsdpSolver );

#ifdef __cplusplus
}
#endif

#endif /* dsdpinitializer_h */
