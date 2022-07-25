#ifndef dsdplanczos_h
#define dsdplanczos_h

#include "structs.h"

// Implementation of Lanczos iteration for SDP stepsize

#ifdef __cplusplus
extern "C" {
#endif

extern DSDP_INT dsdpLanczosInit ( DSDPLanczos *lczSolver );
extern DSDP_INT dsdpLanczosAlloc( DSDPLanczos *lczSolver, DSDP_INT n );
extern DSDP_INT dsdpLanczosStep ( DSDPLanczos *lczSolver, spsMat *S, spsMat *dS, double *lbd, double *delta );
extern void     dsdpLanczosFree ( DSDPLanczos *lczSolver );

#ifdef __cplusplus
}
#endif

#endif /* dsdplanczos_h */
