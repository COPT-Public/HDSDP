#ifndef dsdplanczos_h
#define dsdplanczos_h

#include "dsdphsd.h"
#include "structs.h"
// Implementation of Lanczos iteration for SDP stepsize

#ifdef __cplusplus
extern "C" {
#endif

extern DSDP_INT dsdpLanczos( spsMat *S, spsMat *dS, double *lbd, double *delta );

#ifdef __cplusplus
}
#endif

#endif /* dsdplanczos_h */
