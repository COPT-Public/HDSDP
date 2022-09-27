#ifndef dsdpsort_h
#define dsdpsort_h

#include "dsdphsd.h"

#ifdef __cplusplus
extern "C" {
#endif

extern DSDP_INT checkIsOrdered( DSDP_INT *idx, DSDP_INT n );
extern void dsdpSort ( double *data, DSDP_INT *idxbase, DSDP_INT low, DSDP_INT high );
extern void dsdpSort2( DSDP_INT *data, DSDP_INT *idxbase, DSDP_INT low, DSDP_INT high );

#ifdef __cplusplus
}
#endif

#endif /* dsdpsort_h */
