#ifndef dsdpproxmeasure_h
#define dsdpproxmeasure_h

/* Get the measure of proximity in DSDP */
#include "dsdpsolver.h"

#ifdef __cplusplus
extern "C" {
#endif

extern DSDP_INT dsdpgetPhaseAProxMeasure( HSDSolver *dsdpSolver, double newmu );
extern DSDP_INT dsdpgetPhaseBProxMeasure( HSDSolver *dsdpSolver, double *muub, double *mulb );

#ifdef __cplusplus
}
#endif

#endif /* dsdpproxmeasure_h */
