#ifndef stepheur_h
#define stepheur_h

#include "dsdphsd.h"
#include "dsdpdata.h"
#include "structs.h"
#include "dsdpsolver.h"
#include "hsd.h"

#ifdef __cplusplus
extern "C" {
#endif

extern DSDP_INT getMaxStep( HSDSolver *dsdpSolver, double *maxStep );
extern DSDP_INT takeStep  ( HSDSolver *dsdpSolver, double step );

#ifdef __cplusplus
}
#endif

#endif /* stepheur_h */
