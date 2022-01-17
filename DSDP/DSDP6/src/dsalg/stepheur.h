#ifndef stepheur_h
#define stepheur_h

/* Implement step-size related routines */
#include "dsdphsd.h"
#include "dsdpdata.h"
#include "structs.h"
#include "dsdpsolver.h"
#include "hsd.h"

#ifdef __cplusplus
extern "C" {
#endif

extern DSDP_INT getMaxStep( HSDSolver *dsdpSolver );
extern DSDP_INT takeStep  ( HSDSolver *dsdpSolver );

#ifdef __cplusplus
}
#endif

#endif /* stepheur_h */
