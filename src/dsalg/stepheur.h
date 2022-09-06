#ifndef stepheur_h
#define stepheur_h

/* Implement step-size related routines */
#include "dsdpsolver.h"

#ifdef __cplusplus
extern "C" {
#endif

extern void computeAdaptivedRate( HSDSolver *dsdpSolver );
extern void getMaxStep( HSDSolver *dsdpSolver );
extern void takeStep  ( HSDSolver *dsdpSolver );
extern void selectMu  ( HSDSolver *dsdpSolver, double *newmu );
extern void dualPotentialReduction ( HSDSolver *dsdpSolver );

#ifdef __cplusplus
}
#endif

#endif /* stepheur_h */
