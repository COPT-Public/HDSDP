#ifndef dsdpdinfeas_h
#define dsdpdinfeas_h
/* Implement phase A of dual scaling */

#include "dsdpsolver.h"

#ifdef __cplusplus
extern "C" {
#endif

extern DSDP_INT DSDPDInfeasEliminator( HSDSolver *dsdpSolver );

#ifdef __cplusplus
}
#endif

#endif /* dsdpdinfeas_h */
