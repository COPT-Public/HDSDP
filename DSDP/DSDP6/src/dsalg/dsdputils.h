#ifndef dsdputils_h
#define dsdputils_h
/* Define DSDP utility routines */

#include "dsdpdata.h"
#include "structs.h"
#include "dsdpsolver.h"
#include "hsd.h"

#ifdef __cplusplus
extern "C" {
#endif

extern DSDP_INT getMatnrm       ( HSDSolver *dsdpSolver, DSDP_INT blockid, DSDP_INT constrid, double *nrm );
extern DSDP_INT addMattoS       ( HSDSolver *dsdpSolver, DSDP_INT blockid, DSDP_INT constrid, double alpha );
extern DSDP_INT getDualObj      ( HSDSolver *dsdpSolver );
extern DSDP_INT getSDPPrimalObjB( HSDSolver *dsdpSolver );

#ifdef __cplusplus
}
#endif

#endif /* dsdputils_h */
