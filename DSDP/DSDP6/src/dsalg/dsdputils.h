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

/* Schur matrix set up */
extern DSDP_INT getSinvASinv       ( HSDSolver *dsdpSolver, DSDP_INT blockid, DSDP_INT constrid, void *SinvASinv );
extern DSDP_INT getTraceASinvASinv ( HSDSolver *dsdpSolver, DSDP_INT blockid, DSDP_INT constrid,
                                     DSDP_INT mattype, DSDP_INT constrid2, void *SinvASinv );

/* Check definiteness of matrices */
extern DSDP_INT getPhaseAS         ( HSDSolver *dsdpSolver, double *y, double tau );
extern DSDP_INT dsdpInCone         ( HSDSolver *dsdpSolver, DSDP_INT *ispsd );

/* Objective */
extern DSDP_INT getDualObj         ( HSDSolver *dsdpSolver );
extern DSDP_INT getSDPPrimalObjB   ( HSDSolver *dsdpSolver );

/* Other utilities */
extern DSDP_INT getMatnrm          ( HSDSolver *dsdpSolver, DSDP_INT blockid, DSDP_INT constrid, double *nrm );
extern DSDP_INT addMattoS          ( HSDSolver *dsdpSolver, DSDP_INT blockid, DSDP_INT constrid, double alpha );


#ifdef __cplusplus
}
#endif

#endif /* dsdputils_h */