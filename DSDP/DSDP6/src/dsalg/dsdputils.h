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

#ifndef DSDPConic
#define SDPConic(x)  SDPCone##x
#define LPConic(x)   LPCone##x
#define BConic(x)    BCone##x
#define DSDPConic(x) DSDPCone##x

#define COPS_GET_A_ONE_NORM GetDataOneNorm
#define COPS_GET_C_ONE_NORM GetObjOneNorm
#define COPS_DO_C_SCALE     DoObjScale

#endif


/* Schur matrix set up */
extern void     invertDualVars     ( HSDSolver *dsdpSolver );
extern DSDP_INT getSinvASinv       ( HSDSolver *dsdpSolver, DSDP_INT blockid, DSDP_INT constrid,
                                     void *SinvASinv );
extern DSDP_INT getTraceASinvASinv ( HSDSolver *dsdpSolver, DSDP_INT blockid, DSDP_INT constrid,
                                     DSDP_INT mattype, DSDP_INT constrid2, void *SinvASinv );

/* Check definiteness of matrices */
extern DSDP_INT getPhaseAS         ( HSDSolver *dsdpSolver, double *y,  double tau );
extern DSDP_INT getPhaseACheckerS  ( HSDSolver *dsdpSolver, double *y,  double tau );
extern DSDP_INT getPhaseAdS        ( HSDSolver *dsdpSolver, double drate, double *dy, double dtau );
extern DSDP_INT getPhaseBS         ( HSDSolver *dsdpSolver, double *y );
extern DSDP_INT getPhaseBCheckerS  ( HSDSolver *dsdpSolver, double *y );
extern DSDP_INT getPhaseBdS        ( HSDSolver *dsdpSolver, double alpha, double *dy, double beta );
extern DSDP_INT dsdpCheckerInCone  ( HSDSolver *dsdpSolver, DSDP_INT *ispsd );
extern DSDP_INT dsdpInCone         ( HSDSolver *dsdpSolver, DSDP_INT *ispsd );

/* Objective */
extern DSDP_INT getDualObj              ( HSDSolver *dsdpSolver );
extern DSDP_INT getSDPPrimalObjPhaseB   ( HSDSolver *dsdpSolver );

/* Other utilities */
extern double   getMatOneNorm      ( HSDSolver *dsdpSolver, DSDP_INT blockid, DSDP_INT constrid               );
extern DSDP_INT getMatFnorm        ( HSDSolver *dsdpSolver, DSDP_INT blockid, DSDP_INT constrid, double *nrm  );
extern DSDP_INT matRScale          ( HSDSolver *dsdpSolver, DSDP_INT blockid, DSDP_INT constrid, double scaler);
extern DSDP_INT addMattoS          ( HSDSolver *dsdpSolver, DSDP_INT blockid, DSDP_INT constrid, double alpha );
extern DSDP_INT addMattoChecker    ( HSDSolver *dsdpSolver, DSDP_INT blockid, DSDP_INT constrid, double alpha );
extern DSDP_INT addMattodS         ( HSDSolver *dsdpSolver, DSDP_INT blockid, DSDP_INT constrid, double alpha );

// Conic operations
extern double DSDPConic( COPS_GET_A_ONE_NORM ) ( HSDSolver *dsdpSolver );
extern double DSDPConic( COPS_GET_C_ONE_NORM ) ( HSDSolver *dsdpSolver );
extern void   DSDPConic( COPS_DO_C_SCALE     ) ( HSDSolver *dsdpSolver );

#ifdef __cplusplus
}
#endif

#endif /* dsdputils_h */
