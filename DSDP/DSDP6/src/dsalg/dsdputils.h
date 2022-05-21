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

#define DUALVAR 0
#define CHECKER 1
#define DELTAS  2

#define COPS_GET_A_ONE_NORM GetDataOneNorm
#define COPS_GET_C_ONE_NORM GetObjOneNorm
#define COPS_DO_C_SCALE     DoObjScale
#define COPS_DO_VAR_SCALE   DoVariableScale
#define COPS_GET_DOBJ       GetDualObjective
#define COPS_CHECK_INCONE   CheckInCone
#define COPS_SYMFAC         SymbolicConicFactorize
#define COPS_GET_SCHUR      GetConicSchur
#define COPS_GET_SCHURVEC   GetConicSchurVec
#define COPS_GET_SLACK      GetConicSlack
#define COPS_CONSTR_EXPR    GetConstrLinearExpr
#define COPS_GET_MAXSTEP    GetMaximumStepsize
#define COPS_STEPDIRECTION  SetupConicDirection
#define COPS_GET_POTENTIAL  GetConicPotentialValue
#define COPS_GET_LOGDET     GetCurrentyLogDet
#define COPS_GET_MERIT      GetConicMeritValue

#endif


/* Schur matrix set up */
extern void     invertDualVars     ( HSDSolver *dsdpSolver );
extern DSDP_INT getSinvASinv       ( HSDSolver *dsdpSolver, DSDP_INT blockid, DSDP_INT constrid,
                                     void *SinvASinv );
extern DSDP_INT getTraceASinvASinv ( HSDSolver *dsdpSolver, DSDP_INT blockid, DSDP_INT constrid,
                                     DSDP_INT mattype, DSDP_INT constrid2, void *SinvASinv );

extern void getBslack           ( HSDSolver *dsdpSolver, vec *y, DSDP_INT type );
extern void getPhaseAS         ( HSDSolver *dsdpSolver, vec *y,  double tau );
extern void getPhaseACheckerS  ( HSDSolver *dsdpSolver, vec *y,  double tau );
extern void getPhaseAdS        ( HSDSolver *dsdpSolver, double drate, vec *dy, double dtau );
extern void getPhaseBS         ( HSDSolver *dsdpSolver, vec *y );
extern void getPhaseBCheckerS  ( HSDSolver *dsdpSolver, vec *y );
extern void getPhaseBdS        ( HSDSolver *dsdpSolver, double alpha, vec *dy, double beta );
extern void dsdpCheckerInCone  ( HSDSolver *dsdpSolver, DSDP_INT *ispsd );
extern void dsdpInCone         ( HSDSolver *dsdpSolver, DSDP_INT *ispsd );
extern double getMaxSDPstep    ( HSDSolver *dsdpSolver, DSDP_INT type   );

/* Objective */
extern DSDP_INT getSDPPrimalObjPhaseB   ( HSDSolver *dsdpSolver );

/* Other utilities */
extern double getMatOneNorm      ( HSDSolver *dsdpSolver, DSDP_INT blockid, DSDP_INT constrid               );
extern void   getMatFnorm        ( HSDSolver *dsdpSolver, DSDP_INT blockid, DSDP_INT constrid, double *nrm  );
extern void   matRScale          ( HSDSolver *dsdpSolver, DSDP_INT blockid, DSDP_INT constrid, double scaler);
extern void   addMattoS          ( HSDSolver *dsdpSolver, DSDP_INT blockid, DSDP_INT constrid, double alpha );
extern void   addMattoChecker    ( HSDSolver *dsdpSolver, DSDP_INT blockid, DSDP_INT constrid, double alpha );
extern void   addMattodS         ( HSDSolver *dsdpSolver, DSDP_INT blockid, DSDP_INT constrid, double alpha );

// Conic operations
extern double   DSDPConic( COPS_GET_A_ONE_NORM ) ( HSDSolver *dsdpSolver );
extern double   DSDPConic( COPS_GET_C_ONE_NORM ) ( HSDSolver *dsdpSolver );
extern void     DSDPConic( COPS_DO_C_SCALE     ) ( HSDSolver *dsdpSolver );
extern void     DSDPConic( COPS_DO_VAR_SCALE   ) ( HSDSolver *dsdpSolver, double tau );
extern void     DSDPConic( COPS_GET_DOBJ       ) ( HSDSolver *dsdpSolver );
extern DSDP_INT DSDPConic( COPS_CHECK_INCONE   ) ( HSDSolver *dsdpSolver, DSDP_INT type );
extern void     DSDPConic( COPS_GET_SCHUR      ) ( HSDSolver *dsdpSolver );
extern void     DSDPConic( COPS_GET_SCHURVEC   ) ( HSDSolver *dsdpSolver, DSDP_INT dInfeas );
extern void     DSDPConic( COPS_SYMFAC         ) ( HSDSolver *dsdpSolver );
extern void     DSDPConic( COPS_GET_SLACK      ) ( HSDSolver *dsdpSolver, DSDP_INT type );
extern void     DSDPConic( COPS_CONSTR_EXPR    ) ( HSDSolver *dsdpSolver, DSDP_INT type,
                                                   double ycoef, vec *y, double tau, double r );
extern double   DSDPConic( COPS_GET_MAXSTEP    ) ( HSDSolver *dsdpSolver, DSDP_INT type );
extern void     DSDPConic( COPS_STEPDIRECTION  ) ( HSDSolver *dsdpSolver );
extern double   DSDPConic( COPS_GET_LOGDET     ) ( HSDSolver *dsdpSolver, vec *y, DSDP_INT *inCone );
#ifdef __cplusplus
}
#endif

#endif /* dsdputils_h */
