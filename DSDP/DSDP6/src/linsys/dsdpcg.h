#ifndef dsdpcg_h
#define dsdpcg_h

#include "structs.h"

#define CG_PRECOND_DIAG      101
#define CG_PRECOND_CHOL      102
#define CG_STATUS_SOLVED     103
#define CG_STATUS_MAXITER    104
#define CG_STATUS_FAILED     105
#define CG_STATUS_INDEFINITE 106
#define CG_STATUS_UNKNOWN    107

#ifdef __cplusplus
extern "C" {
#endif

extern void dsdpCGinit            ( CGSolver *cgSolver );
extern DSDP_INT dsdpCGAlloc       ( CGSolver *cgSolver, DSDP_INT m );
extern void dsdpCGFree            ( CGSolver *cgSolver );
extern void dsdpCGSetTol          ( CGSolver *cgSolver, double tol );
extern void dsdpCGSetMaxIter      ( CGSolver *cgSolver, DSDP_INT maxiter );
extern void dsdpCGSetM            ( CGSolver *cgSolver, schurMat *M );
extern void dsdpCGSetDPre         ( CGSolver *cgSolver, vec   *preCond );
extern void dsdpCGSetCholPre      ( CGSolver *cgSolver, schurMat *preCond );
extern void dsdpCGSetPType        ( CGSolver *cgSolver, DSDP_INT pType );
extern void dsdpCGprepareP        ( CGSolver *cgSolver );
extern void dsdpCGGetStatus       ( CGSolver *cgSolver, DSDP_INT *status );
extern void dsdpCGGetSolStatistic ( CGSolver *cgSolver, DSDP_INT *status, double *resinorm );
extern void dsdpCGSetPreReuse     ( CGSolver *cgSolver, DSDP_INT reuse );
extern void dsdpCGStoreRHS        ( CGSolver *cgSolver, vec *bin );
extern void dsdpCGRestoreRHS      ( CGSolver *cgSolver, vec *bout );
extern void dsdpCGStoreLdiag      ( CGSolver *cgSolver );
extern void dsdpCGRestoreLdiag    ( CGSolver *cgSolver );
extern DSDP_INT dsdpCGSolve       ( CGSolver *cgSolver, vec *b, vec *x0 );

#ifdef __cplusplus
}
#endif

#endif /* dsdpcg_h */
