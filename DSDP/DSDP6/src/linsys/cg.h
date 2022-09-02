/** @file speigs.h
 *  @brief Header for constants and routine list
 *
 * Implement the pre-conditioned conjugate gradient method
 *
 *  @author Wenzhi Gao, Shanghai University of Finance and Economics
 *  @date Aug, 27th, 2022
 *
 */

#ifndef cg_h
#define cg_h

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

extern void cgInit                      ( cgsol *cgSolver );
extern DSDP_INT cgAlloc                 ( cgsol *cgSolver, DSDP_INT m );
extern void cgFree                      ( cgsol *cgSolver );
extern void cgSetTol                    ( cgsol *cgSolver, double tol );
extern void cgSetMaxIter                ( cgsol *cgSolver, DSDP_INT maxiter );
extern void cgSetM                      ( cgsol *cgSolver, schurMat *M );
extern void cgSetDiagPreconditioner     ( cgsol *cgSolver, vec   *preCond );
extern void cgSetCholeskyPreconditioner ( cgsol *cgSolver, schurMat *preCond );
extern void cgSetPreconditionerType     ( cgsol *cgSolver, DSDP_INT pType );
extern void cgPreparePreconditioner     ( cgsol *cgSolver );
extern void cgGetStatus                 ( cgsol *cgSolver, DSDP_INT *status );
extern void cgGetSolStats               ( cgsol *cgSolver, DSDP_INT *status, double *resinorm );
extern void cgSetPreconditionerReuse    ( cgsol *cgSolver, DSDP_INT reuse );
extern void cgStoreRHS                  ( cgsol *cgSolver, vec *bin );
extern void cgRestoreRHS                ( cgsol *cgSolver, vec *bout );
extern DSDP_INT cgsolve                      ( cgsol *cgSolver, vec *b, vec *x0 );

#ifdef __cplusplus
}
#endif

#endif /* cg_h */
