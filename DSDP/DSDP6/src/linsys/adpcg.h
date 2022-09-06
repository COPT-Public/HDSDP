/** @file adpcg.h
 *  @brief Header for basic types and routine list
 *
 * Given a set of positive definite linear systems \f$ A^k x = b^k \f$, adpcg solves them adaptively with
 * pre-conditioning conjugate gradient method.
 *
 * The routine is employed in HDSDP.
 *
 * @author Wenzhi Gao, Shanghai University of Finance and Economics
 * @date Aug 29th, 2022
 *
 */

#ifndef adpcg_h
#define adpcg_h

#include "structs.h"
#include "dsdphsd.h"


#define cgerr(x) printf(x); fatal_error;

/* Return code */
#define CG_OK              (DSDP_RETCODE_OK)
#define CG_ERR             (DSDP_RETCODE_FAILED)

/* Boolean */
#define CG_TRUE            (TRUE)
#define CG_FALSE           (FALSE)

/* Pre-conditioner */
#define CG_PRECOND_DIAG     (10)  ///< Use diagonal as the pre-conditioner
#define CG_PRECOND_CHOL     (11)  ///< Use Cholesky factor as the pre-conditioner

/* Solution status */
#define CG_STATUS_SOLVED    (100) ///< System is solved to desired accuracy within maximum iteration
#define CG_STATUS_MAXITER   (101) ///< System reaches maximum iteration
#define CG_STATUS_FAILED    (102) ///< CG fails to converge
#define CG_STATUS_DIRECT    (103) ///< CG serves as a wrapper for direct solver
#define CG_STATUS_UNKNOWN   (104) ///< CG status is not known. Solution not yet started

/* Auxiliary macros */
#ifndef MIN
#define MIN(a, b) (((a)<(b))?(a):(b))
#endif /* MIN */
#ifndef MAX
#define MAX(a, b) (((a)>(b))?(a):(b))
#endif  /* MAX */

/* Adaptive pre-conditioned conjugate gradient solver */
#ifdef __cplusplus
extern "C" {
#endif

extern DSDP_INT cg_iteration ( adpcg *cg, vec *b, DSDP_INT warm );
extern void     cg_init      ( adpcg *cg );
extern DSDP_INT cg_alloc     ( adpcg *cg, DSDP_INT n, DSDP_INT vsize );
extern void     cg_register  ( adpcg *cg, schurMat *A, vec *diag, schurMat *chol );
extern void     cg_free      ( adpcg *cg );
extern void     cg_setparam  ( adpcg *cg, double tol, DSDP_INT reuse, DSDP_INT maxiter, DSDP_INT restart );
extern void     cg_getstats  ( adpcg *cg, DSDP_INT *status, DSDP_INT *niter, double *rnrm,
                               double *avgsvtime, double *avgfctime, DSDP_INT *nused,
                               DSDP_INT *nmaixter, DSDP_INT *nfactors, DSDP_INT *nrounds,
                               DSDP_INT *nsolverd, DSDP_INT *nsolves );
extern DSDP_INT cg_start     ( adpcg *cg );
extern void     cg_finish    ( adpcg *cg );
extern DSDP_INT cg_solve     ( adpcg *cg, vec *b, vec *x0 );

#ifdef __cplusplus
}
#endif

#endif /* adpcg_h */
