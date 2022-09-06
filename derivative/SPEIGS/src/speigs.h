/** @file speigs.h
 *  @brief Header for basic types and routine list
 *
 * Implement the eigen-decomposition algorithm from DSDP5.8 by Steve Benson.
 *
 * Given a real symmetric matrix A, the routine explores special structures
 * within and computes the full eigen-decomposition of the matrix.
 * In the backend the routine calls Lapack dsyev (or Netlib Eispack) to
 * decompose the pre-processed system.
 *
 * This routine is also employed in HDSDP solver for SDP.
 *
 *  @author Wenzhi Gao, Shanghai University of Finance and Economics
 *  @date Aug, 24th, 2022
 *
 */

#ifndef speigs_h
#define speigs_h

#include <stddef.h>

#ifdef MATLAB_MEX_FILE
#include "mex.h"
typedef mwSize spint;
#define sperr mexErrMsgTxt
#define id "%lld"
#else
#ifdef SPEIG_64
typedef int64_t spint;
#define id "%lld"
#else
typedef int32_t spint;
#define id "%d"
#endif
#define sperr(x) printf(x);
#endif

/* Return code */
#define SP_EIGS_OK           (0)
#define SP_EIGS_ERR          (1)

/* Matrix type */
#define MATRIX_TYPE_ZERO     (0)
#define MATRIX_TYPE_SPARSE   (1)
#define MATRIX_TYPE_GENERAL  (2)
#define MATRIX_TYPE_RANKONE  (3)
#define MATRIX_TYPE_DIAG     (4)
#define MATRIX_TYPE_TWOTWO   (5)

#ifdef __cplusplus
extern "C" {
#endif
extern spint speigs_analyze( spint *Ap,    spint *Ai,     double *Ax,   spint  *dim,
                             spint *iwork, spint *liwork, double *work, spint  *lwork,
                             spint *type,  spint *sn,     double tol,   double gthresh );

extern spint speigs_factorize( spint  *Ap,     spint  *Ai,    double *Ax,    spint  *dim,   spint *aiwork,
                               double *awork,  spint  *type,  spint  *sn,    spint  *iwork, spint *liwork,
                               double *work,   spint  *lwork, double *evals, double *evecs,
                               spint  *rank,   double tol );
#ifdef __cplusplus
}
#endif

#define SPEIG_VER  (1) // Version number

#endif /* speigs_h */
