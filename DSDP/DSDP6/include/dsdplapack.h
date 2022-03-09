#ifndef dsdplapack_h
#define dsdplapack_h

#ifdef __APPLE__
#include "Accelerate/Accelerate.h"
#endif

#include "dsdphsd.h"

/*
Declaration of Lapack and Blas routines in the implemenatation of the
first-order algorithms
*/

/*
UNDERBLAS: Whether routines are given by routinename_(params) or routinename(params)
CAPBLAS  : Whether routines are given by ROUTINENAME(params) or routinename(params)
*/

/*
#ifdef UNDERBLAS
#undef UNDERBLAS
#endif

#ifdef CAPBLAS
#undef CAPBLAS
#endif
*/

#ifndef NOBLAS
#if (defined UNDERBLAS) && (defined CAPBLAS)

/* Blas */
#define dot DDOT_
#define axpy DAXPY_
#define copy DCOPY_
#define norm DNRM2_
#define solve DTRSV_
#define matvec DGEMV_
#define symatvec DSYMV_
#define symmatmat DSYMM_
#define packmatvec DSPMV_
#define vecscal DSCAL_
#define vecdiv DRSCL_
#define maxidx IDAMAX_
#define packr1update DSPR_

/* Lapack */
#define chol DPOTRF_
#define ata DSYRK_
#define computescal DGEEQU_
#define matscal DLAQGE_
#define computeols DGELS_
#define eigs DSYEVX_
#define eigs2 DSYEVR_
#define fnorm DLANSP_
#define packchol DPPTRF_
#define packsolve DPPTRS_
#define packldl DSPTRF_
#define ldlsolve DSPTRS_

#endif

#if (! defined UNDERBLAS) && (defined CAPBLAS)

/* Blas */
#define dot DDOT
#define axpy DAXPY
#define copy DCOPY
#define norm DNRM2
#define solve DTRSV
#define matvec DGEMV
#define symmatmat DSYMM
#define symatvec DSYMV
#define packmatvec DSPMV
#define vecscal DSCAL
#define vecdiv DRSCL
#define maxidx IDAMAX
#define packr1update DSPR

/* Lapack */
#define chol DPOTRF
#define ata DSYRK
#define computescal DGEEQU
#define matscal DLAQGE
#define computeols DGELS
#define eigs DSYEVX
#define eigs2 DSYEVR
#define fnorm DLANSP
#define packchol DPPTRF
#define packsolve DPPTRS
#define packldl DSPTRF
#define ldlsolve DSPTRS

#endif

#if (defined UNDERBLAS) && (! defined CAPBLAS)

/* Blas */
#define dot ddot_
#define axpy daxpy_
#define copy dcopy_
#define norm dnrm2_
#define solve dtrsv_
#define matvec dgemv_
#define symatvec dsymv_
#define symmatmat dsymm_
#define packmatvec dspmv_
#define vecscal dscal_
#define vecdiv drscl_
#define maxidx idamax_
#define packr1update dspr_

/* Lapack */
#define chol dpotrf_
#define ata dsyrk_
#define computescal dgeequ_
#define matscal dlaqge_
#define computeols dgels_
#define eigs dsyevx_
#define eigs2 dsyevr_
#define fnorm dlansp_
#define packchol dpptrf_
#define packsolve dpptrs_
#define packldl dsptrf_
#define ldlsolve dsptrs_

#endif

#if (! defined UNDERBLAS) && (! defined CAPBLAS)

/* Blas */
#define dot ddot
#define axpy daxpy
#define copy dcopy
#define norm dnrm2
#define solve dtrsv
#define matvec dgemv
#define symatvec dsymv
#define symmatmat dsymm
#define packmatvec dspmv
#define vecscal dscal
#define vecdiv drscl
#define maxidx idamax
#define packr1update dspr

/* Lapack */
#define chol dpotrf
#define ata dsyrk
#define computescal dgeequ
#define matscal dlaqge
#define computeols dgels
#define eigs dsyevx
#define eigs2 dsyevr
#define fnorm dlansp
#define packchol dpptrf
#define packsolve dpptrs
#define packldl dsptrf
#define ldlsolve dsptrs

#endif

#endif

#ifndef TRANS
#define TRANS
#define mattrans MKL_Dimatcopy
#endif

#ifdef __cplusplus
extern "C" {
#endif

/* Routine Declaration */
extern double dot( const DSDP_INT *n,
                   const double   *x,
                   const DSDP_INT *incx,
                   const double   *y,
                   const DSDP_INT *incy );

/*
 DDOT forms the dot product of two vectors.
 uses unrolled loops for increments equal to one.
*/

extern void copy( const DSDP_INT *n,
                  const double   *x,
                  const DSDP_INT *incx,
                  const double   *y,
                  const DSDP_INT *incy );
/*
 DCOPY copies a vector, x, to a vector, y.
 uses unrolled loops for increments equal to 1.
*/


extern void axpy( const DSDP_INT *n,
                  const double   *alpha,
                  const double   *x,
                  const DSDP_INT *incx,
                  double         *y,
                  const DSDP_INT *incy );

extern void dger( const DSDP_INT *m,
                  const DSDP_INT *n,
                  const double   *alpha,
                  const double   *x,
                  const DSDP_INT *incx,
                  const double   *y,
                  const DSDP_INT *incy,
                  double         *a,
                  const DSDP_INT *lda );

/*
 DAXPY constant times a vector plus a vector.
 uses unrolled loops for increments equal to one.
*/

extern double norm( const DSDP_INT *n,
                    const double   *x,
                    const DSDP_INT *incx );

extern double fnorm(  const char      *nrm,
                      const char      *uplo,
                      const DSDP_INT  *n,
                      const double    *ap,
                      double          *work );

/*
 DNRM2 returns the euclidean norm of a vector via the function
 name, so that

    DNRM2 := sqrt( x'*x )
*/

extern void solve( const char     *uplo,
                   const char     *trans,
                   const char     *diag,
                   const DSDP_INT *n,
                   const double   *a,
                   const DSDP_INT *lda,
                   double         *x,
                   const DSDP_INT *incx );

/*
 DTRSV  solves one of the systems of equations

    A*x = b,   or   A**T*x = b,

 where b and x are n element vectors and A is an n by n unit, or
 non-unit, upper or lower triangular matrix.

 No test for singularity or near-singularity is included in this
 routine. Such tests must be performed before calling this routine.
*/

extern void matvec( const char     *trans,
                    const DSDP_INT *m,
                    const DSDP_INT *n,
                    const double   *alpha,
                    const double   *a,
                    const DSDP_INT *lda,
                    const double   *x,
                    const DSDP_INT *incx,
                    const double   *beta,
                    double         *y,
                    const DSDP_INT *incy );

/*
 DGEMV  performs one of the matrix-vector operations

    y := alpha*A*x + beta*y,   or   y := alpha*A**T*x + beta*y,

 where alpha and beta are scalars, x and y are vectors and A is an
 m by n matrix.
*/

extern void symatvec( const char     *uplo,
                      const DSDP_INT *n,
                      const double   *alpha,
                      const double   *a,
                      const DSDP_INT *lda,
                      const double   *x,
                      const DSDP_INT *incx,
                      const double   *beta,
                      double         *y,
                      const DSDP_INT *incy );
/*
 DSYMV  performs the matrix-vector  operation

    y := alpha*A*x + beta*y,

 where alpha and beta are scalars, x and y are n element vectors and
 A is an n by n symmetric matrix.
*/

extern void packmatvec( const char *uplo,
                        const DSDP_INT *n,
                        const double *alpha,
                        const double *ap,
                        const double *x,
                        const DSDP_INT *incx,
                        const double *beta,
                        double *y,
                        const DSDP_INT *incy );

extern void dsymm( const char     *side,
                   const char     *uplo,
                   const DSDP_INT *m,
                   const DSDP_INT *n,
                   const double   *alpha,
                   const double   *a,
                   const DSDP_INT *lda,
                   const double   *b,
                   const DSDP_INT *ldb,
                   const double   *beta,
                   double         *c,
                   const DSDP_INT *ldc);

extern DSDP_INT maxidx( const DSDP_INT *n, const double *x, const DSDP_INT *incx );
/*
 IDAMAX finds the index of the first element having maximum absolute value.
*/

extern void chol( const char     *uplo,
                  const DSDP_INT *n,
                  double         *a,
                  const DSDP_INT *lda,
                  DSDP_INT       *info );

/*
 DPOTRF computes the Cholesky factorization of a real symmetric
 positive definite matrix A.

 The factorization has the form
    A = U**T * U,  if UPLO = 'U', or
    A = L  * L**T,  if UPLO = 'L',
 where U is an upper triangular matrix and L is lower triangular.

 This is the block version of the algorithm, calling Level 3 BLAS.
*/

extern void ata( const char     *uplo,
                 const char     *trans,
                 const DSDP_INT *n,
                 const DSDP_INT *k,
                 const double   *alpha,
                 const double   *a,
                 const DSDP_INT *lda,
                 const double   *beta,
                 double         *c,
                 const DSDP_INT *ldc );

/*
 DSYRK  performs one of the symmetric rank k operations

     C := alpha*A*A**T + beta*C,

  or

     C := alpha*A**T*A + beta*C,

  where  alpha and beta  are scalars, C is an  n by n  symmetric matrix
  and  A  is an  n by k  matrix in the first case and a  k by n  matrix
  in the second case.
*/

extern void vecdiv( const DSDP_INT *n,
                    const double   *sa,
                    double         *sx,
                    const DSDP_INT *incx );

/*
 DRSCL multiplies an n-element real vector x by the real scalar 1/a.
 This is done without overflow or underflow as long as
 the final result x/a does not overflow or underflow.
*/

extern void vecscal( const DSDP_INT *n,
                     const double   *a,
                     double         *x,
                     const DSDP_INT *incx );
/*
 DSCAL scales a vector by a constant.
 uses unrolled loops for increment equal to 1.
*/

extern void computeols( const char     *trans,
                        const DSDP_INT *m,
                        const DSDP_INT *n,
                        const DSDP_INT *nrhs,
                        double         *a,
                        const DSDP_INT *lda,
                        double         *b,
                        const DSDP_INT *ldb,
                        double         *work,
                        const DSDP_INT *lwork,
                        DSDP_INT       *info );

/*
 DGELS solves overdetermined or underdetermined real linear systems
 involving an M-by-N matrix A, or its transpose, using a QR or LQ
 factorization of A.  It is assumed that A has full rank.
*/

extern void eigs( const char     *jobz,
                  const char     *range,
                  const char     *uplo,
                  const DSDP_INT *n,
                  double         *a,
                  const DSDP_INT *lda,
                  const double   *vl,
                  const double   *vu,
                  const DSDP_INT *il,
                  const DSDP_INT *iu,
                  const double   *abstol,
                  DSDP_INT       *m,
                  double         *w,
                  double         *z,
                  const DSDP_INT *ldz,
                  double         *work,
                  const DSDP_INT *lwork,
                  DSDP_INT       *iwork,
                  DSDP_INT       *ifail,
                  DSDP_INT       *info );
/*
 DSYEVX computes selected eigenvalues and, optionally, eigenvectors
 of a real symmetric matrix A.  Eigenvalues and eigenvectors can be
 selected by specifying either a range of values or a range of indices
 for the desired eigenvalues.
*/

extern void eigs2( const char     *jobz,
                   const char     *range,
                   const char     *uplo,
                   const DSDP_INT *n,
                   double         *a,
                   const DSDP_INT *lda,
                   const double   *vl,
                   const double   *vu,
                   const DSDP_INT *il,
                   const DSDP_INT *iu,
                   const double   *abstol,
                   DSDP_INT       *m,
                   double         *w,
                   double         *z,
                   const DSDP_INT *ldz,
                   DSDP_INT       *isuppz,
                   double         *work,
                   const DSDP_INT *lwork,
                   DSDP_INT       *iwork,
                   const DSDP_INT *liwork,
                   DSDP_INT       *info );
/*
 DSYEVR computes selected eigenvalues and, optionally, eigenvectors
 of a real symmetric matrix A.  Eigenvalues and eigenvectors can be
 selected by specifying either a range of values or a range of
 indices for the desired eigenvalues.

 DSYEVR first reduces the matrix A to tridiagonal form T with a call
 to DSYTRD.  Then, whenever possible, DSYEVR calls DSTEMR to compute
 the eigenspectrum using Relatively Robust Representations.  DSTEMR
 computes eigenvalues by the dqds algorithm, while orthogonal
 eigenvectors are computed from various "good" L D L^T representations
 (also known as Relatively Robust Representations). Gram-Schmidt
 orthogonalization is avoided as far as possible. More specifically,
 the various steps of the algorithm are as follows.

 For each unreduced block (submatrix) of T,
    (a) Compute T - sigma I  = L D L^T, so that L and D
        define all the wanted eigenvalues to high relative accuracy.
        This means that small relative changes in the entries of D and L
        cause only small relative changes in the eigenvalues and
        eigenvectors. The standard (unfactored) representation of the
        tridiagonal matrix T does not have this property in general.
    (b) Compute the eigenvalues to suitable accuracy.
        If the eigenvectors are desired, the algorithm attains full
        accuracy of the computed eigenvalues only right before
        the corresponding vectors have to be computed, see steps c) and d).
    (c) For each cluster of close eigenvalues, select a new
        shift close to the cluster, find a new factorization, and refine
        the shifted eigenvalues to suitable accuracy.
    (d) For each eigenvalue with a large enough relative separation compute
        the corresponding eigenvector by forming a rank revealing twisted
        factorization. Go back to (c) for any clusters that remain.

 The desired accuracy of the output can be specified by the input
 parameter ABSTOL.

 For more details, see DSTEMR's documentation and:
 - Inderjit S. Dhillon and Beresford N. Parlett: "Multiple representations
   to compute orthogonal eigenvectors of symmetric tridiagonal matrices,"
   Linear Algebra and its Applications, 387(1), pp. 1-28, August 2004.
 - Inderjit Dhillon and Beresford Parlett: "Orthogonal Eigenvectors and
   Relative Gaps," SIAM Journal on Matrix Analysis and Applications, Vol. 25,
   2004.  Also LAPACK Working Note 154.
 - Inderjit Dhillon: "A new O(n^2) algorithm for the symmetric
   tridiagonal eigenvalue/eigenvector problem",
   Computer Science Division Technical Report No. UCB/CSD-97-971,
   UC Berkeley, May 1997.


 Note 1 : DSYEVR calls DSTEMR when the full spectrum is requested
 on machines which conform to the ieee-754 floating point standard.
 DSYEVR calls DSTEBZ and DSTEIN on non-ieee machines and
 when partial spectrum requests are made.

 Normal execution of DSTEMR may create NaNs and infinities and
 hence may abort due to a floating point exception in environments
 which do not handle NaNs and infinities in the ieee standard default
 manner.
*/

extern void packchol( const char *uplo, const DSDP_INT *n, double *ap, DSDP_INT *info );
extern void packsolve( const char      *uplo,
                       const DSDP_INT  *n,
                       const DSDP_INT  *nrhs,
                       const double    *ap,
                       const double    *b,
                       const DSDP_INT  *ldb,
                       const DSDP_INT  *info );

extern void packldl( const char      *uplo,
                    const DSDP_INT  *n,
                    double          *ap,
                    DSDP_INT        *ipiv,
                    DSDP_INT        *info );

void ldlsolve( const char      *uplo,
             const DSDP_INT  *n,
             const DSDP_INT  *nrhs,
             const double    *ap,
             const DSDP_INT  *ipiv,
             double          *b,
             const DSDP_INT  *ldb,
             DSDP_INT        *info );

void packr1update( const char *uplo,
           const DSDP_INT *n,
           const double *alpha,
           const double *x,
           const DSDP_INT *incx,
           double *ap );

/* The routines below are depreciated currently */
extern void computescal( const DSDP_INT *m,
                         const DSDP_INT *n,
                         const double   *a,
                         const DSDP_INT *lda,
                         double         *r,
                         double         *c,
                         double         *rowcnd,
                         double         *colcnd,
                         double         *amax,
                         DSDP_INT       *info );

void matscal( const DSDP_INT *m,
              const DSDP_INT *n,
              double         *a,
              const DSDP_INT *lda,
              const double   *r,
              const double   *c,
              const double   *rowcnd,
              const double   *colcnd,
              const double   *amax,
              char           *equed );


void mattrans( const char    ordering,
               const char    trans,
               size_t        rows,
               size_t        cols,
               const double  alpha,
               double        *AB,
               size_t        lda,
               size_t        ldb );

/* Intel MKL matrix transposition routine */


#ifdef __cplusplus
}
#endif

/* Define macros related to parameter */
#define DSDP_MAT_UP ('U')
#define DSDP_MAT_LOW ('L')
#define DSDP_MAT_FNORM ('F')
#define DSDP_MAT_NOTRANSPOSE ('N')
#define DSDP_MAT_TRANSPOSE ('T')

/* Other utilities */
#define nsym(x) ((DSDP_INT) (((x) + 1 ) * (x) / 2))

#endif /* dsdplapack_h */
