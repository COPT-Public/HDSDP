#ifndef dsdplapack_h
#define dsdplapack_h

#ifdef __APPLE__
#include "Accelerate/Accelerate.h"
#endif

#include "dsdphsd.h"

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
#define onenrm DASUM_
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
#define fullchol DPOTRF_
#define packchol DPPTRF_
#define fullsolve DPOTRS_
#define packsolve DPPTRS_
#define fullldl DSYTRF_
#define packldl DSPTRF_
#define ldlsolve DSYTRS_

#endif

#if (! defined UNDERBLAS) && (defined CAPBLAS)

/* Blas */
#define dot DDOT
#define axpy DAXPY
#define copy DCOPY
#define norm DNRM2
#define onenrm DASUM
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
#define fullchol DPOTRF
#define packchol DPPTRF
#define fullsolve DPOTRS
#define packsolve DPPTRS
#define fullldl DSYTRF
#define packldl DSPTRF
#define ldlsolve DSYTRS

#endif

#if (defined UNDERBLAS) && (! defined CAPBLAS)

/* Blas */
#define dot ddot_
#define axpy daxpy_
#define copy dcopy_
#define norm dnrm2_
#define onenrm dasum_
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
#define fullchol dpotrf_
#define packchol dpptrf_
#define fullsolve dpotrs_
#define packsolve dpptrs_
#define fullldl dsytrf_
#define packldl dsptrf_
#define ldlsolve dsytrs_

#endif

#if (! defined UNDERBLAS) && (! defined CAPBLAS)

/* Blas */
#define dot ddot
#define axpy daxpy
#define copy dcopy
#define norm dnrm2
#define onenrm dasum
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
#define fullchol dpotrf
#define packchol dpptrf
#define fullsolve dpotrs
#define packsolve dpptrs
#define fullldl dsytrf
#define packldl dsptrf
#define ldlsolve dsytrs

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

extern void copy( const DSDP_INT *n,
                  const double   *x,
                  const DSDP_INT *incx,
                  const double   *y,
                  const DSDP_INT *incy );

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

extern double norm( const DSDP_INT *n,
                    const double   *x,
                    const DSDP_INT *incx );

extern double onenrm( const DSDP_INT *n,
                      const double *x,
                      const DSDP_INT *incx );

extern double fnorm(  const char      *nrm,
                      const char      *uplo,
                      const DSDP_INT  *n,
                      const double    *ap,
                      double          *work );

double dlansy( const char     *nrm,
               const char     *uplo,
               const DSDP_INT *n,
               const double   *a,
               const DSDP_INT *lda,
               double         *work );

extern void solve( const char     *uplo,
                   const char     *trans,
                   const char     *diag,
                   const DSDP_INT *n,
                   const double   *a,
                   const DSDP_INT *lda,
                   double         *x,
                   const DSDP_INT *incx );

void dtrsm( const char     *side,
            const char     *uplo,
            const char     *transa,
            const char     *diag,
            const DSDP_INT *m,
            const DSDP_INT *n,
            const double   *alpha,
            const double   *a,
            const DSDP_INT *lda,
            double         *b,
            const DSDP_INT *ldb );

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

extern DSDP_INT maxidx( const DSDP_INT *n,
                        const double   *x,
                        const DSDP_INT *incx );

extern void chol( const char     *uplo,
                  const DSDP_INT *n,
                  double         *a,
                  const DSDP_INT *lda,
                  DSDP_INT       *info );

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

extern void vecdiv( const DSDP_INT *n,
                    const double   *sa,
                    double         *sx,
                    const DSDP_INT *incx );

extern void vecscal( const DSDP_INT *n,
                     const double   *a,
                     double         *x,
                     const DSDP_INT *incx );

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

void fullchol ( const char     *uplo,
                const DSDP_INT *n,
                double         *a,
                const DSDP_INT *lda,
                DSDP_INT       *info );

void dpotri( const char     *uplo,
             const DSDP_INT *n,
             double         *a,
             const DSDP_INT *lda,
             DSDP_INT       *info );

void dpotrs( const char     *uplo,
             const DSDP_INT *n,
             const DSDP_INT *nrhs,
             const double   *a,
             const DSDP_INT *lda,
             double         *b,
             const DSDP_INT *ldb,
             DSDP_INT       *info );

extern void packchol( const char     *uplo,
                      const DSDP_INT *n,
                      double         *ap,
                      DSDP_INT       *info );

extern void packsolve( const char      *uplo,
                       const DSDP_INT  *n,
                       const DSDP_INT  *nrhs,
                       const double    *ap,
                       const double    *b,
                       const DSDP_INT  *ldb,
                       const DSDP_INT  *info );

extern void packldl( const char     *uplo,
                    const DSDP_INT  *n,
                    double          *ap,
                    DSDP_INT        *ipiv,
                    DSDP_INT        *info );

void fullldl( const char      *uplo,
              const DSDP_INT  *n,
              double          *a,
              const DSDP_INT  *lda,
              DSDP_INT        *ipiv,
              double          *work,
              const DSDP_INT  *lwork,
              DSDP_INT        *info );

void ldlsolve( const char      *uplo,
               const DSDP_INT  *n,
               const DSDP_INT  *nrhs,
               const double    *a,
               const DSDP_INT  *lda,
               const DSDP_INT  *ipiv,
               double          *b,
               const DSDP_INT  *ldb,
               DSDP_INT        *info );

void dsyr( const char     *uplo,
           const DSDP_INT *n,
           const double   *alpha,
           const double   *x,
           const DSDP_INT *incx,
           double         *a,
           const DSDP_INT *lda );

void packr1update( const char *uplo,
           const DSDP_INT *n,
           const double *alpha,
           const double *x,
           const DSDP_INT *incx,
           double *ap );


#ifdef __cplusplus
}
#endif

/* Macros for min/max. */
#ifndef MIN
#define MIN(a, b) (((a) < (b)) ? (a) : (b))
#endif /* MIN */
#ifndef MAX
#define MAX(a, b) (((a) > (b)) ? (a) : (b))
#endif  /* MAX */

/* Define macros related to parameter */
#define DSDP_MAT_UP ('U')
#define DSDP_MAT_LOW ('L')
#define DSDP_MAT_FNORM ('F')
#define DSDP_MAT_NOTRANSPOSE ('N')
#define DSDP_MAT_TRANSPOSE ('T')

/* Other utilities */
#define nsym(x) ((DSDP_INT) (((x) + 1 ) * (x) / 2))

/* Some constants */
#define LWORK   (30)
#define IWORK   (12)

#endif /* dsdplapack_h */
