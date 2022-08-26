/** @file spinfo.h
 *  @brief Header defining internal constants for speigs
 *
 *  @author Wenzhi Gao, Shanghai University of Finance and Economics
 *  @date Aug, 24th, 2022
 *
 */

#ifndef spinfo_h
#define spinfo_h

/* Boolean */
#ifndef TRUE
#define TRUE                 (1)
#define FALSE                (0)
#endif

/* Some constants */
#define ROOT                 (7.0710678118654757273731092936941422522068e-01) ///< \f$ \frac{\sqrt{2}}{2} \f$
#define LAPACK_IWORK         (12)
#define LAPACK_LWORK         (30)

#ifdef __cplusplus
extern "C" {
#endif

/**
 * @brief Lapack dense eigen routine
 *
 * The Lapack eigen routine
 */
extern void dsyevr( const char   *jobz,   const char   *range,  const char  *uplo,
                    const spint  *n,            double *a,      const spint *lda,
                    const double *vl,     const double *vu,     const spint *il,
                    const spint  *iu,     const double *abstol,       spint *m,
                          double *w,            double *z,      const spint *ldz,
                          spint  *isuppz,       double *work,   const spint *lwork,
                          spint  *iwork,  const spint  *liwork, spint *info         );
#ifdef __cplusplus
}
#endif

#endif /* spinfo_h */
