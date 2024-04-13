#ifdef HEADERPATH
#include "interface/hdsdp_utils.h"
#include "linalg/dense_opts.h"
#include "linalg/vec_opts.h"
#else
#include "hdsdp_utils.h"
#include "dense_opts.h"
#include "vec_opts.h"
#endif

#include <math.h>

extern void dsymv( const char *uplo, const int *n, const double *alpha,
                   const double *a, const int *lda, const double *x,
                   const int *incx, const double *beta, double *y, const int *incy );

extern void dgemv( char *trans, int *m, int *n, double *alpha,
                   double *a, int *lda, double *x, int *incx,
                   double *beta, double *y, int *incy );

extern void dsyevr( const char *jobz, const char *range, const char *uplo,
                    const int  *n, double *a, const int *lda,
                    const double *vl, const double *vu, const int *il,
                    const int *iu, const double *abstol, int *m,
                    double *w, double *z, const int *ldz, int *isuppz,
                    double *work, const int *lwork, int *iwork, const int *liwork,
                    int *info );

extern void dspmv( const char *uplo, const int *n, const double *alpha,
                   const double *ap, const double *x, const int *incx,
                   const double *beta, double *y, const int *incy );

extern void dsymm( const char *side, const char *uplo, const int *m,
                   const int *n, const double *alpha, const double *a,
                   const int *lda, const double *b, const int *ldb,
                   const double *beta, double *c, const int *ldc );

void dsyr( const char *uplo, const int *n, const double *alpha,
           const double *x, const int *incx, double *a, const int *lda );

extern void dger( const int *m, const int *n, const double *alpha,
                  const double *x, const int *incx, const double *y,
                  const int *incy, double *a, const int *lda );

/* Full dense operations */
extern void fds_symv( int n, double alpha, double *A, double *x, double beta, double *y ) {
    
    dsymv(&HCharConstantUploLow, &n, &alpha, A, &n, x, &HIntConstantOne, &beta, y, &HIntConstantOne);
    
    return;
}

extern int fds_syev( int n, double *U, double *d, double *Y, int m,
                      double *work, int *iwork, int lwork, int liwork ) {
    
    int retcode = HDSDP_RETCODE_OK;
    
    char jobz = 'V', range = 'I', uplo = HCharConstantUploUp;
    int isuppz[4] = {0};
    int il = n - m + 1, iu = n;
    int info = 0;
    
    dsyevr(&jobz, &range, &uplo, &n, U, &n,
           &HDblConstantZero, &HDblConstantZero,
           &il, &iu, &HDblConstantZero, &m, d, Y,
           &n, isuppz, work, &lwork, iwork, &liwork, &info);
    
    if ( info != 0 ) {
        retcode = HDSDP_RETCODE_FAILED;
    }
    
    return retcode;
}

extern void fds_gemv( int m, int n, double *M, double *v, double *y ) {
    
    dgemv(&HCharConstantNoTrans, &m, &n, &HDblConstantOne, M, &m,
          v, &HIntConstantOne, &HDblConstantZero, y, &HIntConstantOne);
    
    return;
}

extern void fds_symm( char side, char uplo, int m, int n, double alpha, double *a, int lda,
                      double *b, int ldb, double beta, double *c, int ldc ) {
    
    dsymm(&side, &uplo, &m, &n, &alpha, a, &lda, b, &ldb, &beta, c, &ldc);
    
    return;
}

extern void fds_ger( int m, int n, double alpha, double *x, int incx,
                    double *y, int incy, double *a, int lda) {
    
    dger(&m, &n, &alpha, x, &incx, y, &incy, a, &lda);
    
    return;
}

extern void fds_trimultiply( int n, double *S, double *X, double *aux, double *XSX ) {
    /* Routine for multiplying three dense matrices X * S * X and adding it to buffer
       Check dataMatDenseKKT3ComputeSinvASinvImpl for more details */
    
    double *XCol = NULL;
    double *SXCol = NULL;
    double *SX = aux;
    
    HDSDP_ZERO(SX, double, n * n);
    for ( int i = 0; i < n; ++i ) {
        XCol = X + n * i;
        SXCol = SX + n * i;
        fds_symv(n, 1.0, S, XCol, 0.0, SXCol);
    }
    
    for ( int i = 0; i < n; ++i ) {
        SXCol = SX + n * i;
        for ( int j = 0; j < i; ++j ) {
            XCol = X + n * j;
            double dDotVal = dot(&n, SXCol, &HIntConstantOne, XCol, &HIntConstantOne);
            FULL_ENTRY(XSX, n, i, j) += dDotVal;
            FULL_ENTRY(XSX, n, j, i) += dDotVal;
        }
        
        XCol = X + n * i;
        FULL_ENTRY(XSX, n, i, i) += \
        dot(&n, SXCol, &HIntConstantOne, XCol, &HIntConstantOne);
    }

    return;
}

extern double fds_dot_fds( int n, double *A, double *B ) {
    
    double dAdotB = 0.0;
    
    int j = 0;
    int i = 0;
    
    double *dACol = A;
    double *dBCol = B;
    
    for ( i = 0; i < n; ++i ) {
        /* The first element of each column is diagonal */
        dAdotB += 0.5 * dACol[0] * dBCol[0];
        for ( j = 1; j < n - i; ++j ) {
            dAdotB += dACol[j] * dBCol[j];
        }
        /* Move pointer to the start of each column */
        dACol += n + 1;
        dBCol += n + 1;
    }
    
    return 2.0 * dAdotB;
}

extern void fds_print( int n, double *A ) {
    
    for ( int i = 0; i < n; ++i ) {
        for ( int j = 0; j < n; ++j ) {
            printf("%6.3e ", A[i + n * j]);
        }
        printf("\n");
    }

    return;
}

/** @brief Scale a packed dense matrix
 *
 */
extern void pds_scal( double a, int n, double *A ) {
    
    int nnz = PACK_NNZ(n), incx = 1;
    scal(&nnz, &a, A, &incx);
    
    return;
}

extern double pds_sum_abs( int n, double *A ) {
    
    double nrm = 0.0, *p = A;
    int nrow = n, incx = 1;
    
    for ( int i = 0; i < n; ++i ) {
        nrm -= 0.5 * fabs(p[0]);
        nrm += nrm1(&nrow, p, &incx);
        p += n - i; nrow -= 1;
    }
    
    return 2.0 * nrm;
}

extern double pds_fro_norm( int n, double *A ) {
    
    double nrm = 0.0, colnrm, *p = A;
    int nrow = n, incx = 1;
    
    for ( int i = 0; i < n; ++i ) {
        nrm -= 0.5 * p[0] * p[0];
        colnrm = nrm2(&nrow, p, &incx);
        nrm += colnrm * colnrm;
        p += n - i; nrow -= 1;
    }
    
    return sqrt(2.0 * nrm);
}

extern void pds_dump( int n, double *A, double *v ) {
    
    for ( int i = 0, j; i < n; ++i ) {
        for ( j = i; j < n; ++j ) {
            FULL_ENTRY(v, n, i, j) = FULL_ENTRY(v, n, j, i) = PACK_ENTRY(A, n, i, j);
        }
    }
    
    return;
}

extern void pds_decompress( int nnz, int *Ci, double *Cx, double *A ) {
    
    for ( int k = 0; k < nnz; ++k ) {
        A[Ci[k]] = Cx[k];
    }
    
    return;
}

/** @brief Check is the matrix is rank one
 *
 */
extern int pds_r1_extract( int n, double *A, double *sgn, double *a ) {
    
    /* Get the first nonzero */
    int i, j, k = 0;
    for ( i = 0; i < n; ++i ) {
        if ( A[k] != 0 ) {
            break;
        }
        k += n - i;
    }
    
    if ( i == n ) {
        return 0;
    }
    
    double s = ( A[k] > 0 ) ? 1.0 : -1.0;
    double v = sqrt(fabs(A[k]));
    double eps = 0.0;
    
    /* Extract diagonal */
    for ( k = 0; k < n; ++k ) {
        a[k] = PACK_ENTRY(A, n, k, i) / v;
    }
    
    int id = 0;
    double *pstart = NULL;
    if ( s == 1.0 ) {
        for ( i = 0; i < n; ++i ) {
            pstart = A + id;
            for ( j = 0; j < n - i; ++j ) {
                eps += fabs(pstart[j] - a[i] * a[i + j]);
            }
            id += n - i;
            if ( eps > 1e-10 ) {
                return 0;
            }
        }
    } else {
        for ( i = 0; i < n; ++i ) {
            pstart = A + id;
            for ( j = 0; j < n - i; ++j ) {
                eps += fabs(pstart[j] + a[i] * a[i + j]);
            }
            id += n - i;
            if ( eps > 1e-10 ) {
                return 0;
            }
        }
    }
    
    *sgn = s;
    return 1;
}

extern double pds_quadform( int n, double *A, double *v, double *aux ) {
    
    dspmv(&HCharConstantUploLow, &n, &HDblConstantOne, A, v,
          &HIntConstantOne, &HDblConstantZero, aux, &HIntConstantOne);
    
    return dot(&n, aux, &HIntConstantOne, v, &HIntConstantOne);
}

extern void pds_spmv( char uplo, int n, double alpha, double *ap, double *x, int incx,
                      double beta, double *y, int incy ) {
    
    dspmv(&uplo, &n, &alpha, ap, x, &incx, &beta, y, &incy);
    
    return;
}

extern void pds_syr( char uplo, int n, double alpha, double *x, int incx, double *a, int lda ) {
    
    dsyr(&uplo, &n, &alpha, x, &incx, a, &lda);
    
    return;
}
