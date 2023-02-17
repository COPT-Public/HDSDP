#include "interface/hdsdp_utils.h"

#include "linalg/dense_opts.h"
#include "linalg/vec_opts.h"

#include <math.h>

extern void dsymv( const char *uplo, const int *n, const double *alpha,
                   const double *a, const int *lda, const double *x,
                   const int *incx, const double *beta, double *y, const int *incy );


/* Full dense operations */
extern void fds_symv( int n, double alpha, double *A, double *x, double beta, double *y ) {
    
    char uplo = 'L';
    int one = 1;
    
    dsymv(&uplo, &n, &alpha, A, &n, x, &one, &beta, y, &one);
    
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
    
    double s = ( A[k] > 0 ) ? 1.0 : 0.0;
    double v = sqrt(fabs(A[k]));
    double eps = 0.0;
    
    /* Extract diagonal */
    for ( k = 0; k < n; ++k ) {
        a[i] = PACK_ENTRY(A, n, k, i) / v;
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
    
    return 0;
}
