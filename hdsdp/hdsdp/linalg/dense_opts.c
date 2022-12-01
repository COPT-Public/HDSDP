#include "dense_opts.h"
#include "hdsdp_utils.h"
#include "vec_opts.h"

#include <math.h>

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
