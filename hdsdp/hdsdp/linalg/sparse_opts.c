/** @brief Implement the sparse matrix operations that used in HDSDP
 *  @author Wenzhi Gao
 *
 */

#include "sparse_opts.h"
#include "vec_opts.h"

#include <math.h>

/* Compressed column operations */
extern void csp_Axpby( int n, int *Ap, int *Ai, double *Ax, double a, double *x, double *y ) {
    
    if ( a == 0.0 ) {
        return;
    }
    
    for ( int i = 0, j; i < n; ++i ) {
        for ( j = Ap[i]; j < Ap[i + 1]; ++j ) {
            y[Ai[j]] += a * x[i] * Ax[j];
        }
    }
    
    return;
}


extern void csp_ATxpby( int n, int *Ap, int *Ai, double *Ax, double a, double *x, double *y ) {
    
    if ( a == 0.0 ) {
        return;
    }
    
    double aTy = 0.0;
    
    for ( int i = 0, j; i < n; ++i ) {
        for ( j = Ap[i], aTy = 0.0; j < Ap[i + 1]; ++j ) {
            aTy += x[Ai[j]] * Ax[j];
        }
        y[i] += a * aTy;
    }
    
    return;
}

extern double csp_sum_abs( int n, int *Ap, int *Ai, double *Ax ) {
    
    double sabs = 0.0;
    for ( int i = 0, j; i < Ap[n]; ++i ) {
        for ( j = Ap[i]; j < Ap[i + 1]; ++i ) {
            sabs += ( Ai[i] == i ) ? 0.5 * fabs(Ax[i]) : fabs(Ax[i]);
        }
    }
    
    return sabs;
}

extern double csp_fro_norm( int n, int *Ap, int *Ai, double *Ax ) {
    
    double nrm = 0.0;
    for ( int i = 0, j; i < Ap[n]; ++i ) {
        for ( j = Ap[i]; j < Ap[i + 1]; ++i ) {
            nrm += ( Ai[i] == i ) ? 0.5 * Ax[i] * Ax[i]: Ax[i] * Ax[i];
        }
    }
    
    return sqrt(nrm);
}

/** @brief Column sparse matrix aApB
 *  This function is called in the following contexts:
 *
 *  1. Updating dual variable S + a \* dS
 *
 */
extern void csp_aApB( int n, int nnz, double a, int *Al, double *Ax, double *Bx ) {
    
    if ( Al ) {
        /* In this case B is sparse and Al indicates where each Ax is located in Bx */
        for ( int i = 0; i < nnz; ++i ) {
            Bx[i] += a * Ax[Al[i]];
        }
        
    } else {
        /* In this case B is dense and axpy is sufficient*/
        
    }
    
    return;
}

/** @brief Get the number of nonzero columns in an csc matrix
 *
 */
extern int csp_nnz_cols ( int n, int *Ap ) {
    
    int nzcols = 0;
    
    for ( int i = 0; i < n; ++i ) {
        nzcols += ( Ap[i + 1] - Ap[i] > 0 );
    }
    
    return nzcols;
}

extern void csp_dump( int n, int *Ap, int *Ai, double *Ax, double *v ) {
    
    for ( int i = 0, j; i < n; ++i ) {
        for ( j = Ap[i]; j < Ap[i + 1]; ++j ) {
            v[i * n + Ai[j]] = v[Ai[j] * n + i] = Ax[j];
        }
    }
    
    return;
}


/* Lower triplet operations */
/** @brief Scale a triplet matrix
 *
 */
extern void tsp_scal( double a, int nnz, double *Ax ) {
    
    int incx = 1;
    scal(&nnz, &a, Ax, &incx);
    
    return;
}

/** @brief Compute one norm of a triplet matrix
 *
 */
extern double tsp_sum_abs( int nnz, int *Ai, int *Aj, double *Ax ) {
    
    double nrm = 0.0;
    for ( int i = 0; i < nnz; ++i ) {
        nrm += ( Ai[i] == Aj[i] ) ? 0.5 * fabs(Ax[i]) : fabs(Ax[i]);
    }
    
    return 2.0 * nrm;
}

extern double tsp_fro_norm( int nnz, int *Ai, int *Aj, double *Ax ) {
    
    double nrm = 0.0;
    for ( int i = 0; i < nnz; ++i ) {
        nrm += ( Ai[i] == Aj[i] ) ? 0.5 * Ax[i] * Ax[i]: Ax[i] * Ax[i];
    }
    
    return sqrt(2.0 * nrm);
}

extern void tsp_dump( int n, int nnz, int *Ai, int *Aj, double *Ax, double *v ) {
    
    int i, j;
    for ( int k = 0; k < nnz; ++k ) {
        i = Ai[k]; j = Aj[k];
        v[n * i + j] = v[n * j + i] = Ax[k];
    }
    
    return;
}

