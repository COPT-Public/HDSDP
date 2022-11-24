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

extern double csp_sum_abs( int n, int *Ap, double *Ax ) {
    
    double sabs = 0.0;
    for ( int i = 0; i < Ap[n]; ++i ) {
        sabs += fabs(Ax[i]);
    }
    
    return sabs;
}

extern double csp_fro_norm( int n, int *Ap, double *Ax ) {
    
    double fnrm = 0.0;
    for ( int i = 0; i < Ap[n]; ++i ) {
        fnrm += Ax[i] * Ax[i];
    }
    
    return sqrt(fnrm);
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

