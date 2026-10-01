/** @brief Implement the sparse matrix operations that used in HDSDP
 *  @author Wenzhi Gao
 *
 */

#ifdef HEADERPATH
#include "linalg/sparse_opts.h"
#include "linalg/vec_opts.h"
#include "interface/hdsdp_utils.h"
#else
#include "sparse_opts.h"
#include "vec_opts.h"
#include "hdsdp_utils.h"
#endif

#include <assert.h>
#include <stdio.h>
#include <math.h>

extern void dgemm( const char *transa, const char *transb, const int *m,
                   const int *n, const int *k, const double *alpha,
                   const double *a, const int *lda, const double *b,
                   const int *ldb, const double *beta, double *c,
                   const int *ldc );

/* Compressed column operations */
extern void csp_Axpy( int n, int *Ap, int *Ai, double *Ax, double a, double *x, double *y ) {
    
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

extern void csp_ATxpy( int n, int *Ap, int *Ai, double *Ax, double a, double *x, double *y ) {
    
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
    for ( int i = 0, j; i < n; ++i ) {
        for ( j = Ap[i]; j < Ap[i + 1]; ++j ) {
            sabs += ( Ai[j] == i ) ? 0.5 * fabs(Ax[j]) : fabs(Ax[j]);
        }
    }
    
    return sabs;
}

extern double csp_fro_norm( int n, int *Ap, int *Ai, double *Ax ) {
    
    double nrm = 0.0;
    for ( int i = 0, j; i < n; ++i ) {
        for ( j = Ap[i]; j < Ap[i + 1]; ++j ) {
            nrm += ( Ai[j] == i ) ? 0.5 * Ax[j] * Ax[j]: Ax[j] * Ax[j];
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
        axpy(&n, &a, Ax, &HIntConstantOne, Bx, &HIntConstantOne);
    }
    
    return;
}

/** @brief Get the number of nonzero columns in an csc matrix
 *
 */
extern int csp_nnz_cols( int n, int *Ap ) {
    
    int nzcols = 0;
    
    for ( int i = 0; i < n; ++i ) {
        nzcols += ( Ap[i + 1] - Ap[i] > 0 );
    }
    
    return nzcols;
}

extern void csp_trimultiply( int n, int *Sp, int *Si, double *Sx, double *X, double *aux, double *XSX ) {
    /* Routine for multiplying X * S * X (S is csc sparse) and adding it to buffer
       Check dataMatSparseKKT3ComputeSinvASinvImpl for more details */
    
    int i, j, k;
    double aVal = 0.0;
    double *Sinv = X;
    double *SinvACol = NULL;
    double *SinvARow = NULL;
    double *SinvRow = NULL;
    double *SinvCol = NULL;
    double *SinvA = aux;
    
    HDSDP_ZERO(SinvA, double, n * n);
    
    for ( i = 0; i < n; ++i ) {
        for ( k = Sp[i]; k < Sp[i + 1]; ++k ) {
            j = Si[k];
            aVal = Sx[k];
            SinvACol = SinvA + n * i;
            SinvRow = Sinv + n * j;
            axpy(&n, &aVal, SinvRow, &HIntConstantOne, SinvACol, &HIntConstantOne);
            
            if ( j != i ) {
                SinvACol = SinvA + n * j;
                SinvRow = Sinv + n * i;
                axpy(&n, &aVal, SinvRow, &HIntConstantOne, SinvACol, &HIntConstantOne);
            }
        }
    }
    
    /* XSX += SinvA^T * Sinv  (one dgemm instead of n*(n+1)/2 strided dot calls) */
    dgemm(&HCharConstantTrans, &HCharConstantNoTrans, &n, &n, &n,
          &HDblConstantOne, SinvA, &n, Sinv, &n, &HDblConstantOne, XSX, &n);

    return;
}

extern double csp_dot_fds( int n, int *Ap, int *Ai, double *Ax, double *B ) {
    
    double dAdotB = 0.0;
    
    for ( int i = 0; i < n; ++i ) {
        for ( int k = Ap[i]; k < Ap[i + 1]; ++k ) {
            int j = Ai[k];
            if ( i == j ) {
                dAdotB += 0.5 * Ax[k] * FULL_ENTRY(B, n, j, i);
            } else {
                dAdotB += Ax[k] * FULL_ENTRY(B, n, j, i);
            }
        }
    }
    
    return 2.0 * dAdotB;
}

extern void csp_max_rowabs( int n, int *Ap, int *Ai, double *Ax, double *row ) {
    
    double x = 0.0;
    for ( int i = 0, j; i < n; ++i ) {
        for ( j = Ap[i]; j < Ap[i + 1]; ++j ) {
            x = fabs(Ax[j]);
            row[Ai[j]] = ( x > row[Ai[j]] ) ? x : row[Ai[j]];
        }
    }
    
    return;
}

extern void csp_min_rownzabs( int m, int n, int *Ap, int *Ai, double *Ax, double *row ) {
    
    for ( int i = 0; i < m; ++i ) {
        row[i] = HDSDP_INFINITY;
    }
    
    double x = 0.0;
    for ( int i = 0, j; i < n; ++i ) {
        for ( j = Ap[i]; j < Ap[i + 1]; ++j ) {
            x = fabs(Ax[j]);
            row[Ai[j]] = ( x < row[Ai[j]] && x > 0 ) ? x : row[Ai[j]];
        }
    }
    
    return;
}

extern void csp_max_colabs( int n, int *Ap, int *Ai, double *Ax, double *col ) {
    
    double cmax = 0.0, x = 0.0;
    for ( int i = 0, j; i < n; ++i ) {
        cmax = 0.0;
        for ( j = Ap[i]; j < Ap[i + 1]; ++j ) {
            x = fabs(Ax[j]);
            cmax = ( x > cmax ) ? x : cmax;
        }
        col[i] = ( cmax > col[i] ) ? cmax : col[i];
    }
    
    return;
}

extern void csp_min_colnzabs( int n, int *Ap, int *Ai, double *Ax, double *col ) {
    
    double cmin = 0.0, x = 0.0;
    
    for ( int i = 0; i < n; ++i ) {
        col[i] = HDSDP_INFINITY;
    }
    
    for ( int i = 0, j; i < n; ++i ) {
        cmin = HDSDP_INFINITY;
        for ( j = Ap[i]; j < Ap[i + 1]; ++j ) {
            x = fabs(Ax[j]);
            cmin = ( x < cmin && x > 0 ) ? x : cmin;
        }
        col[i] = ( cmin < col[i] ) ? cmin : col[i];
    }
    
    return;
}

extern void csp_rowscal( int n, int *Ap, int *Ai, double *Ax, double *row ) {
    
    for ( int i = 0, j; i < n; ++i ) {
        for ( j = Ap[i]; j < Ap[i + 1]; ++j ) {
            Ax[j] = Ax[j] / row[Ai[j]];
        }
    }
    
    return;
}

extern void csp_colscal( int n, int *Ap, int *Ai, double *Ax, double *col ) {
    
    for ( int i = 0, j; i < n; ++i ) {
        for ( j = Ap[i]; j < Ap[i + 1]; ++j ) {
            Ax[j] = Ax[j] / col[i];
        }
    }

    return;
}

extern hdsdp_retcode csp_geoscal( int m, int n, int *Ap, int *Ai, double *Ax, double *D, double *E ) {
    
    /* Scale rows */
    hdsdp_retcode retcode = HDSDP_RETCODE_OK;
    double *dGeoScalRow;
    double *dGeoScalCol;
    
    HDSDP_INIT(dGeoScalRow, double, m);
    HDSDP_INIT(dGeoScalCol, double, n);
    
    HDSDP_MEMCHECK(dGeoScalRow);
    HDSDP_MEMCHECK(dGeoScalCol);
    
    /* Scale row */
    csp_max_rowabs(n, Ap, Ai, Ax, dGeoScalRow);
    csp_min_rownzabs(m, n, Ap, Ai, Ax, D);
    
    for ( int i = 0; i < m; ++i ) {
        D[i] = sqrt(D[i] * dGeoScalRow[i]);
    }
    
    csp_rowscal(n, Ap, Ai, Ax, D);
    
    /* Scale column */
    csp_max_colabs(n, Ap, Ai, Ax, dGeoScalCol);
    csp_min_colnzabs(n, Ap, Ai, Ax, E);
    
    for ( int i = 0; i < n; ++i ) {
        E[i] = sqrt(E[i] * dGeoScalCol[i]);
    }
    
    csp_colscal(n, Ap, Ai, Ax, E);
    
exit_cleanup:
    
    HDSDP_FREE(dGeoScalRow);
    HDSDP_FREE(dGeoScalCol);
    
    return retcode;
}

extern hdsdp_retcode csp_ruizscal( int m, int n, int *Ap, int *Ai, double *Ax, double *D, double *E, int maxIter ) {
    
    hdsdp_retcode retcode = HDSDP_RETCODE_OK;
    
    int nRow = m;
    int nCol = n;
    
    double *dRuizScalDiagRow = D;
    double *dRuizScalDiagCol = E;
    double *dRuizWorkDiagRow = NULL;
    double *dRuizWorkDiagCol = NULL;
    
    /* Allocate workspace */
    HDSDP_INIT(dRuizWorkDiagRow, double, nRow);
    HDSDP_INIT(dRuizWorkDiagCol, double, nCol);
    
    if ( !dRuizWorkDiagRow || !dRuizWorkDiagCol ) {
        retcode = HDSDP_RETCODE_FAILED;
        goto exit_cleanup;
    }
    
    for ( int i = 0; i < maxIter; ++i ) {
        
        HDSDP_ZERO(dRuizWorkDiagRow, double, nRow);
        csp_max_rowabs(nCol, Ap, Ai, Ax, dRuizWorkDiagRow);
        HDSDP_ZERO(dRuizWorkDiagCol, double, nCol);
        csp_max_colabs(nCol, Ap, Ai, Ax, dRuizWorkDiagCol);
        
        /* sqrt operation */
        double maxRuizDiagDeviate = 0.0;
        double ruizDiagDeviate = 0.0;
        
        for ( int j = 0; j < nRow; ++j ) {
            dRuizWorkDiagRow[j] = sqrtl(dRuizWorkDiagRow[j]);
            dRuizScalDiagRow[j] = dRuizScalDiagRow[j] * dRuizWorkDiagRow[j];
            ruizDiagDeviate = fabs(dRuizWorkDiagRow[j] - 1.0);
            maxRuizDiagDeviate = HDSDP_MAX(maxRuizDiagDeviate, ruizDiagDeviate);
        }
        
        for ( int j = 0; j < nCol; ++j ) {
            dRuizWorkDiagCol[j] = sqrtl(dRuizWorkDiagCol[j]);
            dRuizScalDiagCol[j] = dRuizScalDiagCol[j] * dRuizWorkDiagCol[j];
            ruizDiagDeviate = fabs(dRuizWorkDiagCol[j] - 1.0);
            maxRuizDiagDeviate = HDSDP_MAX(maxRuizDiagDeviate, ruizDiagDeviate);
        }
        
        if ( maxRuizDiagDeviate < 1e-08 ) {
            break;
        }
        
        /* Scaling */
        csp_rowscal(nCol, Ap, Ai, Ax, dRuizWorkDiagRow);
        csp_colscal(nCol, Ap, Ai, Ax, dRuizWorkDiagCol);
    }
        
    for ( int i = 0; i < nRow; ++i ) {
        dRuizScalDiagRow[i] = 1.0 / dRuizScalDiagRow[i];
    }
    
    for ( int i = 0; i < nCol; ++i ) {
        dRuizScalDiagCol[i] = 1.0 / dRuizScalDiagCol[i];
    }
        
exit_cleanup:
    HDSDP_FREE(dRuizWorkDiagRow);
    HDSDP_FREE(dRuizWorkDiagCol);
    return retcode;
}

extern void csp_l2scal( int m, int n, int *Ap, int *Ai, double *Ax, double *D, double *E ) {
        
    int nRow = m;
    int nCol = n;
    
    double *dL2WorkDiagRow = D;
    double *dL2WorkDiagCol = E;
    
    for ( int i = 0, j; i < nCol; ++i ) {
        for ( j = Ap[i]; j < Ap[i + 1]; ++j ) {
            dL2WorkDiagRow[Ai[j]] += Ax[j] * Ax[j];
        }
    }
    
    for ( int i = 0; i < nRow; ++i ) {
        dL2WorkDiagRow[i] = sqrt(dL2WorkDiagRow[i]);
        if ( dL2WorkDiagRow[i] == 0.0 ) {
            dL2WorkDiagRow[i] = 1.0;
        }
    }
    
    csp_rowscal(nCol, Ap, Ai, Ax, dL2WorkDiagRow);
    
    for ( int i = 0, j; i < nCol; ++i ) {
        for ( j = Ap[i]; j < Ap[i + 1]; ++j ) {
            dL2WorkDiagCol[i] += Ax[j] * Ax[j];
        }
    }
    
    for ( int i = 0; i < nCol; ++i ) {
        dL2WorkDiagCol[i] = sqrt(dL2WorkDiagCol[i]);
        if ( dL2WorkDiagCol[i] == 0.0 ) {
            dL2WorkDiagCol[i] = 1.0;
        }
    }
    
    csp_colscal(nCol, Ap, Ai, Ax, dL2WorkDiagCol);;
    return;
}

extern void csp_dump( int n, int *Ap, int *Ai, double *Ax, double *v ) {
    
    for ( int i = 0, j; i < n; ++i ) {
        for ( j = Ap[i]; j < Ap[i + 1]; ++j ) {
            v[i * n + Ai[j]] = v[Ai[j] * n + i] = Ax[j];
        }
    }
    
    return;
}

/* Decompress a column */
extern void tsp_decompress( int n, int nnz, int *Ci, double *Cx, int *Ai, int *Aj, double *Ax ) {
    
    int j = 0, idthresh = n;
    
    for ( int k = 0; k < nnz; ++k ) {
        while ( Ci[k] >= idthresh ) {
            j += 1;
            idthresh += n - j;
        }
        Ai[k] = Ci[k] - idthresh + n;
        Aj[k] = j;
        Ax[k] = Cx[k];
    }
    
    return;
}

/** @brief Check if the matrix is rank-one
 *
 * a is an n by one auxiliary vector and on exit, the returned integer values tells
 * if the matrix is rank-one and a would contain the rank-one element
 *
 * In a word A = sgn \* a \* a'
 *
 */
extern int tsp_r1_extract( int n, int nnz, int *Ai, int *Aj, double *Ax, double *sgn, double *a ) {
    
    assert( nnz > 0 );
    
    /* Get the first nonzero */
    int i = Ai[0];
    int j = Aj[0];
    double v = Ax[0];
    
    if ( i != j ) {
        return 0;
    }
    
    if ( nnz == 1 ) {
        *sgn = Ax[0];
        a[i] = 1.0;
        return 1;
    }
    
    double s = ( v > 0 ) ? 1.0 : -1.0;
    v = sqrt(fabs(v));
    
    /* Assume that the matrix is rank-one */
    int k = 0;
    int anz = 0;
    for ( k = 0; k < nnz; ++k ) {
        if ( Aj[k] > i ) {
            break;
        }
        a[Ai[k]] = Ax[k] / v;
        anz += 1;
    }
    
    /* The number of nnzs in the submatrix must match */
    if ( nnz != (int) (anz * ( anz + 1 ) / 2) ) {
        return 0;
    }
    
    if ( k == n ) {
        return 0;
    }
    
    /* Now a contains the rank-one components */
    double eps = 0.0;
    
    if ( s == 1.0 ) {
        /* If condition is broken into two cases */
        for ( k = 0; k < nnz; ++k ) {
            eps += fabs(Ax[k] - a[Ai[k]] * a[Aj[k]]);
        }
    } else {
        for ( k = 0; k < nnz; ++k ) {
            eps += fabs(Ax[k] + a[Ai[k]] * a[Aj[k]]);
        }
    }
    
    if ( eps > 1e-10 ) {
        return 0;
    }
    
    *sgn = s;
    
    return 1;
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

extern double tsp_quadform( int n, int nnz, int *Ai, int *Aj, double *Ax, double *v ) {
    
    double quadform = 0.0, tmp = 0.0;
    
    for ( int k = 0; k < nnz; ++k ) {
        tmp = Ax[k] * v[Ai[k]] * v[Aj[k]];
        if ( Ai[k] == Aj[k] ) {
            quadform += 0.5 * tmp;
        } else {
            quadform += tmp;
        }
    }
    
    return 2.0 * quadform;
}
