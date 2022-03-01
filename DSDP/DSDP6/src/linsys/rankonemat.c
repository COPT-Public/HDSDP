#include <stdio.h>
#include "rankonemat.h"
#include "sparsemat.h"
#include "dsdplapack.h"

// Define constants involving Lapack and Blas
static DSDP_INT one = 1;

extern DSDP_INT r1MatInit( r1Mat *x ) {
    // Initialize rank one matrix
    x->x     = NULL;
    x->dim   = 0;
    x->sign  = 0;
    x->nnz   = 0;
    x->nzIdx = NULL;
    return DSDP_RETCODE_OK;
}

extern DSDP_INT r1MatAlloc( r1Mat *x, const DSDP_INT n ) {
    // Allocate memory for vec
    assert( x->dim == 0 );
    x->dim = n;
    x->x = (double *) calloc(n, sizeof(double));
    x->sign = 1.0;
    return DSDP_RETCODE_OK;
}

extern DSDP_INT r1MatSetData( r1Mat *x, double eigval, double *array ) {
    // Set rank-one data
    assert(x->dim > 0);
    memcpy(x->x, array, sizeof(double) * x->dim);
    x->sign = eigval;
    r1MatCountNnz(x);
    return DSDP_RETCODE_OK;
}

extern DSDP_INT r1denseSpsUpdate( spsMat *sAMat, double alpha, r1Mat *r1BMat ) {
    // Compute A = A + alpha * B where B is a rank-one matrix.
    // A is a sparse matrix that is known to be dense actually
    // Also this routine DOES NOT update Ap and Ai
    DSDP_INT retcode = DSDP_RETCODE_OK;
    DSDP_INT n = sAMat->dim;
    
    assert( sAMat->nnz > 0 );
    
    if (fabs(alpha) < 1e-15) {
        return retcode;
    }
    
    double sign = r1BMat->sign;
    alpha = alpha * sign;
    
    double *array = sAMat->x;
    if (r1BMat->nnz >= 0.8 * n) {
        char uplo = DSDP_MAT_LOW;
        packr1update(&uplo, &n, &alpha, r1BMat->x, &one, array);
    } else {
        double *r1x = r1BMat->x;
        DSDP_INT *nzIdx = r1BMat->nzIdx;
        for (DSDP_INT i = 0; i < r1BMat->nnz; ++i) {
            for (DSDP_INT j = 0; j <= i; ++j) {
                packIdx(array, n, nzIdx[i], nzIdx[j]) += sign * r1x[nzIdx[i]] * r1x[nzIdx[j]];
            }
        }
    }
    
    return retcode;
}

extern double r1Matr1Trace( r1Mat *x, r1Mat *y ) {
    
    // Compute the inner product of two rank one matrices
    /* trace(A1 * A2) = trace( d1 * a1 * a1' * d2 * a2 * a2' )
                      = d1 * d2 * (a1' * a2)^2
      
     y is assumed dense
       ******************************************
       *    Computationally critical routine    *
       ******************************************
     */
    
    double res = 0.0;
    DSDP_INT n = x->dim;
    if (x->nnz < 0.5 * n) {
        // Otherwise do sparse computation
        double *xdata = x->x, *ydata = y->x;
        DSDP_INT *dataidx = x->nzIdx, i;
        for (i = 0; i < x->nnz; ++i) {
            res += xdata[*dataidx] * ydata[*dataidx]; ++dataidx;
        }
    } else {
        // Use dense dot product if rank-one matrix is dense
        res = dot(&x->dim, x->x, &one, y->x, &one);
    }
    
    return (res * res * x->sign * y->sign);
}

extern DSDP_INT r1MatdenseTrace( r1Mat *x, dsMat *A, double *trace ) {
    // Compute the inner product of a rank-1 matrix and dense A
    /* trace (A1 * A2) = trace( d * a1 * a1' * A2 )
                       = d * a1' * A2 * a1
     
     ******************************************
     *    Computationally critical routine    *
     ******************************************
     
     */
    DSDP_INT retcode = DSDP_RETCODE_OK;
    DSDP_INT n = x->dim, nnz = x->nnz;
    
    // Two ways of implementation
    double res = 0.0, *Adata = A->array, *xdata = x->x;
    
    if (nnz >= 0.7 * n) {
        double cval = 0.0, csum;
        
        for (DSDP_INT i = 0, j; i < n; ++i) {
            cval = xdata[i]; csum = 0.5 * packIdx(Adata, n, i, i) * cval;
            for (j = i + 1; j < n; ++j) {
                csum += packIdx(Adata, n, j, i) * xdata[j];
            }
            res += csum * cval;
        }
        /*
         TODO: try the following implementation that avoids casting
         for (DSDP_INT i = 0, j = 0, k = 0; i < n; ++i) {
             cval = xdata[i]; csum = 0.5 * Adata[j] * cval; k = j + 1;
             for (j = k; j < k + n - i; ++j) {
                 csum += Adata[j] * xdata[j - k];
             }
             res += csum * cval;
         }
        */
    } else {
        DSDP_INT cidx, ridx, *nzidx = x->nzIdx;
        double cval = 0.0, csum;
        // x' * A * x = \sum_{i, j} a_{i,j} * x_i * x_j
        for (DSDP_INT col = 0, row; col < nnz; ++col) {
            cidx = nzidx[col];
            cval = xdata[cidx];
            csum = 0.5 * packIdx(Adata, n, cidx, cidx) * cval;
            // Necessary to unroll loop here ?
            for (row = col + 1; row < nnz; ++row) {
                ridx = nzidx[row];
                csum += packIdx(Adata, n, ridx, cidx) * xdata[ridx];
            }
            res += csum * cval;
        }
    }
    
    *trace = 2.0 * res * x->sign;
    
    return retcode;
}

extern DSDP_INT r1MatspsTrace( r1Mat *x, spsMat *A, double *trace ) {
    // Compute the inner product of a rank-1 matrix and sparse A
    DSDP_INT retcode = DSDP_RETCODE_OK;
    
    DSDP_INT n = x->dim, r1nnz = x->nnz, spsnnz = A->nnz, idx, idx2, i, j;
    DSDP_INT *Ap = A->p, *Ai = A->i;
    double *Atimesx = (double *) calloc(n, sizeof(double));
    double *Ax = A->x, *datax = x->x, res = 0.0, coeff = 0.0;
    
    if (r1nnz <= 0.7 * n) {
        for (i = 0; i < r1nnz; ++i) {
            idx = x->nzIdx[i];
            coeff = datax[idx];
            
            if (coeff == 0.0 ||
                Ap[idx] == Ap[idx + 1]) {
                
                if (Ap[idx] == spsnnz) {
                    break;
                } else {
                    continue;
                }
            }
            
            idx2 = Ai[Ap[idx]];
            if (idx2 == idx) {
                Atimesx[idx] += 0.5 * coeff * Ax[Ap[idx]];
            } else {
                Atimesx[idx2] += coeff * Ax[Ap[idx]];
            }
            for (j = Ap[idx] + 1; j < Ap[idx + 1]; ++j) {
                Atimesx[Ai[j]] += coeff * Ax[j];
            }
        }
    } else {
        
        for (idx = 0; idx < n; ++idx) {
            coeff = datax[idx];
            
            if (coeff == 0.0 ||
                Ap[idx] == Ap[idx + 1]) {
                if (Ap[idx] == spsnnz) {
                    break;
                } else {
                    continue;
                }
            }
            
            idx2 = Ai[Ap[idx]];
            if (idx2 == idx) {
                Atimesx[idx] += 0.5 * coeff * Ax[Ap[idx]];
            } else {
                Atimesx[idx2] += coeff * Ax[Ap[idx]];
            }
            for (j = Ap[idx] + 1; j < Ap[idx + 1]; ++j) {
                Atimesx[Ai[j]] += coeff * Ax[j];
            }
        }
    }
    
    if (r1nnz <= 0.7 * n) {
        for (i = 0; i < r1nnz; ++i) {
            idx = x->nzIdx[i];
            res += Atimesx[idx] * datax[idx];
        }
    } else {
        res = dot(&n, Atimesx, &one, datax, &one);
    }
    
    *trace = 2 * res * x->sign;
    DSDP_FREE(Atimesx);
    
    return retcode;
}

extern double r1SinvspsPhaseA( spsMat *A, r1Mat *B, double *Sinv, double *Ry, double *asinv ) {
    // For convenience A is sparse and B is rank-one. But <B * Sinv> is computed
    double res = 0.0, res2 = 0.0, res3 = 0.0, tmp, bij, *Ax = A->x, *Bx = B->x;
    DSDP_INT i, j, k, p, q, in, jn, n = A->dim;
    DSDP_INT *Ai = A->i, *Aj = A->nzHash, *Bi = B->nzIdx;
    
    for (p = 0; p < B->nnz; ++p) {
        for (q = 0; q < p; ++q) {
            bij = Bx[Bi[p]] * Bx[Bi[q]]; in = Bi[p] * n; jn = Bi[q] * n;
            tmp = 0.0; res3 += bij * Sinv[in + Bi[q]];
            for (k = 0; k < A->nnz; ++k) {
                if (Ai[k] == Aj[k]) {
                    tmp += 0.5 * Ax[k] * Sinv[in + Ai[k]] * Sinv[jn + Aj[k]];
                } else {
                    tmp += Ax[k] * Sinv[in + Ai[k]] * Sinv[jn + Aj[k]];
                }
            }
            res += bij * tmp;
        }
        
        q = p; bij = Bx[Bi[p]] * Bx[Bi[q]];
        in = Bi[p] * n; jn = Bi[q] * n; tmp = 0.0;
        res3 += 0.5 * bij * Sinv[in + Bi[q]];
        
        for (k = 0; k < A->nnz; ++k) {
            if (Ai[k] == Aj[k]) {
                tmp += 0.5 * Ax[k] * Sinv[in + Ai[k]] * Sinv[jn + Aj[k]];
            } else {
                tmp += Ax[k] * Sinv[in + Ai[k]] * Sinv[jn + Aj[k]];
            }
        }
        res += 0.5 * bij * tmp;
    }
    
    // <Ry * S^-1 * A_j, S^-1>
    for (p = 0; p < n; ++p) {
        in = p * n;
        for (i = 0; i < B->nnz; ++i) {
            for (j = 0; j < i; ++j) {
                res2 += Bx[Bi[i]] * Bx[Bi[j]] * Sinv[Bi[j] + in] * Sinv[Bi[i] + in];
            }
            j = Bi[i]; tmp = Bx[j] * Sinv[j + in]; res2 += 0.5 * tmp * tmp;
        }
    }
    
    res2 *= *Ry; *Ry = (2.0 * res2 * B->sign); *asinv = (2.0 * res3 * B->sign);
    return (4.0 * res * B->sign);
}

extern double r1SinvspsPhaseB( r1Mat *A, spsMat *B, double *Sinv ) {
    // May be better to write a new one since A is more sparse than B
    return spsSinvr1SinvPhaseB(B, A, Sinv);
}

extern DSDP_INT r1MatdiagTrace( r1Mat *x, double diag, double *trace ) {
    // Compute trace(a * a' * diag * I) = diag * norm(a)^2
    DSDP_INT retcode = DSDP_RETCODE_OK;
            
    if (diag == 0.0) {
        *trace = 0.0;
        return retcode;
    }

    double res = 0.0;
    retcode = r1MatFnorm(x, &res);

    if (x->sign >= 0) {
        *trace = diag * res;
    } else {
        *trace = - diag * res;
    }
    
    return retcode;
}

extern DSDP_INT r1MatCountNnz( r1Mat *x ) {
    // Count the number of nonzero elements and setup nzIdx
    assert ( (x->x) && (x->dim) );
    
    DSDP_INT nnz = 0, n = x->dim, i;
    
    if (!x->nzIdx) {
        x->nzIdx = (DSDP_INT *) calloc(sizeof(DSDP_INT), n);
    } else {
        memset(x->nzIdx, 0, sizeof(DSDP_INT) * n);
    }

    for (i = 0; i < n - 7; ++i) {
        if (x->x[i]) {x->nzIdx[nnz] = i; nnz += 1;} ++i;
        if (x->x[i]) {x->nzIdx[nnz] = i; nnz += 1;} ++i;
        if (x->x[i]) {x->nzIdx[nnz] = i; nnz += 1;} ++i;
        if (x->x[i]) {x->nzIdx[nnz] = i; nnz += 1;} ++i;
        if (x->x[i]) {x->nzIdx[nnz] = i; nnz += 1;} ++i;
        if (x->x[i]) {x->nzIdx[nnz] = i; nnz += 1;} ++i;
        if (x->x[i]) {x->nzIdx[nnz] = i; nnz += 1;} ++i;
        if (x->x[i]) {x->nzIdx[nnz] = i; nnz += 1;}
    }
    
    if (i < n - 3) {
        if (x->x[i]) {x->nzIdx[nnz] = i; nnz += 1;} ++i;
        if (x->x[i]) {x->nzIdx[nnz] = i; nnz += 1;} ++i;
        if (x->x[i]) {x->nzIdx[nnz] = i; nnz += 1;} ++i;
        if (x->x[i]) {x->nzIdx[nnz] = i; nnz += 1;} ++i;
    }
    
    if (i < n - 1) {
        if (x->x[i]) {x->nzIdx[nnz] = i; nnz += 1;} ++i;
        if (x->x[i]) {x->nzIdx[nnz] = i; nnz += 1;} ++i;
    }
    
    if (i < n) {
        if (x->x[i]) {x->nzIdx[nnz] = i; nnz += 1;}
    }
    
    x->nnz = nnz;
    return DSDP_RETCODE_OK;
}

extern DSDP_INT r1MatFree( r1Mat *x ) {
    // Free the allocated memory in vec structure
    
    if (x) {
        assert( x->dim );
        x->sign = 0;
        x->dim  = 0;
        DSDP_FREE(x->x);
        DSDP_FREE(x->nzIdx);
        x->nnz = 0;
    }
    
    return DSDP_RETCODE_OK;
}

extern DSDP_INT r1MatNormalize( r1Mat *x ) {
    
    assert( x->dim );
    double nrm = dnrm2(&x->dim, x->x, &one);
    drscl(&x->dim, &nrm, x->x, &one);
    x->sign = x->sign * nrm * nrm;
    
    return DSDP_RETCODE_OK;
}

extern DSDP_INT r1MatFnorm( r1Mat *x, double *fnrm ) {

    if (x->nnz < 0.6 * x->dim) {
        double res = 0.0, *xdata = x->x;
        DSDP_INT idx, i, *nzidx = x->nzIdx;
        for (i = 0; i < x->nnz; ++i) {
            idx = nzidx[i];
            res += xdata[idx] * xdata[idx];
        }
        *fnrm = res * fabs(x->sign);
    } else {
        *fnrm = norm(&x->dim, x->x, &one);
        *fnrm = (*fnrm) * (*fnrm) * fabs(x->sign);
    }
    
    return DSDP_RETCODE_OK;
}

extern DSDP_INT r1MatOneNorm( r1Mat *x, double *onenrm ) {
    
    double res = 0.0, *xdata = x->x;
    DSDP_INT i, j, *nzIdx = x->nzIdx;
    
    if (x->nnz < 0.6 * x->dim) {
        DSDP_INT nnz = x->nnz;
        for (i = 0; i < nnz; ++i) {
            for (j = 0; j < i; ++j) {
                res += fabs(xdata[nzIdx[i]] * xdata[nzIdx[j]]);
            }
            res += 0.5 * fabs(xdata[nzIdx[i]] * xdata[nzIdx[i]]);
        }
    } else {
        DSDP_INT n = x->dim;
        for (i = 0; i < n; ++i) {
            for (j = 0; j < i; ++j) {
                res += fabs(xdata[i] * xdata[j]);
            }
            res += 0.5 * fabs(xdata[i] * xdata[i]);
        }
    }
    
    *onenrm = 2 * fabs(x->sign) * res;
    return DSDP_RETCODE_OK;
}

extern DSDP_INT r1MatRscale( r1Mat *x, double r ) {
    
    assert( (x->dim) && (r != 0.0) );
    x->sign = x->sign / r;
    return DSDP_RETCODE_OK;
}

extern DSDP_INT r1MatView( r1Mat *x ) {
    
    printf("Matrix View: \n");
    printf("Matrix dimension "ID" : \n", x->dim);
    for (DSDP_INT i = 0; i < x->dim; ++i) {
        printf("%10.3g, ", x->x[i] * x->sign);
    }
    
    printf("\n");
    return DSDP_RETCODE_OK;
}
