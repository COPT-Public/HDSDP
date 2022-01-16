#include "test.h"
#include "dsdppresolve.h"

/* Rank 1 data
 
 r1A =
 
 0         0         0         0         0         0         0         0         0         0
 0    9.2108         0   -0.1914         0   -0.6221   -0.3768    4.5211    4.2763    4.3011
 0         0         0         0         0         0         0         0         0         0
 0   -0.1914         0    0.0040         0    0.0129    0.0078   -0.0939   -0.0888   -0.0894
 0         0         0         0         0         0         0         0         0         0
 0   -0.6221         0    0.0129         0    0.0420    0.0254   -0.3053   -0.2888   -0.2905
 0   -0.3768         0    0.0078         0    0.0254    0.0154   -0.1849   -0.1749   -0.1759
 0    4.5211         0   -0.0939         0   -0.3053   -0.1849    2.2192    2.0990    2.1112
 0    4.2763         0   -0.0888         0   -0.2888   -0.1749    2.0990    1.9854    1.9969
 0    4.3011         0   -0.0894         0   -0.2905   -0.1759    2.1112    1.9969    2.0084
 
 
a =
      0
 3.0349
      0
-0.0631
      0
-0.2050
-0.1241
 1.4897
 1.4090
 1.4172
 
*/

DSDP_INT r1adim = 10;
double r1ax[]   = {0, 3.0349, 0, -0.0631, 0, -0.2050, -0.1241, 1.4897, 1.4090, 1.4172};

DSDP_INT spr1Ap[] = {0,  0,  5,  5,  9,  9, 12, 12, 14, 15, 15};
DSDP_INT spr1Ai[] = {1, 3, 5, 7, 8, 3, 5, 7, 8, 5, 7, 8, 7, 8, 8};
double   spr1Ax[] = {-7.12998114e-01, -6.59064673e-01, -5.67035874e-03, -3.26586563e-01,
                     -7.73455822e-01, -6.09210928e-01, -5.24143480e-03, -3.01882518e-01,
                     -7.14949168e-01, -4.50954464e-05, -2.59729014e-03, -6.15116912e-03,
                     -1.49591957e-01, -3.54279028e-01, -8.39039959e-01};



/* Utility functions */
DSDP_INT genDenseMatfromVec1( double *x, dsMat *A );
DSDP_INT genDenseMatfromArray( double *array, dsMat *A );

static DSDP_INT isDenseRank1Acc( dsMat *dataMat, DSDP_INT *isRank1 ) {
    // Detect if a dense matrix is rank one by directly computing the outer product
    // Slower but accurate
    DSDP_INT retcode = DSDP_RETCODE_OK;
    
    double *A    = dataMat->array;
    double *a    = NULL;
    DSDP_INT n   = dataMat->dim;
    DSDP_INT r1  = TRUE;
    DSDP_INT col = 0;
    DSDP_INT isNeg = FALSE;
    
    // Get the first column that contains non-zero elements
    for (DSDP_INT i = 0; i < n; ++i) {
        col = i;
        if (packIdx(A, n, i, i) != 0) {
            break;
        }
    }
    
    assert( col != n - 1 ); // or it is a zero matrix
    a = (double *) calloc(n, sizeof(double));
    
    double adiag = packIdx(A, n, col, col);
    
    if (adiag < 0) {
        isNeg = TRUE;
        adiag = sqrt(- adiag);
    } else {
        adiag = sqrt(adiag);
    }
    
    if (isNeg) {
        for (DSDP_INT i = col; i < n; ++i) {
            a[i] = - packIdx(A, n, i, col) / adiag;
        }
    } else {
        for (DSDP_INT i = col; i < n; ++i) {
            a[i] = packIdx(A, n, i, col) / adiag;
        }
    }
    
    // Check if A = a * a' by computing ||A - a * a'||_F
    double *start = NULL;
    double err    = 0.0;
    double diff   = 0.0;
    DSDP_INT idx  = 0;
    
    if (isNeg) {
        for (DSDP_INT i = 0; i < n; ++i) {
            start = &A[idx];
            for (DSDP_INT j = 0; j < n - i; ++j) {
                diff = start[j] + a[i] * a[i + j];
                err += diff * diff;
            }
            idx += n - i;
            if (err > 1e-08) {
                r1 = FALSE;
                break;
            }
        }
    } else {
        for (DSDP_INT i = 0; i < n; ++i) {
            start = &A[idx];
            for (DSDP_INT j = 0; j < n - i; ++j) {
                diff = start[j] - a[i] * a[i + j];
                err += diff * diff;
            }
            idx += n - i;
            if (err > 1e-08) {
                r1 = FALSE;
                break;
            }
        }
    }
    
    if (r1) {
        *isRank1 = (DSDP_INT) (1 - 2 * isNeg);
    } else {
        *isRank1 = FALSE;
    }
    
    DSDP_FREE(a);
    
    return retcode;
}

static DSDP_INT isSparseRank1( spsMat *dataMat, DSDP_INT *isRank1 ) {
    // Check if a sparse matrix is rank-one
    DSDP_INT retcode = DSDP_RETCODE_OK;
    DSDP_INT isR1 = TRUE;
    
    DSDP_INT *Ap   = dataMat->p;
    DSDP_INT *Ai   = dataMat->i;
    double   *Ax   = dataMat->x;
    DSDP_INT n     = dataMat->dim;
    DSDP_INT col   = 0;
    DSDP_INT isNeg = FALSE;
    double   err   = 0.0;
    double   diff  = 0.0;
    
    // First detect the first column containing nonzeros
    for (DSDP_INT i = 0; i < n; ++i) {
        col = i;
        if (Ap[i + 1] - Ap[i] > 0) {
            break;
        }
    }
    
    assert( col < n - 1 ); // Otherwise the matrix is empty
    
    if (Ai[0] != col) {
        isR1 = FALSE;
        *isRank1 = isR1;
        return retcode;
    }
    
    double *a = NULL;
    double adiag = 0.0;
    a = (double *) calloc(n, sizeof(double));
    
    adiag = Ax[0];
        
    if (adiag < 0) {
        isNeg = TRUE;
        adiag = - sqrt(-adiag);
    } else {
        adiag = sqrt(adiag);
    }
    
    // Get the sparse rank 1 matrix
    for (DSDP_INT j = Ap[col]; j < Ap[col + 1]; ++j) {
        // If the diagonal is zero but other rows contain non-zeros
        a[Ai[j]] = Ax[j] / adiag;
    }
    
    // Ready to check rank-one property
    if (isNeg) {
        for (DSDP_INT i = 0; i < n; ++i) {
            for (DSDP_INT j = Ap[i]; j < Ap[i + 1]; ++j) {
                diff = Ax[j] + a[i] * a[Ai[j]];
                err += diff * diff;
            }
            if (err > 1e-06) {
                isR1 = FALSE;
                break;
            }
        }
    } else {
        
        for (DSDP_INT i = 0; i < n; ++i) {
            for (DSDP_INT j = Ap[i]; j < Ap[i + 1]; ++j) {
                diff = Ax[j] - a[i] * a[Ai[j]];
                err += diff * diff;
            }
            if (err > 1e-06) {
                isR1 = FALSE;
                break;
            }
        }
    }
    
    DSDP_FREE(a);
    
    if (isR1) {
        *isRank1 = (DSDP_INT) (1 - 2 * isNeg);
    } else {
        *isRank1 = FALSE;
    }
    
    return retcode;
}

DSDP_INT test_presolve(void) {
    
    DSDP_INT retcode = DSDP_RETCODE_OK;
    
    dsMat *data = NULL;
    spsMat *spdata = NULL;
    data = (dsMat *) calloc(1, sizeof(dsMat));
    spdata = (spsMat *) calloc(1, sizeof(spsMat));
    
    retcode = denseMatInit(data); checkCodeFree;
    retcode = denseMatAlloc(data, r1adim, FALSE); checkCodeFree;
    retcode = genDenseMatfromVec1(r1ax, data); checkCodeFree;
    // retcode = denseMatView(data);
     
    retcode = spsMatInit(spdata); checkCodeFree;
    retcode = spsMatAllocData(spdata, r1adim, spr1Ap[r1adim]); checkCodeFree;
    
    memcpy(spdata->p, spr1Ap, sizeof(DSDP_INT) * (r1adim + 1));
    memcpy(spdata->i, spr1Ai, sizeof(DSDP_INT) * spr1Ap[r1adim]);
    memcpy(spdata->x, spr1Ax, sizeof(double) * spr1Ap[r1adim]);
    
    DSDP_INT isRank1 = FALSE;
    retcode = isDenseRank1Acc(data, &isRank1); checkCodeFree;
    if (isRank1) {
        passed("Dense rank 1 detection");
    }
    
    isRank1 = FALSE;
    retcode = isSparseRank1(spdata, &isRank1); checkCodeFree;
    
    if (isRank1) {
        passed("Sparse rank 1 detection");
    }
    
clean_up:
    denseMatFree(data);
    DSDP_FREE(data);
    
    spsMatFree(spdata);
    DSDP_FREE(spdata);
    
    return retcode;
}
