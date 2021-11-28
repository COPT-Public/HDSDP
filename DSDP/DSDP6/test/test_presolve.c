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

/* Utility functions */
DSDP_INT genDenseMatfromVec1( double *x, dsMat *A );
DSDP_INT genDenseMatfromArray( double *array, dsMat *A );

static DSDP_INT isDenseRank1Acc( dsMat *dataMat, DSDP_INT *isRank1 ) {
    // Detect if a dense matrix is rank one by directly computing the  outer product
    DSDP_INT retcode = DSDP_RETCODE_OK;
    
    double *A    = dataMat->array;
    double *a    = NULL;
    DSDP_INT n   = dataMat->dim;
    DSDP_INT r1  = TRUE;
    DSDP_INT col = 0;
    
    // Get the first column that contains non-zero elements
    for (DSDP_INT i = 0; i < n; ++i) {
        if (packIdx(A, n, i, i) > 0) {
            col = i;
            break;
        }
    }
    
    assert( col != n - 1 ); // or it is a zero matrix
    a = (double *) calloc(n, sizeof(double));
    double adiag = sqrt(packIdx(A, n, col, col));
    
    for (DSDP_INT i = col; i < n; ++i) {
        a[i] = packIdx(A, n, i, col) / adiag;
    }
    
    // Check if A = a * a' by computing ||A - a * a'||_F
    double *start = NULL;
    double err    = 0.0;
    double diff   = 0.0;
    DSDP_INT idx  = 0;
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
    
    *isRank1 = r1;
    DSDP_FREE(a);
    
    return retcode;
}

DSDP_INT test_presolve(void) {
    
    DSDP_INT retcode = DSDP_RETCODE_OK;
    
    dsMat *data = NULL;
    data = (dsMat *) calloc(1, sizeof(dsMat));
    
    retcode = denseMatInit(data); checkCodeFree;
    retcode = denseMatAlloc(data, r1adim, FALSE); checkCodeFree;
    retcode = genDenseMatfromVec1(r1ax, data); checkCodeFree;
    retcode = denseMatView(data);
    
    DSDP_INT isRank1 = FALSE;
    retcode = isDenseRank1Acc(data, &isRank1); checkCodeFree;
    
    
clean_up:
    denseMatFree(data);
    DSDP_FREE(data);
    return retcode;
}
