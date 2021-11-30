#include "test.h"
DSDP_INT genDenseMatfromVec1( double *x, dsMat *A ) {
    // Generate packed data
    DSDP_INT n = A->dim;
    double *array = A->array;
    for (DSDP_INT i = 0; i < n; ++i) {
        for (DSDP_INT j = 0; j <= i; ++j) {
            packIdx(array, n, i, j) = - x[i] * x[j];
        }
    }
    
    return DSDP_RETCODE_OK;
}

DSDP_INT genDenseMatfromArray( double *data, dsMat *A ) {
    
    double *array = A->array;
    DSDP_INT n = A->dim;
    for (DSDP_INT i = 0; i < nsym(n); ++i) {
        array[i] = data[i];
    }
    
    return DSDP_RETCODE_OK;
}

