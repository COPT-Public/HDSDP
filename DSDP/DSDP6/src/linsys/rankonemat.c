#include "rankonemat.h"

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

extern DSDP_INT r1MatCountNnz( r1Mat *x ) {
    // Count the number of nonzero elements and setup nzIdx
    assert ( (x->x) && (x->dim) && (!x->nzIdx));
    
    DSDP_INT nnz = 0;
    for (DSDP_INT i = 0; i < x->dim; ++i) {
        if (x->x[i] != 0) {
            nnz += 1;
        }
    }
    
    x->nzIdx = (DSDP_INT *) calloc(sizeof(DSDP_INT), nnz);
    
    DSDP_INT idx = 0;
    for (DSDP_INT i = 0; i < x->dim; ++i) {
        if (x->x[i] != 0) {
            x->nzIdx[idx] = i;
            idx ++;
        }
        
        if (idx == nnz) {
            break;
        }
    }
    
    return DSDP_RETCODE_OK;
}

extern DSDP_INT r1MatFree( r1Mat *x ) {
    // Free the allocated memory in vec structure
    
    assert( x->dim );
    x->sign = 0;
    x->dim  = 0;
    DSDP_FREE(x->x);
    DSDP_FREE(x->nzIdx);
    x->nnz = 0;
    
    return DSDP_RETCODE_OK;
}

extern DSDP_INT r1MatFnorm( r1Mat *x, double *fnrm ) {
    
    assert( x->dim );
    norm(&x->dim, x->x, &one);
    assert(fnrm > 0);
    
    return DSDP_RETCODE_OK;
}

extern DSDP_INT r1MatRscale( r1Mat *x, double r ) {
    
    assert( (x->dim) && (r != 0.0));
    vecdiv(&x->dim, &r, x->x, &one);
    
    return DSDP_RETCODE_OK;
}
