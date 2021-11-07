#ifndef vec_h
#define vec_h

/* Implement vector operations for DSDP-HSD */

#include <stdio.h>
#include "dsdphsd.h"
#include "dsdplapack.h"

// Define constants involving Lapack and Blas
static DSDP_INT one = 1;

typedef struct {

    DSDP_INT     dim;  // Dimension of the vectorß
    double        *x;  // Array storing the data
    
} vec;

DSDP_INT vec_copy( vec * src, vec * dst ) {
    // Copy src to dst
    assert( (src->dim == dst->dim) && (src && dst) );
    
    
    
    return DSDP_RETCODE_OK;
}

DSDP_INT vec_axpy( double alpha, vec * x, double beta, vec * y, vec * z ) {
    // Compute y = alpha * x + beta * y
    
    if (z) {
        
    }
    
    
    return DSDP_RETCODE_OK;
}






#endif /* vec_h */
