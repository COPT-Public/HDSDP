#include <math.h>
#include "vec_mat.h"

extern double nrm2( pot_int *n, double *x, pot_int *incx ) {
    
    assert( *incx == 1 );
    
    double nrm = 0.0;
    
    for ( int i = 0; i < *n; ++i ) {
        nrm += x[i] * x[i];
    }
    
    return sqrt(nrm);
}

extern void axpy( pot_int *n, double *a, double *x, pot_int *incx, double *y, pot_int *incy ) {
    
    assert( *incx == 1 && *incy == 1 );
    
    for ( int i = 0; i < *n; ++i ) {
        y[i] += (*a) * x[i];
    }
    
    return;
}

extern double dot( pot_int *n, double *x, pot_int *incx, double *y, pot_int *incy ) {
    
    assert( *incx == 1 && *incy == 1 );
    
    double dres = 0.0;
    
    for ( int i = 0; i < *n; ++i ) {
        dres += x[i] * y[i];
    }
    
    return dres;
}

extern void rscl( pot_int *n, double *sa, double *sx, pot_int *incx ) {
    
    assert( *incx == 1 );
    double a = *sa;
    
    assert( a != 0.0 );
    assert( a > 0.0 );
    
    if ( fabs(a) < 1e-16 ) {
        a = (a > 0) ? 1e-16 : -1e-16;
    }
    
    for ( int i = 0; i < *n; ++i ) {
        sx[i] = sx[i] / a;
    }
    
    return;
}

extern double sumlogdet( pot_int *n, double *x ) {
    
    double logdet = 0.0;
    for ( int i = 0; i < *n; ++i ) {
        logdet += logl(x[i]);
    }
    
    return logdet;
}
