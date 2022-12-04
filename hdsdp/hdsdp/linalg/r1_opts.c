#include "r1_opts.h"
#include "hdsdp_utils.h"
#include "vec_opts.h"

#include <math.h>

extern double dsr1_sum_abs( int n, double sign, double *factor ) {
    
    double nrm = 0.0;
    int incx = 1;
    
    nrm = nrm1(&n, factor, &incx);
    
    return nrm * nrm * fabs(sign);
}

extern double dsr1_fro_norm( int n, double sign, double *factor ) {
    
    double nrm = 0.0;
    int incx = 1;
    
    nrm = nrm2(&n, factor, &incx);
    
    return nrm * nrm * fabs(sign);
}

extern void dsr1_dump( int n, double sign, double *factor, double *v ) {
    
    char uplolow = UPLOLOW;
    int incx = 1;
    syr(&uplolow, &n, &sign, factor, &incx, v, &n);
    HUtilMatSymmetrize(n, v);
    
    return;
}

extern double spr1_sum_abs( double sign, int nnz, double *factornz ) {
    
    double nrm = 0.0;
    int incx = 1;
    
    nrm = nrm1(&nnz, factornz, &incx);
    
    return nrm * nrm * fabs(sign);
}

extern void spr1_dump( int n, double sign, int nnz, int *nzidx, double *factornz ) {
    
    
    
    return;
}
