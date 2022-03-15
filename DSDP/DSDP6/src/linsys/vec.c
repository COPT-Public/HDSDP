#include <stdio.h>
#include "dsdplapack.h"
#include "vec.h"

#define CHECKVEC(x, y) assert( ((x)->dim == (y)->dim) && ((x) && (y)) )

// Define constants involving Lapack and Blas
static DSDP_INT one = 1;

extern DSDP_INT vec_init( vec *x ) {
    // Initialize vec
    x->x = NULL;
    x->dim = 0;
    return DSDP_RETCODE_OK;
}

extern DSDP_INT vec_alloc( vec *x, const DSDP_INT n ) {
    // Allocate memory for vec
    assert( x->dim == 0 );
    x->dim = n;
    x->x = (double *) calloc(n, sizeof(double));
    return DSDP_RETCODE_OK;
}

extern DSDP_INT vec_copy( vec *src, vec *dst ) {
    // Copy src to dst
    CHECKVEC(src, dst);
    copy(&(src->dim), src->x, &one, dst->x, &one);
    return DSDP_RETCODE_OK;
}

extern DSDP_INT vec_axpy( double alpha, vec *x, vec *y ) {
    // Compute y = alpha * x + y
    CHECKVEC(x, y);
    axpy(&(x->dim), &alpha, x->x, &one, y->x, &one);
    return DSDP_RETCODE_OK;
}

extern DSDP_INT vec_axpby( double alpha, vec *x, double beta, vec *y ) {
    // Compute y = alpha * x + beta * y
    CHECKVEC(x, y);
    vecscal(&(x->dim), &beta, y->x, &one);
    axpy(&(x->dim), &alpha, x->x, &one, y->x, &one);
    return DSDP_RETCODE_OK;
}

extern DSDP_INT vec_zaxpby( vec *z, double alpha, vec *x, double beta, vec *y ) {
    // Compute z = alpha * x + beta * y
    CHECKVEC(z, x);
    CHECKVEC(x, y);
    vec_reset(z);
    axpy(&(x->dim), &alpha, x->x, &one, z->x, &one);
    axpy(&(y->dim), &beta,  y->x, &one, z->x, &one);
    return DSDP_RETCODE_OK;
}

extern DSDP_INT vec_scale( vec *x, double a ) {
    if (a == 1.0) return DSDP_RETCODE_OK;
    dscal(&x->dim, &a, x->x, &one);
    return DSDP_RETCODE_OK;
}

extern DSDP_INT vec_rscale( vec *x, double r ) {
    // Compute x = x / r; No over or under flow.
    assert( x->dim && r != 0);
    
    if (r == 1.0) return DSDP_RETCODE_OK;
    
    DSDP_INT i;
    
    for (i = 0; i < x->dim - 7; ++i) {
        x->x[i] /= r; ++i;
        x->x[i] /= r; ++i;
        x->x[i] /= r; ++i;
        x->x[i] /= r; ++i;
        x->x[i] /= r; ++i;
        x->x[i] /= r; ++i;
        x->x[i] /= r; ++i;
        x->x[i] /= r;
    }
    
    if (i < x->dim - 3) {
        x->x[i] /= r; ++i;
        x->x[i] /= r; ++i;
        x->x[i] /= r; ++i;
        x->x[i] /= r; ++i;
    }
    
    if (i < x->dim - 1) {
        x->x[i] /= r; ++i;
        x->x[i] /= r; ++i;
    }
    
    if (i < x->dim) {
        x->x[i] /= r;
    }
    
    // vecdiv(&x->dim, &r, x->x, &one);
    return DSDP_RETCODE_OK;
}

extern void vec_vdiv( vec *x, vec *y ) {
    // x = x ./ y
    for (DSDP_INT i = 0; i < x->dim; ++i) {
        assert( y->x[i] ); x->x[i] /= y->x[i];
    }
    return;
}

extern void vec_lslack( vec *x, vec *sl, double lb ) {
    // x = x - lb
    for (DSDP_INT i = 0; i < x->dim; ++i) {
        sl->x[i] = x->x[i] - lb;
    }
}

extern void vec_uslack( vec *x, vec *su, double ub ) {
    // x = ub - x
    for (DSDP_INT i = 0; i < x->dim; ++i) {
        su->x[i] = ub - x->x[i];
    }
}

extern double vec_step( vec *s, vec *ds, double beta ) {
    // Compute the largest alpha such that s + alpha * beta * ds >= 0
    double alpha = DSDP_INFINITY;
    
    if (beta == 0.0) {
        return DSDP_INFINITY;
    }
    
    if (beta > 0.0) {
        for (DSDP_INT i = 0; i < s->dim; ++i) {
            alpha = MIN(ds->x[i] / s->x[i], alpha);
        }
        alpha = -1.0 / (alpha * beta);
    } else {
        alpha = - DSDP_INFINITY;
        for (DSDP_INT i = 0; i < s->dim; ++i) {
            alpha = MAX(ds->x[i] / s->x[i], alpha);
        }
        alpha = -1.0 / (alpha * beta);
    }

    if (alpha <= 0.0) {
        alpha = DSDP_INFINITY;
    }
    return alpha;
}

extern DSDP_INT vec_set( vec *x, double val) {
    // Set x = val
    for (int i = 0; i < x->dim; ++i) {
        x->x[i] = val;
    }
    return DSDP_RETCODE_OK;
}

extern DSDP_INT vec_inv( vec *xinv, vec *x ) {
    // Compute xinv
    CHECKVEC(xinv, x);
    for (int i = 0; i < x->dim; ++i) {
        xinv->x[i] = 1.0 / x->x[i];
    }
    return DSDP_RETCODE_OK;
}

extern DSDP_INT vec_invsqr( vec *xinvsq, vec *x ) {
    // Set xinvsqr = x.^(-2)
    CHECKVEC(xinvsq, x);
    for (int i = 0; i < x->dim; ++i) {
        xinvsq->x[i] = 1.0 / (x->x[i] * x->x[i]);
    }
    return DSDP_RETCODE_OK;
}

extern DSDP_INT vec_reset( vec *x ) {
    // Set x = 0
    memset(x->x, 0, sizeof(double) * x->dim);
    return DSDP_RETCODE_OK;
}

extern DSDP_INT vec_norm( vec *x, double *nrm ) {
    // Compute norm(x)
    *nrm = norm(&x->dim, x->x, &one);
    return DSDP_RETCODE_OK;
}

extern DSDP_INT vec_dot( vec *x, vec *y, double *xTy ) {
    // Compute x' * y
    assert( x->dim == y->dim );
    *xTy = dot(&x->dim, x->x, &one, y->x, &one);
    return DSDP_RETCODE_OK;
}

extern DSDP_INT vec_print( vec *x ) {
    // Print x
    printf("Vector dimension:"ID" \n", x->dim);
    for (int i = 0; i < x->dim; ++i) {
        if ((i + 1) % 11 == 0) {
            printf("\n");
        }
        printf("%-10.3g ", x->x[i]);
    }
    printf("\n\n");
    return DSDP_RETCODE_OK;
}

extern DSDP_INT vec_free( vec *x ) {
    
    if (!x) {
        return DSDP_RETCODE_OK;
    }
    
    // Free the allocated memory in vec structure
    x->dim = 0;
    DSDP_FREE(x->x);
    return DSDP_RETCODE_OK;
}
