#include <stdio.h>
#include "dsdplapack.h"
#include "vec.h"

#define CHECKVEC(x, y)

// Define constants involving Lapack and Blas
static DSDP_INT one = 1;

extern void vec_init( vec *x ) {
    // Initialize vec
    x->x = NULL; x->dim = 0;
}

extern DSDP_INT vec_alloc( vec *x, const DSDP_INT n ) {
    // Allocate memory for vec
    assert( x->dim == 0 );
    x->dim = n; x->x = (double *) calloc(n, sizeof(double));
    return DSDP_RETCODE_OK;
}

extern void vec_copy( vec *src, vec *dst ) {
    // Copy src to dst
    CHECKVEC(src, dst);
    copy(&(src->dim), src->x, &one, dst->x, &one);
}

extern void vec_axpy( double alpha, vec *x, vec *y ) {
    // Compute y = alpha * x + y
    CHECKVEC(x, y);
    axpy(&(x->dim), &alpha, x->x, &one, y->x, &one);
}

extern void vec_axpby( double alpha, vec *x, double beta, vec *y ) {
    // Compute y = alpha * x + beta * y
    CHECKVEC(x, y);
    vecscal(&(x->dim), &beta, y->x, &one);
    axpy(&(x->dim), &alpha, x->x, &one, y->x, &one);
}

extern void vec_zaxpby( vec *z, double alpha, vec *x, double beta, vec *y ) {
    // Compute z = alpha * x + beta * y
    CHECKVEC(z, x);
    CHECKVEC(x, y);
    vec_reset(z);
    if (alpha) {
        axpy(&(x->dim), &alpha, x->x, &one, z->x, &one);
    }
    if (beta) {
        axpy(&(y->dim), &beta,  y->x, &one, z->x, &one);
    }
}

extern void vec_scale( vec *x, double a ) {
    if (a == 1.0) return;
    dscal(&x->dim, &a, x->x, &one);
}

extern void vec_rscale( vec *x, double r ) {
    // Compute x = x / r; No over or under flow.
    if (r == 1.0) return;
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
    double alpha = DSDP_INFINITY, tmp;
    
    if (beta == 0.0) {
        return DSDP_INFINITY;
    }
    
    if (beta > 0.0) {
        for (DSDP_INT i = 0; i < s->dim; ++i) {
            tmp = ds->x[i] / s->x[i];
            alpha = MIN(tmp, alpha);
        }
        alpha = -1.0 / (alpha * beta);
    } else {
        alpha = - DSDP_INFINITY;
        for (DSDP_INT i = 0; i < s->dim; ++i) {
            tmp = ds->x[i] / s->x[i];
            alpha = MAX(tmp, alpha);
        }
        alpha = -1.0 / (alpha * beta);
    }

    if (alpha <= 0.0) {
        alpha = DSDP_INFINITY;
    }
    return alpha;
}

extern DSDP_INT vec_incone( vec *s ) {
    for (DSDP_INT i = 0; i < s->dim; ++i) {
        if (s->x[i] <= 0.0) { return FALSE; }
    }
    return TRUE;
}

extern void vec_set( vec *x, double val) {
    // Set x = val
    for (int i = 0; i < x->dim; ++i) {
        x->x[i] = val;
    }
}

extern void vec_inv( vec *xinv, vec *x ) {
    // Compute xinv
    CHECKVEC(xinv, x);
    for (int i = 0; i < x->dim; ++i) {
        xinv->x[i] = 1.0 / x->x[i];
    }
}

extern void vec_invsqr( vec *xinvsq, vec *x ) {
    // Set xinvsqr = x.^(-2)
    CHECKVEC(xinvsq, x);
    for (int i = 0; i < x->dim; ++i) {
        xinvsq->x[i] = 1.0 / (x->x[i] * x->x[i]);
    }
}

extern void vec_reset( vec *x ) {
    // Set x = 0
    memset(x->x, 0, sizeof(double) * x->dim);
}

extern void vec_norm( vec *x, double *nrm ) {
    // Compute norm(x)
    *nrm = norm(&x->dim, x->x, &one);
}

extern double vec_onenorm( vec *x ) {
    return dasum(&x->dim, x->x, &one);
}

extern double vec_infnorm( vec *x ) {
    double nrm = 0.0;
    for (DSDP_INT i = 0; i < x->dim; ++i) {
        nrm = MAX(fabs(x->x[i]), nrm);
    }
    return nrm;
}

extern void vec_dot( vec *x, vec *y, double *xTy ) {
    // Compute x' * y
    *xTy = dot(&x->dim, x->x, &one, y->x, &one);
}

extern void vec_print( vec *x ) {
    // Print x
    printf("Vector dimension:"ID" \n", x->dim);
    for (int i = 0; i < x->dim; ++i) {
        if ((i + 1) % 11 == 0) {
            printf("\n");
        }
        printf("%-10.3g ", x->x[i]);
    }
    printf("\n\n");
}

extern void vec_free( vec *x ) {
    if (!x) return;
    x->dim = 0; DSDP_FREE(x->x);
}
