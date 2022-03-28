#ifndef vec_h
#define vec_h

/* Implement vector operations for DSDP-HSD */

#include "dsdphsd.h"
#include "structs.h"

#ifdef __cplusplus
extern "C" {
#endif

extern DSDP_INT vec_init  ( vec *x );
extern DSDP_INT vec_alloc ( vec *x, const DSDP_INT n );
extern DSDP_INT vec_copy  ( vec *src, vec *dst );
extern DSDP_INT vec_axpy  ( double alpha, vec *x, vec *y );
extern DSDP_INT vec_axpby ( double alpha, vec *x, double beta, vec *y );
extern DSDP_INT vec_zaxpby( vec *z, double alpha, vec *x, double beta, vec *y );
extern DSDP_INT vec_scale ( vec *x, double a );
extern DSDP_INT vec_rscale( vec *x, double r );
extern void     vec_vdiv  ( vec *x, vec *y );
extern void     vec_lslack( vec *x, vec *sl, double lb );
extern void     vec_uslack( vec *x, vec *su, double ub );
extern double   vec_step  ( vec *s, vec *ds, double beta );
extern DSDP_INT vec_incone( vec *s );
extern DSDP_INT vec_set   ( vec *x, double val);
extern DSDP_INT vec_inv   ( vec *xinv, vec *x );
extern DSDP_INT vec_invsqr( vec *xinvsq, vec *x );
extern DSDP_INT vec_norm  ( vec *x, double *nrm );
extern DSDP_INT vec_dot   ( vec *x, vec *y, double *xTy );
extern DSDP_INT vec_reset ( vec *x );
extern DSDP_INT vec_print ( vec *x );
extern DSDP_INT vec_free  ( vec *x );

#ifdef __cplusplus
}
#endif

#endif /* vec_h */


