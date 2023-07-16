#ifndef r1_opts_h
#define r1_opts_h

#ifdef __cplusplus
extern "C" {
#endif

extern double dsr1_sum_abs( int n, double sign, double *factor );
extern double dsr1_fro_norm( int n, double sign, double *factor );
extern void dsr1_dump( int n, double sign, double *factor, double *v );
extern double dsr1_quadform( int n, double sign, double *factor, double *v );
extern double spr1_sum_abs( double sign, int nnz, double *factornz );
extern void spr1_dump( int n, double sign, int nnz, int *nzidx, double *factornz );
extern double spr1_quadform( int n, double sign, int nnz, int *nzidx, double *factornzs, double *v );

#ifdef __cplusplus
}
#endif

#endif /* r1_opts_h */
