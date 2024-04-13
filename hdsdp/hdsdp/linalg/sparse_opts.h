#ifndef sparse_opts_h
#define sparse_opts_h

#ifdef __cplusplus
extern "C" {
#endif

extern void csp_Axpby( int n, int *Ap, int *Ai, double *Ax, double a, double *x, double *y );
extern void csp_ATxpby( int n, int *Ap, int *Ai, double *Ax, double a, double *x, double *y );
extern double csp_sum_abs( int n, int *Ap, int *Ai, double *Ax );
extern double csp_fro_norm( int n, int *Ap, int *Ai, double *Ax );
extern void csp_aApB( int n, int nnz, double a, int *Al, double *Ax, double *Bx );
extern int csp_nnz_cols( int n, int *Ap );
extern void csp_trimultiply( int n, int *Sp, int *Si, double *Sx, double *X, double *aux, double *XSX );
extern double csp_dot_fds( int n, int *Ap, int *Ai, double *Ax, double *B );
extern void csp_dump( int n, int *Ap, int *Ai, double *Ax, double *v );

extern void tsp_decompress( int n, int nnz, int *Ci, double *Cx, int *Ai, int *Aj, double *Ax );
extern int tsp_r1_extract( int n, int nnz, int *Ai, int *Aj, double *Ax, double *sgn, double *a );
extern void tsp_scal( double a, int nnz, double *Ax );
extern double tsp_sum_abs( int nnz, int *Ai, int *Aj, double *Ax );
extern double tsp_fro_norm( int nnz, int *Ai, int *Aj, double *Ax );
extern void tsp_dump( int n, int nnz, int *Ai, int *Aj, double *Ax, double *v );
extern double tsp_quadform( int n, int nnz, int *Ai, int *Aj, double *Ax, double *v );

#ifdef __cplusplus
}
#endif

#endif /* sparse_opts_h */
