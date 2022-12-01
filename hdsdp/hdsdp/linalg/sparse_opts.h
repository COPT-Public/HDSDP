#ifndef sparse_opts_h
#define sparse_opts_h

extern void csp_Axpby( int n, int *Ap, int *Ai, double *Ax, double a, double *x, double *y );
extern void csp_Axpby( int n, int *Ap, int *Ai, double *Ax, double a, double *x, double *y );
extern double csp_sum_abs( int n, int *Ap, int *Ai, double *Ax );
extern double csp_fro_norm( int n, int *Ap, int *Ai, double *Ax );
extern void csp_aApB( int n, int nnz, double a, int *Al, double *Ax, double *Bx );
extern int csp_nnz_cols ( int n, int *Ap );
extern void csp_dump( int n, int *Ap, int *Ai, double *Ax, double *v );

extern void tsp_scal( double a, int nnz, double *Ax );
extern double tsp_sum_abs( int nnz, int *Ai, int *Aj, double *Ax );
extern double tsp_fro_norm( int nnz, int *Ai, int *Aj, double *Ax );
extern void tsp_dump( int n, int nnz, int *Ai, int *Aj, double *Ax, double *v );

#endif /* sparse_opts_h */
