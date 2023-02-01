#ifndef dense_opts_h
#define dense_opts_h

#ifdef __cplusplus
extern "C" {
#endif

extern void pds_scal( double a, int n, double *A );
extern double pds_sum_abs( int n, double *A );
extern double pds_fro_norm( int n, double *A );
extern void pds_dump( int n, double *A, double *v );
extern void pds_decompress( int nnz, int *Ci, double *Cx, double *A );

#ifdef __cplusplus
}
#endif

#endif /* dense_opts_h */
