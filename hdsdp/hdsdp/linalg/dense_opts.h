#ifndef dense_opts_h
#define dense_opts_h

extern void pds_scal( double a, int n, double *A );
extern double pds_sum_abs( int n, double *A );
extern double pds_fro_norm( int n, double *A );
extern void pds_dump( int n, double *A, double *v );

#endif /* dense_opts_h */
