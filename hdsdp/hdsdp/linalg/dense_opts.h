#ifndef dense_opts_h
#define dense_opts_h

#ifdef __cplusplus
extern "C" {
#endif

extern void fds_symv( int n, double alpha, double *A, double *x, double beta, double *y );
extern int fds_syev( int n, double *U, double *d, double *Y, int m,
                     double *work, int *iwork, int lwork, int liwork );
extern void fds_gemv( int m, int n, double *M, double *v, double *y );
extern void fds_symm( char side, char uplo, int m, int n, double alpha, double *a, int lda,
                      double *b, int ldb, double beta, double *c, int ldc );
extern void fds_ger( int m, int n, double alpha, double *x, int incx,
                    double *y, int incy, double *a, int lda);
extern void fds_trimultiply( int n, double *S, double *X, double *aux, double *XSX );
extern double fds_dot_fds( int n, double *A, double *B );
extern void fds_print( int n, double *A );
extern void pds_scal( double a, int n, double *A );
extern double pds_sum_abs( int n, double *A );
extern double pds_fro_norm( int n, double *A );
extern void pds_dump( int n, double *A, double *v );
extern void pds_decompress( int nnz, int *Ci, double *Cx, double *A );
extern int pds_r1_extract( int n, double *A, double *sgn, double *a );
extern double pds_quadform( int n, double *A, double *v, double *aux );
extern void pds_spmv( char uplo, int n, double alpha, double *ap, double *x, int incx,
                     double beta, double *y, int incy );
extern void pds_syr( char uplo, int n, double alpha,
                     double *x, int incx, double *a, int lda );

#ifdef __cplusplus
}
#endif

#endif /* dense_opts_h */
