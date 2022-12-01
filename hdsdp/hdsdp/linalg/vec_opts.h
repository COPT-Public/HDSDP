#ifndef vec_opts_h
#define vec_opts_h

#define TRANS   ('T')
#define NOTRANS ('N')
#define UPLOLOW ('L')
#define UPLOUP ('U')

#ifdef __cplusplus
extern "C" {
#endif

extern double nrm1( int *n, double *x, int *incx );
extern double nrm2( int *n, double *x, int *incx );
extern void axpy( int *n, double *a, double *x, int *incx, double *y, int *incy );
extern double dot( int *n, double *x, int *incx, double *y, int *incy );
extern void scal( int *n, double *sa, double *sx, int *incx );
extern void rscl( int *n, double *sa, double *sx, int *incx );
extern void syr( char *uplo, int *n, double *alpha, double *x, int *incx, double *a, int *lda );
extern int idamax( int *n, double *x, int *incx );
extern double sumlogdet( int *n, double *x );
extern void vvscl( int *n, double *s, double *x );
extern void vvrscl( int *n, double *s, double *x );

#ifdef __cplusplus
}
#endif

#endif /* vec_opts_h */
