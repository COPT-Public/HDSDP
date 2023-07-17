#ifndef vec_opts_h
#define vec_opts_h

#define TRANS   ('T')
#define NOTRANS ('N')
#define UPLOLOW ('L')
#define UPLOUP ('U')
#define SIDELEFT ('L')
#define SIDERIGHT ('R')

/* Constants to be used when calling BLAS. Should never be modified */
static char HCharConstantTrans = TRANS;
static char HCharConstantNoTrans = NOTRANS;
static char HCharConstantUploUp = UPLOUP;
static char HCharConstantUploLow = UPLOLOW;
static int  HIntConstantOne = 1;
static double HDblConstantZero = 0.0;
static double HDblConstantOne = 1.0;
static double HDblConstantMinusOne = -1.0;

#ifdef __cplusplus
extern "C" {
#endif

extern double nrm1( int *n, double *x, int *incx );
extern double nrm2( int *n, double *x, int *incx );
extern void axpy( int *n, double *a, double *x, int *incx, double *y, int *incy );
extern void axpby( int *n, double *a, double *x, int *incx, double *b, double *y, int *incy );
extern double dot( int *n, double *x, int *incx, double *y, int *incy );
extern void scal( int *n, double *sa, double *sx, int *incx );
extern void rscl( int *n, double *sa, double *sx, int *incx );
extern void syr( char *uplo, int *n, double *alpha, double *x, int *incx, double *a, int *lda );
extern void symv( char *uplo, int *n, double *alpha, double *a, int *lda, double *x,
                  int *incx, double *beta, double *y, const int *incy );
extern int idamax( int *n, double *x, int *incx );
extern double sumlogdet( int *n, double *x );
extern void vvscl( int *n, double *s, double *x );
extern void vvrscl( int *n, double *s, double *x );
extern double normalize( int *n, double *a );

#ifdef __cplusplus
}
#endif

#endif /* vec_opts_h */
