#ifndef vec_mat_h
#define vec_mat_h

#include "pot_def.h"

static double potDblConstantOne = 1.0;
static double potDblConstantMinusOne = -1.0;
static double potDblConstantZero = 0.0;
static pot_int potIntConstantOne = 1;
static pot_int potIntConstantZero = 0;

extern double nrm2( pot_int *n, double *x, pot_int *incx );
extern void axpy( pot_int *n, double *a, double *x, pot_int *incx, double *y, pot_int *incy );
extern double dot( pot_int *n, double *x, pot_int *incx, double *y, pot_int *incy );
extern void scal( pot_int *n, double *sa, double *sx, pot_int *incx );
extern void rscl( pot_int *n, double *sa, double *sx, pot_int *incx );
extern pot_int idamax( pot_int *n, double *x, pot_int *incx );

extern pot_int psyev( pot_int n, double *U, double *d, double *Y,
                     double *work, pot_int *iwork, pot_int lwork, pot_int liwork );
extern void pgemv( pot_int m, pot_int n, double *M, double *v, double *y );

extern double sumlogdet( pot_int *n, double *x );
extern void vvscl( pot_int *n, double *s, double *x );
extern void vvrscl( pot_int *n, double *s, double *x );

extern void spMatAxpy( int n, int *Ap, int *Ai, double *Ax, double a, double *x, double *y );
extern void spMatATxpy( int n, int *Ap, int *Ai, double *Ax, double a, double *x, double *y );
extern void spMatMaxRowAbs( int n, int *Ap, int *Ai, double *Ax, double *row );
extern void spMatMaxColAbs( int n, int *Ap, int *Ai, double *Ax, double *col );
extern void spMatRowScal( int n, int *Ap, int *Ai, double *Ax, double *row );
extern void spMatColScal( int n, int *Ap, int *Ai, double *Ax, double *col );
extern int spMatBuildQMat( int qm, int qn, int *Qp, int *Qi, double *Qx,
                           int am, int an, int *Ap, int *Ai, double *Ax,
                           double *b, double *c );
extern int spMatRuizScal( int m, int n, int *Ap, int *Ai, double *Ax, double *D, double *E, int maxIter );
extern int spMatL2Scal( int m, int n, int *Ap, int *Ai, double *Ax, double *D, double *E );

#endif /* vec_mat_h */
