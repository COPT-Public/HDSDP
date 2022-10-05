#ifndef vec_mat_h
#define vec_mat_h

#include "pot_solver.h"

static double potDblConstantOne = 1.0;
static double potDblConstantZero = 0.0;
static pot_int potIntConstantOne = 1;
static pot_int potIntConstantZero = 0;

extern double nrm2( pot_int *n, double *x, pot_int *incx );
extern void axpy( pot_int *n, double *a, double *x, pot_int *incx, double *y, pot_int *incy );
extern double dot( pot_int *n, double *x, pot_int *incx, double *y, pot_int *incy );
extern void rscl( pot_int *n, double *sa, double *sx, pot_int *incx );
extern double sumlogdet( pot_int *n, double *x );

#endif /* vec_mat_h */
