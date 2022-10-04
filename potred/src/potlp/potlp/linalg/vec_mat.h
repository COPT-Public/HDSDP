#ifndef vec_mat_h
#define vec_mat_h

#include "potlp.h"

static double potDblConstantOne = 1.0;
static double potDblConstantZero = 0.0;
static pot_int potIntConstantOne = 1;
static pot_int potIntConstantZero = 0;

extern double nrm2( pot_int *n, double *x, pot_int *incx );
extern void rscl( pot_int *n, double *sa, double *sx, pot_int *incx );

#endif /* vec_mat_h */
