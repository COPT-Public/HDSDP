#ifndef rankonemat_h
#define rankonemat_h

/* Implement the rank one matrix */
#include <stdio.h>
#include "dsdphsd.h"
#include "vec.h"


#define r1MatAlloc  vec_alloc
#define r1MatInit   vec_init
#define r1MatFree   vec_free
#define r1MatFnorm  vec_norm
#define r1MatRscale vec_rscale

#endif /* rankonemat_h */

