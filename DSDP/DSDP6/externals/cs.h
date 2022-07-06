/* ========================================================================== */
/* CXSparse/Include/cs.h file */
/* ========================================================================== */

/* This is the CXSparse/Include/cs.h file.  It has the same name (cs.h) as
   the CSparse/Include/cs.h file.  The 'make install' for SuiteSparse installs
   CXSparse, and this file, instead of CSparse.  The two packages have the same
   cs.h include filename, because CXSparse is a superset of CSparse.  Any user
   program that uses CSparse can rely on CXSparse instead, with no change to the
   user code.  The #include "cs.h" line will work for both versions, in user
   code, and the function names and user-visible typedefs from CSparse all
   appear in CXSparse.  For experimenting and changing the package itself, I
   recommend using CSparse since it's simpler and easier to modify.  For
   using the package in production codes, I recommend CXSparse since it has
   more features (support for complex matrices, and both int and long
   versions).
 */

/* ========================================================================== */

#ifndef _CXS_H
#define _CXS_H

#include <stdlib.h>
#include <limits.h>
#include <math.h>
#include <stdio.h>
#include "dsdphsd.h"

#ifdef MATLAB_MEX_FILE
#include "mex.h"
#endif

#define NCOMPLEX

#ifdef DSDP64
#ifndef
#define CS_LONG
#endif
#else
#ifdef CS_LONG
#undef CS_LONG
#endif
#endif

#ifdef __cplusplus
extern "C" {
#endif

#define CS_VER 3                    /* CXSparse Version */
#define CS_SUBVER 2
#define CS_SUBSUB 0
#define CS_DATE "Sept 12, 2017"       /* CSparse release date */
#define CS_COPYRIGHT "Copyright (c) Timothy A. Davis, 2006-2016"
#define CXSPARSE

#include "config.h"
#define cs_long_t       SuiteSparse_long
#define cs_long_t_id    SuiteSparse_long_id
#define cs_long_t_max   SuiteSparse_long_max


#ifndef CS_LONG

typedef struct cs_di_sparse  /* matrix in compressed-column or triplet form */
{
    int nzmax ;     /* maximum number of entries */
    int m ;         /* number of rows */
    int n ;         /* number of columns */
    int *p ;        /* column pointers (size n+1) or col indices (size nzmax) */
    int *i ;        /* row indices, size nzmax */
    double *x ;     /* numerical values, size nzmax */
    int nz ;        /* # of entries in triplet matrix, -1 for compressed-col */
} cs_di ;

int cs_di_entry (cs_di *T, int i, int j, double x) ;
cs_di *cs_di_compress (const cs_di *T) ;
double cs_di_norm (const cs_di *A) ;
int cs_di_print (const cs_di *A, int brief) ;

/* utilities */
void *cs_di_calloc (int n, size_t size) ;
void *cs_di_free (void *p) ;
void *cs_di_realloc (void *p, int n, size_t size, int *ok) ;
cs_di *cs_di_spalloc (int m, int n, int nzmax, int values, int t) ;
cs_di *cs_di_spfree (cs_di *A) ;
int cs_di_sprealloc (cs_di *A, int nzmax) ;
void *cs_di_malloc (int n, size_t size) ;

/* utilities */
double cs_di_cumsum (int *p, int *c, int n) ;
cs_di *cs_di_done (cs_di *C, void *w, void *x, int ok) ;
int *cs_di_idone (int *p, cs_di *C, void *w, int ok) ;
#else

/* --- primary CSparse routines and data structures ------------------------- */

typedef struct cs_dl_sparse  /* matrix in compressed-column or triplet form */
{
    cs_long_t nzmax ; /* maximum number of entries */
    cs_long_t m ;     /* number of rows */
    cs_long_t n ;     /* number of columns */
    cs_long_t *p ;    /* column pointers (size n+1) or col indlces (size nzmax) */
    cs_long_t *i ;    /* row indices, size nzmax */
    double *x ;     /* numerical values, size nzmax */
    cs_long_t nz ;    /* # of entries in triplet matrix, -1 for compressed-col */
} cs_dl ;

cs_long_t cs_dl_entry (cs_dl *T, cs_long_t i, cs_long_t j, double x) ;
cs_dl *cs_dl_compress (const cs_dl *T) ;
double cs_dl_norm (const cs_di *A) ;
cs_long_t cs_dl_print (const cs_dl *A, cs_long_t brief) ;

/* utilities */
void *cs_dl_calloc (cs_long_t n, size_t size) ;
void *cs_dl_free (void *p) ;
void *cs_dl_realloc (void *p, cs_long_t n, size_t size, cs_long_t *ok) ;
cs_dl *cs_dl_spalloc (cs_long_t m, cs_long_t n, cs_long_t nzmax, cs_long_t values, cs_long_t t) ;
cs_dl *cs_dl_spfree (cs_dl *A) ;
cs_long_t cs_dl_sprealloc (cs_dl *A, cs_long_t nzmax) ;
void *cs_dl_malloc (cs_long_t n, size_t size) ;

/* utilities */
double cs_dl_cumsum (cs_long_t *p, cs_long_t *c, cs_long_t n) ;
cs_dl *cs_dl_done (cs_dl *C, void *w, void *x, cs_long_t ok) ;
cs_long_t *cs_dl_idone (cs_long_t *p, cs_dl *C, void *w, cs_long_t ok) ;

#endif
/* -------------------------------------------------------------------------- */
/* Macros for constructing each version of CSparse */
/* -------------------------------------------------------------------------- */

#ifdef CS_LONG

#define CS_INT cs_long_t
#define CS_INT_MAX cs_long_t_max
#define CS_ID cs_long_t_id

#ifdef CS_COMPLEX
#define CS_ENTRY cs_complex_t
#define CS_NAME(nm) cs_cl ## nm
#define cs cs_cl
#else
#define CS_ENTRY double
#define CS_NAME(nm) cs_dl ## nm
#define cs cs_dl
#endif

#else

#define CS_INT int
#define CS_INT_MAX INT_MAX
#define CS_ID "%d"
#ifdef CS_COMPLEX
#define CS_ENTRY cs_complex_t
#define CS_NAME(nm) cs_ci ## nm
#define cs cs_ci
#else
#define CS_ENTRY double
#define CS_NAME(nm) cs_di ## nm
#define cs cs_di
#endif

#endif

#ifdef CS_COMPLEX
#define CS_REAL(x) creal(x)
#define CS_IMAG(x) cimag(x)
#define CS_CONJ(x) conj(x)
#define CS_ABS(x) cabs(x)
#else
#define CS_REAL(x) (x)
#define CS_IMAG(x) (0.)
#define CS_CONJ(x) (x)
#define CS_ABS(x) fabs(x)
#endif

#define CS_MAX(a,b) (((a) > (b)) ? (a) : (b))
#define CS_MIN(a,b) (((a) < (b)) ? (a) : (b))
#define CS_FLIP(i) (-(i)-2)
#define CS_UNFLIP(i) (((i) < 0) ? CS_FLIP(i) : (i))
#define CS_MARKED(w,j) (w [j] < 0)
#define CS_MARK(w,j) { w [j] = CS_FLIP (w [j]) ; }
#define CS_CSC(A) (A && (A->nz == -1))
#define CS_TRIPLET(A) (A && (A->nz >= 0))

/* --- primary CSparse routines and data structures ------------------------- */
#define cs_entry CS_NAME (_entry)
#define cs_compress CS_NAME (_compress)
#define cs_norm CS_NAME (_norm)
#define cs_print CS_NAME (_print)

/* utilities */
#define cs_calloc CS_NAME (_calloc)
#define cs_free CS_NAME (_free)
#define cs_realloc CS_NAME (_realloc)
#define cs_spalloc CS_NAME (_spalloc)
#define cs_spfree CS_NAME (_spfree)
#define cs_sprealloc CS_NAME (_sprealloc)
#define cs_malloc CS_NAME (_malloc)

/* --- secondary CSparse routines and data structures ----------------------- */
#define css CS_NAME (s)
#define csn CS_NAME (n)
#define csd CS_NAME (d)
/* utilities */
#define cs_cumsum CS_NAME (_cumsum)
/* utilities */
#define cs_dalloc CS_NAME (_dalloc)
#define cs_done CS_NAME (_done)

#ifdef __cplusplus
}
#endif

#endif
