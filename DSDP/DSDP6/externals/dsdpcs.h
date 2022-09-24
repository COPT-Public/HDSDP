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

#ifndef _DSDPCS_
#define _DSDPCS_

#include "dsdphsd.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct hdsdp_cs_sparse  /* matrix in compressed-column or triplet form */
{
    DSDP_INT nzmax ;     /* maximum number of entries */
    DSDP_INT m ;         /* number of rows */
    DSDP_INT n ;         /* number of columns */
    DSDP_INT *p ;        /* column pointers (size n+1) or col indices (size nzmax) */
    DSDP_INT *i ;        /* row indices, size nzmax */
    double *x ;          /* numerical values, size nzmax */
    DSDP_INT nz ;        /* # of entries in triplet matrix, -1 for compressed-col */
} dcs ;

DSDP_INT dcs_entry (dcs *T, DSDP_INT i, DSDP_INT j, double x) ;
dcs *dcs_compress (const dcs *T) ;
double dcs_norm (const dcs *A) ;
DSDP_INT dcs_print (const dcs *A, DSDP_INT brief) ;

/* utilities */
void *dcs_calloc (DSDP_INT n, size_t size) ;
void *dcs_free (void *p) ;
void *dcs_realloc (void *p, DSDP_INT n, size_t size, DSDP_INT *ok) ;
dcs *dcs_spalloc (DSDP_INT m, DSDP_INT n, DSDP_INT nzmax, DSDP_INT values, DSDP_INT t) ;
dcs *dcs_spfree (dcs *A) ;
DSDP_INT dcs_sprealloc (dcs *A, DSDP_INT nzmax) ;
void *dcs_malloc (DSDP_INT n, size_t size) ;

/* utilities */
double dcs_cumsum (DSDP_INT *p, DSDP_INT *c, DSDP_INT n) ;
dcs *dcs_done (dcs *C, void *w, void *x, DSDP_INT ok) ;
DSDP_INT *dcs_idone (DSDP_INT *p, dcs *C, void *w, DSDP_INT ok) ;

#define IS_CSC(A) (A && (A->nz == -1))
#define IS_TRIPLET(A) (A && (A->nz >= 0))

#ifdef __cplusplus
}
#endif

#endif
