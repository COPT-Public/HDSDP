#ifndef rank1opts_h
#define rank1opts_h

/* Implement the rank one matrix */
#include "structs.h"

#ifdef __cplusplus
extern "C" {
#endif

extern DSDP_INT r1MatInit       ( r1Mat *x );
extern DSDP_INT r1MatAlloc      ( r1Mat *x, const DSDP_INT n );
extern DSDP_INT r1MatSetData    ( r1Mat *x, double eigval, double *array );

extern double r1Matr1Trace      ( r1Mat *x, r1Mat  *y );
extern void r1MatdenseTrace     ( r1Mat *x, dsMat  *A, double *trace );
extern double r1MatSinvSolve    ( const double *Sinv, r1Mat *x, double *ASinv, double *aux, double *asinv, double Ry );
extern double r1MatSinvASinv    ( const double *Sinv, r1Mat *x, const double *ASinv );
extern double r1Sinvr1          ( r1Mat *A, r1Mat  *B, double *Sinv );
extern double r1RySinv          ( r1Mat *B, double *Sinv, double *asinv, double Ry, double *aux );
extern double r1Sinvsps         ( spsMat *A, r1Mat *B, double *Sinv );

extern DSDP_INT r1MatdiagTrace  ( r1Mat *x, double diag, double *trace );
extern double   r1MatFullTrace  ( r1Mat *x, double *S, double *aux );
extern DSDP_INT r1MatCountNnz   ( r1Mat *x );
extern DSDP_INT r1MatFree       ( r1Mat *x );
extern DSDP_INT r1MatNormalize  ( r1Mat *x );
extern DSDP_INT r1MatFnorm      ( r1Mat *x, double *fnrm );
extern double   r1MatOneNorm    ( r1Mat *x );
extern DSDP_INT r1MatScale      ( r1Mat *x, double a );
extern DSDP_INT r1MatRscale     ( r1Mat *x, double r );
extern void r1MatGetNnzs        ( r1Mat *x );
extern void r1MatCheckSparsity  ( r1Mat *x, DSDP_INT *isdense, double thresh );
extern void r1MatGetSymbolic    ( r1Mat *x, DSDP_INT *hash, DSDP_INT *firstNnz, DSDP_INT *nnzs );
extern DSDP_INT r1MatIsConstant ( r1Mat *x );
extern DSDP_INT r1MatView       ( r1Mat *x );

#ifdef __cplusplus
}
#endif

#endif /* rank1opts_h */

