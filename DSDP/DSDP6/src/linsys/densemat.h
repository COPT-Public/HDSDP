#ifndef densemat_h
#define densemat_h

/* Implement the dense matrix data structure for DSDP */
#include <stdio.h>
#include "dsdphsd.h"
#include "dsdplapack.h"
#include "vec.h"
#include "structs.h"

/*
   In DSDP dense data structure represents the lower triangular part of a symmetric matrix.
   (This is also known as the packed format)
*/

#ifdef __cplusplus
extern "C" {
#endif

/* Structure operations */
extern DSDP_INT denseMatInit     ( dsMat *dMat );
extern DSDP_INT denseMatAlloc    ( dsMat *dMat, DSDP_INT dim, DSDP_INT doFactor );
extern DSDP_INT denseMatFree     ( dsMat *dMat );

/* Basic operations */
extern DSDP_INT denseMataXpbY    ( double alpha, dsMat *dXMat, double beta, dsMat *dYMat );
extern DSDP_INT denseMatFnorm    ( dsMat  *dMat, double *fnrm );
extern DSDP_INT denseMatRscale   ( dsMat *dXMat, double r );

/* Factorization and linear system solver */
extern DSDP_INT denseMatFactorize   ( dsMat * dAMat );
extern DSDP_INT denseVecSolve    ( dsMat *dAMat, vec    *dbVec, double *Ainvb );
extern DSDP_INT denseSpsSolve    ( dsMat *dAMat, spsMat *sBMat, double *AinvB );
extern DSDP_INT denseDsSolve     ( dsMat *dAMat, dsMat  *dBMat, double *AinvB );

/* Schur matrix assembly */
extern DSDP_INT denseSpsTrace    ( dsMat *dAMat, spsMat *sBMat, double *trace );

/* Utilities */
extern DSDP_INT denseMatScatter  ( dsMat *dMat, vec *b, DSDP_INT k );
extern DSDP_INT denseMatFill     ( dsMat *dMat, double *fulldMat );
extern DSDP_INT denseMatGetdiag  ( dsMat *dMat, vec *diag );
extern DSDP_INT denseMatView     ( dsMat *dMat );

#ifdef __cplusplus
}
#endif

#endif /* densemat_h */
