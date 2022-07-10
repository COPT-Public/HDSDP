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
extern void denseMatInit     ( dsMat *dMat                                              );
extern DSDP_INT denseMatAlloc    ( dsMat *dMat, DSDP_INT dim, DSDP_INT doFactor             );
extern void denseMatFree     ( dsMat *dMat                                              );

extern void dsr1check( dsMat *dataMat, DSDP_INT *isRank1 );
extern void dsr1extract( spsMat *dataMat, double *a, DSDP_INT isNeg );

/* Basic operations */
extern void denseMataXpbY    ( double alpha, dsMat *dXMat, double beta, dsMat *dYMat    );
extern void denseMataAxpby   ( dsMat *dAMat, double alpha, vec *x, double beta, vec *Ax );
extern void denseMatAdddiag  ( dsMat *dAMat, double d                                   );
extern void denseMatAdddiagVec( dsMat *dAMat, vec *d                                    );
extern double   denseMatxTAx     ( dsMat *dAMat, double *aux,  double *x                    );
extern void denseMatFnorm    ( dsMat  *dMat, double *fnrm                               );
extern double   denseMatOneNorm  ( dsMat *dMat                                              );
extern void denseMatScale    ( dsMat *dXMat, double a                                   );
extern void denseMatRscale   ( dsMat *dXMat, double r                                   );

/* Factorization and linear system solver */
extern DSDP_INT denseMatFactorize( dsMat * dAMat                                            );
extern void denseArrSolveInp ( dsMat *S, DSDP_INT nrhs, double *B                       );
extern void denseVecSolve    ( dsMat *dAMat, vec    *dbVec, double *Ainvb               );
extern void denseSpsSolve    ( dsMat *dAMat, spsMat *sBMat, double *AinvB               );

/* Schur matrix assembly */
extern void denseSpsTrace    ( dsMat *dAMat, spsMat *sBMat, double *trace               );
extern void denseDsTrace     ( dsMat *dAMat, dsMat *dBMat,  double *trace               );
extern double   denseSinvASinv   ( const double *sv, dsMat *A, const double *asv        );
extern double   denseDiagTrace   ( dsMat *dAMat, double diag                                );
extern double   denseFullTrace   ( dsMat *dMat, double *S                                   );
extern double   SinvDsSinv       ( const double *Sinv, double *aux, dsMat  *A, dsMat *SinvASinv );
extern double   denseSinvSolve   ( const double *Sinv, dsMat *A, double *ASinv, double *asinv, double Ry );
extern double   denseSinvSolve2  ( double *Sinv, dsMat *A, double *asinv, double Ry );

/* Utilities */
extern void denseMatStoreFactor ( dsMat *dMat, rkMat *factor                            );
extern rkMat*   denseMatGetFactor( dsMat *dMat                                              );
extern void denseMatGetRank  ( dsMat *dMat, DSDP_INT *rank                              );
extern void denseMatFillLow  ( dsMat *dMat, double *fulldMat                            );
extern void denseMatFill     ( dsMat *dMat, double *fulldMat                            );
extern void     denseMatReflex   ( dsMat *dMat );
extern void     denseMatGetdiag  ( dsMat *dMat, vec *diag                                   );
extern void denseMatMinEig   ( dsMat *dMat, double *minEig                              );
extern void denseMatView     ( dsMat *dMat                                              );
extern void denseMatReset    ( dsMat *dMat                                              );
extern void denseMatResetFactor( dsMat *dMat                                            );

#ifdef __cplusplus
}
#endif

#endif /* densemat_h */
