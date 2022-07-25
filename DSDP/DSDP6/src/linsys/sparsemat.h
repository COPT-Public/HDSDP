#ifndef sparsemat_h
#define sparsemat_h

/* Implement the sparse matrix data structure for DSDP */
#include "structs.h"

#ifdef __cplusplus
extern "C" {
#endif

/* Structure operations */
extern DSDP_INT spsMatInit     ( spsMat *sMat );
extern DSDP_INT spsMatAlloc    ( spsMat *sMat, DSDP_INT dim );
extern DSDP_INT spsMatAllocData( spsMat *sMat, DSDP_INT dim, DSDP_INT nnz );
extern void spsNominalLinkSinv ( spsMat *sMat, double *Sinv );
extern void spsMatFree         ( spsMat *sMat );

extern void spsMatr1check      ( spsMat *dataMat, DSDP_INT *isRank1 );
extern void spsMatr1extract    ( dsMat *dataMat, double *a, DSDP_INT isNeg );

/* Basic operations */
extern void spsMatAx           ( spsMat *A, vec *x, vec *Ax );
extern void spsMatAx2          ( spsMat *A, vec *x, vec *Ax );
extern double spsMatxTAx       ( spsMat *A, double *x );
extern void spsMataXpbY        ( double alpha,  spsMat *sXMat, double beta, spsMat *sYMat,  DSDP_INT *sumHash );
extern void spsMatAdddiag      ( spsMat *sMat,  double d, DSDP_INT *sumHash );
extern void spsMatAddds        ( spsMat *sXMat, double alpha, dsMat *sYMat );
extern void spsMatAddr1        ( spsMat *sXMat, double alpha, r1Mat *r1YMat, DSDP_INT *sumHash );
extern void spsMatAddrk        ( spsMat *sXMat, double alpha, rkMat *rkYMat, DSDP_INT *sumHash );
extern void spsMatFnorm        ( spsMat *sMat,  double *fnrm );
extern double spsMatOneNorm    ( spsMat *sMat );
extern void spsMatScale        ( spsMat *sXMat, double alpha );
extern void spsMatRscale       ( spsMat *sXMat, double r );

/* Factorization and linear system solver */
extern void spsMatSymbolic     ( spsMat *sAMat );
extern void spsMatFactorize    ( spsMat *sAMat );
extern void spsMatLDLFactorize ( spsMat *sAMat );
extern void spsArrSolveInp     ( spsMat *sAMat, DSDP_INT nrhs, double *B, double *aux );
extern void spsMatVecSolve     ( spsMat *sAMat, vec    *sbVec, double *Ainvb );
extern void spsMatVecFSolve    ( spsMat *sAmat, vec *sbVec, vec *Ainvb );
extern void spsMatVecBSolve    ( spsMat *sAmat, vec *sbVec, vec *Ainvb );
extern void spsMatGetX         ( spsMat *S, spsMat *dS, double *LinvSLTinv   );
extern double spsMatGetAlpha   ( DSDPLanczos *lczSolver, spsMat *S, spsMat *dS, double *alpha );

/* Schur matrix assembly */
extern double SinvSpSinv       ( const double *Sinv, double *aux, spsMat *A, dsMat *SinvASinv );
extern double SinvRkSinv       ( spsMat *S, rkMat *A, rkMat *SinvASinv );
extern double SinvR1Sinv       ( spsMat *S, r1Mat *A, r1Mat *SinvASinv );
extern double spsSinvSolve     ( const double *Sinv, spsMat *A, double *ASinv, double *asinv, double Ry );
extern double spsSinvASinv     ( const double *Sinv, spsMat *A, const double *ASinv );
extern double spsRySinv        ( spsMat *A, double *Sinv, double *asinv, double Ry );
extern double spsSinvspsSinv   ( spsMat *A, spsMat *B, double *Sinv );
extern double spsSinvr1Sinv    ( spsMat *A, r1Mat  *B, double *Sinv );
extern double spsFullTrace     ( spsMat *A, double *S );

/* Utilities */
extern void spsMatIspd         ( spsMat *sMat, DSDP_INT *ispd );
extern void dsdpGetAlphaLS     ( spsMat *S, spsMat *dS, spsMat *Scker, double alphamax, double *alpha, DSDP_INT *sumHash );
extern double spsMatGetlogdet  ( spsMat *sMat, double *aux );
extern void spsMatScatter      ( spsMat *sMat, vec *b, DSDP_INT k );
extern void spsMatInverse      ( spsMat *sMat, double *Sinv, double *aux );
extern void spsMatStoreFactor  ( spsMat *sMat, rkMat *factor );
extern rkMat* spsMatGetFactor  ( spsMat *sMat );
extern DSDP_INT spsMatGetRank  ( spsMat *sMat );
extern void spsMatFillLower    ( spsMat *sMat, double *lowFullMat );
extern void spsMatFillLower2   ( spsMat *sMat, dsMat *lowMat    );
extern void spsMatFill         ( spsMat *sMat, double *fulldMat );
extern void spsMatReset        ( spsMat *sMat );
extern void spsMatGetSymbolic  ( spsMat *sMat, DSDP_INT *hash, DSDP_INT *firstNnz, DSDP_INT *nnzs );
extern DSDP_INT spsMatIsDiag   ( spsMat *sMat );
extern double spsMatGetXbound  ( spsMat *sMat, vec *b );
extern void spsMatView         ( spsMat *sMat );
extern void     spsMatLinvView ( spsMat *S );
extern void     spsMatInvView  ( spsMat *S );
extern void     spsMatExport   ( spsMat *A );
/* Eigen value routines */
#ifdef MKL_FEAST
extern void spsMatMaxEig       ( spsMat *sMat, double *maxEig );
extern void spsMatMinEig       ( spsMat *sMat, double *minEig );
#endif

#ifdef __cplusplus
}
#endif

#endif /* sparsemat_h */
