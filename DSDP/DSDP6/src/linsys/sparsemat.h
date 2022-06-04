#ifndef sparsemat_h
#define sparsemat_h

/* Implement the sparse matrix data structure for DSDP */
#include "cs.h"
#include "dsdphsd.h"
#include "dsdppardiso.h"
#include "dsdplapack.h"
#include "rankonemat.h"
#include "vec.h"
#include "structs.h"
#include "dsdplanczos.h"


/* Apply CXSparse to store sparse data
 
 Note that in DSDP, we use cs_dl CSC matrix format as the data structure
 to hold sparse matrices.
 
 In practice, the sparse data structure is used for both LP and SDP
 
 For LP, SpsMat stores the complete matrix A
 For SDP, SpsMat holds the lower triangular part of the symmetric matrix S, A or C
 so that the representation is also the CSR format of the upper triangular part
 which is acceptable by Pardiso (or potentially Cholmod)
 
*/

#ifdef __cplusplus
extern "C" {
#endif

/* Structure operations */
extern DSDP_INT spsMatInit         ( spsMat *sMat );
extern DSDP_INT spsMatAlloc        ( spsMat *sMat, DSDP_INT dim );
extern DSDP_INT spsMatAllocData    ( spsMat *sMat, DSDP_INT dim, DSDP_INT nnz );
extern void     spsNominalLinkSinv ( spsMat *sMat, double *Sinv );
extern void spsMatFree         ( spsMat *sMat );

/* Basic operations */
extern void spsMatAx           ( spsMat *A, vec *x, vec *Ax );
extern void     spsMatAx2          ( spsMat *A, vec *x, vec *Ax );
extern double   spsMatxTAx         ( spsMat *A, double *x       );
extern void spsMataXpbY        ( double alpha,  spsMat *sXMat, double beta,
                                    spsMat *sYMat,  DSDP_INT *sumHash );
extern void spsMatAdddiag      ( spsMat *sMat,  double d, DSDP_INT *sumHash );
extern void spsMatAddds        ( spsMat *sXMat, double alpha, dsMat *sYMat );
extern void spsMatAddr1        ( spsMat *sXMat, double alpha, r1Mat *r1YMat, DSDP_INT *sumHash );
extern void spsMatAddrk        ( spsMat *sXMat, double alpha, rkMat *rkYMat, DSDP_INT *sumHash );
extern void spsMatFnorm        ( spsMat *sMat,  double *fnrm );
extern double   spsMatOneNorm      ( spsMat *sMat );
extern void spsMatScale        ( spsMat *sXMat, double alpha );
extern void spsMatRscale       ( spsMat *sXMat, double r );

/* Factorization and linear system solver */
extern void spsMatSymbolic     ( spsMat *sAMat );
extern void spsMatFactorize    ( spsMat *sAMat );
extern void spsArrSolveInp     ( spsMat *sAMat, DSDP_INT nrhs, double *B, double *aux );
extern void spsMatVecSolve     ( spsMat *sAMat, vec    *sbVec, double *Ainvb );
extern void spsMatVecFSolve    ( spsMat *sAmat, vec *sbVec, vec *Ainvb );
extern void spsMatVecBSolve    ( spsMat *sAmat, vec *sbVec, vec *Ainvb );
extern void spsMatGetX         ( spsMat *S, spsMat *dS, double *LinvSLTinv   );
extern double   dsdpGetAlpha       ( DSDPLanczos *lczSolver, spsMat *S, spsMat *dS, double *alpha );

/* Schur matrix assembly */
extern double SinvSpSinv   ( const double *Sinv, double *aux, spsMat *A, dsMat *SinvASinv );
extern double SinvRkSinv   ( spsMat *S, rkMat *A, rkMat *SinvASinv );
extern double SinvR1Sinv   ( spsMat *S, r1Mat *A, r1Mat *SinvASinv );

extern double spsSinvSolve ( const double *Sinv, spsMat *A, double *ASinv, double *asinv, double Ry );
extern double spsSinvASinv ( const double *Sinv, spsMat *A, const double *ASinv );

extern double spsRySinv            ( spsMat *A, double *Sinv, double *asinv, double Ry );
extern double spsSinvspsSinv       ( spsMat *A, spsMat *B, double *Sinv );
extern double spsSinvr1Sinv        ( spsMat *A, r1Mat  *B, double *Sinv );
extern double spsFullTrace         ( spsMat *A, double *S               );

/* Eigen value routines */
extern void spsMatMaxEig       ( spsMat *sMat, double *maxEig );
extern void spsMatMinEig       ( spsMat *sMat, double *minEig );

/* Utilities */
extern void spsMatIspd         ( spsMat *sMat, DSDP_INT *ispd );
extern void dsdpGetAlphaLS     ( spsMat *S, spsMat *dS, spsMat *Scker, double alphamax,
                                     double *alpha, DSDP_INT *sumHash );
extern double   spsMatGetlogdet    ( spsMat *sMat, double *aux );
extern void spsMatScatter      ( spsMat *sMat, vec *b, DSDP_INT k );
extern void     spsMatInverse      ( spsMat *sMat, double *Sinv, double *aux );
extern void spsMatStoreFactor  ( spsMat *sMat, rkMat *factor );
extern rkMat*   spsMatGetFactor    ( spsMat *sMat );
extern DSDP_INT spsMatGetRank      ( spsMat *sMat );
extern void spsMatFillLower    ( spsMat *sMat, double *lowFullMat );
extern void spsMatFillLower2   ( spsMat *sMat, dsMat *lowMat    );
extern void spsMatFill         ( spsMat *sMat, double *fulldMat );
extern void spsMatReset        ( spsMat *sMat );
extern void spsMatView         ( spsMat *sMat );
extern void     spsMatLinvView     ( spsMat *S );
extern void     spsMatInvView      ( spsMat *S );
extern void     spsMatExport       ( spsMat *A );

#ifdef __cplusplus
}
#endif

#endif /* sparsemat_h */
