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
extern DSDP_INT spsMatFree         ( spsMat *sMat );

/* Basic operations */
extern DSDP_INT spsMatAx           ( spsMat *A, vec *x, vec *Ax );
extern double   spsMatxTAx         ( spsMat *A, double *x       );
extern DSDP_INT spsMataXpbY        ( double alpha,  spsMat *sXMat, double beta,
                                    spsMat *sYMat,  DSDP_INT *sumHash );
extern DSDP_INT spsMatAdddiag      ( spsMat *sMat,  double d, DSDP_INT *sumHash );
extern DSDP_INT spsMatAddds        ( spsMat *sXMat, double alpha, dsMat *sYMat );
extern DSDP_INT spsMatAddr1        ( spsMat *sXMat, double alpha, r1Mat *r1YMat, DSDP_INT *sumHash );
extern DSDP_INT spsMatAddrk        ( spsMat *sXMat, double alpha, rkMat *rkYMat, DSDP_INT *sumHash );
extern DSDP_INT spsMatFnorm        ( spsMat *sMat,  double *fnrm );
extern DSDP_INT spsMatOneNorm      ( spsMat *sMat,  double *onenrm );
extern DSDP_INT spsMatScale        ( spsMat *sXMat, double alpha );
extern DSDP_INT spsMatRscale       ( spsMat *sXMat, double r );

/* Factorization and linear system solver */
extern DSDP_INT spsMatSymbolic     ( spsMat *sAMat );
extern DSDP_INT spsMatFactorize    ( spsMat *sAMat );
extern DSDP_INT spsMatVecSolve     ( spsMat *sAMat, vec    *sbVec, double *Ainvb );
extern DSDP_INT spsMatVecFSolve    ( spsMat *sAmat, vec *sbVec, vec *Ainvb );
extern DSDP_INT spsMatVecBSolve    ( spsMat *sAmat, vec *sbVec, vec *Ainvb );
extern DSDP_INT spsMatLspLSolve    ( spsMat *S,     spsMat *dS,    spsMat *spaux );
extern DSDP_INT spsMatGetX         ( spsMat *S, spsMat *dS, double *LinvSLTinv   );
extern DSDP_INT dsdpGetAlpha       ( spsMat *S, spsMat *dS, spsMat *spaux, double *alpha );
extern DSDP_INT dsdpGetAlphaLS     ( spsMat *S, spsMat *dS, spsMat *Scker,
                                     double alphamax, double *alpha, DSDP_INT *sumHash );

/* Schur matrix assembly */
extern double spsSinvSpSinvSolve   ( const double *Sinv, double *aux, spsMat *A, dsMat *SinvASinv );
extern double spsSinvDsSinvSolve   ( const double *Sinv, double *aux, dsMat  *A, dsMat *SinvASinv );
extern double spsSinvRkSinvSolve   ( spsMat *S, rkMat *A, rkMat *SinvASinv );
extern double spsSinvR1SinvSolve   ( spsMat *S, r1Mat *A, r1Mat *SinvASinv );
extern double spsRySinv            ( spsMat *A, double *Sinv, double *asinv, double Ry );
extern double spsSinvspsSinv       ( spsMat *A, spsMat *B, double *Sinv );
extern double spsSinvr1Sinv        ( spsMat *A, r1Mat  *B, double *Sinv );

/* Eigen value routines */
extern DSDP_INT spsMatMaxEig       ( spsMat *sMat, double *maxEig );
extern DSDP_INT spsMatMinEig       ( spsMat *sMat, double *minEig );

/* Utilities */
extern DSDP_INT spsMatIspd         ( spsMat *sMat, DSDP_INT *ispd );
extern DSDP_INT dsdpGetAlphaLS     ( spsMat *S, spsMat *dS, spsMat *Scker, double alphamax,
                                     double *alpha, DSDP_INT *sumHash );
extern double   spsMatGetlogdet    ( spsMat *sMat, double *aux );
extern DSDP_INT spsMatScatter      ( spsMat *sMat, vec *b, DSDP_INT k );
extern void     spsMatInverse      ( spsMat *sMat, double *Sinv, double *aux );
extern DSDP_INT spsMatStoreFactor  ( spsMat *sMat, rkMat *factor );
extern rkMat*   spsMatGetFactor    ( spsMat *sMat );
extern DSDP_INT spsMatGetRank      ( spsMat *sMat );
extern DSDP_INT spsMatFillLower    ( spsMat *sMat, double *lowFullMat );
extern DSDP_INT spsMatFill         ( spsMat *sMat, double *fulldMat );
extern DSDP_INT spsMatReset        ( spsMat *sMat );
extern DSDP_INT spsMatView         ( spsMat *sMat );
extern void     spsMatLinvView     ( spsMat *S );

#ifdef __cplusplus
}
#endif

#endif /* sparsemat_h */
