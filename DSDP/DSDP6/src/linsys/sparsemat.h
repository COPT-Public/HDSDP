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
extern DSDP_INT spsMataXpbY        ( double alpha, spsMat *sXMat, double beta,
                                    spsMat *sYMat, DSDP_INT *sumHash );
extern DSDP_INT spsMatAdddiag      ( spsMat *sMat, double d, DSDP_INT *sumHash );
extern DSDP_INT spsMatAddds        ( spsMat *sXMat, double alpha, dsMat *sYMat );
extern DSDP_INT spsMatAddr1        ( spsMat *sXMat, double alpha, r1Mat *r1YMat, DSDP_INT *sumHash );
extern DSDP_INT spsMatFnorm        ( spsMat *sMat, double *fnrm );
extern DSDP_INT spsMatScale        ( spsMat *sXMat, double alpha );
extern DSDP_INT spsMatRscale       ( spsMat *sXMat, double r );

/* Factorization and linear system solver */
extern DSDP_INT spsMatSymbolic     ( spsMat *sAMat );
extern DSDP_INT spsMatFactorize    ( spsMat *sAMat );
extern DSDP_INT spsMatVecSolve     ( spsMat *sAMat, vec    *sbVec, double *Ainvb );
extern DSDP_INT spsMatVecFSolve    ( spsMat *sAmat, vec *sbVec, vec *Ainvb );
extern DSDP_INT spsMatVecBSolve    ( spsMat *sAmat, vec *sbVec, vec *Ainvb );
extern DSDP_INT spsMatSpSolve      ( spsMat *sAMat, spsMat *sBMat, double *AinvB );
extern DSDP_INT spsMatDsSolve      ( spsMat *sAMat, dsMat  *sBMat, double *AinvB );
extern DSDP_INT spsMatLspLSolve    ( spsMat *S,     spsMat *dS,    spsMat *spaux );
extern DSDP_INT dsdpGetAlpha       ( spsMat *S, spsMat *dS, spsMat *spaux, double *alpha );
extern DSDP_INT dsdpGetAlphaLS     ( spsMat *S, spsMat *dS, spsMat *Scker,
                                     double alphamax, double *alpha, DSDP_INT *sumHash );

/* Schur matrix assembly */
extern DSDP_INT spsSinvSpSinvSolve ( spsMat *S, spsMat *A, dsMat *SinvASinv, double *asinv );
extern DSDP_INT spsSinvDsSinvSolve ( spsMat *S, dsMat  *A, dsMat *SinvASinv, double *asinv );
extern DSDP_INT spsSinvR1SinvSolve ( spsMat *S, r1Mat  *A, r1Mat *SinvASinv, double *asinv );

/* Eigen value routines */
extern DSDP_INT spsMatMaxEig       ( spsMat *sMat, double *maxEig );
extern DSDP_INT spsMatMinEig       ( spsMat *sMat, double *minEig );

/* Utilities */
extern DSDP_INT spsMatIspd         ( spsMat *sMat, DSDP_INT *ispd );
extern DSDP_INT dsdpGetAlphaLS     ( spsMat *S, spsMat *dS, spsMat *Scker, double alphamax,
                                     double *alpha, DSDP_INT *sumHash );
extern DSDP_INT spsMatGetlogdet    ( spsMat *sMat, double *logdet );
extern DSDP_INT spsMatScatter      ( spsMat *sMat, vec *b, DSDP_INT k );
extern DSDP_INT spsMatFill         ( spsMat *sMat, double *fulldMat );
extern DSDP_INT spsMatReset        ( spsMat *sMat );
extern DSDP_INT spsMatView         ( spsMat *sMat );

#ifdef __cplusplus
}
#endif

#endif /* sparsemat_h */
