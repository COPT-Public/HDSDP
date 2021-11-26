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
extern DSDP_INT spsMatAllocData    ( spsMat *sMat, DSDP_INT dim, DSDP_INT nna );
extern DSDP_INT spsMatFree         ( spsMat *sMat );

/* Basic operations */
extern DSDP_INT spsMataXpbY        ( double alpha, spsMat *sXMat, double beta, spsMat *sYMat );
extern DSDP_INT spsMatFnorm        ( spsMat *sMat, double *fnrm );

/* Factorization and linear system solver */
extern DSDP_INT spFactorize        ( spsMat *sAMat );
extern DSDP_INT spsVecSolve        ( spsMat *sAMat, vec    *sbVec, double *Ainvb );
extern DSDP_INT spsSpSolve         ( spsMat *sAMat, spsMat *sBMat, double *AinvB );
extern DSDP_INT spsDsSolve         ( spsMat *sAMat, dsMat  *sBMat, double *AinvB );
extern DSDP_INT spsR1Solve         ( spsMat *sAMat, r1Mat  *sBMat, double *AinvB );

/* Schur matrix assembly */
extern DSDP_INT spsSinvSpSinvSolve ( spsMat *S, spsMat *A, dsMat *SinvASinv, double *asinv );
extern DSDP_INT spsSinvDsSinvSolve ( spsMat *S, dsMat  *A, dsMat *SinvASinv, double *asinv );

/* Utilities */
extern DSDP_INT spsMatScatter      ( spsMat *sMat, vec *b, DSDP_INT k );
extern DSDP_INT spsMatFill         ( spsMat *sMat, double *fulldMat );
extern DSDP_INT spsMatView         ( spsMat *sMat );

#ifdef __cplusplus
}
#endif

#endif /* sparsemat_h */
