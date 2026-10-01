/** @file hdsdp\_mineig.c
 *  @brief Header for computing the minimum perturbation applied to a symmetric matrix to make it positive definitego
 *
 * Given a symmetric matrix A, this routine computes the smallest multiple of identity \lambda so that
 *  A + \lambda \cdot I
 * becomes positive semidefinite. Two strategies are implemented to ensure scalability/robustness
 *
 *  1. Apply Lanczos iteration to estimate the smallest eigenvalue of A
 *  2. (Robust option) Apply bisection to the perturbation \lamba until Cholesky decomposition succeeds
 *
 *  Method 1 requires one or two Cholesky decompositions and multiple matrix-vector products.
 *  Method 2 is robust and requires multiple Cholesky decompositions.
 *
 *  The implementation of the algorithm accepts the input matrix in CSC/CSR format. Also a Cholesky decomposition routine is required.
 *
 */

#ifndef hdsdp_mineig_h
#define hdsdp_mineig_h

#ifdef HEADERPATH
#include "linalg/def_hdsdp_mineig.h"
#include "linalg/hdsdp_linsolver.h"
#include "interface/hdsdp_utils.h"
#else
#include "def_hdsdp_mineig.h"
#include "hdsdp_linsolver.h"
#include "hdsdp_utils.h"
#endif

#ifdef __cplusplus
extern "C" {
#endif

/** @brief Create the data structure and allocate necessary memory */
extern hdsdp_retcode HMinEigCreate( hdsdp_mineig **peig );

/** @brief Allocate the internal working array
 *  @param[in] nCol Size of the symmetric input matrix
 *  @param[in] nKrylovDim Dimension of the Krylov subspace
 */
extern hdsdp_retcode HMinEigInit( hdsdp_mineig *eig, int nCol, int nKrylovDim );

/** @brief Driver routine that (approximately) computes the smallest multiple of identity \lambda \cdot I  that ensures symmatric matrix
 *  A + \lambda \cdot I
 *  can successfully pass Cholesky decomposition
 *
 *  The driver routine will NOT modify or copy the input data.
 *
 *  @param[in] eig Pointer to the data structure
 *  @param[in] colMatBeg CSC matrix beg. The input matrix should have nonzero diagonal
 *  @param[in] colMatIdx CSC matrix idx
 *  @param[in] colMatElem CSC matrix elem
 *  @param[in] chol Cholesky certificate. The routine takes a Cholesky decomposition data structure and the CSC representation a symmetric matrix and outputs whether Cholesky is successful.
 *  @param[in] algo Min eigenvalue solution routine. MIN_EIG_LANCZOS applies Lanczos iteration and MIN_EIG_BISECTION applies bisection
 *  @param[out] dPerturb The approximate minimum perturbation
 *
 */
extern hdsdp_retcode HMinEigComputeMinEig( hdsdp_mineig *eig, int *colMatBeg, int *colMatIdx, double *colMatElem,
                                          void *chol, int (*ispsd) ( void *, int *, int *, double *, int * ), mineig_algo algo, double *dPerturb );

/** @brief Clear internal memory */
extern void HMinEigClear( hdsdp_mineig *eig );

/** @brief Clear the internal memory and destroy the data structure */
extern void HMinEigDestroy( hdsdp_mineig **peig );

#ifdef __cplusplus
}
#endif


#endif /* hdsdp_mineig_h */
