/** @file speigs.c
 *  @brief The implementation of sparse eigen decomposition routine for HDSDP
 *
 * A set of routines that factorize very sparse matrices that typically arise from semi-definite programming problems.
 * The routines detect the special structures of the matrix and accelerate the factorization procedure.
 *
 * This implementation is less regirous since it is called consecutively many times in HDSDP.
 * We therefore refer to https://github.com/leavesgrp/SPEIGS for a more careful implementation of the SPEIG routine library
 *
 *  @author Wenzhi Gao, Shanghai University of Finance and Economics
 *  @date Aug, 27th, 2022
 *
 */

/* Include headers */
#include <string.h>
#include <math.h>
#include "speigs.h"
#include "dsdplapack.h"

#ifndef ROOT
#define ROOT (7.0710678118654757273731092936941422522068e-01) ///< \f$ \frac{\sqrt{2}}{2} \f$
#endif

/* Define Lapack-related constants */
static char jobz = 'V', range = 'A', uplolow = 'L';
static double abstol = 0.0, eps = 1e-10;

/** @brief Reset the speig struct
 *  @param[in] ef  The SPEIG struct
 *
 *  The function resets the internal memory of SPEIGS for the next factorizatiin
 */
static void speigReset( speigfac *ef ) {
    
    memset(ef->perm,     0, sizeof(DSDP_INT) * ef->nmax);
    memset(ef->pinv,     0, sizeof(DSDP_INT) * ef->nmax);
    memset(ef->dworkmat, 0, sizeof(double)   * ef->nmax * ef->nmax);
    
    return;
}

/** @brief Extract the column nonzero statistics of SPEIG
 *  @param[in] ef  The SPEIG struct
 *  @param[in] n Size of matrix
 *  @param[in] nnz Number of nonzeros
 *  @param[in] Ai Row indices of triplet format
 *  @param[in] Aj Column indices of the triplet format
 *  @return Type of the matrix. 0: Nothing found 1: diagonal 2: two-two
 *  The function resets the internal memory of SPEIGS for the next factorizatiin
 */
static DSDP_INT speigGetColNnzStats( speigfac *ef, DSDP_INT n, DSDP_INT nnz,
                                     DSDP_INT *Ai, DSDP_INT *Aj ) {
    
    memset(ef->cnz, 0, sizeof(DSDP_INT) * n);
    DSDP_INT *c = ef->cnz, i, j, k, isDiag = TRUE, isTwo = TRUE;
    
    for (k = 0; k < nnz; ++k) {
        i = Ai[k]; j = Aj[k]; c[i] += 1;
        if (i != j) {
            c[j] += 1; isDiag = FALSE;
        }
    }
    
    for (i = 0; i < n; ++i) {
        if (c[i] > 1) {
            isTwo = FALSE; break;
        }
    }
    
    if (isDiag)  { return 1; }
    if (isTwo) { return 2; }
    
    return 0;
}

/** @brief Factorize a diagonal or two-two matrix
 *  @param[in] ef  The SPEIG struct
 *  @param[in] n Size of matrix
 *  @param[in] nnz Number of nonzeros
 *  @param[in] stype The special type of the matrix.
 *  @param[in] Ai Row indices of triplet format
 *  @param[in] Aj Column indices of the triplet format
 *  @param[in] Ax Data in the triplet format
 *  @param[out] evals Eigen-values computed
 *  @param[out] evecs Eigen-vectors computed
 *
 *  The function factorizes diagonal matrix or two-two matrix using Givens' rotations
 */
static void speigGetSpecialFactor( speigfac *ef, DSDP_INT n, DSDP_INT nnz,
                                   DSDP_INT stype, DSDP_INT *Ai, DSDP_INT *Aj,
                                   double *Ax, double *evals, double *evecs ) {
    
    DSDP_INT i, j, k, l; double *v = evecs, *e = evals;
    
    if (stype == 1) { // Diagonal
        for (k = 0; k < nnz; ++k) {
            v = &evecs[k * n];
            v[Ai[k]] = 1.0; e[k] = Ax[k];
        }
    } else { // Two-two
        for (k = l = 0; k < nnz; ++k) {
            i = Ai[k]; j = Aj[k];
            if (i == j) {
                v[i] = 1.0;
                e[l] = Ax[k]; v += n; ++l;
            } else {
                v[j] = ROOT; v[i] = ROOT;
                e[l] = Ax[k]; v += n; ++l;
                v[i] = ROOT; v[j] = -ROOT;
                e[l] = -Ax[k]; v += n; ++l;
            }
        }
    }
    
    return;
}

/** @brief Compute the permutation that gives the submatrix
 *  @param[in] ef  The SPEIG struct
 *  @param[in] n Size of matrix
 *  @param[in] nnz Number of nonzeros
 *  @param[in] Ai Row indices of triplet format
 *  @param[in] Aj Column indices of the triplet format
 *  @param[in] Ax Data in the triplet format
 *  @return Size of the submatrix
 *
 *  The function constructs permutation and the submatrix for further factorization
 */
static DSDP_INT speigGetSubMat( speigfac *ef, DSDP_INT n, DSDP_INT nnz,
                              DSDP_INT *Ai, DSDP_INT *Aj, double *Ax ) {
    // Prepare sub matrix, permutation and inverse permutation
    DSDP_INT sn = 0, i, j, k, *perm = ef->perm;
    double *a = ef->dworkmat;
    
    for (i = 0; i < n; ++i) {
        if (ef->cnz[i]) {
            ef->perm[i] = sn; ef->pinv[sn] = i; ++sn;
        }
    }
    
    for (k = 0; k < nnz; ++k) {
        i = Ai[k]; j = Aj[k];
        a[sn * perm[j] + perm[i]] = Ax[k];
    }
    
    return sn;
}

/** @brief Compute dense eigen decomposition
 *  @param[in] ef  The SPEIG struct
 *  @param[in] n Size of matrix
 *  @param[in] a double array that contains the dense matrix to be factorized
 *  @param[out] evals Eigen-values computed
 *  @param[out] evecs Eigen-vectors computed
 *
 *  The function should not be invoked in the current parameter set up.
 */
static void speigdenseFactor( speigfac *ef, DSDP_INT n,
                              double *a, double *evals, double *evecs) {
    
    double *dwork = ef->dwork, *dworkmat = ef->dworkmat;
    memset(evals, 0, sizeof(double) * n);
    memset(evecs, 0, sizeof(double) * n * n);
    memcpy(dworkmat, a, sizeof(double) * n * n);
    
    DSDP_INT *iwork = ef->iwork, *iworkup = ef->iworkup, \
              lwork = ef->lwork, liwork = ef->liwork, m, info;
    
    dsyevr(&jobz, &range, &uplolow, &n, dworkmat,
           &n, NULL, NULL, NULL, NULL, &abstol, &m,
           evals, evecs, &n, iworkup, dwork,
           &lwork, iwork, &liwork, &info);
    
    return;
}

/** @brief Compute sparse eigen-value decomposition using submatrix
 *  @param[in] ef  The SPEIG struct
 *  @param[in] n Size of matrix
 *  @param[in] sn Size of the submatrix
 *  @param[out] evals Eigen-values computed
 *  @param[out] evecs Eigen-vectors computed
 *
 *  The function should not be invoked in the current parameter set up.
 */
static void speigSubMatFactor( speigfac *ef, DSDP_INT n, DSDP_INT sn,
                               double *evals, double *evecs ) {
    // Do eigen decomposition
    double *dwork = ef->dwork, *dwmat = ef->dworkmat,\
                    *dweval = ef->dworkevl, *dwevec = ef->dworkevc;
    double *sev = NULL, *ev = NULL;
    DSDP_INT lwork = ef->lwork, liwork = ef->liwork, m, info;
    DSDP_INT *iwork = ef->iwork, *iworkup = ef->iworkup;
    DSDP_INT *iperm = ef->pinv, i, j, r;
    
    memset(evals, 0, sizeof(double) * n);
    memset(evecs, 0, sizeof(double) * n * n);
    
    dsyevr(&jobz, &range, &uplolow, &sn, dwmat,
           &sn, NULL, NULL, NULL, NULL, &abstol, &m,
           dweval, dwevec, &sn, iworkup, dwork,
           &lwork, iwork, &liwork, &info);
    
    for (i = r = 0; i < sn; ++i) {
        if (fabs(dweval[i]) > eps) {
            evals[r] = dweval[i];
            sev = dwevec + sn * i; ev = evecs + n * r;
            for (j = 0; j < sn; ++j) {
                ev[iperm[j]] = sev[j];
            }
            ++r;
        }
    }
    
    return;
}

/** @brief Initialize SPEIG struct
 *  @param[in] ef  The SPEIG struct
 *
 *  All the pointers are set NULL and lengths are initialized 0
 */
extern void speigInit( speigfac *ef ) {
    
    ef->nmax  = 0; ef->lwork = 0; ef->liwork = 0;
    ef->perm  = NULL; ef->pinv = NULL; ef->cnz = NULL;
    ef->dwork = NULL; ef->dworkmat = NULL;
    ef->dworkevc = NULL; ef->dworkevl = NULL;
    ef->iwork  = NULL; ef->iworkup = NULL;
    
    return;
}

/** @brief Allocate internal memory for SPEIG
 *  @param[in] ef  The SPEIG struct
 *  @param[in] nmax The maximum size of matrices to be processed
 *  @return DSDP_RETCODE_OK if memory is successfully allocated
 * Allocate memory for SPEIG
 */
extern DSDP_INT speigAlloc( speigfac *ef, DSDP_INT nmax ) {
    
    DSDP_INT retcode = DSDP_RETCODE_OK;
    ef->nmax     = nmax;  ///< The largest among the matrix flow
    ef->lwork    = LWORK * nmax; ef->liwork = IWORK * nmax; ///< Lapack lwork. >= nmax * 26. liwork >= nmax * 10
    ef->perm     = (DSDP_INT *) calloc(nmax, sizeof(DSDP_INT)); ///< Submatrix permutation
    ef->pinv     = (DSDP_INT *) calloc(nmax, sizeof(DSDP_INT)); ///< Inverse submatrix permutation
    ef->cnz      = (DSDP_INT *) calloc(nmax, sizeof(DSDP_INT)); ///< Number of nonzeros in each column
    ef->dwork    = (double   *) calloc(ef->lwork, sizeof(double)); ///< Lapack double working space
    ef->dworkmat = (double   *) calloc(nmax * nmax, sizeof(double)); ///< Holder for submatrix
    ef->dworkevl = (double   *) calloc(nmax, sizeof(double)); ///< Holder for eigen-value
    ef->dworkevc = (double   *) calloc(nmax * nmax, sizeof(double)); ///< Holder for eigen-vector
    ef->iwork    = (DSDP_INT *) calloc(ef->liwork, sizeof(DSDP_INT)); ///< Lapack iwork array
    ef->iworkup  = (DSDP_INT *) calloc(2 * nmax, sizeof(DSDP_INT)); ///< Lapack isuppz array
    
    if (!ef->perm || !ef->pinv || !ef->cnz || !ef->dwork ||
        !ef->dworkmat || !ef->dworkevl || !ef->dworkevc ||
        !ef->iwork || !ef->iworkup) {
        retcode = DSDP_RETCODE_FAILED;
    }
    
    return retcode;
}

/** @brief Apply SPEIGS to a sparse matrix
 *  @param[in] ef  The SPEIG struct
 *  @param[in] A  Maximum size of matrices to be processed
 *  @param[out] evals Eigen-values
 *  @param[out] evecs Eigen-vectors
 *  @return TRUE is factorization is finished.
 */
extern DSDP_INT speigSpFactor( speigfac *ef, spsMat *A,
                       double *evals, double *evecs ) {
    
    if (A->nnz < 10 || A->nnz >= A->dim * 3 || A->nnz >= 15000) {
        return FALSE;
    }
    DSDP_INT stype = 0, nsub = 0;
    stype = speigGetColNnzStats(ef, A->dim, A->nnz, A->i, A->cidx);
    
    memset(evals, 0, sizeof(double) * A->dim);
    memset(evecs, 0, sizeof(double) * A->dim * A->dim);
    
    if (stype) {
        speigGetSpecialFactor(ef, A->dim, A->nnz, stype,
                              A->i, A->cidx, A->x, evals, evecs);
    } else {
        speigReset(ef);
        nsub = speigGetSubMat(ef, A->dim, A->nnz, A->i, A->cidx, A->x);
        speigSubMatFactor(ef, A->dim, nsub, evals, evecs);
    }
    
    return TRUE;
}

/** @brief Compute dense eigen decomposition
 *  @param[in] ef  The SPEIG struct
 *  @param[in] A Size of matrix
 *  @param[out] evals Eigen-values computed
 *  @param[out] evecs Eigen-vectors computed
 *
 *  The function should not be invoked in the current parameter set up.
 */
extern void speigDfac( speigfac *ef, dsMat *A, double *evals, double *evecs ) {
    
    speigdenseFactor(ef, A->dim, A->array, evals, evecs);
    return;
}

/** @brief Free the internal memory of SPEIG struct
 *  @param[in] ef  The SPEIG struct
 *
 *  The function frees all the internal memory allocated by SPEIG. The SPEIG pointer itself is not freed
 */
extern void speigFree( speigfac *ef ) {
    
    DSDP_FREE(ef->perm);     DSDP_FREE(ef->pinv);
    DSDP_FREE(ef->cnz);      DSDP_FREE(ef->dwork);
    DSDP_FREE(ef->dworkmat); DSDP_FREE(ef->dworkevl);
    DSDP_FREE(ef->dworkevl); DSDP_FREE(ef->iwork);
    DSDP_FREE(ef->iworkup);
    
    ef->nmax = 0; ef->lwork = 0; ef->liwork = 0;
    return;
}
