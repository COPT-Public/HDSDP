#include <string.h>
#include "speigs.h"
#include "dsdpfeast.h"
#include "densemat.h"
#include "sparsemat.h"

static char etype[] = "Eigen Interface";
static char jobz = 'V';
static char range = 'A';
static char uplolow = 'L';
static double abstol = 0.0;
static double eps = 1e-10;

/*
 Implement a highly-efficient eigen-decomposition interface using LAPACK or Intel FEAST
 as the backend. The routine is derived from DSDP5.8 by Steve J. Benson.
 */

static void speigReset( speigfac *eigfac ) {
    memset(eigfac->colnnz,   0, sizeof(DSDP_INT) * eigfac->nmax);
    memset(eigfac->perm,     0, sizeof(DSDP_INT) * eigfac->nmax);
    memset(eigfac->pinv,     0, sizeof(DSDP_INT) * eigfac->nmax);
    memset(eigfac->dworkmat, 0, sizeof(double)   * eigfac->nmax * eigfac->nmax);
}

static DSDP_INT speigGetCStats( speigfac *eigfac, DSDP_INT n,
                                DSDP_INT nnz, DSDP_INT *Ai, DSDP_INT *Aj ) {
    // Count the number of nonzeros for each column
    DSDP_INT *colnnz = eigfac->colnnz, i, j, k;
    DSDP_INT diag = TRUE, elem1 = TRUE;
    
    for (k = 0; k < nnz; ++k) {
        i = Ai[k]; j = Aj[k]; colnnz[i] += 1;
        if (i != j) {
            colnnz[j] += 1; diag = FALSE;
        }
    }
    
    for (i = 0; i < n; ++i) {
        if (colnnz[i] > 1) {
            elem1 = FALSE;
        }
    }
    
    if (diag)  { return 1; }
    if (elem1) { return 2; }
    
    return 0;
}

static DSDP_INT denser1check( speigfac *eigfac, DSDP_INT n, double *A ) {
    // Detect if a dense matrix is rank one by directly computing the outer product
    // Slower but accurate
    
    double *a = eigfac->dworkevl; memset(a, 0, sizeof(double) * n);
    DSDP_INT i, j, r1  = TRUE, col = 0, isNeg = FALSE;
    
    // Get the first column that contains non-zero elements
    for (i = 0, j = 0; i < n; ++i) {
        if (A[j] != 0) { break; }
        j += n - i;
    }
    
    if (i == n - 1 && !A[j]) {
        return FALSE;
    }
    
    col = i; double adiag = A[j];
    
    if (adiag < 0) {
        isNeg = TRUE; adiag = sqrt(- adiag);
    } else {
        adiag = sqrt(adiag);
    }
    
    for (i = col; i < n; ++i) {
        a[j] = A[j] / adiag;
        j += n - i;
    }
    
    // Check if A = a * a' by computing ||A - a * a'||_F
    double *start = NULL, err = 0.0, diff   = 0.0;
    DSDP_INT idx  = 0;
    
    if (isNeg) {
        for (i = 0; i < n; ++i) {
            start = &A[idx];
            for (j = 0; j < n - i; ++j) {
                diff = start[j] + a[i] * a[i + j]; err += diff * diff;
            }
            idx += n - i;
            if (err > eps) { r1 = FALSE; break; }
        }
    } else {
        for (i = 0; i < n; ++i) {
            start = &A[idx];
            for (j = 0; j < n - i; ++j) {
                diff = start[j] - a[i] * a[i + j]; err += diff * diff;
            }
            idx += n - i;
            if (err > eps) { r1 = FALSE; break; }
        }
    }
    
    if (r1) {
        return (1 - 2 * isNeg);
    }
    
    return FALSE;
}

#define ROOT 7.0710678118654757273731092936941422522068e-01
static void speigGetSpecialFactor( speigfac *eigfac, DSDP_INT nnz, DSDP_INT n,
                                   DSDP_INT special, DSDP_INT *Ai, DSDP_INT *Aj,
                                   double *Ax, double *eigvals, double *eigvecs ) {
    
    DSDP_INT k; double *eigvec;
    // Check if special factorization is available
    if (special == 1) { // Diagonal
        for (k = 0; k < nnz; ++k) {
            eigvec = &eigvecs[k * n];
            eigvec[Ai[k]] = 1.0; eigvals[k] = Ax[k];
        }
    } else {
        DSDP_INT i, j, idx; double tmp = ROOT;
        for (k = idx = 0; k < nnz; ++k) {
            i = Ai[k]; j = Aj[k];
            eigvec = eigvecs + idx * n;
            if (i == j) {
                eigvals[idx] = Ax[k]; idx += 1;
                eigvec[i] = 1.0; eigvec += n;
            } else {
                eigvec[j] = eigvec[i] = tmp;
                eigvals[idx] = Ax[k]; idx += 1;
                eigvec += n;
                eigvec[i] = tmp; eigvec[j] = -tmp;
                eigvals[idx] = - Ax[k]; idx += 1;
            }
        }
    }
}

static DSDP_INT speigGetSMat( speigfac *eigfac, DSDP_INT n, DSDP_INT nnz,
                              DSDP_INT *Ai, DSDP_INT *Aj, double *Ax ) {
    // Prepare sub matrix, permutation and inverse permutation
    DSDP_INT nsub = 0, i, j, k;
    for (i = 0; i < n; ++i) {
        if (eigfac->colnnz[i] > 0) {
            eigfac->perm[i] = nsub;
            eigfac->pinv[nsub] = i;
            ++nsub;
        }
    }
    
    DSDP_INT *perm = eigfac->perm;
    double *dmat = eigfac->dworkmat;
    // Construct sub-matrix
    for (k = 0; k < nnz; ++k) {
        i = Ai[k]; j = Aj[k];
        dmat[nsub * perm[i] + perm[j]] = Ax[k];
        if (i != j) {
            dmat[nsub * perm[j] + perm[i]] = Ax[k];
        }
    }
    
    return nsub;
}

static void speigdenseFactor( speigfac *eigfac, DSDP_INT n,
                              double *array, double *eigvals, double *eigvecs) {
    
    // Do eigen decomposition
    double *dwork = eigfac->dwork, *dworkmat = eigfac->dworkmat;
    
    memset(eigvals, 0, sizeof(double) * n);
    memset(eigvecs, 0, sizeof(double) * n * n);
    memcpy(dworkmat, array, sizeof(double) * n * n);
    
    DSDP_INT *iwork = eigfac->iwork, *iworkup = eigfac->iworkup;
    DSDP_INT lwork = eigfac->lwork, liwork = eigfac->liwork, m, info;
    
    dsyevr(&jobz, &range, &uplolow, &n, dworkmat,
           &n, NULL, NULL, NULL, NULL, &abstol, &m,
           eigvals, eigvecs, &n, iworkup, dwork,
           &lwork, iwork, &liwork, &info);
    
    assert( info == 0 );
}

static void speigFactor( speigfac *eigfac, DSDP_INT n, DSDP_INT nsub,
                             double *eigvals, double *eigvecs ) {
    // Do eigen decomposition
    double *dwork = eigfac->dwork, \
           *dworkmat = eigfac->dworkmat,\
           *dworkeval = eigfac->dworkevl, \
           *dworkevec = eigfac->dworkevc;
    
    memset(eigvals, 0, sizeof(double) * n);
    memset(eigvecs, 0, sizeof(double) * n * n);
    
    DSDP_INT *iwork = eigfac->iwork, *iworkup = eigfac->iworkup;
    DSDP_INT lwork = eigfac->lwork, liwork = eigfac->liwork, m, info;
    
    dsyevr(&jobz, &range, &uplolow, &nsub, dworkmat,
           &nsub, NULL, NULL, NULL, NULL, &abstol, &m,
           dworkeval, dworkevec, &nsub, iworkup, dwork,
           &lwork, iwork, &liwork, &info);
    
    assert( info == 0 );
    
    // Extract eigenvalues and transform back into the old array
    DSDP_INT *iperm = eigfac->pinv, i, j, rank = 0;
    double *subvec = NULL, *vec = NULL;
    for (i = 0; i < nsub; ++i) {
        if (fabs(dworkeval[i]) > eps) {
            eigvals[rank] = dworkeval[i];
            subvec = &dworkevec[nsub * i]; vec = &eigvecs[n * rank];
            for (j = 0; j < nsub; ++j) {
                vec[iperm[j]] = subvec[j];
            }
            ++rank;
        }
    }
}

extern void speigInit( speigfac *eigfac ) {
    
    eigfac->nmax   = 0;
    eigfac->lwork  = 0; eigfac->liwork = 0;
    eigfac->perm   = NULL; eigfac->pinv    = NULL;
    eigfac->colnnz = NULL; eigfac->dwork   = NULL;
    eigfac->dworkmat = NULL; eigfac->dworkevc = NULL;
    eigfac->dworkevl = NULL;
    eigfac->iwork  = NULL; eigfac->iworkup = NULL;
    eigfac->factorMethod = EIG_FACTOR_FEAST;
}

extern DSDP_INT speigAlloc( speigfac *eigfac, DSDP_INT nmax ) {
    
    DSDP_INT retcode = DSDP_RETCODE_OK;
    
    eigfac->nmax     = nmax;
    eigfac->lwork    = LWORK * nmax; eigfac->liwork = IWORK * nmax;
    eigfac->perm     = (DSDP_INT *) calloc(nmax, sizeof(DSDP_INT));
    eigfac->pinv     = (DSDP_INT *) calloc(nmax, sizeof(DSDP_INT));
    eigfac->colnnz   = (DSDP_INT *) calloc(nmax, sizeof(DSDP_INT));
    eigfac->dwork    = (double   *) calloc(eigfac->lwork, sizeof(double));
    eigfac->dworkmat = (double   *) calloc(nmax * nmax, sizeof(double));
    eigfac->dworkevl = (double   *) calloc(nmax, sizeof(double));
    eigfac->dworkevc = (double   *) calloc(nmax * nmax, sizeof(double));
    eigfac->iwork    = (DSDP_INT *) calloc(eigfac->liwork, sizeof(DSDP_INT));
    eigfac->iworkup  = (DSDP_INT *) calloc(2 * nmax, sizeof(DSDP_INT));
    
    return retcode;
}

extern DSDP_INT speigSfac( speigfac *eigfac, spsMat *A,
                       double *eigvals, double *eigvecs ) {
    if (A->nnz < 5 || A->nnz >= A->dim * 3 || A->nnz >= 15000) {
        return FALSE;
    }
    DSDP_INT special = 0, nsub = 0; speigReset(eigfac);
    special = speigGetCStats(eigfac, A->dim, A->nnz, A->i, A->cidx);
    
    memset(eigvals, 0, sizeof(double) * A->dim);
    memset(eigvecs, 0, sizeof(double) * A->dim * A->dim);
    
    if (special) {
        speigGetSpecialFactor(eigfac, A->nnz, A->dim, special,
                              A->i, A->cidx, A->x, eigvals, eigvecs);
    } else {
        nsub = speigGetSMat(eigfac, A->dim, A->nnz, A->i, A->cidx, A->x);
        speigFactor(eigfac, A->dim, nsub, eigvals, eigvecs);
    }
    return TRUE;
}

extern void speigDfac( speigfac *eigfac, dsMat *A,
                       double *eigvals, double *eigvecs ) {
    speigdenseFactor(eigfac, A->dim, A->array, eigvals, eigvecs);
}

extern void speigFree( speigfac *eigfac ) {
    
    DSDP_FREE(eigfac->perm);     DSDP_FREE(eigfac->pinv);
    DSDP_FREE(eigfac->colnnz);   DSDP_FREE(eigfac->dwork);
    DSDP_FREE(eigfac->dworkmat); DSDP_FREE(eigfac->dworkevl);
    DSDP_FREE(eigfac->dworkevl); DSDP_FREE(eigfac->iwork);
    DSDP_FREE(eigfac->iworkup);
    
    eigfac->nmax = 0; eigfac->lwork = 0; eigfac->liwork = 0;
    eigfac->factorMethod = 0;
}
