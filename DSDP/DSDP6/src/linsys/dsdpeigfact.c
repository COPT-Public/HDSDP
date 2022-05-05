#include <string.h>
#include "dsdpfeast.h"
#include "densemat.h"
#include "sparsemat.h"
/* Implement eigenvalue decomposition using MKL FEAST algorithm */

static char etype[] = "Feast algorithm";
static char uploup = 'U';
static char uplolow = 'L';

static DSDP_INT C2FIndex( spsMat *A ) {
    // Convert C 0-based index into Fortran 1-based index
    // Somehow stupid and looking for FEAST to add one option
    
    DSDP_INT *Ai = A->i, *Ap = A->p;
    register DSDP_INT i;
    
    for (i = 0; i < A->nnz - 7; ++i) {
        Ai[i] += 1; ++i;
        Ai[i] += 1; ++i;
        Ai[i] += 1; ++i;
        Ai[i] += 1; ++i;
        Ai[i] += 1; ++i;
        Ai[i] += 1; ++i;
        Ai[i] += 1; ++i;
        Ai[i] += 1;
    }
    
    if (i < A->nnz - 3) {
        Ai[i] += 1; ++i;
        Ai[i] += 1; ++i;
        Ai[i] += 1; ++i;
        Ai[i] += 1; ++i;
    }
    
    if (i < A->nnz - 1) {
        Ai[i] += 1; ++i;
        Ai[i] += 1; ++i;
    }
    
    if (i < A->nnz) {
        Ai[i] += 1;
    }
    
    for (i = 0; i < A->dim - 6; ++i) {
        Ap[i] += 1; ++i;
        Ap[i] += 1; ++i;
        Ap[i] += 1; ++i;
        Ap[i] += 1; ++i;
        Ap[i] += 1; ++i;
        Ap[i] += 1; ++i;
        Ap[i] += 1; ++i;
        Ap[i] += 1;
    }
    
    if (i < A->dim - 2) {
        Ap[i] += 1; ++i;
        Ap[i] += 1; ++i;
        Ap[i] += 1; ++i;
        Ap[i] += 1; ++i;
    }
    
    if (i < A->dim) {
        Ap[i] += 1; ++i;
        Ap[i] += 1; ++i;
    }
    
    if (i < A->dim + 1) {
        Ap[i] += 1;
    }
    
    return DSDP_RETCODE_OK;
}

static DSDP_INT F2CIndex( spsMat *A ) {
    // Convert C 0-based index into Fortran 1-based index
    // Somehow stupid and looking for FEAST to add one option
    
    DSDP_INT *Ai = A->i, *Ap = A->p;
    register DSDP_INT i;
    
    for (i = 0; i < A->nnz - 7; ++i) {
        Ai[i] -= 1; ++i;
        Ai[i] -= 1; ++i;
        Ai[i] -= 1; ++i;
        Ai[i] -= 1; ++i;
        Ai[i] -= 1; ++i;
        Ai[i] -= 1; ++i;
        Ai[i] -= 1; ++i;
        Ai[i] -= 1;
    }
    
    if (i < A->nnz - 3) {
        Ai[i] -= 1; ++i;
        Ai[i] -= 1; ++i;
        Ai[i] -= 1; ++i;
        Ai[i] -= 1; ++i;
    }
    
    if (i < A->nnz - 1) {
        Ai[i] -= 1; ++i;
        Ai[i] -= 1; ++i;
    }
    
    if (i < A->nnz) {
        Ai[i] -= 1;
    }
    
    for (i = 0; i < A->dim - 6; ++i) {
        Ap[i] -= 1; ++i;
        Ap[i] -= 1; ++i;
        Ap[i] -= 1; ++i;
        Ap[i] -= 1; ++i;
        Ap[i] -= 1; ++i;
        Ap[i] -= 1; ++i;
        Ap[i] -= 1; ++i;
        Ap[i] -= 1;
    }
    
    if (i < A->dim - 2) {
        Ap[i] -= 1; ++i;
        Ap[i] -= 1; ++i;
        Ap[i] -= 1; ++i;
        Ap[i] -= 1; ++i;
    }
    
    if (i < A->dim) {
        Ap[i] -= 1; ++i;
        Ap[i] -= 1; ++i;
    }
    
    if (i < A->dim + 1) {
        Ap[i] -= 1;
    }
    
    return DSDP_RETCODE_OK;
}

extern DSDP_INT factorizeSpecial( spsMat *A, double *eigvals, double *eigvecs, DSDP_INT *isSpecial ) {
    // Check if the sparse matrix follows special pattern
    // Currently detect diagonal matrix and matrices with at most one entry per column
    DSDP_INT retcode = DSDP_RETCODE_OK;
    DSDP_INT isDiag = TRUE, isElem1 = TRUE, nnz = A->nnz, n = A->dim;
    
    if (nnz > n) {
        *isSpecial = FALSE; return retcode;
    }
    
    memset(eigvals, 0, sizeof(double) * n);
    DSDP_INT *Ai = A->i, *Aj = A->nzHash, idx, i, j, k;
    double *Ax = A->x;
    
    for (k = 0; k < nnz; ++k) {
        i = Ai[k]; j = Aj[k];
        if (eigvals[i] || eigvals[j]) {
            *isSpecial = FALSE; return retcode;
        }
        if (i != j) {isDiag = FALSE;}
        eigvals[i] = eigvals[j] = 1;
    }
    
    *isSpecial = TRUE;
    double *eigvec = NULL, tmp = sqrt(2.0) / 2;
    memset(eigvals, 0, sizeof(double) * n);
    memset(eigvecs, 0, sizeof(double) * n * n);
    
    if (isDiag) {
        for (k = 0; k < A->nnz; ++k) {
            eigvec = &eigvecs[k * n];
            eigvec[Ai[k]] = 1.0; eigvals[k] = Ax[k];
        }
    } else if (isElem1) {
        for (k = idx = 0; k < A->nnz; ++k) {
            i = Ai[k]; j = Aj[k];
            eigvec = eigvecs + idx * n;
            if (i == j) {
                eigvals[idx] = Ax[k]; idx += 1;
            } else {
                eigvec[j] = eigvec[i] = tmp;
                eigvals[idx] = Ax[k]; idx += 1;
                eigvec += n;
                eigvec[i] = tmp; eigvec[j] = -tmp;
                eigvals[idx] = - Ax[k]; idx += 1;
            }
        }
    }
    
    return retcode;
}

extern DSDP_INT factorizeSparseData2( spsMat *A, double elow, double *eigvals, double *eigvecs ) {
    
    DSDP_INT retcode = DSDP_RETCODE_OK;
    DSDP_INT n = A->dim;
    assert( n );
    DSDP_INT loop, neigsfound, info;
    
    memset(eigvals, 0, sizeof(double) * n);
    memset(eigvecs, 0, sizeof(double) * n * n);
    double *Ax = (double *) calloc(n * n, sizeof(double));
    double epsout = 0.0, eup = - elow;
    double *res = (double *) calloc(n, sizeof(double));
    
    spsMatFillLower(A, Ax);
    
    dfeast_syev(&uplolow, &n, Ax, &n, fpm,
                &epsout, &loop, &elow, &eup,
                &n, eigvals, eigvecs, &neigsfound,
                res, &info);
    
    if (info) {
        DSDP_FREE(res);
        DSDP_FREE(Ax);
        error(etype, "Feast algorithm fails on dense data. \n");
    }
    
    DSDP_FREE(res);
    DSDP_FREE(Ax);
    return retcode;
}

extern DSDP_INT factorizeSparseData( spsMat *A, double elow, double *eigvals, double *eigvecs ) {
    
    /*
       Apply FEAST algorithm to sparse data
       Lower triangular in CSC => Upper triangular in CSR
    */
    DSDP_INT retcode = DSDP_RETCODE_OK;
    DSDP_INT n = A->dim;
    assert( n );
    memset(eigvals, 0, sizeof(double) * n);
    memset(eigvecs, 0, sizeof(double) * n * n);
    
    C2FIndex(A);
    DSDP_INT *Ap = A->p, *Ai = A->i, loop, neigsfound, info;
    double *Ax = A->x, epsout = 0.0, eup = - elow;
    double *res = (double *) calloc(n, sizeof(double));
    
    dfeast_scsrev(&uploup, &n, Ax, Ap, Ai, fpm,
                  &epsout, &loop, &elow, &eup,
                  &n, eigvals, eigvecs, &neigsfound,
                  res, &info);
    
    if (info) {
        DSDP_FREE(res);
        error(etype, "Feast algorithm fails on sparse data. \n");
    }
    
    DSDP_FREE(res);
    F2CIndex(A);
    return retcode;
}

extern DSDP_INT factorizeDenseData( dsMat *A, double elow, double *eigvals, double *eigvecs ) {
    
    /*
       Apply FEAST algorithm to dense data
       Lower triangular is supplied
    */
    
    DSDP_INT retcode = DSDP_RETCODE_OK;
    DSDP_INT n = A->dim, loop, neigsfound, info;
    assert( n );
    double *Ax = (double *) calloc(n * n, sizeof(double));
    double epsout = 0.0, eup = - elow;
    double *res = (double *) calloc(n, sizeof(double));
    
    denseMatFillLow(A, Ax);
    dfeast_syev(&uplolow, &n, Ax, &n, fpm,
                &epsout, &loop, &elow, &eup,
                &n, eigvals, eigvecs, &neigsfound,
                res, &info);
    
    if (info) {
        DSDP_FREE(res);
        error(etype, "Feast algorithm fails on dense data. \n");
    }
    
    DSDP_FREE(res);
    DSDP_FREE(Ax);
    return retcode;
}
