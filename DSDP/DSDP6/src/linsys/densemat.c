#include "densemat.h"
#include "sparsemat.h"
#include "rankkmat.h"
#include "dsdpeigfact.h"

// Error type
static char etype[] = "Dense Operation Error";
static DSDP_INT one = 1;
static double dzero = 0.0;
static double done = 1.0;
static char uplolow = DSDP_MAT_LOW;

/* Internal Lapack Wrapper */
static DSDP_INT packFactorize( dsMat *S ) {
    
    /* Factorize the dsMat matrix */
    DSDP_INT retcode = DSDP_RETCODE_OK;
    
    assert( (!S->isFactorized) && (S->dim > 0) );
    if (S->isFactorized) {
        error(etype, "Matrix is already factorized. \n");
    }
    
    DSDP_INT n    = S->dim;
    char uplo     = DSDP_MAT_LOW;
    DSDP_INT info = 0;
    
    memcpy(S->lfactor, S->array, sizeof(double) * nsym(n));
    
    if (!S->isillCond) {
        packchol(&uplo, &n, S->lfactor, &info);
        if (info > 0) {
            S->isillCond = TRUE;
            packldl(&uplo, &n, S->lfactor, S->ipiv, &info);
        }
    } else {
        packldl(&uplo, &n, S->lfactor, S->ipiv, &info);
    }
    
    if (info < 0) {
        error(etype, "Illegal value detected in packed dense format. \n");
    }
    
    S->isFactorized = TRUE;
    return retcode;
}

static DSDP_INT packSolve( dsMat *S, DSDP_INT nrhs, double *B, double *X ) {
    
    /* Solve the linear system S * X = B using Lapack packed format */
    DSDP_INT retcode = DSDP_RETCODE_OK;
    assert( (S->isFactorized) && (S->dim > 0) && (nrhs > 0) );
    
    if (!S->isFactorized) {
        error(etype, "Matrix is not yet factorized. \n");
    }
    
    char uplo     = DSDP_MAT_LOW;
    DSDP_INT n    = S->dim;
    DSDP_INT ldb  = S->dim;
    DSDP_INT info = 0;
    
    // Copy solution data
    memcpy(X, B, sizeof(double) * nrhs * n);
    
    if (S->isillCond) {
        ldlsolve(&uplo, &n, &nrhs, S->lfactor, S->ipiv, X, &ldb, &info);
    } else {
        packsolve(&uplo, &n, &nrhs, S->lfactor, X, &ldb, &info);
    }
    
    if (info < 0) {
        error(etype, "Packed linear system solution failed. \n");
    }
    
    return retcode;
}

/* Structure operations */
extern DSDP_INT denseMatInit( dsMat *dMat ) {
    
    // Initialize dense matrix
    DSDP_INT retcode = DSDP_RETCODE_OK;
    dMat->dim     = 0;
    dMat->array   = NULL;
    dMat->lfactor = NULL;
    dMat->ipiv    = NULL;
    dMat->isFactorized = FALSE;
    dMat->isillCond = FALSE;
    
    return retcode;
}

extern DSDP_INT denseMatAlloc( dsMat *dMat, DSDP_INT dim, DSDP_INT doFactor ) {
    
    // Allocate memory for dense matrix data
    DSDP_INT retcode = DSDP_RETCODE_OK;
    assert( dMat->dim == 0 );
    
    dMat->dim   = dim;
    if (FALSE) {
        dMat->array = (double *) calloc(dim * dim, sizeof(double)); // For dsaux
    } else {
        dMat->array = (double *) calloc((DSDP_INT) nsym(dim), sizeof(double));
    }
    
    if (doFactor) {
        dMat->lfactor = (double *) calloc((DSDP_INT) nsym(dim), sizeof(double));
        dMat->ipiv = (DSDP_INT *) calloc(dim, sizeof(DSDP_INT));
    }
    
    return retcode;
}

extern DSDP_INT denseMatFree( dsMat *dMat ) {
    
    // Free memory allocated
    DSDP_INT retcode = DSDP_RETCODE_OK;
    
    if (dMat) {
        dMat->dim = 0;
        dMat->isillCond = FALSE;
        dMat->isFactorized = FALSE;
        DSDP_FREE(dMat->array);
        DSDP_FREE(dMat->lfactor);
        DSDP_FREE(dMat->ipiv);
    }
    
    if (dMat->factor) {
        rkMatFree(dMat->factor);
    }
    
    return retcode;
}

/* Basic operations */
extern DSDP_INT denseMataXpbY( double alpha, dsMat *dXMat, double beta, dsMat *dYMat ) {
    
    // Matrix operation. Let dYMat = alpha * dXMat + beta * dYMat
    DSDP_INT retcode = DSDP_RETCODE_OK;
    
    assert( dXMat->dim == dYMat->dim );
    assert((!dXMat->isFactorized) && (!dYMat->isFactorized));
    
    if (dXMat->isFactorized || dYMat->isFactorized) {
        error(etype, "Adding a factorized matrix. \n");
    }
    
    DSDP_INT dim = nsym(dXMat->dim);
    
    if (beta == 0.0) {
        if (alpha == 0.0) {
            memset(dYMat->array, 0, sizeof(double) * dim);
        } else {
            memcpy(dYMat->array, dXMat->array, sizeof(double) * dim);
            if (alpha != 1.0) {
                vecscal(&dim, &alpha, dYMat->array, &one);
            }
        }
    } else {
        
        if (beta != 1.0) {
            vecscal(&dim, &beta, dYMat->array, &one);
        }
        if (alpha != 0.0) {
            axpy(&dim, &alpha, dXMat->array, &one, dYMat->array, &one);
        }
    }
    
    return retcode;
}

extern DSDP_INT denseMataAxpby( dsMat *dAMat, double alpha, vec *x, double beta, vec *Ax ) {
    // Compute Ax = A * x
    DSDP_INT retcode = DSDP_RETCODE_OK;
    char uplo = 'L';
    dspmv(&uplo, &x->dim, &alpha, dAMat->array, x->x, &one, &beta, Ax->x, &one);
    return retcode;
}

extern DSDP_INT denseMatAdddiag( dsMat *dAMat, double d ) {
    // A = A * d * I
    DSDP_INT retcode = DSDP_RETCODE_OK;
    double *array = dAMat->array;
    for (DSDP_INT i = 0, idx = 0, n = dAMat->dim; i < n; ++i) {
        array[idx] += d;
        idx += n - i;
    }
    return retcode;
}

extern double denseMatxTAx( dsMat *dAMat, double *aux, double *x ) {
    // Compute quadratic form x' * A * x
    packmatvec(&uplolow, &dAMat->dim, &done, dAMat->array, x,
               &one, &dzero, aux, &one);
    return dot(&dAMat->dim, x, &one, aux, &one);
}

extern DSDP_INT denseMatRscale( dsMat *dXMat, double r ) {
    // Scale a matrix by reciprocical without over/under flow
    
    DSDP_INT retcode = DSDP_RETCODE_OK;
    assert( dXMat->dim > 0 );
    DSDP_INT n = nsym(dXMat->dim);
    vecdiv(&n, &r, dXMat->array, &one);
    
    if (dXMat->factor) {
        rkMatRscale(dXMat->factor, r);
    }
    
    return retcode;
}

extern DSDP_INT denseMatFnorm( dsMat *dMat, double *fnrm ) {
    
    DSDP_INT retcode = DSDP_RETCODE_OK;
    assert( dMat->dim > 0);
    
    char nrm = DSDP_MAT_FNORM, uplo = DSDP_MAT_LOW;
    double *work = NULL; DSDP_INT rank;
    denseMatGetRank(dMat, &rank);
    
    if (rank < dMat->dim * 0.3) {
        rkMatFnorm(dMat->factor, fnrm);
    } else {
        *fnrm = fnorm(&nrm, &uplo, &dMat->dim, dMat->array, work);
    }
    
    return retcode;
}

extern DSDP_INT denseMatOneNorm( dsMat *dMat, double *onenrm ) {
    
    DSDP_INT retcode = DSDP_RETCODE_OK;
    DSDP_INT i, n = dMat->dim;
    double nrm = 0.0;
    
    for (i = 0; i < nsym(n); ++i) { nrm += fabs(dMat->array[i]); }
    for (i = 0; i < n; ++i) { nrm -= 0.5 * fabs(packIdx(dMat->array, n, i, i)); }
    *onenrm = 2 * nrm;
    
    return retcode;
}

/* Factorization and linear system solver */
extern DSDP_INT denseMatFactorize( dsMat *dAMat ) {
    
    // Dense packed matrix cholesky factorization
    DSDP_INT retcode = DSDP_RETCODE_OK;
    retcode = packFactorize(dAMat);
    return retcode;
}

extern DSDP_INT denseVecSolve( dsMat *dAMat, vec *dbVec, double *Ainvb ) {
    
    // Solve a dense system A * x = b
    DSDP_INT retcode = DSDP_RETCODE_OK;
    assert(dAMat->dim == dbVec->dim);
    retcode = packSolve(dAMat, 1, dbVec->x, Ainvb);
    return retcode;
}

extern DSDP_INT denseArrSolveInp( dsMat *S, DSDP_INT nrhs, double *B ) {
    
    /* Solve the linear system S * X = B using Lapack packed format */
    DSDP_INT retcode = DSDP_RETCODE_OK;
    assert( (S->isFactorized) && (S->dim > 0) && (nrhs > 0));
    
    if (!S->isFactorized) {
        error(etype, "Matrix is not yet factorized. \n");
    }
    
    char uplo     = DSDP_MAT_LOW;
    DSDP_INT n    = S->dim;
    DSDP_INT ldb  = S->dim;
    DSDP_INT info = 0;
    
    if (S->isillCond) {
        ldlsolve(&uplo, &n, &nrhs, S->lfactor, S->ipiv, B, &ldb, &info);
    } else {
        packsolve(&uplo, &n, &nrhs, S->lfactor, B, &ldb, &info);
    }
    
    if (info < 0) {
        error(etype, "Packed linear system solution failed. \n");
        retcode = DSDP_RETCODE_FAILED;
    }
    
    return retcode;
}

extern DSDP_INT denseSpsSolve( dsMat *dAMat, spsMat *sBMat, double *AinvB ) {
    
    // Solve a dense system A * X = B for sparse B
    DSDP_INT retcode = DSDP_RETCODE_OK;
    
    assert( dAMat->dim == sBMat->dim );
    DSDP_INT n = dAMat->dim;
    
    if (n > DSDP_MEMORY_THRESHOLD) {
        /* Memory friendly strategy as in sparse case */
        vec b;
        vec *pb = &b;
        vec_init(pb);
        vec_alloc(pb, n);
        
        for (DSDP_INT k = 0; k < n; ++k) {
            // Fill b by the k th column of sBMat
            retcode = spsMatScatter(sBMat, pb, k);
            retcode = denseVecSolve(dAMat, pb, &AinvB[k * n]);
        }
        
        vec_free(pb);
        
    } else {
        /* Violent solve: not recommended */
        double *B = NULL;
        B = (double *) calloc(n * n, sizeof(double));
        retcode = spsMatFill(sBMat, B);
        retcode = packSolve(dAMat, n, B, AinvB);
    }
    
    return retcode;
}

extern DSDP_INT denseDsSolve( dsMat *dAMat, dsMat *dBMat, double *AinvB ) {
    
    // Solve a dense system A * X = B for dense A and B
    DSDP_INT retcode = DSDP_RETCODE_OK;
    
    assert( dAMat->dim == dBMat->dim );
    DSDP_INT n = dAMat->dim;
    
    if (n > DSDP_MEMORY_THRESHOLD) {
        /* Memory friendly strategy as in sparse case */
        vec b;
        vec *pb = &b;
        vec_init(pb);
        vec_alloc(pb, n);
        
        for (DSDP_INT k = 0; k < n; ++k) {
            // Fill b by the k th column of sBMat
            retcode = denseMatScatter(dBMat, pb, k);
            retcode = denseVecSolve(dAMat, pb, &AinvB[k * n]);
        }
        
        vec_free(pb);
        
    } else {
        /* Violent solve: not recommended */
        double *B = NULL;
        B = (double *) calloc(n * n, sizeof(double));
        retcode = denseMatFill(dBMat, B);
        retcode = packSolve(dAMat, n, B, AinvB);
    }
    
    return retcode;
}

/* Schur matrix assembly */
extern DSDP_INT denseSpsTrace( dsMat *dAMat, spsMat *sBMat, double *trace ) {
    
    // Compute trace (A * B) for dense A and sparse B.
    /*
     trace (A1 * A2) = \sum_{i, j} a_{i,j} * b_{i, j}
     
     ******************************************
     *    Computationally critical routine    *
     ******************************************
     */
    DSDP_INT retcode = DSDP_RETCODE_OK;
    assert( dAMat->dim == sBMat->dim );
    
    DSDP_INT n = dAMat->dim;
    double *A = dAMat->array, *Bx = sBMat->x, t = 0.0, tmp = 0.0;
    DSDP_INT i, j, k, nnz = sBMat->nnz, *Bp = sBMat->p, *Bi = sBMat->i;
    *trace = 0.0;
    
    if (nnz == nsym(n)) {
        *trace = ddot(&nnz, A, &one, Bx, &one);
        for (DSDP_INT i = 0; i < n; ++i) {
            if (Bi[Bp[i]] == i) {
                *trace -= 0.5 * packIdx(A, n, i, i) * Bx[Bp[i]];
            }
        }
        *trace *= 2;
        return retcode;
    }
    
    if (sBMat->nzHash) {
        for (i = 0; i < nnz; ++i) {
            j = sBMat->nzHash[i]; k = Bi[i];
            tmp = Bx[i] * packIdx(A, n, k, j);
            if (k == j) {
                *trace += 0.5 * tmp;
            } else {
                *trace += tmp;
            }
        }
        *trace *= 2;
        return retcode;
    }
    
    for (i = 0; i < n; ++i) {
        k = Bp[i];
        // The sparse column is empty
        if (k == Bp[i + 1]) {
            if (k == nnz) {
                break;
            } else {
                continue;
            }
        }
        if (Bi[k] == i) {
            tmp = 0.5 * Bx[k] * packIdx(A, n, Bi[k], i);
        } else {
            tmp = Bx[k] * packIdx(A, n, Bi[k], i);
        }
        for (k = Bp[i] + 1; k < Bp[i + 1]; ++k) {
            tmp += Bx[k] * packIdx(A, n, Bi[k], i);
        }
        t += tmp;
    }
    
    *trace = 2.0 * t;
    return retcode;
}

extern DSDP_INT denseDsTrace( dsMat *dAMat, dsMat *dBMat, double *trace ) {
    // Compute trace (A * B) for dense A and dense B
    
    /*
     ******************************************
     *    Computationally critical routine    *
     ******************************************
    */
    
    DSDP_INT retcode = DSDP_RETCODE_OK;
    assert( dAMat->dim == dBMat->dim );
    assert((!dAMat->isFactorized) && (!dBMat->isFactorized));
    
    DSDP_INT n = dAMat->dim;
    double *A = dAMat->array, *B = dBMat->array, res = 0.0;
    
    for (DSDP_INT i = 0, j = 0, k = 0; i < n; ++i) {
        res += 0.5 * A[j] * B[j]; k = j;
        for (j = k + 1; j < k + n - i; ++j) {
            res += A[j] * B[j];
        }
    }
    
    *trace = 2.0 * res;
    return retcode;
}

extern DSDP_INT denseDiagTrace( dsMat *dAMat, double diag, double *trace ) {
    
    // Compute trace( A * diag * I ) = diag * trace( A ). Used for Ry
    DSDP_INT retcode = DSDP_RETCODE_OK;
    DSDP_INT n = dAMat->dim;
    double *array = dAMat->array, mattrace = 0.0;
    register DSDP_INT i, idx = 0;
    
    if (n > 64) {
        
        for (i = 0; i < n - 7; ++i) {
            mattrace += array[idx]; idx += n - i; ++i;
            mattrace += array[idx]; idx += n - i; ++i;
            mattrace += array[idx]; idx += n - i; ++i;
            mattrace += array[idx]; idx += n - i; ++i;
            mattrace += array[idx]; idx += n - i; ++i;
            mattrace += array[idx]; idx += n - i; ++i;
            mattrace += array[idx]; idx += n - i; ++i;
            mattrace += array[idx]; idx += n - i;
        }
        
        if (i < n - 3) {
            mattrace += array[idx]; idx += n - i; ++i;
            mattrace += array[idx]; idx += n - i; ++i;
            mattrace += array[idx]; idx += n - i; ++i;
            mattrace += array[idx]; idx += n - i; ++i;
        }
        
        if (i < n - 1) {
            mattrace += array[idx]; idx += n - i; ++i;
            mattrace += array[idx]; idx += n - i; ++i;
        }
        
        if (i < n) {
            mattrace += array[idx]; idx += n - i;
        }
        
    } else {
        for (i = 0; i < n; ++i) {
            mattrace += packIdx(array, n, i, i);
        }
    }
    
    *trace = diag * mattrace;
    
    return retcode;
}

/* Utilities */
extern DSDP_INT denseMatScatter( dsMat *dMat, vec *b, DSDP_INT k ) {
    
    // Scatter b = dMat(:, k)
    DSDP_INT retcode = DSDP_RETCODE_OK;
    assert( dMat->dim == b->dim );
    
    // Index conversion
    DSDP_INT n   = dMat->dim;
    DSDP_INT idx = (DSDP_INT) (2 * n - k + 1) * k / 2; // simplified packIdx
    memcpy(&b->x[k], &dMat->array[idx], sizeof(double) * (n - k));
    idx = k;
    
    for (DSDP_INT i = 0; i < k; ++i) {
        b->x[i] = dMat->array[idx];
        idx += n - i - 1;
    }
    
    return retcode;
}

extern DSDP_INT denseMatStoreFactor( dsMat *dMat, rkMat *factor ) {
    
    DSDP_INT retcode = DSDP_RETCODE_OK;
    dMat->factor = factor;
    return retcode;
}

extern rkMat* denseMatGetFactor( dsMat *dMat ) {
    return dMat->factor;
}

extern DSDP_INT denseMatGetRank( dsMat *dMat, DSDP_INT *rank ) {
    
    if (dMat->factor) {
        *rank = dMat->factor->rank;
    } else {
        *rank = DSDP_INFINITY;
    }
    
    return DSDP_RETCODE_OK;
}

extern DSDP_INT denseMatFillLow( dsMat *dMat, double *fulldMat ) {
    
    // Fill packed matrix into a square array
    DSDP_INT retcode = DSDP_RETCODE_OK;
    DSDP_INT n = dMat->dim, idx;
    double *x = NULL;
    
    for (DSDP_INT k = 0; k < n; ++k) {
        idx = (DSDP_INT) (2 * n - k + 1) * k / 2;
        x = &(fulldMat[k * n + k]);
        memcpy(x, &dMat->array[idx], sizeof(double) * (n - k));
    }
    
    return retcode;
}

extern DSDP_INT denseMatFill( dsMat *dMat, double *fulldMat ) {
    
    // Fill packed matrix to full (there is no structure for symmetric full dense matrix)
    DSDP_INT retcode = DSDP_RETCODE_OK;
    DSDP_INT n = dMat->dim, idx;
    double *x = NULL;
    
    for (DSDP_INT k = 0; k < n; ++k) {
        idx = (DSDP_INT) (2 * n - k + 1) * k / 2;
        x = &(fulldMat[k * n + k]);
        memcpy(x, &dMat->array[idx], sizeof(double) * (n - k));
        idx = k;
        
        x = &(fulldMat[k * n]);
        for (DSDP_INT i = 0; i < k; ++i) {
            x[i] = dMat->array[idx];
            idx += n - i - 1;
        }
    }
    
    return retcode;
}

extern DSDP_INT denseMatGetdiag( dsMat *dMat, vec *diag ) {
    
    // diag = diag(dMat)
    DSDP_INT retcode = DSDP_RETCODE_OK;
    assert( diag->dim == dMat->dim );
    
    DSDP_INT n    = dMat->dim;
    double *x = diag->x, *array = dMat->array;
    register DSDP_INT i, idx = 0;
    
    if (n >= 10000) {
        
        for (i = 0; i < n - 7; ++i) {
            x[i] = array[idx]; idx += n - i; ++i;
            x[i] = array[idx]; idx += n - i; ++i;
            x[i] = array[idx]; idx += n - i; ++i;
            x[i] = array[idx]; idx += n - i; ++i;
            x[i] = array[idx]; idx += n - i; ++i;
            x[i] = array[idx]; idx += n - i; ++i;
            x[i] = array[idx]; idx += n - i; ++i;
            x[i] = array[idx]; idx += n - i;
        }
        
        if (i < n - 3) {
            x[i] = array[idx]; idx += n - i; ++i;
            x[i] = array[idx]; idx += n - i; ++i;
            x[i] = array[idx]; idx += n - i; ++i;
            x[i] = array[idx]; idx += n - i; ++i;
        }
        
        if (i < n - 1) {
            x[i] = array[idx]; idx += n - i; ++i;
            x[i] = array[idx]; idx += n - i; ++i;
        }
        
        if (i < n) {
            x[i] = array[idx]; idx += n - i; ++i;
        }
        
    } else {
        for (i = 0; i < n; ++i) {
            x[i] = array[idx];
            idx += n - i;
        }
    }
    
    return retcode;
}

extern DSDP_INT denseMatReset( dsMat *dMat ) {
    
    // Reset M to be 0
    DSDP_INT retcode = DSDP_RETCODE_OK;
    
    DSDP_INT n = dMat->dim;
    assert( n > 0 );
    memset(dMat->array, 0, sizeof(double) * nsym(n));
    
    return retcode;
}

extern DSDP_INT denseMatResetFactor( dsMat *dMat ) {
    
    // Reset Cholesky factor
    DSDP_INT retcode = DSDP_RETCODE_OK;
    DSDP_INT n = dMat->dim;
    assert( n > 0 );
    
    if (dMat->isFactorized) {
        dMat->isFactorized = FALSE;
        memset(dMat->lfactor, 0, sizeof(double) * nsym(n));
    }
    
    return retcode;
}

extern DSDP_INT denseMatMinEig( dsMat *dMat, double *minEig ) {
    // Compute the minimum eigenvalue of matrix X (for DIMACS Error)
    DSDP_INT retcode = DSDP_RETCODE_OK;
    DSDP_INT dim = dMat->dim, il = 1, iu = 1, neigs, info, lwork = 30, iwork = 12;
    
    double *X = (double *) calloc(dim * dim, sizeof(double));
    denseMatFill(dMat, X);
    double *eigvals = (double *) calloc(dim, sizeof(double));
    double *eigvecs = (double *) calloc(dim, sizeof(double));
    double *eigaux = (double *) calloc(dim * lwork, sizeof(double));
    double *d = (double *) calloc(dim, sizeof(double));

    DSDP_INT *eigintaux = (DSDP_INT *) calloc(dim * iwork, sizeof(DSDP_INT));
    DSDP_INT isuppz[2] = {0};
    
    char jobz = 'V', range = 'I', uplo = DSDP_MAT_UP;
    
    double alpha = -1.0;
    *minEig = DSDP_INFINITY;
    
    lwork *= dim;
    iwork *= dim;
    
    dsyevr(&jobz, &range, &uplo, &dim, X, &dim,
           NULL, NULL, &il, &iu, &alpha, &neigs, d,
           eigvecs, &dim, isuppz, eigaux, &lwork, eigintaux,
           &iwork, &info);
    
    *minEig = d[0];
    
    DSDP_FREE(X);
    DSDP_FREE(d);
    DSDP_FREE(eigvals);
    DSDP_FREE(eigvecs);
    DSDP_FREE(eigaux);
    DSDP_FREE(eigintaux);
    
//    denseMatFnorm(dMat, &nrm);
//    factorizeDenseData(dMat, -nrm, eigvals, eigvecs);
    return retcode;
}

extern DSDP_INT denseMatView( dsMat *dMat ) {
    
    // Print the upper triangular part of the matrix
    DSDP_INT retcode = DSDP_RETCODE_OK;
    DSDP_INT n = dMat->dim;
    printf("Matrix view: \n");
    
    for (DSDP_INT i = 0; i < n; ++i) {
        printf("R"ID": " , i);
        for (DSDP_INT j = 0; j < i; ++j) {
            printf("%-10.3g ", 0.0);
        }
        for (DSDP_INT j = i; j < n; ++j) {
            printf("%-10.3g ", packIdx(dMat->array, n, j, i));
        }
        printf("\n");
    }
    
    return retcode;
}
