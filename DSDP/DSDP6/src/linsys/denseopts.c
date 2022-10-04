#include "denseopts.h"
#include "dsdplapack.h"
#include "structs.h"
#include "rankkopts.h"

// Error type
static char etype[] = "Dense Operation Error";
static DSDP_INT one = 1;
static double dzero = 0.0;
static double done = 1.0;
static char uplolow = DSDP_MAT_LOW;

#ifndef NB
#define NB 16
#endif

/* Internal Lapack Wrapper */
static DSDP_INT fullFactorize( dsMat *S ) {
    // Factorize the dense Schur matrix
    DSDP_INT retcode = DSDP_RETCODE_OK;
    
    DSDP_INT n = S->dim, info = 0; char uplo = DSDP_MAT_LOW;
    memcpy(S->lfactor, S->array, sizeof(double) * n * n);
    
    if (!S->isillCond) {
        dpotrf(&uplo, &n, S->lfactor, &n, &info);
        if (info > 0) {
            S->isillCond = TRUE;
            memcpy(S->lfactor, S->array, sizeof(double) * n * n);
            dsytrf(&uplo, &n, S->lfactor, &n, S->ipiv, S->work, &S->lwork, &info);
        }
    } else {
        dsytrf(&uplo, &n, S->lfactor, &n, S->ipiv, S->work, &S->lwork, &info);
    }
    
    if (info < 0) {
        printf("Illegal value detected in full dense format. \n");
        fatal_error;
    }
    
    S->isFactorized = TRUE;
    return retcode;
}

static DSDP_INT fullSolve( dsMat *S, DSDP_INT nrhs, double *B, double *X ) {
    // Solve the linear system S * X = B using Lapack full format
    DSDP_INT retcode = DSDP_RETCODE_OK;
    char uplo = DSDP_MAT_LOW; DSDP_INT n = S->dim, info = 0;
    // Copy solution data
    memcpy(X, B, sizeof(double) * nrhs * n);
    
    if (S->isillCond) {
        dsytrs(&uplo, &n, &nrhs, S->lfactor, &n, S->ipiv, X, &n, &info);
    } else {
        dpotrs(&uplo, &n, &nrhs, S->lfactor, &n, X, &n, &info);
    }
    
    if (info < 0) {
        error(etype, "Full linear system solution failed. \n");
    }

    return retcode;
}

static DSDP_INT fullSolveInplace( dsMat *S, DSDP_INT nrhs, double *B ) {
    // Solve the linear system S * X = B inplace
    DSDP_INT retcode = DSDP_RETCODE_OK;
    char uplo = DSDP_MAT_LOW; DSDP_INT n = S->dim, info = 0;
    
    if (S->isillCond) {
        dsytrs(&uplo, &n, &nrhs, S->lfactor, &n, S->ipiv, B, &n, &info);
    } else {
        dpotrs(&uplo, &n, &nrhs, S->lfactor, &n, B, &n, &info);
    }
    
    if (info < 0) {
        error(etype, "Full linear system solution failed. \n");
    }
    
    return retcode;
}

/* Structure operations */
extern void denseMatInit( dsMat *dMat ) {
    // Initialize dense matrix
    dMat->dim = 0; dMat->work = NULL; dMat->array = NULL;
    dMat->lfactor = NULL; dMat->ipiv = NULL; dMat->factor = NULL;
    dMat->lwork = 0; dMat->isFactorized = FALSE; dMat->isillCond = FALSE;
}

extern DSDP_INT denseMatAlloc( dsMat *dMat, DSDP_INT dim, DSDP_INT doFactor ) {
    
    // Allocate memory for dense matrix data
    DSDP_INT retcode = DSDP_RETCODE_OK; dMat->dim = dim;
    
    if (doFactor) {
        dMat->array = (double *) calloc(dim * dim, sizeof(double));
        dMat->lfactor = (double *) calloc(dim * dim, sizeof(double));
        dMat->ipiv  = (DSDP_INT *) calloc(dim, sizeof(DSDP_INT));
        dMat->lwork = NB * dim;
        dMat->work  = (double *) calloc(dMat->lwork, sizeof(double));
        
        if (!dMat->array || !dMat->lfactor || !dMat->ipiv || !dMat->lwork || !dMat->work) {
            printf("| Failed to allocate memory. \n");
            fatal_error_msg(etype);
        }
    } else {
        dMat->array = (double *) calloc((DSDP_INT) nsym(dim), sizeof(double));
        dMat->lwork = 0;
    }
    
    return retcode;
}

extern void denseMatFree( dsMat *dMat ) {
    // Free memory allocated
    if (!dMat) { return; }
    if (!dMat->factor) {
        dMat->dim = 0; dMat->isillCond = FALSE;
        dMat->isFactorized = FALSE; dMat->lwork = 0;
        DSDP_FREE(dMat->array); DSDP_FREE(dMat->lfactor);
        DSDP_FREE(dMat->ipiv); DSDP_FREE(dMat->work);
    } else {
        rkMatFree(dMat->factor);
    }
}

extern void dsr1check( dsMat *dataMat, DSDP_INT *isRank1 ) {
    // Is a dense packed array rank-one? (using matrix product)
    double *A = dataMat->array, *start = NULL, *a = NULL, adiag, err = 0.0, diff = 0.0;
    DSDP_INT n = dataMat->dim, i, j, r1 = TRUE, col = 0, isNeg = FALSE, idx = 0;

    // Get the first column that contains non-zero elements
    for (i = 0; i < n; ++i) {
        if (packIdx(A, n, i, i) != 0) { break; }
    }
    if (i == n) {*isRank1 = FALSE; return;}
    if (i >= n - 1 && !packIdx(A, n, i, i)) {
        *isRank1 = FALSE; return;
    }
    col = i; a = (double *) calloc(n, sizeof(double));
    adiag = packIdx(A, n, col, col);
    if (adiag < 0) {
        isNeg = TRUE; adiag = sqrt(-adiag);
    } else {
        adiag = sqrt(adiag);
    }
    for (i = col; i < n; ++i) {
        a[i] = packIdx(A, n, i, col) / adiag;
    }
    // Check if A = a * a' by computing ||A - a * a'||_F
    if (isNeg) {
        for (i = 0; i < n; ++i) {
            start = &A[idx];
            for (j = 0; j < n - i; ++j) {
                diff = start[j] + a[i] * a[i + j];
                err += diff * diff;
            }
            idx += n - i;
            if (err > 1e-08) { r1 = FALSE; break; }
        }
    } else {
        for (i = 0; i < n; ++i) {
            start = &A[idx];
            for (j = 0; j < n - i; ++j) {
                diff = start[j] - a[i] * a[i + j];
                err += diff * diff;
            }
            idx += n - i;
            if (err > 1e-08) { r1 = FALSE; break; }
        }
    }
    *isRank1 = (r1) ? (DSDP_INT) (1 - 2 * isNeg) : FALSE;
    DSDP_FREE(a);
}

extern void dsr1extract( spsMat *dataMat, double *a, DSDP_INT isNeg ) {
    // Extract the rank 1 data from sparse data structure
    double *Ax  = dataMat->x, adiag = 0.0;
    DSDP_INT n = dataMat->dim, col = 0, *Ap = dataMat->p, *Ai = dataMat->i, i;
    memset(a, 0, sizeof(double) * n);
    for (i = 0; i < n; ++i) { col = i; if (Ap[i + 1] - Ap[i] > 0) { break; } }
    adiag = (isNeg == -1) ? sqrt(-Ax[0]) : sqrt(Ax[0]);
    if (adiag != adiag) { fatal_error; }
    for (i = Ap[col]; i < Ap[col + 1]; ++i) { a[Ai[i]] = Ax[i] / adiag; }
}

/* Basic operations */
extern void denseMataXpbY( double alpha, dsMat *dXMat, double beta, dsMat *dYMat ) {
    // Matrix operation. Let dYMat = alpha * dXMat + beta * dYMat
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
}

extern void denseMataAxpby( dsMat *dAMat, double alpha, vec *x, double beta, vec *Ax ) {
    // Compute Ax = A * x
    char uplo = DSDP_MAT_LOW;
    dsymv(&uplo, &x->dim, &alpha, dAMat->array, &x->dim, x->x, &one, &beta, Ax->x, &one);
    // dspmv(&uplo, &x->dim, &alpha, dAMat->array, x->x, &one, &beta, Ax->x, &one);
}

extern void denseMatAdddiag( dsMat *dAMat, double d ) {
    // A = A * d * I
    double *array = dAMat->array;
    for (DSDP_INT i = 0, idx = 0, n = dAMat->dim; i < n; ++i) {
        array[idx] += d;
        idx += n - i;
    }
}

extern void denseMatAdddiagVec( dsMat *dAMat, vec *d ) {
    // A = A * d * diag(d)
    double *array = dAMat->array;
    for (DSDP_INT i = 0, idx = 0, n = dAMat->dim; i < n; ++i) {
        array[idx] += d->x[i];
        idx += n - i;
    }
}

extern double denseMatxTAx( dsMat *dAMat, double *aux, double *x ) {
    // Compute quadratic form x' * A * x
    packmatvec(&uplolow, &dAMat->dim, &done, dAMat->array, x,
               &one, &dzero, aux, &one);
    return dot(&dAMat->dim, x, &one, aux, &one);
}

extern void denseMatScale( dsMat *dXMat, double a ) {
    if (a == 1.0) { return; }
    DSDP_INT n = nsym(dXMat->dim);
    dscal(&n, &a, dXMat->array, &one);
    if (dXMat->factor) {
        rkMatScale(dXMat->factor, a);
    }
}

extern void denseMatRscale( dsMat *dXMat, double r ) {
    // Scale a matrix by reciprocical without over/under flow
    DSDP_INT n = nsym(dXMat->dim);
    vecdiv(&n, &r, dXMat->array, &one);
    if (dXMat->factor) {
        rkMatRscale(dXMat->factor, r);
    }
}

extern void denseMatFnorm( dsMat *dMat, double *fnrm ) {
    // Note that we only compute the F-norm of data
    char nrm = DSDP_MAT_FNORM, uplo = DSDP_MAT_LOW;
    double *work = NULL; DSDP_INT rank;
    denseMatGetRank(dMat, &rank);
    *fnrm = fnorm(&nrm, &uplo, &dMat->dim, dMat->array, work);
}

extern double denseMatOneNorm( dsMat *dMat ) {
    
    DSDP_INT i, n = dMat->dim;
    double nrm = 0.0;
    for (i = 0; i < nsym(n); ++i) { nrm += fabs(dMat->array[i]); }
    for (i = 0; i < n; ++i) { nrm -= 0.5 * fabs(packIdx(dMat->array, n, i, i)); }
    return 2 * nrm;
}

/* Factorization and linear system solver */
extern DSDP_INT denseMatFactorize( dsMat *dAMat ) {
    // Dense full matrix cholesky factorization
    return fullFactorize(dAMat);;
}

extern void denseVecSolve( dsMat *dAMat, vec *dbVec, double *Ainvb ) {
    // Solve a dense system A * x = b
    fullSolve(dAMat, 1, dbVec->x, Ainvb);
}

extern void denseArrSolveInp( dsMat *S, DSDP_INT nrhs, double *B ) {
    // Solve the linear system S * X = B inplace using Lapack packed format
    fullSolveInplace(S, nrhs, B);
}

/* Schur matrix assembly */
extern void denseSpsTrace( dsMat *dAMat, spsMat *sBMat, double *trace ) {
    // Compute trace (A * B) for dense A and sparse B.
    DSDP_INT n = dAMat->dim;
    double *A = dAMat->array, *Bx = sBMat->x, tmp = 0.0;
    DSDP_INT i, j, k, nnz = sBMat->nnz, *Bi = sBMat->i;
    *trace = 0.0;
    
    if (sBMat->nominalsps) {
        for (i = 0; i < n; ++i) {
            tmp += 0.5 * packIdx(A, n, i, i) * Bx[i * n + i];
            for (j = 0; j < i; ++j) {
                tmp += packIdx(A, n, i, j) * Bx[j * n + i];
            }
        }
        *trace = tmp * 2; return;
    }
    
    if (sBMat->cidx) {
        for (i = 0; i < nnz; ++i) {
            j = sBMat->cidx[i]; k = Bi[i];
            tmp = Bx[i] * packIdx(A, n, k, j);
            *trace += (k == j) ? 0.5 * tmp : tmp;
        }
        *trace *= 2; return;
    }
    fatal_error
}

extern void denseDsTrace( dsMat *dAMat, dsMat *dBMat, double *trace ) {
    // Compute trace (A * B) for dense A and dense B
    DSDP_INT n = dAMat->dim;
    double *A = dAMat->array, *B = dBMat->array, res = 0.0;
    for (DSDP_INT i = 0, j = 0, k = 0; i < n; ++i) {
        res += 0.5 * A[j] * B[j]; k = j;
        for (j = k + 1; j < k + n - i; ++j) {
            res += A[j] * B[j];
        }
    }
    *trace = 2.0 * res;
}

extern double SinvDsSinv( const double *Sinv, double *aux, dsMat *A, dsMat *SinvASinv ) {
    // Routine for setting up the Schur matrix
    DSDP_INT n = A->dim, i, j;
    double *Adata = A->array, *ASinv = aux, *ASinvi,
           *SinvASinvdata = SinvASinv->array, res = 0.0;
    const double *Sinvi;
    
    // Get A * S^-1 (A is dense and in packed form). Calling level-2 blas
    memset(ASinv, 0, sizeof(double) * n * n);
    for (i = 0; i < n; ++i) {
        Sinvi = Sinv + n * i; ASinvi = ASinv + n * i;
        dspmv(&uplolow, &n, &done, Adata, Sinvi, &one, &dzero, ASinvi, &one);
    }
    // Set up inv(S) * A * inv(S) n^3
    for (j = 0; j < n; ++j) {
        ASinvi = ASinv + n * j; res += ASinvi[j];
        for (i = 0; i <= j; ++i) {
            Sinvi = Sinv + n * i;
            packIdx(SinvASinvdata, n, j, i) = ddot(&n, ASinvi, &one, Sinvi, &one);
        }
    }
    return res;
}

extern double denseSinvSolve( const double *Sinv, dsMat *A, double *ASinv, double *asinv, double Ry ) {
    // Routine for setting up the Schur matrix. Seems rarely invoked...
    DSDP_INT n = A->dim, i, j, k;
    double *Ax = A->array, *ASinvi, coeff = 0.0, res = 0.0, res2 = 0.0;
    const double *Sinvj; memset(ASinv, 0, sizeof(double) * n * n);
    
    for (k = 0, j = 0; k < n; ++k) {
        coeff = Ax[j];
        if (fabs(coeff) >= 1e-15) {
            Sinvj = Sinv + k * n; ASinvi = ASinv + k;
            daxpy(&n, &coeff, Sinvj, &one, ASinvi, &n);
        }
        for (i = 1; i < n - k; ++i) {
            coeff = Ax[i + j]; if (fabs(coeff) < 1e-20) { continue; }
            Sinvj = Sinv + k * n; ASinvi = ASinv + (k + i);
            daxpy(&n, &coeff, Sinvj, &one, ASinvi, &n);
            Sinvj = Sinv + (k + i) * n; ASinvi = ASinv + k;
            daxpy(&n, &coeff, Sinvj, &one, ASinvi, &n);
        }
        j += n - k;
    }
    
    for (i = 0; i < n; ++i) { res += ASinv[i * n + i]; }
    *asinv = res; if (Ry == 0.0) { return 0.0; }
    for (k = 0; k < n; ++k) {
        Sinvj = Sinv + n * k; ASinvi = ASinv + n * k;
        res2 += ddot(&n, Sinvj, &one, ASinvi, &one);
    }
    
    return (res2 * Ry);
}

extern double denseSinvSolve2( double *Sinv, dsMat *A, double *asinv, double Ry ) {
    // Routine for more efficiently setting up asinvrysinv in corrector step ASinv is not set up
    DSDP_INT n = A->dim, i, j;
    double asinvasinv = 0.0, asinvres = 0.0, aij = 0.0, *Sinvi, *Sinvj, *Adata = A->array;
    for (i = 0; i < n; ++i) {
        Sinvi = Sinv + n * i; Sinvj = Sinv;
        for (j = 0; j < i; ++j) {
            aij = packIdx(Adata, n, i, j);
            asinvasinv += aij * ddot(&n, Sinvi, &one, Sinvj, &one);
            asinvres += aij * Sinvi[j]; Sinvj += n;
        }
        aij = packIdx(Adata, n, i, i);
        asinvasinv += 0.5 * aij * ddot(&n, Sinvi, &one, Sinvj, &one);
        asinvres += 0.5 * aij * Sinvi[i];
    }
    *asinv = asinvres * 2.0; return (2.0 * asinvasinv * Ry);
}

extern double denseSinvASinv( const double *sv, dsMat *A, const double *asv ) {
    // Compute A_j * S^-1 * (A_i S^-1)
    DSDP_INT n = A->dim, i, j, k;
    double *ax = A->array, c = 0.0, t = 0.0; const double *sk, *ask;
    for (k = 0, j = 0; k < n; ++k) {
        c = ax[j];
        if (fabs(c) >= 1e-20) {
            sk = sv + n * k; ask = asv + n * k;
            t += c * ddot(&n, sk, &one, ask, &one);
        }
        for (i = 1; i < n - k; ++i) {
            c = ax[i + j]; if (fabs(c) < 1e-20) { continue; }
            sk = sv + n * k; ask = asv + n * (k + i);
            t += c * ddot(&n, sk, &one, ask, &one);
            sk = sv + n * (k + i); ask = asv + n * k;
            t += c * ddot(&n, sk, &one, ask, &one);
        }
        j += n - k;
    }
    return t;
}

extern double denseDiagTrace( dsMat *dAMat, double diag ) {
    // Compute trace( A * diag * I ) = diag * trace( A ). Used for Ry
    DSDP_INT n = dAMat->dim, i, k = 0;
    double *x = dAMat->array, t = 0.0;
    for (i = 0; i < n; ++i) { t += x[k]; k += n - i; }
    return (diag * t);
}

extern double denseFullTrace( dsMat *dMat, double *S ) {
    // Compute the inner product between a dense and a full dense matrix
    DSDP_INT n = dMat->dim, i, k;
    double t = 0.0, *x = dMat->array, *s = S;
    for (i = 0; i < n; ++i) {
        t -= 0.5 * x[0] * s[0]; k = n - i;
        t += ddot(&k, x, &one, s, &one);
        x += k; s += n + 1;
    }
    return (2.0 * t);
}

/* Utilities */
extern void denseMatStoreFactor( dsMat *dMat, rkMat *factor ) {
    
    dMat->factor = factor;
}

extern rkMat* denseMatGetFactor( dsMat *dMat ) {
    
    return dMat->factor;
}

extern void denseMatGetRank( dsMat *dMat, DSDP_INT *rank ) {
    
    *rank = (dMat->factor) ? dMat->factor->rank : (dMat->dim * 1000);
}

extern void denseMatFillLow( dsMat *dMat, double *fulldMat ) {
    // Fill packed matrix into a square array
    DSDP_INT n = dMat->dim, idx, k; double *x = NULL;
    
    for (k = 0; k < n; ++k) {
        idx = (DSDP_INT) (2 * n - k + 1) * k / 2;
        x = &(fulldMat[k * n + k]);
        memcpy(x, &dMat->array[idx], sizeof(double) * (n - k));
    }
}

extern void denseMatFillLow2( dsMat *dMat, double *lowData ) {
    
    if (!lowData || !dMat || !dMat->array) {
        return;
    }
    
    memcpy(lowData, dMat->array, sizeof(double) * nsym(dMat->dim));
}

extern void denseMatFill( dsMat *dMat, double *fulldMat ) {
    // Fill packed matrix to full (there is no structure for symmetric full dense matrix)
    DSDP_INT n = dMat->dim, idx, i, k; double *x = NULL;
    
    for (k = 0; k < n; ++k) {
        idx = (DSDP_INT) (2 * n - k + 1) * k / 2;
        x = &(fulldMat[k * n + k]);
        memcpy(x, &dMat->array[idx], sizeof(double) * (n - k));
        idx = k;
        
        x = &(fulldMat[k * n]);
        for (i = 0; i < k; ++i) {
            x[i] = dMat->array[idx]; idx += n - i - 1;
        }
    }
}

extern void denseMatReflex( dsMat *dMat ) {
    // Fill the upper triangular of the dense matrix by the corresponding lower triangular elements
    DSDP_INT n = dMat->dim, i, j;
    double *p1, *p2;
    for (i = 0; i < n; ++i) {
        p1 = dMat->array + i * n + (i + 1);
        p2 = dMat->array + (i + 1) * n + i;
        for (j = 0; j < n - i - 1; ++j) {
            *p2 = *p1; ++p1; p2 += n;
        }
    }
}

extern void denseMatGetdiag( dsMat *dMat, vec *diag ) {
    // diag = diag(dMat).
    DSDP_INT n = dMat->dim, i, idx = 0;
    double *x = diag->x, *array = dMat->array;
    for (i = 0; i < n; ++i) {
        x[i] = array[idx]; idx += n - i;
    }
}

extern void denseMatReset( dsMat *dMat ) {
    // Reset M or a dense matrix to 0
    if (dMat->lwork) {
        memset(dMat->array, 0, sizeof(double) * dMat->dim * dMat->dim);
    } else {
        memset(dMat->array, 0, sizeof(double) * nsym(dMat->dim));
    }
}

extern void denseMatResetFactor( dsMat *dMat ) {
    // Reset Cholesky factor
    if (dMat->isFactorized) {
        memset(dMat->lfactor, 0, sizeof(double) * dMat->dim * dMat->dim);
        dMat->isFactorized = FALSE;
    }
}

extern void denseMatMinEig( dsMat *dMat, double *minEig ) {
    // Compute the minimum eigenvalue of matrix X (for DIMACS Error)
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
    double alpha = -1.0; *minEig = DSDP_INFINITY; lwork *= dim; iwork *= dim;
    dsyevr(&jobz, &range, &uplo, &dim, X, &dim,
           NULL, NULL, &il, &iu, &alpha, &neigs, d,
           eigvecs, &dim, isuppz, eigaux, &lwork, eigintaux,
           &iwork, &info);
    *minEig = d[0];
    DSDP_FREE(X); DSDP_FREE(d); DSDP_FREE(eigvals); DSDP_FREE(eigvecs);
    DSDP_FREE(eigaux); DSDP_FREE(eigintaux);
}

extern void denseMatView( dsMat *dMat ) {
    // Print the upper triangular part of the matrix
    DSDP_INT n = dMat->dim, i, j;
    printf("Matrix view: \n");
    
    for (i = 0; i < n; ++i) {
        printf("R"ID": " , i);
        for (j = 0; j < i; ++j) {
            printf("%-10.3g ", 0.0);
        }
        for (j = i; j < n; ++j) {
            printf("%-10.3g ", packIdx(dMat->array, n, j, i));
        }
        printf("\n");
    }
}
