#include "sparsemat.h"
#include "densemat.h"
#include "dsdpfeast.h"
#include "dsdpdata.h"

// Enable hash sum check
#ifdef VERIFY_HASH
#undef VERIFY_HASH
#endif

static char etype[] = "Sparse/Dual matrix";
static DSDP_INT one = 1;
static double done = 1.0;
static double dzero = 0.0;
static char uplolow = 'L';
static char trans = 'T';
static char notrans = 'N';
static DSDP_INT phaseSolve = PARDISO_SOLVE;
static DSDP_INT errorSolve;

/* Internal Pardiso Wrapper */
static void pardisoSymFactorize( spsMat *S ) {
    /* Factorize the spsMat matrix */
    DSDP_INT phase = PARDISO_SYM;
    pardiso(S->pdsWorker, &maxfct, &mnum, &mtype, &phase, &S->dim,
            S->x, S->p, S->i, &idummy, &idummy, PARDISO_PARAMS_CHOLESKY,
            &msglvl, NULL, NULL, &errorSolve);
    if (errorSolve) {
        printf("[Pardiso Error]: Matrix factorization failed."
               " Error code: "ID" \n", errorSolve);
    }
}

/* Internal Pardiso Wrapper */
static void pardisoNumFactorize( spsMat *S ) {
    DSDP_INT phase = PARDISO_FAC;
    // Invoke pardiso to do symbolic analysis and Cholesky factorization
    pardiso(S->pdsWorker, &maxfct, &mnum, &mtype, &phase, &S->dim,
            S->x, S->p, S->i, &idummy, &idummy, PARDISO_PARAMS_CHOLESKY,
            &msglvl, NULL, NULL, &errorSolve) ;
    if (errorSolve) {
        printf("| [Pardiso Error]: Matrix factorization failed with"
               " code: "ID". \n", errorSolve);
    }
}

static void pardisoIndefiniteNumFactorize( spsMat *S ) {
    DSDP_INT phase = PARDISO_FAC;
    // Invoke pardiso to do symbolic analysis and Cholesky factorization
    pardiso(S->pdsWorker, &maxfct, &mnum, &mtype, &phase, &S->dim,
            S->x, S->p, S->i, &idummy, &idummy, PARDISO_PARAMS_LDL,
            &msglvl, NULL, NULL, &errorSolve) ;
    if (errorSolve) {
        printf("| [Pardiso Error]: Matrix factorization failed with"
               " code: "ID". \n", errorSolve);
    }
}

static void pardisoForwardSolve( spsMat *S, DSDP_INT nrhs, double *B, double *aux, DSDP_INT overwrite ) {
    // pardiso forward solve
    DSDP_INT phase = PARDISO_FORWARD, *param;
    param = (overwrite) ? PARDISO_PARAMS_FORWARD_BACKWORD : PARDISO_PARAMS_FORWARD_BACKWORD_LANCZOS;
    pardiso(S->pdsWorker, &maxfct, &mnum, &mtype, &phase, &S->dim,
            NULL, S->p, S->i, &idummy, &nrhs, param, &msglvl,
            B, aux, &errorSolve);;
    if (errorSolve) {
        printf("[Pardiso Error]: Matrix backward solve failed."
               " Error code: "ID" \n", errorSolve);
    }
}

static void pardisoBackwardSolve( spsMat *S, DSDP_INT nrhs, double *B, double *aux, DSDP_INT overwrite ) {
    // pardiso backward solve
    DSDP_INT phase = PARDISO_BACKWARD, *param;
    param = (overwrite) ? PARDISO_PARAMS_FORWARD_BACKWORD : PARDISO_PARAMS_FORWARD_BACKWORD_LANCZOS;
    pardiso(S->pdsWorker, &maxfct, &mnum, &mtype, &phase, &S->dim,
            NULL, S->p, S->i, &idummy, &nrhs, param, &msglvl,
            B, aux, &errorSolve);
    if (errorSolve) {
        printf("[Pardiso Error]: Matrix backward solve failed."
               " Error code: "ID" \n", errorSolve);
    }
}

static void pardisoSolve( spsMat *S, DSDP_INT nrhs, double *B, double *X ) {
    /* Solve the linear system S * X = B using Pardiso */
    // Invoke pardiso to perform solution
    pardiso(S->pdsWorker, &maxfct, &mnum, &mtype, &phaseSolve, &S->dim,
            S->x, S->p, S->i, &idummy,
            &nrhs, PARDISO_PARAMS_CHOLESKY, &msglvl,
            B, X, &errorSolve);
}

static void pardisoSolveInplace( spsMat *S, DSDP_INT nrhs, double *B, double *aux ) {
    /* Solve the linear system S * X = B using Pardiso */
    // Invoke pardiso to perform solution
    pardiso(S->pdsWorker, &maxfct, &mnum, &mtype, &phaseSolve, &S->dim,
            S->x, S->p, S->i, &idummy,
            &nrhs, PARDISO_PARAMS_CHOLESKY_INPLACE, &msglvl,
            B, aux, &errorSolve);
}

static void pardisoSolvespsRHS( spsMat *S, DSDP_INT nnz, DSDP_INT *nzidx,
                                   double *b, double *x ) {
    
    DSDP_INT n = S->dim, i;
    if (nnz > n / 4) {
        pardisoSolve(S, 1, b, x);
    } else {
        double *Sinv = S->Sinv, a;
        memset(x, 0, sizeof(double) * n);
        for (i = 0; i < nnz; ++i) {
            a = b[nzidx[i]];
            daxpy(&n, &a, &Sinv[nzidx[i] * n], &one, x, &one);
        }
    }
}

static void pardisoFree( spsMat *S ) {
    /* Free the internal structure of pardiso */
    DSDP_INT phase = PARDISO_FREE;
    pardiso(S->pdsWorker, &maxfct, &mnum, &mtype, &phase, &S->dim,
            S->x, S->p, S->i, &idummy, &one,
            PARDISO_PARAMS_CHOLESKY, &msglvl, &done, &done, &errorSolve);
}

/* Internal Lapack Wrapper*/
static void lapackNumFactorize( spsMat *S ) {
    // memcpy(S->Sinv, S->x, sizeof(double) * S->nnz);
    dpotrf(&uplolow, &S->dim, S->x, &S->dim, &errorSolve);
}

static void lapackForwardSolveBatch( spsMat *S, DSDP_INT nrhs, double *B ) {
    dtrsm(&uplolow, &uplolow, &notrans, &notrans,
          &S->dim, &nrhs, &done, S->x, &S->dim, B, &S->dim);
}

static void lapackBackwardSolveBatch( spsMat *S, DSDP_INT nrhs, double *B ) {
    dtrsm(&uplolow, &uplolow, &trans, &notrans,
          &S->dim, &nrhs, &done, S->x, &S->dim, B, &S->dim);
}


static void lapackForwardSolve( spsMat *S, DSDP_INT nrhs, double *B, double *aux, DSDP_INT overwrite) {
    memcpy(aux, B, sizeof(double) * S->dim * nrhs);
    for (DSDP_INT i = 0; i < nrhs; ++i) {
        dtrsv(&uplolow, &notrans, &notrans, &S->dim,
              S->x, &S->dim, &aux[i * S->dim], &one);
    }
}

static void lapackBackwardSolve( spsMat *S, DSDP_INT nrhs, double *B, double *aux, DSDP_INT overwrite) {
    memcpy(aux, B, sizeof(double) * S->dim * nrhs);
    for (DSDP_INT i = 0; i < nrhs; ++i) {
        dtrsv(&uplolow, &trans, &notrans, &S->dim,
              S->x, &S->dim, &aux[i * S->dim], &one);
    }
}

static void lapackSolve( spsMat *S, DSDP_INT nrhs, double *B, double *X ) {
    // memcpy(X, B, sizeof(double) * nrhs * S->dim);
    dsymm(&uplolow, &uplolow, &S->dim, &nrhs, &done,
          S->Sinv, &S->dim, B, &S->dim, &dzero, X, &S->dim);
    // dpotrs(&uplolow, &S->dim, &nrhs, S->x, &S->dim, X, &S->dim, &errorSolve);
}

static void lapackSolveInplace( spsMat *S, DSDP_INT nrhs, double *B, double *aux ) {
    memcpy(aux, B, sizeof(double) * nrhs * S->dim);
    dsymm(&uplolow, &uplolow, &S->dim, &nrhs, &done,
          S->Sinv, &S->dim, aux, &S->dim, &dzero, B, &S->dim);
    assert( errorSolve == 0 );
}

static void lapackSolvespsRHS( spsMat *S, DSDP_INT nnz, DSDP_INT *nzidx,
                               double *b, double *x ) {
    // Solve linear system by X = S^-1 * B
    DSDP_INT n = S->dim, i;
    if (nnz > n / 4) {
        lapackSolve(S, 1, b, x);
    } else {
        double *Sinv = S->Sinv, a;
        memset(x, 0, sizeof(double) * n);
        for (i = 0; i < nnz; ++i) {
            a = b[nzidx[i]];
            daxpy(&n, &a, &Sinv[nzidx[i] * n], &one, x, &one);
        }
    }
}

static void arrayTranspose( double *A, DSDP_INT n ) {
    double tmp = 0.0; DSDP_INT i, j;
    for (i = 0; i < n; ++i) {
        for (j = i + 1; j < n; ++j) {
            tmp = A[i * n + j];
            A[i * n + j] = A[j * n + i];
            A[j * n + i] = tmp;
        }
    }
}

static void arraySymmetrize( double *A, DSDP_INT n ) {
    // Symmetrize lower triangular matrix A from Lapack
    DSDP_INT i;
    register double *pc, *pr; pc = pr = A;
    for (i = 0; i < n; ++i) {
        dcopy(&i, pr, &n, pc, &one);
        pc += n; pr += 1;
    }
}

/* Structure operations */
extern DSDP_INT spsMatInit( spsMat *sMat ) {
    // Initialize sparse matrix structure
    DSDP_INT retcode   = DSDP_RETCODE_OK;
    sMat->dim = 0; sMat->factor = NULL;
    sMat->p = NULL; sMat->i = NULL; sMat->x = NULL;
    sMat->isFactorized = FALSE; sMat->cidx = NULL;
    sMat->nominalsps = FALSE; sMat->Sinv = NULL;
    memset(sMat->pdsWorker, 0, PARDISOINDEX * sizeof(void *));
    return retcode;
}

extern DSDP_INT spsMatAlloc( spsMat *sMat, DSDP_INT dim ) {
    // Allocate memory for sparse matrix variable
    DSDP_INT retcode = DSDP_RETCODE_OK;
    assert( sMat->dim == 0 );
    sMat->dim = dim;
    sMat->p = (DSDP_INT *) calloc(dim + 1, sizeof(DSDP_INT));
    sMat->i = (DSDP_INT *) calloc(nsym(dim), sizeof(DSDP_INT));
    sMat->x = (double *) calloc(nsym(dim), sizeof(double));
    return retcode;
}

extern DSDP_INT spsMatAllocData( spsMat *sMat, DSDP_INT dim, DSDP_INT nnz ) {
    // Allocate memory for sparse matrix
    DSDP_INT retcode = DSDP_RETCODE_OK;
    assert( sMat->dim == 0 ); assert( nnz <= nsym(dim));
    sMat->dim = dim; sMat->nnz = nnz;
    if (nnz < nsym(dim)) {
        sMat->p = (DSDP_INT *) calloc(dim + 1, sizeof(DSDP_INT));
        sMat->i = (DSDP_INT *) calloc(nnz, sizeof(DSDP_INT));
        sMat->cidx = (DSDP_INT *) calloc(nnz, sizeof(DSDP_INT));
        sMat->x = (double *) calloc(nnz, sizeof(double));
    } else {
        sMat->nominalsps = TRUE;
        sMat->x = (double *) calloc(dim * dim, sizeof(double));
        sMat->nnz = dim * dim;
    }
    
    return retcode;
}

extern void spsNominalLinkSinv( spsMat *sMat, double *Sinv ) {
    sMat->Sinv = Sinv;
}

extern void spsMatFree( spsMat *sMat ) {
    // Free memory allocated
    sMat->dim = 0; sMat->nnz = 0;
    // Note that we first free p, i and x before calling pardiso to destroy the working array
    if (sMat->isFactorized) {
        pardisoFree(sMat); sMat->isFactorized = FALSE;
    }
    if (sMat->factor) {
        rkMatFree(sMat->factor); DSDP_FREE(sMat->factor);
    }
    DSDP_FREE(sMat->p); DSDP_FREE(sMat->i);
    DSDP_FREE(sMat->x); DSDP_FREE(sMat->cidx);
    sMat->nominalsps = FALSE;
}

extern void spsMatr1check( spsMat *dataMat, DSDP_INT *isRank1 ) {
    // Check if a sparse matrix is rank-one
    DSDP_INT *Ap = dataMat->p, *Ai = dataMat->i, r1 = TRUE;
    double *Ax = dataMat->x, *a = NULL, err = 0.0, diff = 0.0, adiag = 0.0;
    DSDP_INT n = dataMat->dim, i, j, col = 0, isNeg = FALSE, nnz = 0;
    // First detect the first column containing nonzeros
    for (i = 0; i < n; ++i) {
        col = i; if (Ap[i + 1] - Ap[i] > 0) { break; }
    }
    if (Ai[0] != col) {
        r1 = FALSE; *isRank1 = r1; return;
    }
    a = (double *) calloc(n, sizeof(double));
    adiag = Ax[0];
    if (adiag < 0) {
        isNeg = TRUE; adiag = - sqrt(-adiag);
    } else {
        adiag = sqrt(adiag);
    }
    // Get the sparse rank 1 matrix
    for (i = Ap[col]; i < Ap[col + 1]; ++i) {
        // If the diagonal is zero but other rows contain non-zeros
        a[Ai[i]] = Ax[i] / adiag; nnz += (Ax[i] != 0);
    }
    if (dataMat->nnz != nsym(nnz)) { r1 = FALSE; }
    if (r1) {
        for (i = col + 1; i < n; ++i) {
            if (Ap[i] > Ap[col + 1] && Ap[i] < dataMat->nnz) {
                if (Ai[Ap[i]] < col) { r1 = FALSE; break; }
            }
        }
    }
    if (r1) {
        if (isNeg) {
            for (i = 0; i < n; ++i) {
                for (j = Ap[i]; j < Ap[i + 1]; ++j) {
                    diff = Ax[j] + a[i] * a[Ai[j]]; err += fabs(diff);
                }
                if (err > 1e-10) { r1 = FALSE; break; }
            }
        } else {
            for (i = 0; i < n; ++i) {
                for (j = Ap[i]; j < Ap[i + 1]; ++j) {
                    diff = Ax[j] - a[i] * a[Ai[j]]; err += fabs(diff);
                }
                if (err > 1e-10) { r1 = FALSE; break; }
            }
        }
    }
    *isRank1 = (r1) ? (DSDP_INT) (1 - 2 * isNeg) : FALSE;
    DSDP_FREE(a);
}

extern void spsMatr1extract( dsMat *dataMat, double *a, DSDP_INT isNeg ) {
    // Extract the rank 1 data from dense data structure
    double *A = dataMat->array, adiag;
    DSDP_INT n = dataMat->dim, col = 0, i;
    memset(a, 0, sizeof(double) * n);
    for (i = 0; i < n; ++i) { col = i; if (packIdx(A, n, i, i) != 0) { break; } }
    adiag = packIdx(A, n, col, col);
    adiag = (isNeg == -1) ? sqrt(-adiag) : sqrt(adiag);
    for (i = col; i < n; ++i) { a[i] = packIdx(A, n, i, col) / adiag; }
}

/* Basic operations */
extern void spsMatAx( spsMat *A, vec *x, vec *Ax ) {
    /* Sparse matrix multiplication for Lanczos Method
      Warning: This method is only for Lanczos since we compute Ax = - A * x
     */
    DSDP_INT *cidx = A->cidx, *Ap = A->p, *Ai = A->i;
    DSDP_INT n = A->dim, nnz = A->nnz, idx, i, j, k;
    double *Axdata = A->x, *Axres = Ax->x, *xdata = x->x, coeff = -1.0;
    vec_reset(Ax);
    
    if (nnz == n * n) {
        assert( A->nominalsps );
        dsymv(&uplolow, &n, &coeff, Axdata, &n,
              xdata, &one, &dzero, Axres, &one);
    } else if (cidx) {
        for (k = 0; k < nnz; ++k) {
            i = Ai[k]; j = cidx[k]; Axres[i] -= Axdata[k] * xdata[j];
            if (i != j) { Axres[j] -= Axdata[k] * xdata[i]; }
        }
    } else {
        for (i = 0; i < n ; ++i) {
            coeff = xdata[i]; if (coeff == 0.0) { continue; }
            idx = Ap[i];
            
            if (Ai[idx] == i) {
                Axres[i] -= coeff * Axdata[idx];
            } else {
                Axres[Ai[idx]] -= coeff * Axdata[idx];
                Axres[i] -= xdata[Ai[idx]] * Axdata[idx];
            }
            
            for (j = Ap[i] + 1; j < Ap[i + 1]; ++j) {
                idx = Ai[j];
                Axres[idx] -= coeff * Axdata[j];
                Axres[i] -= xdata[idx] * Axdata[j];
            }
        }
    }
}

extern void spsMatAx2( spsMat *A, vec *x, vec *Ax ) {
    // Sparse matrix multiplication for Schur matrix in CG
    assert( A->dim == x->dim && x->dim == Ax->dim );
    DSDP_INT *Ap = A->p, *Ai = A->i, n = A->dim, idx;
    double *Axdata = A->x, *Axres = Ax->x, *xdata = x->x, coeff = -1.0;
    vec_reset(Ax);
    
    for (DSDP_INT i = 0, j; i < n ; ++i) {
        coeff = xdata[i]; if (coeff == 0.0) { continue; }
        idx = Ap[i];
        
        if (Ai[idx] == i) {
            Axres[i] += coeff * Axdata[idx];
        } else {
            Axres[Ai[idx]] += coeff * Axdata[idx];
            Axres[i] += xdata[Ai[idx]] * Axdata[idx];
        }
        
        for (j = Ap[i] + 1; j < Ap[i + 1]; ++j) {
            idx = Ai[j];
            Axres[idx] += coeff * Axdata[j];
            Axres[i] += xdata[idx] * Axdata[j];
        }
    }
}

extern double spsMatxTAx( spsMat *A, double *x ) {
    // Compute quadratic form x' * A * x
    register DSDP_INT *Ai = A->i, *cidx = A->cidx, i, j, k;
    register double res = 0.0, tmp = 0.0, *Ax = A->x;
    for (k = 0; k < A->nnz; ++k) {
        i = cidx[k]; j = Ai[k]; tmp = Ax[k] * x[i] * x[j];
        res += (i == j) ? (0.5 * tmp) : tmp;
    }
    return (2.0 * res);
}

extern void spsMataXpbY( double alpha, spsMat *sXMat, double beta,
                             spsMat *sYMat, DSDP_INT *sumHash ) {
    
    // Matrix axpy operation: let sYMat = alpha * sXMat + beta * sYMat
    // Note that sYMat must have a hash table (unless nomilally sparse)
    if (beta != 1.0) { spsMatScale(sYMat, beta); }
    if (alpha == 0.0) { return; }
    DSDP_INT dim = sXMat->dim, *Ap = sXMat->p, *Ai = sXMat->i;
    double *Ax = sXMat->x, *Bx = sYMat->x;
    if (sumHash) {
        for (DSDP_INT i = 0, j; i < dim; ++i) {
            for (j = Ap[i]; j < Ap[i + 1]; ++j) {
                Bx[packIdx(sumHash, dim, Ai[j], i)] += alpha * Ax[j];
            }
        }
    } else {
        // sYMat is actually dense and Bx is an (n + 1) * n / 2 array
        assert( sYMat->nnz == dim * dim && sYMat->nominalsps );
        // If adding sparse data e.g. S += A * y[i]
        if (sXMat->cidx) {
            for (DSDP_INT i, j, k = 0; k < sXMat->nnz; ++k) {
                i = Ai[k]; j = sXMat->cidx[k];
                fullIdx(Bx, dim, i, j) += alpha * Ax[k];
            }
        // If adding iteration e.g. S += dS
        } else {
            assert( sXMat->nnz == sYMat->nnz );
            daxpy(&sXMat->nnz, &alpha, Ax, &one, Bx, &one);
        }
    }
}

extern void spsMatAdddiag( spsMat *sMat, double d, DSDP_INT *sumHash ) {
    // Add a diagonal element to a sparse matrix with hash table
    if (d == 0.0) { return; }
    DSDP_INT dim = sMat->dim;
    
    if (sumHash) {
        for (DSDP_INT i = 0; i < dim; ++i) {
            sMat->x[packIdx(sumHash, dim, i, i)] += d;
        }
    } else {
        assert( sMat->nominalsps );
        for (DSDP_INT i = 0; i < dim; ++i) {
            sMat->x[i * dim + i] += d;
        }
    }
}

extern void spsMatAddds( spsMat *sXMat, double alpha, dsMat *dsYMat ) {
    // Add a dense matrix to a nominally sparse matrix
    assert( sXMat->nominalsps );
    if (alpha == 0.0) { return; }
    DSDP_INT dim = sXMat->dim, i, j;
//    axpy(&dim, &alpha, dsYMat->array, &one, sXMat->x, &one);
    for (i = 0; i < dim; ++i) {
        for (j = 0; j <= i; ++j) {
            fullIdx(sXMat->x, dim, i, j) += \
                alpha * packIdx(dsYMat->array, dim, i, j);
        }
    }
}

extern void spsMatAddr1( spsMat *sXMat, double alpha, r1Mat *r1YMat, DSDP_INT *sumHash ) {
    // Add a rank 1 matrix to a sparse matrix
    if (alpha == 0.0) { return; }
    DSDP_INT dim = sXMat->dim, *hash = sumHash, *nzIdx = r1YMat->nzIdx;
    double sign = r1YMat->sign * alpha, *rx = r1YMat->x;

    if (sumHash) {
        for (DSDP_INT i = 0; i < r1YMat->nnz; ++i) {
            for (DSDP_INT j = 0; j <= i; ++j) {
                sXMat->x[packIdx(hash, dim, nzIdx[i], nzIdx[j])] +=
                sign * rx[nzIdx[i]] * rx[nzIdx[j]];
            }
        }
    } else {
        DSDP_INT i, j; assert( sXMat->nominalsps );
        if (r1YMat->nnz < 0.5 * dim) {
            for (i = 0; i < r1YMat->nnz; ++i) {
                for (j = 0; j <= i; ++j) {
                    fullIdx(sXMat->x, dim, nzIdx[i], nzIdx[j]) += \
                    sign * rx[nzIdx[i]] * rx[nzIdx[j]];
                }
            }
        } else {
            dsyr(&uplolow, &dim, &sign, r1YMat->x, &one, sXMat->x, &dim);
            // packr1update(&uplo, &dim, &sign, r1YMat->x, &one, sXMat->x);
        }
    }
}

extern void spsMatAddrk( spsMat *sXMat, double alpha, rkMat *rkYMat, DSDP_INT *sumHash ) {
    // Add a rank k matrix to a sparse matrix
    if (alpha == 0.0) { return; }
    for (DSDP_INT i = 0; i < rkYMat->rank; ++i) {
        spsMatAddr1(sXMat, alpha, rkYMat->data[i], sumHash);
    }
}

extern void spsMatScale( spsMat *sXMat, double alpha ) {
    // Scale a sparse matrix by some number.
    if (alpha == 1.0) { return; }
    vecscal(&sXMat->nnz, &alpha, sXMat->x, &one);
}

extern void spsMatRscale( spsMat *sXMat, double r ) {
    // Scale a sparse matrix by the reciprocical of some number. No overflow or under flow
    if (r == 1.0) return;
    if (fabs(r) < 1e-25) {
        printf("Dividing a matrix by 0. \n");
        assert( FALSE );
    }
    if (sXMat->factor) { rkMatRscale(sXMat->factor, r); }
    vecdiv(&sXMat->nnz, &r, sXMat->x, &one);
}

extern void spsMatFnorm( spsMat *sMat, double *fnrm ) {
    
    // Matrix Fronenius norm
    DSDP_INT i, n = sMat->dim;
    double nrm = 0.0; i = spsMatGetRank(sMat);
    if (i < 0.1 * n) {
        rkMatFnorm(sMat->factor, fnrm); return;
    }
    
    if (sMat->cidx) {
        DSDP_INT j, k;
        for (i = 0; i < sMat->nnz; ++i) {
            j = sMat->i[i]; k = sMat->cidx[i];
            if (j == k) {
                nrm += sMat->x[i] * sMat->x[i];
            } else {
                nrm += 2 * sMat->x[i] * sMat->x[i];
            }
        }
        *fnrm = sqrt(nrm);
    } else {
        char ntype = 'F', low = DSDP_MAT_LOW;
        *fnrm = dlansy(&ntype, &low, &n, sMat->x, &n, NULL);
        
    }
}

extern double spsMatOneNorm( spsMat *sMat ) {
    
    // Element-wise sum of absolute values
    double nrm = 0.0;
    DSDP_INT i, j, k;
    if (sMat->cidx) {
        for (i = 0; i < sMat->nnz; ++i) {
            j = sMat->i[i]; k = sMat->cidx[i];
            nrm += (j == k) ? 0.5 * fabs(sMat->x[i]) : fabs(sMat->x[i]);
        }
        nrm *= 2;
    } else {
        assert( FALSE );
        for (i = 0; i < sMat->nnz; ++i) {
            for (j = 0; j <= i; ++j) {
                nrm += fullIdx(sMat->x, sMat->dim, i, j);
            }
        }
    }
    return nrm;
}

/* Factorization and linear system solver */
extern void spsMatSymbolic( spsMat *sAMat ) {
    // Symbolic analysis
    if (sAMat->nominalsps) {
        sAMat->isFactorized = TRUE;
    } else {
        pardisoSymFactorize(sAMat);
        sAMat->isFactorized = TRUE;
    }
}

extern void spsMatFactorize( spsMat *sAMat ) {
    // Sparse matrix Cholesky decomposition
    if (sAMat->nominalsps) {
        lapackNumFactorize(sAMat);
    } else {
        pardisoNumFactorize(sAMat);
    }
}

extern void spsMatIndefiniteFactorize( spsMat *sAMat ) {
    // Sparse matrix Cholesky decomposition
    assert( !sAMat->nominalsps );
    pardisoIndefiniteNumFactorize(sAMat);
}

extern void spsArrSolveInp( spsMat *sAMat, DSDP_INT nrhs, double *B, double *aux ) {
    // Sparse matrix operation X = A \ b
    if (sAMat->nominalsps) {
        lapackSolveInplace(sAMat, nrhs, B, aux);
    } else {
        pardisoSolveInplace(sAMat, nrhs, B, aux);
    }
}

extern void spsMatVecSolve( spsMat *sAMat, vec *sbVec, double *Ainvb ) {
    // Sparse matrix operation X = A \ b
    if (sAMat->nominalsps) {
        lapackSolve(sAMat, 1, sbVec->x, Ainvb);
    } else {
        pardisoSolve(sAMat, 1, sbVec->x, Ainvb);
    }
}

extern void spsMatVecFSolve( spsMat *sAmat, vec *sbVec, vec *Ainvb ) {
    /* Forward solve L * x = b */
    if (sAmat->nominalsps) {
        lapackForwardSolve(sAmat, 1, sbVec->x, Ainvb->x, FALSE);
    } else {
        pardisoForwardSolve(sAmat, 1, sbVec->x, Ainvb->x, FALSE);
    }
    
}

extern void spsMatVecBSolve( spsMat *sAmat, vec *sbVec, vec *Ainvb ) {
    /* Backward solve L' * x = b */
    if (sAmat->nominalsps) {
        lapackBackwardSolve(sAmat, 1, sbVec->x, Ainvb->x, FALSE);
    } else {
        pardisoBackwardSolve(sAmat, 1, sbVec->x, Ainvb->x, FALSE);
    }
}

extern void spsMatGetX( spsMat *S, spsMat *dS, double *LinvSLTinv ) {
    // Routine for retrieving X
    DSDP_INT n = S->dim, i, j;
    assert( dS->dim == n );
    double *fulldS = LinvSLTinv, *aux = (double *) calloc(n * n, sizeof(double)), tmp = 0.0;
    spsMatFill(dS, fulldS);
    
    if (S->nominalsps) {
        lapackForwardSolveBatch(S, n, fulldS);
    } else {
        pardisoForwardSolve(S, n, fulldS, aux, TRUE);
    }
    // Transpose TODO: Avoid transpose using dtrsm
    arrayTranspose(fulldS, n);
    if (S->nominalsps) {
        lapackForwardSolveBatch(S, n, fulldS);
    } else {
        pardisoForwardSolve(S, n, fulldS, aux, TRUE);
    }
    // I + L^-1 * dS * LT^-1
    for (i = 0; i < n; ++i) {
        fulldS[i * n + i] += 1.0;
        for (j = i + 1; j < n; ++j) {
            tmp = (fulldS[i * n + j] + fulldS[j * n + i]) * 0.5;
            fulldS[i * n + j] = fulldS[j * n + i] = tmp;
        }
    }
    
    if (S->nominalsps) {
        lapackBackwardSolveBatch(S, n, fulldS);
    } else {
        pardisoBackwardSolve(S, n, fulldS, aux, TRUE);
    }
    
    arrayTranspose(fulldS, n);
    
    if (S->nominalsps) {
        lapackBackwardSolveBatch(S, n, fulldS);
    } else {
        pardisoBackwardSolve(S, n, fulldS, aux, TRUE);
    }
    
    // Fix numerical instability
    for (i = 0; i < n; ++i) {
        for (j = i + 1; j < n; ++j) {
            tmp = (fulldS[i * n + j] + fulldS[j * n + i]) * 0.5;
            fulldS[i * n + j] = fulldS[j * n + i] = tmp;
        }
    }
    DSDP_FREE(aux);
}

/* DSDP routine for computing the stepsize in the SDP cone */
extern double dsdpGetAlpha( DSDPLanczos *lczSolver, spsMat *S, spsMat *dS, double *alpha ) {
    // Get the maximum alpha such that S + alpha * dS is PSD
    double lbd = 0.0, delta = 0.0;
    if (S->dim == 1) {
        *alpha = (dS->x[0] > 0) ? DSDP_INFINITY : -(S->x[0] * S->x[0]) / dS->x[0];
        return *alpha;
    }
    dsdpLanczosStep(lczSolver, S, dS, &lbd, &delta);
    *alpha = (lbd + delta <= 0) ? DSDP_INFINITY : 1.0 / (lbd + delta);
    return *alpha;
}

/* M3 Technique */
extern double SinvSpSinv( const double *Sinv, double *aux, spsMat *A, dsMat *SinvASinv ) {
    
    // Routine for setting up the Schur matrix
    /*
        Compute inv(S) * A * inv(S) for sparse A
    */
    
    DSDP_INT n = A->dim, *Ai = A->i, *cidx = A->cidx, i, j, k;
    double *SinvASinvdata = SinvASinv->array;
    double *SinvA = aux, *Ax = A->x, *SinvAj, coeff = 0.0, res = 0.0;
    const double *Sinvi; memset(SinvA, 0, sizeof(double) * n * n);
    
    // n f_sigma
    for (k = 0; k < A->nnz; ++k) {
        i = Ai[k]; j = cidx[k]; coeff = Ax[k];
        SinvAj = SinvA + n * j; Sinvi = Sinv + n * i;
        daxpy(&n, &coeff, Sinvi, &one, SinvAj, &one);
        if (i != j) {
            SinvAj = SinvA + n * i; Sinvi = Sinv + n * j;
            daxpy(&n, &coeff, Sinvi, &one, SinvAj, &one);
        }
    }
    
    // Set up inv(S) * A * inv(S) n^3
    for (j = 0; j < n; ++j) {
        SinvAj = SinvA + j; res += SinvAj[j * n];
        for (i = 0; i <= j; ++i) {
            Sinvi = Sinv + n * i;
            packIdx(SinvASinvdata, n, j, i) = ddot(&n, SinvAj, &n, Sinvi, &one);
        }
    }
    
    return res;
}

extern double SinvRkSinv( spsMat *S, rkMat *A, rkMat *SinvASinv ) {
    double res = 0.0; SinvASinv->rank = A->rank;
    for (DSDP_INT i = 0; i < A->rank; ++i) {
        res += SinvR1Sinv(S, A->data[i], SinvASinv->data[i]);
    }
    return res;
}

extern double SinvR1Sinv( spsMat *S, r1Mat *A, r1Mat *SinvASinv ) {
    
    // Routine for setting up the Schur matrix
    DSDP_INT i;
    double *xSinvASinv = SinvASinv->x, *xA = A->x, res = 0.0;
    SinvASinv->sign = A->sign; SinvASinv->nnz = S->dim;
    if (S->nominalsps) {
//        lapackSolve(S, 1, xA, xSinvASinv);
        lapackSolvespsRHS(S, A->nnz, A->nzIdx, xA, xSinvASinv);
    } else {
        pardisoSolvespsRHS(S, A->nnz, A->nzIdx, xA, xSinvASinv);
        // pardisoSolve(S, 1, xA, xSinvASinv);
    }
    
    for (i = 0; i < S->dim; ++i) {
        res += xA[i] * xSinvASinv[i];
    }
    
    return (res * A->sign);
}

/* M4 Technique */
extern double spsSinvSolve( const double *Sinv, spsMat *A, double *ASinv, double *asinv, double Ry ) {
    /*
        Compute A * inv(S) for A (with only lower triangular part)
    */
    DSDP_INT n = A->dim, *Ai = A->i, *cidx = A->cidx, i, j, k;
    double *Ax = A->x, *ASinvi, coeff = 0.0, res = 0.0, res2 = 0.0;
    const double *Sinvj; memset(ASinv, 0, sizeof(double) * n * n);
    
    for (k = 0; k < A->nnz; ++k) {
        i = Ai[k]; j = cidx[k]; coeff = Ax[k];
        // Add coeff * j-th column of Sinv to the i-th row of ASinv
        Sinvj = Sinv + j * n; ASinvi = ASinv + i;
        daxpy(&n, &coeff, Sinvj, &one, ASinvi, &n); // Potential stride optimization
        if (i != j) {
            Sinvj = Sinv + i * n; ASinvi = ASinv + j;
            daxpy(&n, &coeff, Sinvj, &one, ASinvi, &n);
        }
    }
    
    // Compute <A * Sinv>
    for (i = 0; i < n; ++i) {
        res += ASinv[i * n + i];
    }
    
    *asinv = res; if (Ry == 0.0) { return 0.0; }
    
    for (k = 0; k < n; ++k) {
        Sinvj = Sinv + n * k; ASinvi = ASinv + n * k;
        res2 += ddot(&n, Sinvj, &one, ASinvi, &one);
    }
    
    return (res2 * Ry);
}

extern double spsSinvASinv( const double *Sinv, spsMat *A, const double *ASinv ) {
    /*
        Compute A_j * S^-1 * (A_i S^-1)
    */
    DSDP_INT n = A->dim, *Ai = A->i, *cidx = A->cidx, i, j, k;
    double *Ax = A->x, coeff = 0.0, res = 0.0;
    const double *Sinvi, *ASinvj;
    
    for (k = 0; k < A->nnz; ++k) {
        i = Ai[k]; j = cidx[k]; coeff = Ax[k];
        Sinvi = Sinv + n * i; ASinvj = ASinv + n * j;
        res += coeff * ddot(&n, Sinvi, &one, ASinvj, &one);
        if (i != j) {
            Sinvi = Sinv + n * j; ASinvj = ASinv + n * i;
            res += coeff * ddot(&n, Sinvi, &one, ASinvj, &one);
        }
    }
    
    return res;
}

/* M5 Technique */
extern double spsRySinv( spsMat *A, double *Sinv, double *asinv, double Ry ) {
    // Set up <A * S^-1> and <S^-1 * A * S^-1, Ry> by direct computation.
    // This version is used in Phase A
    double res = 0.0, res2 = 0.0, *Ax = A->x;
    DSDP_INT i, p, q, in, n = A->dim, *Ai = A->i, *Aj = A->cidx;
    
    // <A_i, S^-1>
    for (p = 0; p < A->nnz; ++p) {
        i = Ai[p] * n + Aj[p];
        res += (Ai[p] == Aj[p]) ? (0.5 * Ax[p] * Sinv[i]) : (Ax[p] * Sinv[i]);
    }
    // <Ry * S^-1 * A_j, S^-1>
    if (Ry == 0.0) { *asinv = 2.0 * res; return 0.0; }
    
    for (p = 0; p < n; ++p) {
        in = p * n;
        for (q = 0; q < A->nnz; ++q) {
            if (Ai[q] == Aj[q]) {
                res2 += 0.5 * Ax[q] * Sinv[Ai[q] + in] * Sinv[Aj[q] + in];
            } else {
                res2 += Ax[q] * Sinv[Ai[q] + in] * Sinv[Aj[q] + in];
            }
        }
    }
    *asinv = 2.0 * res; return (2.0 * Ry * res2);
}

extern double spsSinvspsSinv( spsMat *A, spsMat *B, double *Sinv ) {
    // Set up <A_i * S^-1 * A_j, S^-1> by direct computation.
    // This simple version is used in Phase B
    double res = 0.0, tmp, aij, *Ax = A->x, *Bx = B->x;
    DSDP_INT i, j, p, q, in, jn, n = A->dim;
    DSDP_INT *Ai = A->i, *Aj = A->cidx, *Bi = B->i, *Bj = B->cidx;
    
    for (p = 0; p < A->nnz; ++p) {
        aij = Ax[p]; i = Ai[p]; j = Aj[p];
        in = i * n; jn = j * n; tmp = 0.0;
        for (q = 0; q < B->nnz; ++q) {
            if (Bi[q] == Bj[q]) {
                tmp += Bx[q] * Sinv[in + Bi[q]] * Sinv[jn + Bj[q]];
            } else {
                tmp += Bx[q] * Sinv[in + Bi[q]] * Sinv[jn + Bj[q]];
                tmp += Bx[q] * Sinv[in + Bj[q]] * Sinv[jn + Bi[q]];
            }
        }
        
        res += (i == j) ? (0.5 * aij * tmp) : (aij * tmp);
    }
    
    return 2.0 * res;
}

extern double spsSinvr1Sinv( spsMat *A, r1Mat *B, double *Sinv ) {
    double res = 0.0, tmp, aij, bij, *Ax = A->x, *Bx = B->x;
    DSDP_INT i, j, k, p, q, in, jn, n = A->dim;
    DSDP_INT *Ai = A->i, *Aj = A->cidx, *Bi = B->nzIdx;
    
    for (p = 0; p < A->nnz; ++p) {
        aij = Ax[p]; i = Ai[p]; j = Aj[p];
        in = i * n; jn = j * n; tmp = 0.0;
        for (q = 0; q < B->nnz; ++q) {
            for (k = 0; k < q; ++k) {
                bij = Bx[Bi[q]] * Bx[Bi[k]];
                tmp += bij * Sinv[in + Bi[q]] * Sinv[jn + Bi[k]];
                tmp += bij * Sinv[jn + Bi[q]] * Sinv[in + Bi[k]];
            }
            k = Bi[q]; tmp += Bx[k] * Bx[k] * Sinv[in + k] * Sinv[jn + k];
        }
        res += (i == j) ? (0.5 * aij * tmp) : (aij * tmp);
    }
    return (2.0 * res * B->sign);
}

extern double spsFullTrace( spsMat *A, double *S ) {
    
    register double res = 0.0, *Ax = A->x;
    DSDP_INT i, j, k, *Ai= A->i, *Aj = A->cidx, n = A->dim;
    for (k = 0; k < A->nnz; ++k) {
        i = Ai[k]; j = Aj[k];
        res += (i == j) ? 0.5 * Ax[k] * S[i * n + j] : Ax[k] * S[i * n + j];
    }
    return 2.0 * res;
}

/* Eigen value routines */
extern void spsMatMaxEig( spsMat *sMat, double *maxEig ) {
    // Eigen value utility: compute the maximum eigenvalue of a matrix
    DSDP_INT n = sMat->dim, info = 0;
    
    double *eigvec = (double *) calloc(n, sizeof(double));
    sparse_matrix_t A = NULL;
    mkl_sparse_d_create_csr(&A, SPARSE_INDEX_BASE_ZERO,
                            n, n, sMat->p, sMat->p + 1,
                            sMat->i, sMat->x);
    info = mkl_sparse_d_ev(&MAX_EIG, pm, A, dsdp_descr,
                           k0, &k, maxEig, eigvec, &resi);
    if (info != SPARSE_STATUS_SUCCESS) {
        printf("Maximum eigen value computation failed. \n");
    }
    mkl_sparse_destroy( A ); // eigvec is freed here
}

extern void spsMatMinEig( spsMat *sMat, double *minEig ) {
    // Eigen value utility: compute the minimum eigenvalue of a matrix
    DSDP_INT n = sMat->dim, info = 0;
    double *eigvec = (double *) calloc(n, sizeof(double));
    sparse_matrix_t A = NULL;
    mkl_sparse_d_create_csr(&A, SPARSE_INDEX_BASE_ZERO,
                            n, n, sMat->p, sMat->p + 1,
                            sMat->i, sMat->x);
    info = mkl_sparse_d_ev(&MIN_EIG, pm, A, dsdp_descr,
                           k0, &k, minEig, eigvec, &resi);
    if (info != SPARSE_STATUS_SUCCESS) {
        printf("Minimum eigen value computation failed. \n");
    }
    mkl_sparse_destroy( A );
}

/* Other utilities */
extern void spsMatIspd( spsMat *sMat, DSDP_INT *ispd ) {
    // A critical routine that determines whether a matrix is positive definite
    if (sMat->nominalsps) {
        dpotrf(&uplolow, &sMat->dim, sMat->x, &sMat->dim, &errorSolve);
        assert( errorSolve >= 0 );
        *ispd = (errorSolve == 0) ? TRUE : FALSE;
    } else {
        DSDP_INT phase = PARDISO_FAC;
        // Invoke pardiso to do symbolic analysis and Cholesky factorization
        pardiso(sMat->pdsWorker, &maxfct, &mnum, &mtype, &phase, &sMat->dim,
                sMat->x, sMat->p, sMat->i, &idummy, &idummy, PARDISO_PARAMS_PSD_CHECK,
                &msglvl, NULL, NULL, &errorSolve);
            
        if (errorSolve == 0) {
            sMat->isFactorized = TRUE; *ispd = TRUE;
        } else if (errorSolve == -4) {
            *ispd = FALSE;
        } else {
           printf("Pardiso fails for some reason. \n");
        }
    }
}

extern double spsMatGetlogdet( spsMat *sMat, double *aux ) {
    // Compute log det S
    double res = 0.0; DSDP_INT i;
    if (sMat->nominalsps) {
        for (i = 0; i < sMat->dim; ++i) {
            res += log(fullIdx(sMat->x, sMat->dim, i, i));
        }
    } else {
        pardiso_getdiag((const void **) sMat->pdsWorker,
                        aux, &aux[sMat->dim], &one, &errorSolve);
        for (i = 0; i < sMat->dim; ++i) { res += log(aux[i]); }
    }
    return (sMat->nominalsps) ? 2.0 * res : res;
}

extern void spsMatInverse( spsMat *sMat, double *Sinv, double *aux ) {
    if (sMat->nominalsps) {
        memcpy(sMat->Sinv, sMat->x, sizeof(double) * sMat->nnz);
        dpotri(&uplolow, &sMat->dim, sMat->Sinv, &sMat->dim, &errorSolve);
        assert( errorSolve == 0 );
        arraySymmetrize(Sinv, sMat->dim);
    } else {
        pardisoSolveInplace(sMat, sMat->dim, Sinv, aux);
    }
    return;
}

extern void spsMatStoreFactor( spsMat *sMat, rkMat *factor ) {
    // Save factorized data
    sMat->factor = factor;
}

extern rkMat* spsMatGetFactor( spsMat *sMat ) {
    return sMat->factor;
}

extern DSDP_INT spsMatGetRank( spsMat *sMat ) {
    return (sMat->factor) ? sMat->factor->rank : (sMat->dim) * 1000;
}

extern void spsMatFillLower( spsMat *sMat, double *lowFullMat ) {
    DSDP_INT n = sMat->dim, ni = 0, *Ap = sMat->p, *Ai = sMat->i, i, j;
    double *Ax = sMat->x;
    for (i = 0; i < n; ++i) {
        ni = n * i;
        for (j = Ap[i]; j < Ap[i + 1]; ++j) {
            lowFullMat[ni + Ai[j]] = Ax[j];
        }
    }
}

extern void spsMatFillLower2( spsMat *sMat, dsMat *lowMat ) {
    DSDP_INT n = sMat->dim, *Ap = sMat->p, *Ai = sMat->i, i, j;
    double *Ax = sMat->x;
    if (sMat->nominalsps) {
        double *p1 = Ax, *p2 = lowMat->array;
        for (i = 0; i < n; ++i) {
            memcpy(p2, p1, sizeof(double) * (n - i));
            p1 += n + 1; p2 += n - i;
        }
    } else {
        for (i = 0; i < n; ++i) {
            for (j = Ap[i]; j < Ap[i + 1]; ++j) {
                packIdx(lowMat->array, n, Ai[j], i) = Ax[j];
            }
        }
    }
}

extern void spsMatFill( spsMat *sMat, double *fulldMat ) {
    // Fill sparse matrix to full (there is no structure for symmetric full dense matrix)
    DSDP_INT n = sMat->dim, *Ap = sMat->p, *Ai = sMat->i, i, j;
    double *Ax = sMat->x;
    
    if (sMat->nominalsps) {
        memcpy(fulldMat, Ax, sizeof(double) * sMat->nnz);
        arraySymmetrize(fulldMat, n);
    } else {
        for (i = 0; i < n; ++i) {
            for (j = Ap[i]; j < Ap[i + 1]; ++j) {
                fulldMat[i * n + Ai[j]] = Ax[j];
                fulldMat[Ai[j] * n + i] = Ax[j];
            }
        }
    }
}

extern void spsMatReset( spsMat *sMat ) {
    assert( sMat->dim > 0 && !sMat->factor ); // Never reset an (eigen) factorized sparse matrix
    memset(sMat->x, 0, sizeof(double) * sMat->nnz);
}

extern void spsMatGetSymbolic( spsMat *sMat, DSDP_INT *hash, DSDP_INT *firstNnz, DSDP_INT *nnzs ) {
    DSDP_INT i, j, k, n = sMat->dim;
    if (sMat->i[0] == 0 && sMat->p[1] > 0) { *firstNnz = TRUE; }
    for (k = 0; k < sMat->nnz; ++k) {
        i = sMat->i[k]; j = sMat->cidx[k];
        if (packIdx(hash, n, i, j) == 0) {
            packIdx(hash, n, i, j) = 1; *nnzs = *nnzs + 1;
        }
    }
}

extern DSDP_INT spsMatIsDiagonal( spsMat *sMat ) {
    if (sMat->nnz != sMat->dim) {
        return FALSE;
    }
    for (DSDP_INT i = 0; i < sMat->dim; ++i) {
        if (sMat->cidx[i] != sMat->i[i]) {
            return FALSE;
        }
    }
    return TRUE;
}

extern double spsMatGetXbound( spsMat *sMat, vec *b ) {
    double trX = 0.0;
    for (DSDP_INT i = 0; i < sMat->dim; ++i) {
        trX += b->x[i] / sMat->x[i];
    }
    return trX;
}

// Debugging
extern void spsMatView( spsMat *sMat ) {
    // View a sparse matrix by calling CXSparse cs_print
    if (sMat->nominalsps) {
        printf("Matrix is dense. \n"); return;
    }
    cs mat; mat.p = sMat->p; mat.i = sMat->i; mat.x = sMat->x;
    mat.nz = -1; mat.nzmax = sMat->nnz; mat.m = sMat->dim;
    mat.n = sMat->dim; cs_print(&mat, FALSE);
}

extern void spsMatLinvView( spsMat *S ) {
    // Lanczos debugging routine. Print P^-1 L^-1 to the screen
    assert( S->isFactorized );
    DSDP_INT n = S->dim, i, j;
    
    double *eye = (double *) calloc(n * n, sizeof(double));
    double *Linv = (double *) calloc(n * n, sizeof(double));
    for (i = 0; i < n; ++i) { eye[n * i + i] = 1.0; }
    pardisoForwardSolve(S, n, eye, Linv, FALSE);
    printf("L^-1 Matrix view: \n");
    for (i = 0; i < n; ++i) {
        for (j = 0; j < n; ++j) {
            printf("%50.50e, ", Linv[n * j + i]);
        }
        printf("\n");
    }
    DSDP_FREE(eye); DSDP_FREE(Linv);
}

extern void spsMatInvView( spsMat *S ) {
    // Lanczos debugging routine. Print P^-1 L^-1 to the screen
    assert( S->isFactorized );
    DSDP_INT n = S->dim, i, j;
    
    double *eye = (double *) calloc(n * n, sizeof(double));
    double *Sinv = (double *) calloc(n * n, sizeof(double));
    for (i = 0; i < n; ++i) { eye[n * i + i] = 1.0; }
    pardisoSolve(S, n, eye, Sinv);
    printf("Inverse Matrix view: \n");
    for (i = 0; i < n; ++i) {
        for (j = 0; j < n; ++j) {
            printf("%20.20e, ", Sinv[i * n + j]);
        }
        printf("\n");
    }
    DSDP_FREE(eye); DSDP_FREE(Sinv);
}

extern void spsMatExport( spsMat *A ) {
    DSDP_INT *Ai = A->i, *Aj = A->cidx, nnz = A->nnz;
    double *Ax = A->x;
    printf("i: \n");
    for (DSDP_INT i = 0; i < nnz; ++i) { printf("%d, ", Ai[i]); }
    printf("\nj: \n");
    for (DSDP_INT i = 0; i < nnz; ++i) { printf("%d, ", Aj[i]); }
    printf("\nx: \n");
    for (DSDP_INT i = 0; i < nnz; ++i) { printf("%20.20e, ", Ax[i]); }
    printf("\n");
    return;
}
