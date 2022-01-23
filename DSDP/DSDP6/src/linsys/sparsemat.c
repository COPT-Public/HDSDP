#include "sparsemat.h"
#include "densemat.h"
#include "dsdpfeast.h"

// Enable hash sum check
#ifndef VERIFY_HASH
#define VERIFY_HASH
#endif

static char *etype = "Sparse matrix";
static DSDP_INT one = 1;

/* Internal Pardiso Wrapper */
static DSDP_INT pardisoSymFactorize( spsMat *S ) {
    /* Factorize the spsMat matrix */
    DSDP_INT retcode = DSDP_RETCODE_OK;
    
    // Extract the lower triangular part from CSparse structure
    DSDP_INT *Sp = S->p;
    DSDP_INT *Si = S->i;
    double   *Sx = S->x;
    
    // Get the pardiso parameter
    DSDP_INT phase = PARDISO_SYM;
    DSDP_INT error = 0;
    DSDP_INT n     = S->dim;
    
    // Invoke pardiso to do symbolic analysis and Cholesky factorization
    pardiso(S->pdsWorker, &maxfct, &mnum, &mtype, &phase, &n,
            Sx, Sp, Si, &idummy, &idummy, PARDISO_PARAMS_CHOLESKY,
            &msglvl, NULL, NULL, &error);
    
    assert( error == PARDISO_OK );
    
    if (error) {
        printf("[Pardiso Error]: Matrix factorization failed."
               " Error code: "ID" \n", error);
        error(etype, "Pardiso failes to factorize. \n");
    }
    
    // Complete the factorization
    S->isFactorized = TRUE;
    return retcode;
}

static DSDP_INT pardisoFactorize( spsMat *S ) {
    
    /* Factorize the spsMat matrix */
    DSDP_INT retcode = DSDP_RETCODE_OK;
    
    // Extract the lower triangular part from CSparse structure
    DSDP_INT *Sp = S->p;
    DSDP_INT *Si = S->i;
    double   *Sx = S->x;
    
    // Get the pardiso parameter
    DSDP_INT phase = PARDISO_SYM_FAC;
    DSDP_INT error = 0;
    DSDP_INT n     = S->dim;
    
    // Invoke pardiso to do symbolic analysis and Cholesky factorization
    pardiso(S->pdsWorker, &maxfct, &mnum, &mtype, &phase, &n,
            Sx, Sp, Si, &idummy, &idummy, PARDISO_PARAMS_CHOLESKY,
            &msglvl, NULL, NULL, &error);
    
    assert( error == PARDISO_OK );
    
    if (error) {
        printf("[Pardiso Error]: Matrix factorization failed."
               " Error code: "ID" \n", error);
        error(etype, "Pardiso failes to factorize. \n");
    }
    
    // Complete the factorization
    S->isFactorized = TRUE;
    return retcode;
}

/* Internal Pardiso Wrapper */
static DSDP_INT pardisoNumFactorize( spsMat *S ) {
    
    /*
     Numerically factorize the spsMat matrix
     Reuse symbolic ordering after the first iteration
     */
    
    DSDP_INT retcode = DSDP_RETCODE_OK;
    assert( S->isFactorized == TRUE );
    
    // Extract the lower triangular part from CSparse structure
    DSDP_INT *Sp = S->p;
    DSDP_INT *Si = S->i;
    double   *Sx = S->x;
    
    // Get the pardiso parameter
    DSDP_INT phase = PARDISO_FAC;
    DSDP_INT error = 0;
    DSDP_INT n     = S->dim;
    
    // Invoke pardiso to do symbolic analysis and Cholesky factorization
    pardiso(S->pdsWorker, &maxfct, &mnum, &mtype, &phase, &n,
            Sx, Sp, Si, &idummy, &idummy, PARDISO_PARAMS_CHOLESKY,
            &msglvl, NULL, NULL, &error);
    
    assert( error == PARDISO_OK );
    
    if (error) {
        printf("[Pardiso Error]: Matrix factorization failed."
               " Error code: "ID" \n", error);
        error(etype, "Pardiso failes to factorize. \n");
    }
    
    return retcode;
}

static DSDP_INT pardisoForwardSolve( spsMat *S, DSDP_INT nrhs, double *B, double *aux ) {
    // pardiso forward solve
    DSDP_INT retcode = DSDP_RETCODE_OK;
    
    assert( S->isFactorized == TRUE );
    
    // Extract the lower triangular part from CSparse structure
    DSDP_INT *Sp = S->p;
    DSDP_INT *Si = S->i;
    double   *Sx = S->x;
    
    // Get the pardiso parameter
    DSDP_INT phase = PARDISO_FORWARD;
    DSDP_INT error = 0;
    DSDP_INT n     = S->dim;
    
    // Invoke pardiso to do symbolic analysis and Cholesky factorization
    pardiso(S->pdsWorker, &maxfct, &mnum, &mtype, &phase, &n,
            Sx, Sp, Si, &idummy, &nrhs, PARDISO_PARAMS_FORWARD_BACKWORD,
            &msglvl, B, aux, &error);
    
    assert( error == PARDISO_OK );
    
    if (error) {
        printf("[Pardiso Error]: Matrix forward solve failed."
               " Error code: "ID" \n", error);
        error(etype, "Pardiso failes to factorize. \n");
    }
    
    return retcode;
}

static DSDP_INT pardisoBackwardSolve( spsMat *S, DSDP_INT nrhs, double *B, double *aux ) {
    // pardiso forward solve
    DSDP_INT retcode = DSDP_RETCODE_OK;
    
    assert( S->isFactorized == TRUE );
    
    // Extract the lower triangular part from CSparse structure
    DSDP_INT *Sp = S->p;
    DSDP_INT *Si = S->i;
    double   *Sx = S->x;
    
    // Get the pardiso parameter
    DSDP_INT phase = PARDISO_BACKWARD;
    DSDP_INT error = 0;
    DSDP_INT n     = S->dim;
    
    // Invoke pardiso to do symbolic analysis and Cholesky factorization
    pardiso(S->pdsWorker, &maxfct, &mnum, &mtype, &phase, &n,
            Sx, Sp, Si, &idummy, &nrhs, PARDISO_PARAMS_FORWARD_BACKWORD,
            &msglvl, B, aux, &error);
    
    assert( error == PARDISO_OK );
    
    if (error) {
        printf("[Pardiso Error]: Matrix backward solve failed."
               " Error code: "ID" \n", error);
        error(etype, "Pardiso failes to factorize. \n");
    }
    
    DSDP_FREE(aux);
    return retcode;
}

static DSDP_INT pardisoSolve( spsMat *S, DSDP_INT nrhs, double *B, double *X ) {
    
    /* Solve the linear system S * X = B using Pardiso */
    DSDP_INT retcode = DSDP_RETCODE_OK;
    assert( S->isFactorized == TRUE );
    DSDP_INT phase = PARDISO_SOLVE;
    DSDP_INT error = 0;
    DSDP_INT n     = S->dim;
    
    // Invoke pardiso to perform solution
    pardiso(S->pdsWorker, &maxfct, &mnum, &mtype, &phase, &n,
            S->x, S->p, S->i, &idummy,
            &nrhs, PARDISO_PARAMS_CHOLESKY, &msglvl,
            B, X, &error);
    
    assert( error == PARDISO_OK );
    
    if (error) {
        printf("[Pardiso Error]: Matrix solution failed."
               " Error code: "ID" \n", error);
        error(etype, "Pardiso failes to solve system. \n");
    }
    
    return retcode;
}

static DSDP_INT pardisoFree( spsMat *S ) {
    
    /* Free the internal structure of pardiso */
    DSDP_INT retcode = DSDP_RETCODE_OK;
    assert( S->isFactorized == TRUE );
    
    DSDP_INT phase = PARDISO_FREE;
    DSDP_INT error = 0;
    DSDP_INT n     = S->dim;
    
    pardiso(S->pdsWorker, &maxfct, &mnum, &mtype, &phase, &n,
            S->x, S->p, S->i, &idummy, &idummy,
            PARDISO_PARAMS_CHOLESKY, &msglvl, NULL, NULL, &error);
    
    assert( error == PARDISO_OK );
    
    if (error) {
        printf("[Pardiso Error]: Pardiso free failed."
               " Error code: "ID" \n", error);
        error(etype, "Pardiso failes to free. \n");
    }
    
    return retcode;
}

static void fastTranspose( double *A, double *AT, DSDP_INT row,
                           DSDP_INT col, DSDP_INT m, DSDP_INT n ) {
    
    // Cache oblivious transposition of A matrix
    // Reference: https://github.com/iainkfraser/cache_transpose
    if (m > 4 || n > 4) {
        if (n >= m) {
            DSDP_INT half = n / 2;
            fastTranspose(A, AT, row, col, m, half);
            fastTranspose(A, AT, row, col + half, m, half);
        } else {
            DSDP_INT half = m / 2;
            fastTranspose(A, AT, row, col, half, n);
            fastTranspose(A, AT, row + half, col, half, n);
        }
    } else {
        DSDP_INT r = row + n;
        DSDP_INT c = col + m;
        for (DSDP_INT i = row; i < r; ++i) {
            for (DSDP_INT j = col; j < c; ++j) {
                AT[n * j + i] = A[n * i + j];
            }
        }
    }
}

/* Structure operations */
extern DSDP_INT spsMatInit( spsMat *sMat ) {
    
    // Initialize sparse matrix structure
    DSDP_INT retcode   = DSDP_RETCODE_OK;
    
    sMat->dim          = 0;
    sMat->p            = NULL;
    sMat->i            = NULL;
    sMat->x            = NULL;
    sMat->isFactorized = FALSE;
    
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
    
    // Allocate memory for sparse matrix data
    DSDP_INT retcode = DSDP_RETCODE_OK;
    assert( sMat->dim == 0 );
    assert( nnz <= nsym(dim));
    
    sMat->dim = dim;
    sMat->nnz = nnz;
    
    // Currently do nothing here. TODO: Do Optimization if S is dense
    if (nnz < nsym(dim)) {
        sMat->p = (DSDP_INT *) calloc(dim + 1, sizeof(DSDP_INT));
        sMat->i = (DSDP_INT *) calloc(nnz, sizeof(DSDP_INT));
    } else {
        sMat->p = (DSDP_INT *) calloc(dim + 1, sizeof(DSDP_INT));
        sMat->i = (DSDP_INT *) calloc(nsym(dim), sizeof(DSDP_INT));
    }
    
    sMat->x = (double *) calloc(nnz, sizeof(double));
    
    return retcode;
}

extern DSDP_INT spsMatFree( spsMat *sMat ) {
    
    // Free memory allocated
    DSDP_INT retcode = DSDP_RETCODE_OK;
    sMat->dim = 0;
    sMat->nnz = 0;
    
    // Note that we first free p, i and x before calling pardiso to destroy the working array
    DSDP_FREE(sMat->p);
    DSDP_FREE(sMat->i);
    DSDP_FREE(sMat->x);
    
    if (sMat->isFactorized) {
        retcode = pardisoFree(sMat);
        sMat->isFactorized = FALSE;
    }
        
    return retcode;
}

/* Basic operations */
extern DSDP_INT spsMataXpbY( double alpha, spsMat *sXMat, double beta,
                             spsMat *sYMat, DSDP_INT *sumHash ) {
    
    // Matrix axpy operation: let sYMat = alpha * sXMat + beta * sYMat
    // Note that sYMat must have a hash table
    DSDP_INT retcode = DSDP_RETCODE_OK;

    // The matrices added is NEVER factorized since they purely serve as constant data
    assert ( sXMat->dim == sYMat->dim );
    assert ((!sXMat->isFactorized));
    if (sXMat->isFactorized) {
        error(etype, "Adding a factorized matrix. \n");
    }
    
    if (beta != 1.0) {
        retcode = spsMatScale(sYMat, beta);
    }
    
    // Case if alpha = 0.0
    if (alpha == 0.0) {
        return retcode;
    }
    
    DSDP_INT dim = sXMat->dim;
    DSDP_INT *Ap = sXMat->p;
    DSDP_INT *Ai = sXMat->i;
    double   *Ax = sXMat->x;
    double   *Bx = sYMat->x;
    DSDP_INT hash = 0;
    
    for (DSDP_INT i = 0; i < dim; ++i) {
        for (DSDP_INT j = Ap[i]; j < Ap[i + 1]; ++j) {
#ifdef VERIFY_HASH
            hash = packIdx(sumHash, dim, Ai[j], i);
            assert( hash > 0 || j == 0 );
            Bx[hash] += alpha * Ax[j];
#else
            Bx[packIdx(sumHash, dim, Ai[j], i)] += alpha * Ax[j];
#endif
        }
    }

    return retcode;
}

extern DSDP_INT spsMatAdddiag( spsMat *sMat, double d, DSDP_INT *sumHash ) {
    
    // Add a diagonal element to a sparse matrix with hash table
    DSDP_INT retcode = DSDP_RETCODE_OK;
    
    assert( sMat->dim );
    
    if (d == 0.0) {
        return retcode;
    }
    
    DSDP_INT *hash = sumHash;
    DSDP_INT dim = sMat->dim;
    DSDP_INT idx = 0;
    
    for (DSDP_INT i = 0; i < dim; ++i) {
#ifdef VERIFY_HASH
        idx = packIdx(hash, dim, i, i);
        assert( idx > 0 || i == 0 );
        sMat->x[idx] += d;
#else
        sMat->x[packIdx(hash, dim, i, i)] += d;
#endif
    }
    
    return retcode;
}

extern DSDP_INT spsMatAddds( spsMat *sXMat, double alpha, dsMat *dsYMat ) {
    
    // Add a dense matrix to a nominally sparse matrix
    DSDP_INT retcode = DSDP_RETCODE_OK;
    
    assert( sXMat->nnz == nsym(sXMat->dim) );
    assert( sXMat->dim == dsYMat->dim );
    
    if (alpha == 0.0) {
        return retcode;
    }
    
    DSDP_INT dim = sXMat->dim;
    axpy(&dim, &alpha, dsYMat->array, &one, sXMat->x, &one);
    
    return retcode;
}

extern DSDP_INT spsMatAddr1( spsMat *sXMat, double alpha, r1Mat *r1YMat, DSDP_INT *sumHash ) {
    
    // Add a rank 1 matrix to a sparse matrix
    DSDP_INT retcode = DSDP_RETCODE_OK;
    
    assert( sXMat->dim == r1YMat->dim);
    
    if (alpha == 0.0) {
        return retcode;
    }
    
    DSDP_INT dim = sXMat->dim;
    DSDP_INT idx = 0;
    DSDP_INT *hash = sumHash;
    DSDP_INT *nzIdx = r1YMat->nzIdx;
    
    double sign = (double) r1YMat->sign;
    double *rx = r1YMat->x;
    
    for (DSDP_INT i = 0; i < r1YMat->nnz; ++i) {
        for (DSDP_INT j = 0; j <= i; ++j) {
#ifdef VERIFY_HASH
            idx = packIdx(hash, dim, nzIdx[i], nzIdx[j]);
            assert( idx > 0 || i == 0 );
            sXMat->x[idx] += alpha * sign * rx[nzIdx[i]] * rx[nzIdx[j]];
#else
            sXMat->x[packIdx(hash, dim, nzIdx[i], nzIdx[j])] +=
            alpha * sign * rx[nzIdx[i]] * rx[nzIdx[j]];
#endif
        }
    }
    
    return retcode;
}

extern DSDP_INT spsMatScale( spsMat *sXMat, double alpha ) {
    // Scale a sparse matrix by some number.
    DSDP_INT retcode = DSDP_RETCODE_OK;
    assert( sXMat->dim );
    double *x = sXMat->x;
    vecscal(&sXMat->nnz, &alpha, x, &one);
    return retcode;
}

extern DSDP_INT spsMatRscale( spsMat *sXMat, double r ) {
    // Scale a sparse matrix by the reciprocical of some number. No overflow or under flow
    DSDP_INT retcode = DSDP_RETCODE_OK;
    assert( sXMat->dim );
    
    if (fabs(r) < 1e-10) {
        error(etype, "Dividing a matrix by 0. \n");
    }
    
    double *x = sXMat->x;
    vecdiv(&sXMat->p[sXMat->dim], &r, x, &one);
    return retcode;
}

extern DSDP_INT spsMatFnorm( spsMat *sMat, double *fnrm ) {
    
    // Matrix (Approximate) Fronenius norm
    DSDP_INT retcode = DSDP_RETCODE_OK;
    assert( sMat->dim > 0 );
    
    double nrm = 0.0;
    nrm = norm(&sMat->p[sMat->dim], sMat->x, &one);
    nrm = 2 * nrm * nrm;
    
    // Retrieve the diagonal elements
    for (DSDP_INT i = 0; i < sMat->dim; ++i) {
        if (sMat->i[sMat->p[i]] == i) {
            nrm -= sMat->x[sMat->p[i]] * sMat->x[sMat->p[i]];
        }
    }
    
    *fnrm = sqrt(nrm);
    return retcode;
}

/* Factorization and linear system solver */
extern DSDP_INT spsMatSymbolic( spsMat *sAMat ) {
    
    // Symbolic analysis 
    DSDP_INT retcode = DSDP_RETCODE_OK;
    retcode = pardisoSymFactorize(sAMat);
    return retcode;
}

extern DSDP_INT spsMatFactorize( spsMat *sAMat ) {
    
    // Sparse matrix Cholesky decomposition
    // We allow consecutive factorizations involving iteration variables
    DSDP_INT retcode = DSDP_RETCODE_OK;
    retcode = pardisoNumFactorize(sAMat);
    return retcode;
}

extern DSDP_INT spsMatVecSolve( spsMat *sAMat, vec *sbVec, double *Ainvb ) {
    
    // Sparse matrix operation X = A \ B = inv(A) * B for sparse A and B
    // A full dense matrix inv(A) * B is returned
    DSDP_INT retcode = DSDP_RETCODE_OK;
    // A is factorized and its size must agree with b
    assert( sAMat->dim == sbVec->dim );
    assert( sAMat->isFactorized );
    
    retcode = pardisoSolve(sAMat, 1, sbVec->x, Ainvb);
    
    return retcode;
}

extern DSDP_INT spsMatSpSolve( spsMat *sAMat, spsMat *sBMat, double *AinvB ) {
    
    // Sparse matrix operation X = A \ B = inv(A) * B for sparse A and B
    // A full dense matrix inv(A) * B is returned
    DSDP_INT retcode = DSDP_RETCODE_OK;
    
    // A is factorized, B is not and sizes must agree
    assert( sAMat->dim == sAMat->dim );
    assert( sAMat->isFactorized );
    
    DSDP_INT n = sAMat->dim;
    
    if (n >= DSDP_MEMORY_THRESHOLD) {
        /* Memory friendly strategy, [S \ a1, ..., S \ an] sequentially */
        vec b;
        vec *pb = &b;
        vec_init(pb);
        vec_alloc(pb, n);
        
        for (DSDP_INT k = 0; k < n; ++k) {
            // Fill b by the k th column of sBMat
            retcode = spsMatScatter(sAMat, pb, k);
            retcode = spsMatVecSolve(sAMat, pb, &AinvB[k * n]);
        }
        
        vec_free(pb);

    } else {
        /* Violent parallel solve */
        double *B = NULL;
        B = (double *) calloc(n * n, sizeof(double));
        retcode = spsMatFill(sBMat, B);
        retcode = pardisoSolve(sAMat, n, B, AinvB);
        DSDP_FREE(B);
    }
    
    return retcode;
}

extern DSDP_INT spsMatDsSolve( spsMat *sAMat, dsMat *sBMat, double *AinvB ) {
    
    // Matrix operation X = A \ B = inv(A) * B for sparse A and dense B
    // A full dense matrix inv(A) * B is returned
    DSDP_INT retcode = DSDP_RETCODE_OK;
    
    // A is factorized, B is not and sizes must agree
    assert( sAMat->dim == sBMat->dim );
    assert( sAMat->isFactorized );
    
    DSDP_INT n = sAMat->dim;
    
    if (n >= DSDP_MEMORY_THRESHOLD) {
        /* Memory-friendly strategy, [S \ a1,..., S \ an] sequentially */
        vec b;
        vec *pb = &b;
        vec_init(pb);
        vec_alloc(pb, n);
        
        for (DSDP_INT k = 0; k < n; ++k) {
            // Fill b by the k th column of sBMat
            retcode = denseMatScatter(sBMat, pb, k);
            retcode = spsMatVecSolve(sAMat, pb, &AinvB[k * n]);
        }
        
        vec_free(pb);
        
    } else {
        /* Violent parallel solve */
        double *B = NULL;
        B = (double *) calloc(n * n, sizeof(double));
        retcode = denseMatFill(sBMat, B); checkCode;
        retcode = pardisoSolve(sAMat, n, B, AinvB); checkCode;
        DSDP_FREE(B);
    }
    
    return retcode;
}

extern DSDP_INT spsMatLspLSolve( spsMat *S, spsMat *dS, spsMat *spaux ) {
    // Routine for computing SDP cone maximum stepsize
    DSDP_INT retcode = DSDP_RETCODE_OK;
    
    DSDP_INT n = S->dim;
    assert( dS->dim == n );
    assert( spaux->dim == n);
    
    double *fulldS = (double *) calloc(n * n, sizeof(double));
    double *aux    = (double *) calloc(n * n, sizeof(double));
    
    // L^-1 dS
    retcode = spsMatFill(dS, fulldS);
    retcode = pardisoForwardSolve(S, n, fulldS, aux);
    
    // Transpose
    // fastTranspose(fulldS, aux, 0, 0, n, n);
    double tmp = 0.0;
    for (DSDP_INT i = 0; i < n; ++i) {
        for (DSDP_INT j = 0; j <= i; ++j) {
            tmp = fulldS[i * n + j];
            fulldS[i * n + j] = fulldS[j * n + i];
            fulldS[j * n + i] = tmp;
        }
    }
    
    // L^-T (L^-1 dS)
    // TODO: write a Lanczos algorithm to compute the maximum eigenvalue
    // Currently FEAST is used by transferring the matrix into dense format and is INEFFICIENT
    DSDP_INT *Ap = spaux->p;
    DSDP_INT *Ai = spaux->i;
    double   *Ax = spaux->x;
    
    retcode = pardisoForwardSolve(S, n, fulldS, aux);
    
    DSDP_INT nnz = 0;
    DSDP_INT start = 0;
    for (DSDP_INT i = 0; i < n; ++i) {
        start = i * n;
        for (DSDP_INT j = i; j < n; ++j) {
            tmp = fulldS[start + j];
            if (fabs(tmp) > 1e-12) {
                Ax[nnz] = tmp;
                Ai[nnz] = j;
                nnz += 1;
            }
        }
        Ap[i + 1] = nnz;
    }
    
    spaux->nnz = nnz;
    DSDP_FREE(fulldS);
    DSDP_FREE(aux);
    
    return retcode;
}

/* DSDP routine for computing the stepsize in the SDP cone */
extern DSDP_INT dsdpGetAlpha( spsMat *S, spsMat *dS, spsMat *spaux, double *alpha ) {
    // Get the maximum alpha such that S + alpha * dS is PSD
    DSDP_INT retcode = DSDP_RETCODE_OK;
    double mineig = 0.0;
    
    DSDP_INT n = S->dim;
    assert( n == dS->dim && n == spaux->dim );
    // Two ways for computation
    // The first way computes l1(L^-1 dS L^-T)
    retcode = spsMatLspLSolve(S, dS, spaux); checkCode;
    retcode = spsMatMinEig(spaux, &mineig); checkCode;
    
    if (mineig > 0) {
        *alpha = DSDP_INFINITY;
    } else {
        *alpha = - 1.0 / mineig;
    }
    
    return retcode;
}

extern DSDP_INT dsdpGetAlphaLS( spsMat *S, spsMat *dS, spsMat *Scker,
                                double alphamax, double *alpha, DSDP_INT *sumHash ) {
    // Get the maximum alpha such that S + alpha * dS is PSD by line-search
    DSDP_INT retcode = DSDP_RETCODE_OK;
    
    double step = 1 / alphamax;
    DSDP_INT ispsd = FALSE;
    
    for (DSDP_INT i = 0; ; ++i) {
        if (step <= 1e-08) {
            *alpha = 0.0;
            break;
        }
        
        memcpy(Scker->x, S->x, sizeof(double) * S->nnz);
        spsMataXpbY(step, dS, 1.0, Scker, sumHash);
        spsMatIspd(Scker, &ispsd);
        
        if (ispsd) {
            *alpha = step;
            break;
        }
        step *= 0.8;
    }

    return retcode;
}

extern DSDP_INT spsSinvSpSinvSolve( spsMat *S, spsMat *A, dsMat *SinvASinv, double *asinv ) {
    
    // Routine for setting up the Schur matrix
    /*
        Compute inv(S) * A * inv(S) for sparse A by inv(S) * transpose(inv(S) * A)
        The routine performs Cholesky - Transpose - Cholesky sequentially
     
    */
    DSDP_INT retcode = DSDP_RETCODE_OK;
    assert( S->dim == A->dim );
    assert( S->dim == SinvASinv->dim );
    assert( S->isFactorized );
    
    DSDP_INT n         = S->dim;
    double *SinvA      = NULL;
    double *fSinvASinv = NULL;
    SinvA      = (double *) calloc(n * n, sizeof(double));
    fSinvASinv = (double *) calloc(n * n, sizeof(double));
    
    // Solve to get inv(S) * A
    retcode    = spsMatSpSolve(S, A, SinvA);
    double tmp = 0.0;
    
    // Transpose
    *asinv = 0.0;
    for (DSDP_INT i = 0; i < n; ++i) {
        for (DSDP_INT j = 0; j <= i; ++j) {
            if (i == j) {
                *asinv += SinvA[n * i + i];
            } else {
                tmp = SinvA[i * n + j];
                SinvA[i * n + j] = SinvA[j * n + i];
                SinvA[j * n + i] = tmp;
            }
        }
    }
    
    // Solve
    retcode      = pardisoSolve(S, n, SinvA, fSinvASinv);
    DSDP_INT idx = 0;
    
    // Extract solution
    for (DSDP_INT k = 0; k < n; ++k) {
        memcpy(&(SinvASinv->array[idx]), &(fSinvASinv[k * n + k]),
               sizeof(double) * (n - k));
        idx += n - k;
    }
    
    DSDP_FREE(SinvA);
    DSDP_FREE(fSinvASinv);
    return retcode;
}

extern DSDP_INT spsSinvDsSinvSolve( spsMat *S, dsMat *A, dsMat *SinvASinv, double *asinv ) {
    
    // Routine for setting up the Schur matrix
    /*
        Compute inv(S) * A * inv(S) for packed A by inv(S) * transpose(inv(S) * A)
        The routine performs Cholesky - Transpose - Cholesky sequentially
    
    */
    
    DSDP_INT retcode = DSDP_RETCODE_OK;
    
    assert( S->dim == A->dim );
    assert( S->dim == SinvASinv->dim );
    assert( S->isFactorized );
    
    DSDP_INT n          = S->dim;
    double *SinvA       = NULL;
    double *fSinvASinv  = NULL;
    SinvA               = (double *) calloc(n * n, sizeof(double));
    fSinvASinv          = (double *) calloc(n * n, sizeof(double));
    retcode             = spsMatDsSolve(S, A, SinvA);
    double tmp          = 0.0;
    *asinv              = 0.0;
    
    // TODO: Replace the transpose by a cache-aware version
    for (DSDP_INT i = 0; i < n; ++i) {
        for (DSDP_INT j = 0; j <= i; ++j) {
            if (i == j) {
                *asinv += SinvA[n * i + i];
            } else {
                tmp = SinvA[i * n + j];
                SinvA[i * n + j] = SinvA[j * n + i];
                SinvA[j * n + i] = tmp;
            }
        }
    }
    
    retcode = pardisoSolve(S, n, SinvA, fSinvASinv);
    DSDP_INT idx = 0;
    for (DSDP_INT k = 0; k < n; ++k) {
        memcpy(&(SinvASinv->array[idx]), &(fSinvASinv[k * n + k]),
               sizeof(double) * (n - k));
        idx += n - k;
    }
    DSDP_FREE(SinvA);
    DSDP_FREE(fSinvASinv);
        
    return retcode;
}

extern DSDP_INT spsSinvR1SinvSolve( spsMat *S, r1Mat *A, r1Mat *SinvASinv, double *asinv ) {
    
    // Routine for setting up the Schur matrix
    DSDP_INT retcode = DSDP_RETCODE_OK;
    
    assert( S->dim == A->dim );
    assert( SinvASinv->dim == S->dim );
    
    DSDP_INT n = S->dim;
    
    if (S->dim != A->dim) {
        error(etype, "Matrix size mismatch. \n");
    }
    
    double *xSinvASinv = SinvASinv->x;
    double *xA = A->x;
    
    retcode = pardisoSolve(S, 1, xA, xSinvASinv); checkCode;
    
    double res = 0.0;
    
    if (n > 64) {
        for (DSDP_INT i = 0; i < n - n % 4; i+=4) {
            res += xA[i    ] * xSinvASinv[i    ];
            res += xA[i + 1] * xSinvASinv[i + 1];
            res += xA[i + 2] * xSinvASinv[i + 2];
            res += xA[i + 3] * xSinvASinv[i + 3];
        }
        
        for (DSDP_INT i = n - n % 4; i < n; ++i) {
            res += xA[i] * xSinvASinv[i];
        }
    } else {
        for (DSDP_INT i = 0; i < n; ++i) {
            res += xA[i] * xSinvASinv[i];
        }

    }
    
    *asinv = res;
    
    return retcode;
}

/* Eigen value routines */
extern DSDP_INT spsMatMaxEig( spsMat *sMat, double *maxEig ) {
    // Eigen value utility: compute the maximum eigenvalue of a matrix
    DSDP_INT retcode = DSDP_RETCODE_OK;
    DSDP_INT n = sMat->dim;
    DSDP_INT info = 0;
    
    double *eigvec = (double *) calloc(n, sizeof(double));
    
    sparse_matrix_t A = NULL;
    
    mkl_sparse_d_create_csr(&A, SPARSE_INDEX_BASE_ZERO,
                            n, n, sMat->p, sMat->p + 1,
                            sMat->i, sMat->x);
    
    info = mkl_sparse_d_ev(&MAX_EIG, pm, A, dsdp_descr,
                           k0, &k, maxEig, eigvec, &resi);
    
    if (info != SPARSE_SUCCESS) {
        error(etype, "Maximum eigen value computation failed. \n");
    }
    
    mkl_sparse_destroy( A );
    
    return retcode;
}

extern DSDP_INT spsMatMinEig( spsMat *sMat, double *minEig ) {
    // Eigen value utility: compute the minimum eigenvalue of a matrix
    DSDP_INT retcode = DSDP_RETCODE_OK;
    DSDP_INT n = sMat->dim;
    DSDP_INT info = 0;
    
    double *eigvec = (double *) calloc(n, sizeof(double));
    
    sparse_matrix_t A = NULL;
    
    mkl_sparse_d_create_csr(&A, SPARSE_INDEX_BASE_ZERO,
                            n, n, sMat->p, sMat->p + 1,
                            sMat->i, sMat->x);
    
    info = mkl_sparse_d_ev(&MIN_EIG, pm, A, dsdp_descr,
                           k0, &k, minEig, eigvec, &resi);
    
    if (info != SPARSE_SUCCESS) {
        error(etype, "Minimum eigen value computation failed. \n");
    }
    
    mkl_sparse_destroy( A );
    
    return retcode;
}

/* Other utilities */
extern DSDP_INT spsMatIspd( spsMat *sMat, DSDP_INT *ispd ) {
    // A critical routine that determines whether a matrix is positive definite
    DSDP_INT retcode = DSDP_RETCODE_OK;
    
    // Get the pardiso parameter
    DSDP_INT phase = PARDISO_FAC;
    DSDP_INT error = 0;
    
    // Invoke pardiso to do symbolic analysis and Cholesky factorization
    pardiso(sMat->pdsWorker, &maxfct, &mnum, &mtype, &phase, &sMat->dim,
            sMat->x, sMat->p, sMat->i, &idummy, &idummy, PARDISO_PARAMS_CHOLESKY,
            &msglvl, NULL, NULL, &error);
        
    if (error == 0) {
        sMat->isFactorized = TRUE;
        *ispd = TRUE;
    } else if (error == -4) {
        *ispd = FALSE;
    } else {
        error(etype, "Pardiso failes for some reason. \n");
    }
    
    return retcode;
}

extern DSDP_INT spsMatGetlogdet( spsMat *sMat, double *logdet ) {
    // Compute log det S
    DSDP_INT retcode = DSDP_RETCODE_OK;
    
    assert( sMat->isFactorized );
    double *array = (double *) calloc(2 * sMat->dim, sizeof(double));
    double res = 0.0;
    DSDP_INT ecode = PARDISO_OK;
    
    pardiso_getdiag((const void **) sMat->pdsWorker, array,
                    &array[sMat->dim], &one, &ecode);
    
    if (ecode != PARDISO_OK) {
        DSDP_FREE(array);
        error(etype, "Pardiso Failed to get diagonal entries. \n");
    }
    
    for (DSDP_INT i = 0; i < sMat->dim; ++i) {
        res += log(array[i]);
    }
    
    *logdet = res;
    DSDP_FREE(array);
    return retcode;
}

extern DSDP_INT spsMatScatter( spsMat *sMat, vec *b, DSDP_INT k ) {
    
    // Let b = sMat(:, k)
    DSDP_INT retcode = DSDP_RETCODE_OK;
    DSDP_INT dim = sMat->dim;
    
    assert( dim == b->dim );
    
    DSDP_INT *Ap = sMat->p;
    DSDP_INT *Ai = sMat->i;
    double   *Ax = sMat->x;
    
    memset(b->x, 0, sizeof(double) * dim);
    
    // Copy the lower triangular part of sMat(:, k) to x
    for (DSDP_INT i = Ap[k]; i < Ap[k + 1]; ++i) {
        b->x[Ai[i]] = Ax[i];
    }
    
    // Search for the rest
    for (DSDP_INT j = 0; j < k; ++j) {
        for (DSDP_INT i = Ap[j]; i < Ap[j + 1]; ++i) {
            if (Ai[i] > k) {
                break;
            }
            if (Ai[i] == k) {
                b->x[j] = Ax[i];
                break;
            }
        }
    }
    
    return retcode;
}

extern DSDP_INT spsMatFill( spsMat *sMat, double *fulldMat ) {
    
    // Fill sparse matrix to full (there is no structure for symmetric full dense matrix)
    DSDP_INT retcode = DSDP_RETCODE_OK;
    assert( sMat->dim > 0 );
    
    DSDP_INT n = sMat->dim;
    vec b;
    vec *pb = &b;
    vec_init(pb);
    vec_alloc(pb, n);
    
    for (DSDP_INT k = 0; k < n; ++k) {
        retcode = spsMatScatter(sMat, pb, k); checkCode;
        memcpy(&(fulldMat[k * n]), pb->x, sizeof(double) * n);
    }
    
    vec_free(pb);
    return retcode;
}

extern DSDP_INT spsMatReset( spsMat *sMat ) {
    
    // Reset a sparse matrix to be 0
    DSDP_INT retcode = DSDP_RETCODE_OK;
    assert( sMat->dim > 0 );
    memset(sMat->x, 0, sizeof(double) * sMat->nnz);
    return retcode;
}

extern DSDP_INT spsMatView( spsMat *sMat ) {
    
    // View a sparse matrix by calling CXSparse cs_print
    DSDP_INT retcode = DSDP_RETCODE_OK;
    assert( sMat->dim > 0 );
    
    cs mat;
    mat.p = sMat->p;
    mat.i = sMat->i;
    mat.x = sMat->x;
    mat.nz = -1;
    mat.nzmax = sMat->nnz;
    mat.m = sMat->dim;
    mat.n = sMat->dim;
    cs_print(&mat, FALSE);
    
    return retcode;
}
