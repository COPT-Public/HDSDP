#include "sparsemat.h"
#include "densemat.h"

/* Internal Pardiso Wrapper */
static DSDP_INT pardisoFactorize( spsMat *S ) {
    
    /* Factorize the spsMat matrix */
    DSDP_INT retcode = DSDP_RETCODE_OK;
    
    // Extract the lower triangular part from CSparse structure
    DSDP_INT *Sp = S->cscMat->p;
    DSDP_INT *Si = S->cscMat->i;
    double   *Sx = S->cscMat->x;
    
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
        retcode = DSDP_RETCODE_FAILED;
        return retcode;
    }
    
    // Complete the factorization
    S->isFactorized = TRUE;
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
            S->cscMat->x, S->cscMat->p, S->cscMat->i, &idummy,
            &idummy, PARDISO_PARAMS_CHOLESKY, &msglvl,
            B, X, &error);
    
    assert( error == PARDISO_OK );
    
    if (error) {
        printf("[Pardiso Error]: Matrix solution failed."
               " Error code: "ID" \n", error);
        retcode = DSDP_RETCODE_FAILED;
        return retcode;
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
            S->cscMat->x, S->cscMat->p, S->cscMat->i, &idummy, &idummy,
            PARDISO_PARAMS_CHOLESKY, &msglvl, NULL, NULL, &error);
    
    assert( error == PARDISO_OK );
    
    if (error) {
        printf("[Pardiso Error]: Pardiso free failed."
               " Error code: "ID" \n", error);
        retcode = DSDP_RETCODE_FAILED;
        return retcode;
    }
    
    return retcode;
}

/* Structure operations */
extern DSDP_INT spsMatInit( spsMat *sMat ) {
    
    // Initialize sparse matrix structure
    DSDP_INT retcode   = DSDP_RETCODE_OK;
    sMat->dim          = 0;
    sMat->cscMat       = NULL;
    sMat->isFactorized = FALSE;
    memset(sMat->pdsWorker, 0, PARDISOINDEX * sizeof(void *));
    
    return retcode;
}

extern DSDP_INT spsMatAlloc( spsMat *sMat, DSDP_INT dim ) {
    
    // Allocate memory for sparse matrix data
    DSDP_INT retcode = DSDP_RETCODE_OK;
    assert( sMat->dim == 0 );
    
    sMat->dim = dim;
    sMat->cscMat = cs_spalloc(dim, dim, nsym(dim), TRUE, TRUE);
    
    return retcode;
}

extern DSDP_INT spsMatFree( spsMat *sMat ) {
    
    // Free memory allocated
    DSDP_INT retcode = DSDP_RETCODE_OK;
    sMat->dim = 0;
    
    // Note that we first free p, i and x before calling pardiso to destroy the working array
    cs_spfree(sMat->cscMat);
    pardisoFree(sMat);
    
    if (sMat->isFactorized) {
        sMat->isFactorized = FALSE;
    }
    
    return retcode;
}

/* Basic operations */
extern DSDP_INT spsMataXpbY( double alpha, spsMat *sXMat, double beta, spsMat *sYMat ) {
    
    // Matrix axpy operation: let sYMat = alpha * sXMat + sYMat
    DSDP_INT retcode = DSDP_RETCODE_OK;

    // Two matrices for adding are NEVER factorized since they purely serve as constant data
    assert ( sXMat->dim == sYMat->dim );
    assert ((!sXMat->isFactorized) && (!sYMat->isFactorized));
    if (sXMat->isFactorized || sYMat->isFactorized) {
        printf("[Sparse Operation Error]: Adding a factorized matrix. \n");
    }
    
    // Case if alpha = 0.0
    if (alpha == 0.0) {
        return retcode;
    }
    
    // Declare a place holder
    cs *aXpY     = NULL;
    DSDP_INT n   = sXMat->dim;
    cs *Xmat     = sXMat->cscMat;
    DSDP_INT nnz = sXMat->cscMat->nz;
    
    // Case if beta = 0.0
    if (beta == 0.0) {
        aXpY = cs_spalloc(n, n, sXMat->cscMat->nz, TRUE, FALSE);
        memcpy(aXpY->i, Xmat->i, sizeof(DSDP_INT) * nnz);
        memcpy(aXpY->p, Xmat->p, sizeof(DSDP_INT) * (n + 1));
        for (DSDP_INT i = 0; i < nnz; ++i) {
            aXpY->x[i] = alpha * aXpY->x[i];
        }
    } else {
        // General case with modified cs_add
        aXpY = cs_sadd(sXMat->cscMat, sYMat->cscMat, alpha, beta);
    }
    
    cs_spfree(sYMat->cscMat);
    sYMat->cscMat = aXpY;
    return retcode;
}

extern DSDP_INT spsMatFnorm( spsMat *sMat, double *fnorm ) {
    
    // Matrix Fronenius norm
    DSDP_INT retcode = DSDP_RETCODE_OK;
    assert( sMat->dim > 0 );
    DSDP_INT incx = 1;
    *fnorm = norm(&sMat->cscMat->nz, sMat->cscMat->x, &incx);
    
    return retcode;
}

/* Factorization and linear system solver */
extern DSDP_INT spsFactorize( spsMat *sAMat ) {
    
    // Sparse matrix Cholesky decomposition
    DSDP_INT retcode = DSDP_RETCODE_OK;
    
    // Call the internal method
    retcode = pardisoFactorize(sAMat);
    
    return retcode;
}

extern DSDP_INT spsVecSolve( spsMat *sAMat, vec *sbVec, double *Ainvb ) {
    
    // Sparse matrix operation X = A \ B = inv(A) * B for sparse A and B
    // A full dense matrix inv(A) * B is returned
    DSDP_INT retcode = DSDP_RETCODE_OK;
    // A is factorized and its size must agree with b
    assert( (sAMat->dim == sbVec->dim) );
    assert( (sAMat->isFactorized) );
    
    retcode = pardisoSolve(sAMat, 1, sbVec->x, Ainvb);
    
    return retcode;
}

extern DSDP_INT spsSpSolve( spsMat *sAMat, spsMat *sBMat, double *AinvB ) {
    
    // Sparse matrix operation X = A \ B = inv(A) * B for sparse A and B
    // A full dense matrix inv(A) * B is returned
    DSDP_INT retcode = DSDP_RETCODE_OK;
    
    // A is factorized, B is not and sizes must agree
    assert( (sAMat->dim == sAMat->dim) );
    assert( (sAMat->isFactorized) );
    
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
            retcode = spsVecSolve(sAMat, pb, &AinvB[k * n]);
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

extern DSDP_INT spsDsSolve( spsMat *sAMat, dsMat *sBMat, double *AinvB ) {
    
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
            retcode = spsVecSolve(sAMat, pb, &AinvB[k * n]);
        }
        
        vec_free(pb);
        
    } else {
        /* Violent parallel solve */
        double *B = NULL;
        B = (double *) calloc(n * n, sizeof(double));
        retcode = denseMatFill(sBMat, B);
        retcode = pardisoSolve(sAMat, n, B, AinvB);
        DSDP_FREE(B);
    }
    
    return retcode;
}

extern DSDP_INT spsR1Solve( spsMat *sAMat, r1Mat *sBMat, double *AinvB ) {
    
    // Matrix operation X = A \ B = inv(A) * B for sparse A and dense B
    // A full dense matrix inv(A) * B is returned
    DSDP_INT retcode = DSDP_RETCODE_OK;
    
    // A is factorized, B is not and sizes must agree
    assert( sAMat->dim == sBMat->dim );
    assert( sAMat->isFactorized );
    
    // TODO: Implement the sparse \ rank 1 operation
    
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
    retcode    = spsSpSolve(S, A, SinvA);
    double tmp = 0.0;
    
    // Transpose
    *asinv = 0.0;
    for (DSDP_INT i = 0; i < n; ++i) {
        for (DSDP_INT j = 0; j < n; ++j) {
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
    
    // TODO: Implement the SinvASinv operation for packed dense A
    
        
    return retcode;
}


/* Other utilities */
extern DSDP_INT spsMatScatter( spsMat *sMat, vec *b, DSDP_INT k ) {
    
    // Let b = sMat(:, k)
    DSDP_INT retcode = DSDP_RETCODE_OK;
    assert( sMat->dim == b->dim );
    
    DSDP_INT *Ap = sMat->cscMat->p;
    DSDP_INT *Ai = sMat->cscMat->i;
    double   *Ax = sMat->cscMat->x;
    
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
        spsMatScatter(sMat, pb, k);
        memcpy(&(fulldMat[k * n]), pb->x, sizeof(double) * n);
    }
    
    vec_free(pb);
    return retcode;
}

extern DSDP_INT spsMatView( spsMat *sMat ) {
    
    // View a sparse matrix by calling CXSparse cs_print
    DSDP_INT retcode = DSDP_RETCODE_OK;
    assert( sMat->dim > 0 );
    cs_print(sMat->cscMat, TRUE);
    
    return retcode;
}

