#include "densemat.h"
#include "sparsemat.h"

// Error type
static char etype[] = "Dense Operation Error";
static DSDP_INT one = 1;

/* Internal Lapack Wrapper */
static DSDP_INT packFactorize( dsMat *S ) {
    
    /* Factorize the dsMat matrix */
    DSDP_INT retcode = DSDP_RETCODE_OK;
    
    assert( (!S->isFactorized) && (S->dim > 0));
    if (S->isFactorized) {
        error(etype, "Matrix is already factorized. \n");
    }
    
    DSDP_INT n    = S->dim;
    char uplo     = DSDP_MAT_LOW;
    DSDP_INT info = 0;
    
    memcpy(S->lfactor, S->array, sizeof(double) * nsym(n));
    packchol(&uplo, &n, S->lfactor, &info);
    
    if (info < 0) {
        error(etype, "Illegal value detected in packed dense format. \n");
    } else if (info > 0) {
        error(etype, "The matrix is non-PSD. \n");
    }
    
    S->isFactorized = TRUE;
    return retcode;
}

static DSDP_INT packSolve( dsMat *S, DSDP_INT nrhs, double *B, double *X ) {
    
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
    
    // Copy solution data
    memcpy(X, B, sizeof(double) * nrhs * n);
    
    dpptrs(&uplo, &n, &nrhs, S->lfactor, X, &ldb, &info);
    
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
    
    return retcode;
}

extern DSDP_INT denseMatAlloc( dsMat *dMat, DSDP_INT dim, DSDP_INT doFactor ) {
    
    // Allocate memory for dense matrix data
    DSDP_INT retcode = DSDP_RETCODE_OK;
    assert( dMat->dim == 0 );
    
    dMat->dim   = dim;
    dMat->array = (double *) calloc((DSDP_INT) nsym(dim), sizeof(double));
    
    if (doFactor) {
        dMat->lfactor = (double *) calloc((DSDP_INT) nsym(dim), sizeof(double));
    }
    
    return retcode;
}

extern DSDP_INT denseMatFree( dsMat *dMat ) {
    
    // Free memory allocated
    DSDP_INT retcode = DSDP_RETCODE_OK;
    
    if (dMat) {
        dMat->dim = 0;
        DSDP_FREE(dMat->array);
        DSDP_FREE(dMat->lfactor);
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

extern DSDP_INT denseMatxTAx( dsMat *dAMat, vec *x, double *xTAx ) {
    
    // Compute quadratic for m x' * A * x 
    DSDP_INT retcode = DSDP_RETCODE_OK;
    
    assert( dAMat->dim == x->dim );
    
    char uplo    = DSDP_MAT_LOW;
    DSDP_INT dim = dAMat->dim;
    double alpha = 1.0;
    double beta  = 0.0;
    
    double *y = (double *) calloc(dim, sizeof(double));
    packmatvec(&uplo, &dim, &alpha, dAMat->array, x->x,
               &one, &beta, y, &one);
    
    *xTAx = ddot(&dim, x->x, &one, y, &one);
    
    DSDP_FREE(y);
    return retcode;
}


extern DSDP_INT denseMatRscale( dsMat *dXMat, double r ) {
    // Scale a matrix by reciprocical without over/under flow
    
    DSDP_INT retcode = DSDP_RETCODE_OK;
    
    assert( dXMat->dim );
    double *array = dXMat->array;
    DSDP_INT n = nsym(dXMat->dim);
    
    vecdiv(&n, &r, array, &one);
    
    return retcode;
}

extern DSDP_INT denseMatFnorm( dsMat *dMat, double *fnrm ) {
    
    DSDP_INT retcode = DSDP_RETCODE_OK;
    assert(dMat->dim > 0);
    
    char nrm     = DSDP_MAT_FNORM;
    char uplo    = DSDP_MAT_LOW;
    double *work = NULL;
    
    *fnrm = fnorm(&nrm, &uplo, &dMat->dim, dMat->array, work);
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
    
    dpptrs(&uplo, &n, &nrhs, S->lfactor, B, &ldb, &info);
    
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
    
    // Compute trace (A * B) for dense A and sparse B. Used for trace(SinvASinv * A)
    DSDP_INT retcode = DSDP_RETCODE_OK;
    assert( dAMat->dim == sBMat->dim );
    assert((!dAMat->isFactorized) && (!sBMat->isFactorized));
    
    DSDP_INT n = dAMat->dim;
    
    // t stores trace, tmp stores the element-wise product
    double   t = 0.0;
    double tmp = 0.0;
    
    // Get sparse data structure
    double   *A  = dAMat->array;
    DSDP_INT *Bp = sBMat->p;
    DSDP_INT *Bi = sBMat->i;
    double   *Bx = sBMat->x;
    DSDP_INT   s = 0;
    
    for (DSDP_INT i = 0; i < n; ++i) {
        for (DSDP_INT k = Bp[i]; k < Bp[i + 1]; ++k) {
            s   = Bi[k];
            tmp = Bx[k] * A[s + (DSDP_INT) (i * (2 * n - i - 1) / 2)];
            if (Bi[k] == i) {
                tmp = tmp / 2;
            }
            t += tmp;
        }
    }
    
    *trace = 2.0 * t;
    return retcode;
}

extern DSDP_INT denseDsTrace( dsMat *dAMat, dsMat *dBMat, double *trace ) {
    
    // Compute trace (A * B) for dense A and dense B
    DSDP_INT retcode = DSDP_RETCODE_OK;
    assert( dAMat->dim == dBMat->dim );
    assert((!dAMat->isFactorized) && (!dBMat->isFactorized));
    
    DSDP_INT n = dAMat->dim;
    double *A = dAMat->array;
    double *B = dBMat->array;
    double res = 0.0;
    DSDP_INT nmat = nsym(n);
    
    // Take care of the diagonal elements
    res = ddot(&nmat, A, &one, B, &one);
    
    // Subtract half of diagonal off
    DSDP_INT idx = 0;
    for (DSDP_INT i = 0; i < n; ++i) {
        res -= 0.5 * A[idx] * B[idx];
        idx += n - i;
    }
    
    *trace = 2.0 * res;
    
    return retcode;
}

extern DSDP_INT denseDiagTrace( dsMat *dAMat, double diag, double *trace ) {
    
    // Compute trace( A * diag * I ) = diag * trace( A ). Used for Ry
    DSDP_INT retcode = DSDP_RETCODE_OK;
        
    DSDP_INT n = dAMat->dim;
    double *array = dAMat->array;
    double mattrace = 0.0;
    
    if (n > 64) {
        for (DSDP_INT i = 0; i < n - n % 4; i+=4) {
            mattrace += array[(DSDP_INT) (2 * n - i + 1) * (i    ) / 2];
            mattrace += array[(DSDP_INT) (2 * n - i    ) * (i + 1) / 2];
            mattrace += array[(DSDP_INT) (2 * n - i - 1) * (i + 2) / 2];
            mattrace += array[(DSDP_INT) (2 * n - i - 2) * (i + 3) / 2];
        }
        
        for (DSDP_INT i = n - n % 4; i < n; ++i) {
            mattrace += array[(DSDP_INT) (2 * n - i + 1) * i / 2];
        }
        
    } else {
        for (DSDP_INT i = 0; i < n; ++i) {
            mattrace += array[(DSDP_INT) (2 * n - i + 1) * i / 2];
        }
    }
    
    *trace = diag * mattrace;
    
    return retcode;
}

/* Utilities */
extern DSDP_INT denseMatScatter( dsMat *dMat, vec *b, DSDP_INT k ) {
    
    // Scatter dMat(:, k) into b
    DSDP_INT retcode = DSDP_RETCODE_OK;
    assert( dMat->dim == b->dim );
    
    // Index conversion
    DSDP_INT n   = dMat->dim;
    DSDP_INT idx = (DSDP_INT) (2 * n - k + 1) * k / 2;
    memcpy(&b->x[k], &dMat->array[idx], sizeof(double) * (n - k));
    idx = k;
    
    for (DSDP_INT i = 0; i < k; ++i) {
        b->x[i] = dMat->array[idx];
        idx += n - i - 1;
    }
    
    return retcode;
}

extern DSDP_INT denseMatFill( dsMat *dMat, double *fulldMat ) {
    
    // Fill packed matrix to full (there is no structure for symmetric full dense matrix)
    DSDP_INT retcode = DSDP_RETCODE_OK;
    DSDP_INT n       = dMat->dim;
    DSDP_INT idx     = 0;
    double *x        = NULL;
    
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
    double *x     = diag->x;
    double *array = dMat->array;
        
    if (n > 64) {
        for (DSDP_INT i = 0; i < n - n % 4; i+=4) {
            x[i    ] = array[(DSDP_INT) (2 * n - i + 1) * (i    ) / 2];
            x[i + 1] = array[(DSDP_INT) (2 * n - i    ) * (i + 1) / 2];
            x[i + 2] = array[(DSDP_INT) (2 * n - i - 1) * (i + 2) / 2];
            x[i + 3] = array[(DSDP_INT) (2 * n - i - 2) * (i + 3) / 2];
        }
        
        for (DSDP_INT i = n - n % 4; i < n; ++i) {
            x[i] = array[(DSDP_INT) (2 * n - i + 1) * i / 2];
        }
    } else {
        for (DSDP_INT i = 0; i < n; ++i) {
            x[i] = array[(DSDP_INT) (2 * n - i + 1) * i / 2];
        }
    }
    
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
