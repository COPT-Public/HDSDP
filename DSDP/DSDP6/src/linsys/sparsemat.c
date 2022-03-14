#include "sparsemat.h"
#include "densemat.h"
#include "dsdpfeast.h"
#include "dsdpdata.h"
#include "dsdplanczos.h"

// Enable hash sum check
#ifdef VERIFY_HASH
#undef VERIFY_HASH
#endif

static char *etype = "Sparse matrix";
static DSDP_INT one = 1;
static double dzero = 0.0;
static char uplolow = 'L';

/* Internal Pardiso Wrapper */
static DSDP_INT pardisoSymFactorize( spsMat *S ) {
    /* Factorize the spsMat matrix */
    DSDP_INT retcode = DSDP_RETCODE_OK;
    // Extract the lower triangular part from CSparse structure
    DSDP_INT *Sp = S->p, *Si = S->i;
    double   *Sx = S->x;
    // Get the pardiso parameter
    DSDP_INT phase = PARDISO_SYM, error = 0, n = S->dim;
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
    
    DSDP_INT retcode = DSDP_RETCODE_OK; assert( S->isFactorized == TRUE );
    DSDP_INT *Sp = S->p, *Si = S->i, phase = PARDISO_FAC, error = 0, n = S->dim;
    double   *Sx = S->x;
    
    // Invoke pardiso to do symbolic analysis and Cholesky factorization
    pardiso(S->pdsWorker, &maxfct, &mnum, &mtype, &phase, &n,
            Sx, Sp, Si, &idummy, &idummy, PARDISO_PARAMS_CHOLESKY,
            &msglvl, NULL, NULL, &error);
    
    if (error == -4) {
        double eig = 0.0; spsMatMinEig(S, &eig);
        printf("Minimum Eigenvalue: %g \n", eig);
    }
    
    assert( error == PARDISO_OK );
    
    if (error) {
        printf("[Pardiso Error]: Matrix factorization failed."
               " Error code: "ID" \n", error);
        error(etype, "Pardiso failes to factorize. \n");
    }
    
    return retcode;
}

static DSDP_INT pardisoForwardSolve( spsMat *S, DSDP_INT nrhs, double *B, double *aux, DSDP_INT overwrite ) {
    // pardiso forward solve
    DSDP_INT retcode = DSDP_RETCODE_OK; assert( S->isFactorized == TRUE );
    DSDP_INT *Sp = S->p, *Si = S->i;
    DSDP_INT phase = PARDISO_FORWARD, error = 0, n = S->dim, *param;
    param = (overwrite) ? PARDISO_PARAMS_FORWARD_BACKWORD : PARDISO_PARAMS_FORWARD_BACKWORD_LANCZOS;
    pardiso(S->pdsWorker, &maxfct, &mnum, &mtype, &phase, &n,
            NULL, Sp, Si, &idummy, &nrhs, param, &msglvl,
            B, aux, &error); assert( error == PARDISO_OK );
    if (error) {
        printf("[Pardiso Error]: Matrix backward solve failed."
               " Error code: "ID" \n", error);
        error(etype, "Pardiso failes to factorize. \n");
    }
    return retcode;
}

static DSDP_INT pardisoBackwardSolve( spsMat *S, DSDP_INT nrhs, double *B, double *aux, DSDP_INT overwrite ) {
    
    // pardiso forward solve
    DSDP_INT retcode = DSDP_RETCODE_OK; assert( S->isFactorized == TRUE );
    DSDP_INT *Sp = S->p, *Si = S->i;
    DSDP_INT phase = PARDISO_BACKWARD, error = 0, n = S->dim, *param;
    param = (overwrite) ? PARDISO_PARAMS_FORWARD_BACKWORD : PARDISO_PARAMS_FORWARD_BACKWORD_LANCZOS;
    pardiso(S->pdsWorker, &maxfct, &mnum, &mtype, &phase, &n,
            NULL, Sp, Si, &idummy, &nrhs, param, &msglvl,
            B, aux, &error); assert( error == PARDISO_OK );
    if (error) {
        printf("[Pardiso Error]: Matrix backward solve failed."
               " Error code: "ID" \n", error);
        error(etype, "Pardiso failes to factorize. \n");
    }
    return retcode;
}

static DSDP_INT pardisoSolve( spsMat *S, DSDP_INT nrhs, double *B, double *X ) {
    
    /* Solve the linear system S * X = B using Pardiso */
    DSDP_INT retcode = DSDP_RETCODE_OK;
    assert( S->isFactorized == TRUE );
    DSDP_INT phase = PARDISO_SOLVE, error = 0, n = S->dim;
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

static DSDP_INT pardisoSolveInplace( spsMat *S, DSDP_INT nrhs, double *B, double *aux ) {
    
    /* Solve the linear system S * X = B using Pardiso */
    DSDP_INT retcode = DSDP_RETCODE_OK;
    assert( S->isFactorized == TRUE );
    DSDP_INT phase = PARDISO_SOLVE, error = 0, n = S->dim;
    // Invoke pardiso to perform solution
    pardiso(S->pdsWorker, &maxfct, &mnum, &mtype, &phase, &n,
            S->x, S->p, S->i, &idummy,
            &nrhs, PARDISO_PARAMS_CHOLESKY_INPLACE, &msglvl,
            B, aux, &error);
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
    DSDP_INT phase = PARDISO_FREE, error = 0, n = S->dim;
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

static void arrayTranspose( double *A, DSDP_INT n ) {
    double tmp = 0.0;
    for (DSDP_INT i = 0; i < n; ++i) {
        for (DSDP_INT j = i + 1; j < n; ++j) {
            tmp = A[i * n + j];
            A[i * n + j] = A[j * n + i];
            A[j * n + i] = tmp;
        }
    }
}

/* Structure operations */
extern DSDP_INT spsMatInit( spsMat *sMat ) {
    
    // Initialize sparse matrix structure
    DSDP_INT retcode   = DSDP_RETCODE_OK;
    sMat->dim = 0; sMat->factor = NULL;
    sMat->p = NULL; sMat->i = NULL; sMat->x = NULL;
    sMat->isFactorized = FALSE; sMat->nzHash = NULL;
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
    assert( sMat->dim == 0 ); assert( nnz <= nsym(dim));
    
    sMat->dim = dim; sMat->nnz = nnz;
    
    // Currently do nothing here. TODO: Do Optimization if S is dense
    if (nnz < nsym(dim)) {
        sMat->p = (DSDP_INT *) calloc(dim + 1, sizeof(DSDP_INT));
        sMat->i = (DSDP_INT *) calloc(nnz, sizeof(DSDP_INT));
    } else {
        sMat->p = (DSDP_INT *) calloc(dim + 1, sizeof(DSDP_INT));
        sMat->i = (DSDP_INT *) calloc(nsym(dim), sizeof(DSDP_INT));
    }
    
    sMat->x = (double *) calloc(nnz, sizeof(double));
    sMat->nzHash = (DSDP_INT *) calloc(nnz, sizeof(DSDP_INT));
    
    return retcode;
}

extern DSDP_INT spsMatFree( spsMat *sMat ) {
    
    // Free memory allocated
    DSDP_INT retcode = DSDP_RETCODE_OK;
    sMat->dim = 0; sMat->nnz = 0;
    
    // Note that we first free p, i and x before calling pardiso to destroy the working array
    DSDP_FREE(sMat->p); DSDP_FREE(sMat->i); DSDP_FREE(sMat->x);
    DSDP_FREE(sMat->nzHash);
    
    if (sMat->isFactorized) {
        retcode = pardisoFree(sMat); sMat->isFactorized = FALSE;
    }
    if (sMat->factor) {
        rkMatFree(sMat->factor); DSDP_FREE(sMat->factor);
    }
        
    return retcode;
}

/* Basic operations */
extern DSDP_INT spsMatAx( spsMat *A, vec *x, vec *Ax ) {
    /* Sparse matrix multiplication for Lanczos Method
      Warning: This method is only for Lanczos since we compute Ax = - A * x
     */
    DSDP_INT retcode = DSDP_RETCODE_OK;
    
    assert( A->dim == x->dim && x->dim == Ax->dim );
    DSDP_INT *cidx = A->nzHash, *Ap = A->p, *Ai = A->i;
    DSDP_INT n = A->dim, nnz = A->nnz, idx;
    double *Axdata = A->x, *Axres = Ax->x, *xdata = x->x, coeff = -1.0;
    vec_reset(Ax);
    
    if (nnz == nsym(n)) {
        dspmv(&uplolow, &n, &coeff, Axdata, xdata, &one, &dzero, Axres, &one);
    } else if (nnz <= 0.8 * n * n && cidx) {
        for (DSDP_INT i, j, k = 0; k < nnz; ++k) {
            i = Ai[k]; j = cidx[k]; Axres[i] -= Axdata[k] * xdata[j];
            if (i != j) { Axres[j] -= Axdata[k] * xdata[i]; }
        }
    } else {
        for (DSDP_INT i = 0; i < n ; ++i) {
            coeff = xdata[i]; if (coeff == 0.0) { continue; }
            idx = Ap[i];
            
            if (Ai[idx] == i) {
                Axres[i] -= coeff * Axdata[idx];
            } else {
                Axres[Ai[idx]] -= coeff * Axdata[idx];
                Axres[i] -= xdata[Ai[idx]] * Axdata[idx];
            }
            
            for (DSDP_INT j = Ap[i] + 1; j < Ap[i + 1]; ++j) {
                idx = Ai[j];
                Axres[idx] -= coeff * Axdata[j];
                Axres[i] -= xdata[idx] * Axdata[j];
            }
        }
    }
    
    return retcode;
}

extern double spsMatxTAx( spsMat *A, double *x ) {
    // Compute quadratic form x' * A * x
    register DSDP_INT *Ai = A->i, *nzHash = A->nzHash, i, j, k;
    register double res = 0.0, tmp = 0.0, *Ax = A->x;
    for (k = 0; k < A->nnz; ++k) {
        i = nzHash[k]; j = Ai[k]; tmp = Ax[k] * x[i] * x[j];
        res += (i == j) ? (0.5 * tmp) : tmp;
    }
    return (2.0 * res);
}

extern DSDP_INT spsMataXpbY( double alpha, spsMat *sXMat, double beta,
                             spsMat *sYMat, DSDP_INT *sumHash ) {
    
    // Matrix axpy operation: let sYMat = alpha * sXMat + beta * sYMat
    // Note that sYMat must have a hash table
    DSDP_INT retcode = DSDP_RETCODE_OK;
    if (beta != 1.0) { retcode = spsMatScale(sYMat, beta); }
    if (alpha == 0.0) { return retcode; }
    DSDP_INT dim = sXMat->dim, *Ap = sXMat->p, *Ai = sXMat->i;
    double   *Ax = sXMat->x, *Bx = sYMat->x;
    if (sumHash) {
        for (DSDP_INT i = 0, j, hash; i < dim; ++i) {
            for (j = Ap[i]; j < Ap[i + 1]; ++j) {
#ifdef VERIFY_HASH
                hash = packIdx(sumHash, dim, Ai[j], i);
                assert( hash > 0 || j == 0 );
                Bx[hash] += alpha * Ax[j];
#else
                Bx[packIdx(sumHash, dim, Ai[j], i)] += alpha * Ax[j];
#endif
            }
        }
    } else {
        // sYMat is actually dense and Bx is an (n + 1) * n / 2 array
        assert( sYMat->nnz == nsym(dim) );
        // If adding sparse data e.g. S += A * y[i]
        if (sXMat->nzHash) {
            for (DSDP_INT i, j, k = 0; k < sXMat->nnz; ++k) {
                i = Ai[k]; j = sXMat->nzHash[k];
                packIdx(Bx, dim, i, j) += alpha * Ax[k];
            }
        // If adding iteration e.g. S += dS
        } else {
            assert( sXMat->nnz == sYMat->nnz );
            daxpy(&sXMat->nnz, &alpha, Ax, &one, Bx, &one);
        }
    }

    return retcode;
}

extern DSDP_INT spsMatAdddiag( spsMat *sMat, double d, DSDP_INT *sumHash ) {
    
    // Add a diagonal element to a sparse matrix with hash table
    DSDP_INT retcode = DSDP_RETCODE_OK;
    if (d == 0.0) { return retcode; }
    DSDP_INT dim = sMat->dim, idx = 0;
    
    if (sumHash) {
        for (DSDP_INT i = 0; i < dim; ++i) {
#ifdef VERIFY_HASH
            idx = packIdx(sumHash, dim, i, i);
            assert( idx > 0 || i == 0 );
            sMat->x[idx] += d;
#else
            sMat->x[packIdx(sumHash, dim, i, i)] += d;
#endif
        }
    } else {
        for (DSDP_INT i = 0; i < dim; ++i) {
            sMat->x[idx] += d; idx += dim - i;
        }
    }
    
    return retcode;
}

extern DSDP_INT spsMatAddds( spsMat *sXMat, double alpha, dsMat *dsYMat ) {
    
    // Add a dense matrix to a nominally sparse matrix
    DSDP_INT retcode = DSDP_RETCODE_OK;
    assert( sXMat->nnz == nsym(sXMat->dim) );
    if (alpha == 0.0) { return retcode; }
    DSDP_INT dim = nsym(sXMat->dim);
    axpy(&dim, &alpha, dsYMat->array, &one, sXMat->x, &one);
    
    return retcode;
}

extern DSDP_INT spsMatAddr1( spsMat *sXMat, double alpha, r1Mat *r1YMat, DSDP_INT *sumHash ) {
    
    // Add a rank 1 matrix to a sparse matrix
    DSDP_INT retcode = DSDP_RETCODE_OK;
    if (alpha == 0.0) { return retcode; }
    DSDP_INT dim = sXMat->dim, idx = 0, *hash = sumHash, *nzIdx = r1YMat->nzIdx;
    double sign = r1YMat->sign * alpha, *rx = r1YMat->x;
    
    if (sumHash) {
        for (DSDP_INT i = 0; i < r1YMat->nnz; ++i) {
            for (DSDP_INT j = 0; j <= i; ++j) {
#ifdef VERIFY_HASH
                idx = packIdx(hash, dim, nzIdx[i], nzIdx[j]);
                assert( idx > 0 || i == 0 );
                sXMat->x[idx] += sign * rx[nzIdx[i]] * rx[nzIdx[j]];
#else
                sXMat->x[packIdx(hash, dim, nzIdx[i], nzIdx[j])] +=
                sign * rx[nzIdx[i]] * rx[nzIdx[j]];
#endif
            }
        }
    } else {
        DSDP_INT i, j;
        if (r1YMat->nnz < 0.5 * dim) {
            for (i = 0; i < r1YMat->nnz; ++i) {
                for (j = 0; j <= i; ++j) {
                    packIdx(sXMat->x, dim, nzIdx[i], nzIdx[j]) += \
                    sign * rx[nzIdx[i]] * rx[nzIdx[j]];
                }
            }
        } else {
            char uplo = DSDP_MAT_LOW;
            packr1update(&uplo, &dim, &sign, r1YMat->x, &one, sXMat->x);
        }
    }
    
    return retcode;
}

extern DSDP_INT spsMatAddrk( spsMat *sXMat, double alpha, rkMat *rkYMat, DSDP_INT *sumHash ) {

    // Add a rank k matrix to a sparse matrix
    DSDP_INT retcode = DSDP_RETCODE_OK;
    if (alpha == 0.0) { return retcode; }
    for (DSDP_INT i = 0; i < rkYMat->rank; ++i) {
        spsMatAddr1(sXMat, alpha, rkYMat->data[i], sumHash);
    }
    return retcode;
}

extern DSDP_INT spsMatScale( spsMat *sXMat, double alpha ) {
    // Scale a sparse matrix by some number.
    vecscal(&sXMat->nnz, &alpha, sXMat->x, &one);
    return DSDP_RETCODE_OK;
}

extern DSDP_INT spsMatRscale( spsMat *sXMat, double r ) {
    // Scale a sparse matrix by the reciprocical of some number. No overflow or under flow
    
    DSDP_INT retcode = DSDP_RETCODE_OK;
    assert( sXMat->dim );
    if (fabs(r) < 1e-15) {
        error(etype, "Dividing a matrix by 0. \n");
    }
    if (sXMat->factor) { rkMatRscale(sXMat->factor, r); }
    vecdiv(&sXMat->nnz, &r, sXMat->x, &one);
    
    return retcode;
}

extern DSDP_INT spsMatFnorm( spsMat *sMat, double *fnrm ) {
    
    // Matrix Fronenius norm
    DSDP_INT retcode = DSDP_RETCODE_OK;
    assert( sMat->dim > 0 );
    
    DSDP_INT i, n = sMat->dim;
    double nrm = 0.0; i = spsMatGetRank(sMat);
    if (i < 0.1 * n) {
        rkMatFnorm(sMat->factor, fnrm); return retcode;
    }
    
    if (sMat->nzHash) {
        DSDP_INT j, k;
        for (i = 0; i < sMat->nnz; ++i) {
            j = sMat->i[i]; k = sMat->nzHash[i];
            if (j == k) {
                nrm += sMat->x[i] * sMat->x[i];
            } else {
                nrm += 2 * sMat->x[i] * sMat->x[i];
            }
        }
        *fnrm = sqrt(nrm);
    } else if (sMat->nnz == nsym(n)) {
        char ntype = 'F', low = DSDP_MAT_LOW;
        *fnrm = dlansp(&ntype, &low, &n, sMat->x, NULL);
    }
    
    return retcode;
}

extern DSDP_INT spsMatOneNorm( spsMat *sMat, double *onenrm ) {
    
    // Element-wise sum of absolute values
    DSDP_INT retcode = DSDP_RETCODE_OK;
    
    double nrm = 0.0;
    DSDP_INT i, j, k;
    if (sMat->nzHash) {
        for (i = 0; i < sMat->nnz; ++i) {
            j = sMat->i[i]; k = sMat->nzHash[i];
            nrm += (j == k) ? 0.5 * fabs(sMat->x[i]) : fabs(sMat->x[i]);
        }
        nrm *= 2;
    } else {
        assert( sMat->nnz == nsym(sMat->dim) );
        for (i = 0; i < sMat->nnz; ++i) { nrm += fabs(sMat->x[i]); }
        j = 0; k = sMat->dim;
        for (i = 0; i < k; ++i) {
            nrm -= 0.5 * fabs(sMat->x[j]); j += k - i;
        }
        nrm *= 2;
    }
    *onenrm = nrm; return retcode;
}

/* Factorization and linear system solver */
extern DSDP_INT spsMatSymbolic( spsMat *sAMat ) {
    // Symbolic analysis
    return pardisoSymFactorize(sAMat);
}

extern DSDP_INT spsMatFactorize( spsMat *sAMat ) {
    // Sparse matrix Cholesky decomposition
    return pardisoNumFactorize(sAMat);
}

extern DSDP_INT spsMatVecSolve( spsMat *sAMat, vec *sbVec, double *Ainvb ) {
    // Sparse matrix operation X = A \ b
    return pardisoSolve(sAMat, 1, sbVec->x, Ainvb);;
}

extern DSDP_INT spsMatVecFSolve( spsMat *sAmat, vec *sbVec, vec *Ainvb ) {
    /* Forward solve L * x = b */
    return pardisoForwardSolve(sAmat, 1, sbVec->x, Ainvb->x, FALSE);
}

extern DSDP_INT spsMatVecBSolve( spsMat *sAmat, vec *sbVec, vec *Ainvb ) {
    /* Backward solve L' * x = b */
    return pardisoBackwardSolve(sAmat, 1, sbVec->x, Ainvb->x, FALSE);
}

extern DSDP_INT spsMatLspLSolve( spsMat *S, spsMat *dS, spsMat *spaux ) {
    // Routine for computing SDP cone maximum stepsize
    // TODO: Accelerate the routine using Sinv (if available)
    DSDP_INT retcode = DSDP_RETCODE_OK;
    DSDP_INT n = S->dim, *Ap = spaux->p, *Ai = spaux->i, i, j;
    double *fulldS = (double *) calloc(n * n, sizeof(double));
    double *aux    = (double *) calloc(n * n, sizeof(double));
    double *Ax     = spaux->x, tmp = 0.0;
    
    // L^-1 dS
    retcode = spsMatFill(dS, fulldS);
    retcode = pardisoForwardSolve(S, n, fulldS, aux, TRUE);

    for (i = 0; i < n; ++i) {
        for (j = 0; j <= i; ++j) {
            tmp = fulldS[i * n + j];
            fulldS[i * n + j] = fulldS[j * n + i];
            fulldS[j * n + i] = tmp;
        }
    }
    // L^-T (L^-1 dS)
    retcode = pardisoForwardSolve(S, n, fulldS, aux, TRUE);
    DSDP_INT nnz = 0, start = 0;
    for (i = 0; i < n; ++i) {
        start = i * n;
        for (j = i; j < n; ++j) {
            tmp = fulldS[start + j];
            Ax[nnz] = tmp; Ai[nnz] = j; ++nnz;
        }
        Ap[i + 1] = nnz;
    }
    
    spaux->nnz = nnz; DSDP_FREE(fulldS); DSDP_FREE(aux);
    return retcode;
}

extern DSDP_INT spsMatGetX( spsMat *S, spsMat *dS, double *LinvSLTinv ) {
    // Routine for retrieving X
    DSDP_INT retcode = DSDP_RETCODE_OK;
    
    DSDP_INT n = S->dim, i, j;
    assert( dS->dim == n );
    double *fulldS = LinvSLTinv, *aux = (double *) calloc(n * n, sizeof(double)), tmp = 0.0;
    spsMatFill(dS, fulldS); pardisoForwardSolve(S, n, fulldS, aux, TRUE);
    // Transpose
    arrayTranspose(fulldS, n); pardisoForwardSolve(S, n, fulldS, aux, TRUE);
    // I + L^-1 * dS * LT^-1
    for (i = 0; i < n; ++i) {
        fulldS[i * n + i] += 1.0;
        for (j = i + 1; j < n; ++j) {
            tmp = (fulldS[i * n + j] + fulldS[j * n + i]) * 0.5;
            fulldS[i * n + j] = fulldS[j * n + i] = tmp;
        }
    }
    pardisoBackwardSolve(S, n, fulldS, aux, TRUE);
    arrayTranspose(fulldS, n);
    pardisoBackwardSolve(S, n, fulldS, aux, TRUE);
    
    // Fix numerical instability
    for (i = 0; i < n; ++i) {
        for (j = i + 1; j < n; ++j) {
            tmp = (fulldS[i * n + j] + fulldS[j * n + i]) * 0.5;
            fulldS[i * n + j] = fulldS[j * n + i] = tmp;
        }
    }
    DSDP_FREE(aux); return retcode;
}

/* DSDP routine for computing the stepsize in the SDP cone */
extern DSDP_INT dsdpGetAlpha( spsMat *S, spsMat *dS, spsMat *spaux, double *alpha ) {
    // Get the maximum alpha such that S + alpha * dS is PSD
    DSDP_INT retcode = DSDP_RETCODE_OK;
    DSDP_INT n = S->dim;
    double mineig = 0.0, lbd = 0.0, delta = 0.0;
    
    // Check block size of 1
    if (n == 1) {
        *alpha = (dS->x[0] > 0) ? DSDP_INFINITY : -S->x[0] / dS->x[0];
        return retcode;
    }
    // Lanczos iteration
//    retcode = dsdpLanczos(S, dS, &lbd, &delta);
    
    if (1 || lbd != lbd || delta != delta || (lbd + delta > 1e+10)) {
        // MKL extremal routine
        retcode = spsMatLspLSolve(S, dS, spaux); checkCode;
        retcode = spsMatMinEig(spaux, &mineig); checkCode;
        *alpha = (mineig >= 0) ? DSDP_INFINITY : - 1.0 / mineig;
    } else {
        *alpha = (lbd + delta < 0) ? DSDP_INFINITY : 1.0 / (lbd + delta);
    }
    return retcode;
}

extern DSDP_INT dsdpGetAlphaLS( spsMat *S, spsMat *dS, spsMat *Scker,
                                double alphamax, double *alpha, DSDP_INT *sumHash ) {
    // Get the maximum alpha such that S + alpha * dS is PSD by line-search
    DSDP_INT retcode = DSDP_RETCODE_OK;
    
    double step = 1 / alphamax, *src = NULL;
    DSDP_INT ispsd = FALSE;
    spsMat *buffer = NULL;
    
    if (Scker) {
        src = S->x; buffer = Scker;
    } else {
        assert( FALSE );
        src = (double *) calloc(S->nnz, sizeof(double));
        memcpy(src, S->x, sizeof(double) * S->nnz); buffer = S;
    }
    
    for (DSDP_INT i = 0; ; ++i) {
        if (step <= 1e-05) { *alpha = 0.0; break; }
        memcpy(src, S->x, sizeof(double) * S->nnz);
        spsMataXpbY(step, dS, 1.0, buffer, sumHash); spsMatIspd(buffer, &ispsd);
        if (ispsd) { *alpha = step; break; }
        step *= 0.8;
    }

    if (!Scker) { DSDP_FREE(src); }
    return retcode;
}

/* M3 Technique */
extern double SinvSpSinv( const double *Sinv, double *aux, spsMat *A, dsMat *SinvASinv ) {
    
    // Routine for setting up the Schur matrix
    /*
        Compute inv(S) * A * inv(S) for sparse A
    */
    
    DSDP_INT n = A->dim, *Ai = A->i, *nzHash = A->nzHash, i, j, k;
    double *SinvASinvdata = SinvASinv->array;
    double *SinvA = aux, *Ax = A->x, *SinvAj, coeff = 0.0, res = 0.0;
    const double *Sinvi; memset(SinvA, 0, sizeof(double) * n * n);
    
    // n f_sigma
    for (k = 0; k < A->nnz; ++k) {
        i = Ai[k]; j = nzHash[k]; coeff = Ax[k];
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
    pardisoSolve(S, 1, xA, xSinvASinv);
    
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
    DSDP_INT n = A->dim, *Ai = A->i, *nzHash = A->nzHash, i, j, k;
    double *Ax = A->x, *ASinvi, coeff = 0.0, res = 0.0, res2 = 0.0;
    const double *Sinvj; memset(ASinv, 0, sizeof(double) * n * n);
    
    for (k = 0; k < A->nnz; ++k) {
        i = Ai[k]; j = nzHash[k]; coeff = Ax[k];
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
    DSDP_INT n = A->dim, *Ai = A->i, *nzHash = A->nzHash, i, j, k;
    double *Ax = A->x, coeff = 0.0, res = 0.0;
    const double *Sinvi, *ASinvj;
    
    for (k = 0; k < A->nnz; ++k) {
        i = Ai[k]; j = nzHash[k]; coeff = Ax[k];
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
    DSDP_INT i, p, q, in, n = A->dim, *Ai = A->i, *Aj = A->nzHash;
    
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
    DSDP_INT *Ai = A->i, *Aj = A->nzHash, *Bi = B->i, *Bj = B->nzHash;
    
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
    DSDP_INT *Ai = A->i, *Aj = A->nzHash, *Bi = B->nzIdx;
    
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

/* Eigen value routines */
extern DSDP_INT spsMatMaxEig( spsMat *sMat, double *maxEig ) {
    // Eigen value utility: compute the maximum eigenvalue of a matrix
    DSDP_INT retcode = DSDP_RETCODE_OK;
    DSDP_INT n = sMat->dim, info = 0;
    
    double *eigvec = (double *) calloc(n, sizeof(double));
    sparse_matrix_t A = NULL;
    mkl_sparse_d_create_csr(&A, SPARSE_INDEX_BASE_ZERO,
                            n, n, sMat->p, sMat->p + 1,
                            sMat->i, sMat->x);
    info = mkl_sparse_d_ev(&MAX_EIG, pm, A, dsdp_descr,
                           k0, &k, maxEig, eigvec, &resi);
    if (info != SPARSE_STATUS_SUCCESS) {
        error(etype, "Maximum eigen value computation failed. \n");
    }
    mkl_sparse_destroy( A ); // eigvec is freed here
    return retcode;
}

extern DSDP_INT spsMatMinEig( spsMat *sMat, double *minEig ) {
    // Eigen value utility: compute the minimum eigenvalue of a matrix
    DSDP_INT retcode = DSDP_RETCODE_OK;
    DSDP_INT n = sMat->dim, info = 0;
    double *eigvec = (double *) calloc(n, sizeof(double));
    sparse_matrix_t A = NULL;
    mkl_sparse_d_create_csr(&A, SPARSE_INDEX_BASE_ZERO,
                            n, n, sMat->p, sMat->p + 1,
                            sMat->i, sMat->x);
    info = mkl_sparse_d_ev(&MIN_EIG, pm, A, dsdp_descr,
                           k0, &k, minEig, eigvec, &resi);
    if (info != SPARSE_STATUS_SUCCESS) {
        error(etype, "Minimum eigen value computation failed. \n");
    }
    mkl_sparse_destroy( A );
    return retcode;
}

/* Other utilities */
extern DSDP_INT spsMatIspd( spsMat *sMat, DSDP_INT *ispd ) {
    // A critical routine that determines whether a matrix is positive definite
    DSDP_INT retcode = DSDP_RETCODE_OK;
    DSDP_INT error = 0, phase = PARDISO_FAC;
    // Invoke pardiso to do symbolic analysis and Cholesky factorization
    pardiso(sMat->pdsWorker, &maxfct, &mnum, &mtype, &phase, &sMat->dim,
            sMat->x, sMat->p, sMat->i, &idummy, &idummy, PARDISO_PARAMS_PSD_CHECK,
            &msglvl, NULL, NULL, &error);
        
    if (error == 0) {
        sMat->isFactorized = TRUE;
        *ispd = TRUE;
    } else if (error == -4) {
        *ispd = FALSE;
    } else {
        error(etype, "Pardiso fails for some reason. \n");
    }
    
    return retcode;
}

extern double spsMatGetlogdet( spsMat *sMat, double *aux ) {
    // Compute log det S
    DSDP_INT ecode = PARDISO_OK;
    double res = 0.0;
    pardiso_getdiag((const void **) sMat->pdsWorker,
                    aux, &aux[sMat->dim], &one, &ecode);
    for (DSDP_INT i = 0; i < sMat->dim; ++i) { res += log(aux[i]); }
    return res;
}

extern void spsMatInverse( spsMat *sMat, double *Sinv, double *aux ) {
    pardisoSolveInplace(sMat, sMat->dim, Sinv, aux);
    /* Numerical fix*/
    if (FALSE) {
        double tmp;
        for (DSDP_INT i = 0, n = sMat->dim, j; i < n; ++i) {
            for (j = 0; j < n; ++j) {
                tmp = 0.5 * (Sinv[i * n + j] + Sinv[j * n + i]);
                Sinv[j * n + i] = Sinv[i * n + j] = tmp;
            }
        }
    }
    return;
}

extern DSDP_INT spsMatStoreFactor( spsMat *sMat, rkMat *factor ) {
    // Save factorized data
    DSDP_INT retcode = DSDP_RETCODE_OK;
    sMat->factor = factor;
    return retcode;
}

extern rkMat* spsMatGetFactor( spsMat *sMat ) {
    return sMat->factor;
}

extern DSDP_INT spsMatGetRank( spsMat *sMat ) {
    return (sMat->factor) ? sMat->factor->rank : (sMat->dim) * 1000;
}

extern DSDP_INT spsMatFillLower( spsMat *sMat, double *lowFullMat ) {
    
    DSDP_INT retcode = DSDP_RETCODE_OK;
    DSDP_INT n = sMat->dim, ni = 0, *Ap = sMat->p, *Ai = sMat->i;
    double *Ax = sMat->x;
    for (DSDP_INT i = 0; i < n; ++i) {
        ni = n * i;
        for (DSDP_INT j = Ap[i]; j < Ap[i + 1]; ++j) {
            lowFullMat[ni + Ai[j]] = Ax[j];
        }
    }
    return retcode;
}

extern DSDP_INT spsMatFill( spsMat *sMat, double *fulldMat ) {
    
    // Fill sparse matrix to full (there is no structure for symmetric full dense matrix)
    DSDP_INT retcode = DSDP_RETCODE_OK;
    DSDP_INT n = sMat->dim, *Ap = sMat->p, *Ai = sMat->i, i, j;
    double *Ax = sMat->x;
    for (i = 0; i < n; ++i) {
        for (j = Ap[i]; j < Ap[i + 1]; ++j) {
            fulldMat[i * n + Ai[j]] = Ax[j];
            fulldMat[Ai[j] * n + i] = Ax[j];
        }
    }
    return retcode;
}

extern DSDP_INT spsMatReset( spsMat *sMat ) {
    DSDP_INT retcode = DSDP_RETCODE_OK;
    assert( sMat->dim > 0 && !sMat->factor ); // Never reset an (eigen) factorized sparse matrix
    memset(sMat->x, 0, sizeof(double) * sMat->nnz);
    return retcode;
}

extern DSDP_INT spsMatView( spsMat *sMat ) {
    
    // View a sparse matrix by calling CXSparse cs_print
    DSDP_INT retcode = DSDP_RETCODE_OK;
    cs mat; mat.p = sMat->p; mat.i = sMat->i; mat.x = sMat->x;
    mat.nz = -1; mat.nzmax = sMat->nnz; mat.m = sMat->dim;
    mat.n = sMat->dim; cs_print(&mat, FALSE);
    return retcode;
}

extern void spsMatLinvView( spsMat *S ) {
    // Lanczos debugging routine. Print P^-1 L^-1 to the screen
    assert( S->isFactorized );
    DSDP_INT n = S->dim, i, j;
    
    double *eye = (double *) calloc(n * n, sizeof(double));
    double *Linv = (double *) calloc(n * n, sizeof(double));
    for (i = 0; i < n; ++i) { eye[n * i + i] = 1.0; }
    pardisoForwardSolve(S, n, eye, Linv, FALSE);
    printf("Matrix view: \n");
    for (i = 0; i < n; ++i) {
        for (j = 0; j < n; ++j) {
            printf("%20.12e, ", Linv[i * n + j]);
        }
        printf("\n");
    }
    DSDP_FREE(eye); DSDP_FREE(Linv);
}
