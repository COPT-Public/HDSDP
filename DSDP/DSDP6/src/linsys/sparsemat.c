#include "sparsemat.h"
#include "densemat.h"
#include "dsdpfeast.h"
#include "dsdpdata.h"
#include "dsdplanczos.h"

// Enable hash sum check
#ifndef VERIFY_HASH
#define VERIFY_HASH
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
    
    if (error == -4) {
        double eig = 0.0;
        spsMatMinEig(S, &eig);
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
    DSDP_INT *param;
    
    
    // Invoke pardiso to do symbolic analysis and Cholesky factorization
    if (overwrite) {
        param = PARDISO_PARAMS_FORWARD_BACKWORD;
    } else {
        param = PARDISO_PARAMS_FORWARD_BACKWORD_LANCZOS;
    }
    
    pardiso(S->pdsWorker, &maxfct, &mnum, &mtype, &phase, &n,
            NULL, Sp, Si, &idummy, &nrhs, param,
            &msglvl, B, aux, &error);
    
    assert( error == PARDISO_OK );
    
    if (error) {
        printf("[Pardiso Error]: Matrix forward solve failed."
               " Error code: "ID" \n", error);
        error(etype, "Pardiso failes to factorize. \n");
    }
    
    return retcode;
}

static DSDP_INT pardisoBackwardSolve( spsMat *S, DSDP_INT nrhs, double *B, double *aux, DSDP_INT overwrite ) {
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
    DSDP_INT *param;
    
    // Invoke pardiso to do symbolic analysis and Cholesky factorization
    if (overwrite) {
        param = PARDISO_PARAMS_FORWARD_BACKWORD;
    } else {
        param = PARDISO_PARAMS_FORWARD_BACKWORD_LANCZOS;
    }
    pardiso(S->pdsWorker, &maxfct, &mnum, &mtype, &phase, &n,
            NULL, Sp, Si, &idummy, &nrhs, param, &msglvl,
            B, aux, &error);
    
    assert( error == PARDISO_OK );
    
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

static DSDP_INT pardisoPartialSolve( spsMat *S, DSDP_INT *colNnz, double *xIn, double *xOut ) {
    
    /* Apply pardiso partial solve strategy */
    DSDP_INT retcode = DSDP_RETCODE_OK;
    
    assert( S->isFactorized == TRUE );
    DSDP_INT phase = PARDISO_SOLVE;
    DSDP_INT error = 0;
    DSDP_INT n     = S->dim;
    
    // Invoke pardiso to perform solution
    pardiso(S->pdsWorker, &maxfct, &mnum, &mtype, &phase, &n,
            S->x, S->p, S->i, colNnz, &one, PARDISO_PARAMS_PARTIAL_SOLVE,
            &msglvl, xIn, xOut, &error);
    
    if (error) {
        printf("[Pardiso Error]: Matrix partial solution failed."
               " Error code: "ID" \n", error);
        error(etype, "Pardiso failes to solve partial system. \n");
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
    
    sMat->dim          = 0;
    sMat->p            = NULL;
    sMat->i            = NULL;
    sMat->x            = NULL;
    sMat->isFactorized = FALSE;
    sMat->nzHash       = NULL;
    sMat->factor       = NULL;
    
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
    sMat->nzHash = (DSDP_INT *) calloc(nnz, sizeof(DSDP_INT));
    
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
    DSDP_FREE(sMat->nzHash);
    
    if (sMat->isFactorized) {
        retcode = pardisoFree(sMat);
        sMat->isFactorized = FALSE;
    }
    
    if (sMat->factor) {
        rkMatFree(sMat->factor);
        DSDP_FREE(sMat->factor);
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
    } else if (nnz <= 1.0 * n * n && cidx) {
        for (DSDP_INT i, j, k = 0; k < nnz; ++k) {
            i = Ai[k];
            j = cidx[k];
            Axres[i] -= Axdata[k] * xdata[j];
            if (i != j) {
                Axres[j] -= Axdata[k] * xdata[i];
            }
        }
    } else {
        for (DSDP_INT i = 0; i < n ; ++i) {
            coeff = xdata[i];
            
            if (coeff == 0.0) {
                continue;
            }
            
            idx = Ap[i];
            if (idx == Ap[i + 1]) {
                if (idx == nnz) {
                    break;
                } else {
                    continue;
                }
            }
            
            if (Ai[idx] == i) {
                Axres[i] -= coeff * Axdata[idx];
            } else {
                Axres[Ai[idx]] -= coeff * Axdata[idx];
                Axres[i] -= xdata[Ai[idx]] * Axdata[idx];
            }
            
            for (DSDP_INT j = Ap[i] + 1; j < Ap[i + 1]; ++j) {
                idx = Ai[j];
                Axres[idx] -= coeff * Axdata[j];
                Axres[i]   -= xdata[idx] * Axdata[j];
            }
        }
    }
    
    return retcode;
}

extern double spsMatxTAx( spsMat *A, double *x ) {
    // Compute quadratic form x' * A * x
    /*
     Computationally Critical Routine
     x is dense and A is sparse
    */
    
    register DSDP_INT *Ai = A->i, *nzHash = A->nzHash, i, j, k;
    register double res = 0.0, tmp = 0.0, *Ax = A->x;
    
    for (k = 0; k < A->nnz; ++k) {
        i = nzHash[k]; j = Ai[k];
        tmp = Ax[k] * x[i] * x[j];
        res += (i == j) ? (0.5 * tmp) : tmp;
    }
    
    return 2.0 * res;
}

extern DSDP_INT spsMataXpbY( double alpha, spsMat *sXMat, double beta,
                             spsMat *sYMat, DSDP_INT *sumHash ) {
    
    // Matrix axpy operation: let sYMat = alpha * sXMat + beta * sYMat
    // Note that sYMat must have a hash table
    DSDP_INT retcode = DSDP_RETCODE_OK;

    // The matrices added is NEVER factorized since they purely serve as constant data
    assert ( sXMat->dim == sYMat->dim );
    assert ((!sXMat->isFactorized));
    
    if (beta != 1.0) {
        retcode = spsMatScale(sYMat, beta);
    }
    
    // Case if alpha = 0.0
    if (alpha == 0.0) {
        return retcode;
    }
    
    DSDP_INT dim = sXMat->dim, *Ap = sXMat->p, *Ai = sXMat->i;
    double   *Ax = sXMat->x, *Bx = sYMat->x;
    
    if (sumHash) {
        for (DSDP_INT i = 0, hash; i < dim; ++i) {
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
    } else {
        // sYMat is actually dense and Bx is an (n + 1) * n / 2 array
        assert( sYMat->nnz == nsym(dim) );
        // If adding sparse data e.g. S += A * y[i]
        if (sXMat->nzHash) {
            for (DSDP_INT i, j, k = 0; k < sXMat->nnz; ++k) {
                i = Ai[k];
                j = sXMat->nzHash[k];
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
    
    assert( sMat->dim );
    
    if (d == 0.0) {
        return retcode;
    }
    
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
            sMat->x[idx] += d;
            idx += dim - i;
        }
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
    
    DSDP_INT dim = nsym(sXMat->dim);
    axpy(&dim, &alpha, dsYMat->array, &one, sXMat->x, &one);
    
    return retcode;
}

extern DSDP_INT spsMatAddr1( spsMat *sXMat, double alpha, r1Mat *r1YMat, DSDP_INT *sumHash ) {
    
    // Add a rank 1 matrix to a sparse matrix
    DSDP_INT retcode = DSDP_RETCODE_OK;
    
    assert( sXMat->dim == r1YMat->dim );
    
    if (alpha == 0.0) {
        return retcode;
    }
    
    DSDP_INT dim = sXMat->dim, idx = 0;
    DSDP_INT *hash = sumHash, *nzIdx = r1YMat->nzIdx;
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
    assert( sXMat->dim == rkYMat->dim);
    
    if (alpha == 0.0) {
        return retcode;
    }
    
    for (DSDP_INT i = 0; i < rkYMat->rank; ++i) {
        retcode = spsMatAddr1(sXMat, alpha, rkYMat->data[i], sumHash);
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
    
    if (fabs(r) < 1e-15) {
        error(etype, "Dividing a matrix by 0. \n");
    }
    
    if (sXMat->factor) {
        rkMatRscale(sXMat->factor, r);
    }
    
    double *x = sXMat->x;
    vecdiv(&sXMat->nnz, &r, x, &one);
    return retcode;
}

extern DSDP_INT spsMatFnorm( spsMat *sMat, double *fnrm ) {
    
    // Matrix Fronenius norm
    DSDP_INT retcode = DSDP_RETCODE_OK;
    assert( sMat->dim > 0 );
    
    DSDP_INT idx, i, n = sMat->dim;
    double nrm = 0.0, tmp;
    
    spsMatGetRank(sMat, &i);
    
    if (i < 0.1 * n) {
        rkMatFnorm(sMat->factor, fnrm);
    }
    
    if (sMat->nzHash) {
        DSDP_INT j, k;
        for (i = 0; i < sMat->nnz; ++i) {
            j = sMat->i[i];
            k = sMat->nzHash[i];
            if (j == k) {
                nrm += sMat->x[i] * sMat->x[i];
            } else {
                nrm += 2 * sMat->x[i] * sMat->x[i];
            }
        }
        *fnrm = sqrt(nrm);
    } else if (sMat->nnz == nsym(n)) {
        char ntype = 'F';
        char low = DSDP_MAT_LOW;
        *fnrm = dlansp(&ntype, &low, &n, sMat->x, NULL);
    } else {
        assert( FALSE );
        nrm = norm(&sMat->nnz, sMat->x, &one);
        nrm = 2 * nrm * nrm;
        for (i = 0; i < sMat->dim; ++i) {
            idx = sMat->p[i];
            if (idx == sMat->nnz) { break;}
            if (idx < sMat->p[i + 1] && sMat->i[idx] == i) {
                tmp = sMat->x[idx];
                nrm -= tmp * tmp;
            }
        }
        *fnrm = sqrt(nrm);
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
            j = sMat->i[i];
            k = sMat->nzHash[i];
            if (j == k) {
                nrm += 0.5 * fabs(sMat->x[i]);
            } else {
                nrm += fabs(sMat->x[i]);
            }
        }
        nrm *= 2;
    } else {
        assert( sMat->nnz == nsym(sMat->dim) );
        for (i = 0; i < sMat->nnz; ++i) {
            nrm += fabs(sMat->x[i]);
        }
        
        j = 0; k = sMat->dim;
        for (i = 0; i < k; ++i) {
            nrm -= 0.5 * fabs(sMat->x[j]);
            j += k - i;
        }
        nrm *= 2;
    }

    *onenrm = nrm;
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

extern DSDP_INT spsMatVecFSolve( spsMat *sAmat, vec *sbVec, vec *Ainvb ) {
    /* Forward solve L * x = b */
    return pardisoForwardSolve(sAmat, 1, sbVec->x, Ainvb->x, FALSE);
}

extern DSDP_INT spsMatVecBSolve( spsMat *sAmat, vec *sbVec, vec *Ainvb ) {
    /* Backward solve L' * x = b */
    return pardisoBackwardSolve(sAmat, 1, sbVec->x, Ainvb->x, FALSE);
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
    // A full dense matrix inv(A) * B is overwritten
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
    retcode = pardisoForwardSolve(S, n, fulldS, aux, TRUE);
    
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
    // Currently FEAST is used by transferring the matrix into dense format and is INEFFICIENT
    // But it may still be used at the end of the solver to recover highly accurate X
    DSDP_INT *Ap = spaux->p;
    DSDP_INT *Ai = spaux->i;
    double   *Ax = spaux->x;
    
    retcode = pardisoForwardSolve(S, n, fulldS, aux, TRUE);
    
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

extern DSDP_INT spsMatGetX( spsMat *S, spsMat *dS, double *LinvSLTinv ) {
    // Routine for retrieving X
    DSDP_INT retcode = DSDP_RETCODE_OK;
    
    DSDP_INT n = S->dim;
    assert( dS->dim == n );
    double *fulldS = LinvSLTinv;
    double *aux    = (double *) calloc(n * n, sizeof(double));
    // L^-1 dS
    double tmp = 0.0;
    retcode = spsMatFill(dS, fulldS);
    retcode = pardisoForwardSolve(S, n, fulldS, aux, TRUE);
    // Transpose
    arrayTranspose(fulldS, n);
    retcode = pardisoForwardSolve(S, n, fulldS, aux, TRUE);
    
    // I + L^-1 * dS * LT^-1
    for (DSDP_INT i = 0; i < n; ++i) {
        fulldS[i * n + i] += 1.0;
        for (DSDP_INT j = i + 1; j < n; ++j) {
            tmp = (fulldS[i * n + j] + fulldS[j * n + i]) * 0.5;
            fulldS[i * n + j] = fulldS[j * n + i] = tmp;
        }
    }
    
    pardisoBackwardSolve(S, n, fulldS, aux, TRUE);
    arrayTranspose(fulldS, n);
    pardisoBackwardSolve(S, n, fulldS, aux, TRUE);
    
    // Fix numerical instability
    for (DSDP_INT i = 0; i < n; ++i) {
        for (DSDP_INT j = i + 1; j < n; ++j) {
            tmp = (fulldS[i * n + j] + fulldS[j * n + i]) * 0.5;
            fulldS[i * n + j] = fulldS[j * n + i] = tmp;
        }
    }
    
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
    
    // Check block size of 1
    if (n == 1) {
        if (dS->x[0] > 0) {
            *alpha = DSDP_INFINITY;
        } else {
            *alpha = - S->x[0] / dS->x[0];
        }
        
        return retcode;
    }
    
    // Lanczos iteration
    double lbd = 0.0, delta = 0.0;
    // retcode = dsdpLanczos(S, dS, &lbd, &delta);
    
    if (1 || lbd != lbd || delta != delta || (lbd + delta > 1e+10)) {
        // MKL extremal routine
        retcode = spsMatLspLSolve(S, dS, spaux); checkCode;
        retcode = spsMatMinEig(spaux, &mineig); checkCode;
        if (mineig > 0) {
            *alpha = DSDP_INFINITY;
        } else {
            *alpha = - 1.0 / mineig;
        }
    } else {
        if (lbd + delta < 0) {
            *alpha = DSDP_INFINITY;
        } else {
            *alpha = 1 / (lbd + delta);
        }
    }
    return retcode;
}

extern DSDP_INT dsdpGetAlphaLS( spsMat *S, spsMat *dS, spsMat *Scker,
                                double alphamax, double *alpha, DSDP_INT *sumHash ) {
    // Get the maximum alpha such that S + alpha * dS is PSD by line-search
    DSDP_INT retcode = DSDP_RETCODE_OK;
    
    double step = 1 / alphamax;
    double *src = NULL;
    DSDP_INT ispsd = FALSE;
    spsMat *buffer = NULL;
    
    if (Scker) {
        src = S->x;
        buffer = Scker;
    } else {
        src = (double *) calloc(S->nnz, sizeof(double));
        memcpy(src, S->x, sizeof(double) * S->nnz);
        buffer = S;
    }
    
    for (DSDP_INT i = 0; ; ++i) {
        if (step <= 1e-04) {
            *alpha = 0.0;
            break;
        }
        memcpy(src, S->x, sizeof(double) * S->nnz);
        spsMataXpbY(step, dS, 1.0, buffer, sumHash);
        spsMatIspd(buffer, &ispsd);
        
        if (ispsd) {
            *alpha = step;
            break;
        }
        step *= 0.8;
    }

    if (!Scker) {
        DSDP_FREE(src);
    }
    
    return retcode;
}

extern DSDP_INT spsSinvSpSinvSolve( spsMat *S, spsMat *A, dsMat *SinvASinv, double *asinv ) {
    
    // Routine for setting up the Schur matrix
    /*
        Compute inv(S) * A * inv(S) for sparse A
     
    */
    DSDP_INT retcode = DSDP_RETCODE_OK;
    
    DSDP_INT n = S->dim;
    double *Sinv  = (double *) calloc(n * n, sizeof(double)), *rhs, *x;
    double *ASinv = (double *) calloc(n * n, sizeof(double));
    double *tmpSinvASinv = (double *) calloc(n * n, sizeof(double));
    double *sol   = (double *) calloc(n, sizeof(double));
    DSDP_INT *nzIdx = (DSDP_INT *) calloc(n, sizeof(DSDP_INT));
    
    // Sinv
    rhs = &Sinv[0];
    rhs[0] = 1.0;
    nzIdx[0] = 1;
    pardisoSolve(S, 1, rhs, sol);
    memcpy(rhs, sol, sizeof(double) * n);
    
    for (DSDP_INT i = 1; i < n; ++i) {
        rhs = &Sinv[n * i];
        rhs[i] = 1.0;
        nzIdx[i - 1] = 0;
        nzIdx[i] = 1;
        pardisoSolve(S, 1, rhs, sol);
        memcpy(rhs, sol, sizeof(double) * n);
    }
    
    DSDP_FREE(sol);
    
    // A * Sinv
    DSDP_INT *Ap = A->p, *Ai = A->i, nnz = A->nnz, idx;
    double *Ax = A->x, coeff = 0.0, res = 0.0;
    
    for (DSDP_INT i = 0; i < n; ++i) {
        // sol = ASinv(:, i) = A * Sinv(:, i) = A * x
        x = &Sinv[n * i];
        sol = &ASinv[n * i];
        for (DSDP_INT col = 0; col < n; ++col) {
            coeff = x[col];
            if (fabs(coeff) < 1e-20) {
                continue;
            }
            
            idx = Ap[col];
            if (idx == Ap[col + 1]) {
                if (idx == nnz) {
                    break;
                } else {
                    continue;
                }
            }
            
            if (Ai[idx] == col) {
                sol[col] += 0.5 * coeff * Ax[idx];
            } else {
                sol[Ai[idx]] += coeff * Ax[idx];
            }
            
            for (DSDP_INT k = idx + 1; k < Ap[col + 1]; ++k) {
                sol[Ai[k]] += coeff * Ax[k];
            }
        }
    }
    
    for (DSDP_INT i = 0; i < n; ++i) {
        res += ASinv[n * i + i];
    }
    
    *asinv = res * 2;

    // Sinv * Lower(A) * Sinv
    char side = 'L';
    char uplo = DSDP_MAT_LOW;
    double alpha = 1.0, beta = 0.0, *X = SinvASinv->array;
    
    dsymm(&side, &uplo, &n, &n, &alpha, Sinv, &n,
          ASinv, &n, &beta, tmpSinvASinv, &n);
    
    // Sinv * Lower(A) * Sinv + (Sinv * Lower(A) * Sinv)'
    for (DSDP_INT i = 0; i < n; ++i) {
        for (DSDP_INT j = 0; j <= i; ++j) {
            packIdx(X, n, i, j) = tmpSinvASinv[j * n + i] + tmpSinvASinv[i * n + j];
        }
    }
    
    DSDP_FREE(Sinv);
    DSDP_FREE(ASinv);
    DSDP_FREE(tmpSinvASinv);
    DSDP_FREE(nzIdx);
    return retcode;
}

extern DSDP_INT spsSinvSpSinvSolveSlow( spsMat *S, spsMat *A, dsMat *SinvASinv, double *asinv ) {
    
    // Routine for setting up the Schur matrix
    /*
        Compute inv(S) * A * inv(S) for sparse A
     
    */
    DSDP_INT retcode = DSDP_RETCODE_OK;
    assert( S->dim == A->dim );
    assert( S->dim == SinvASinv->dim );
    assert( S->isFactorized );
    
    DSDP_INT n         = S->dim;
    DSDP_INT inc       = n + 1;
    DSDP_INT idx       = 0;
    double *SinvA      = NULL;
    double *fSinvASinv = NULL;
    SinvA      = (double *) calloc(n * n, sizeof(double));
    fSinvASinv = (double *) calloc(n * n, sizeof(double));
    
    // Solve to get inv(S) * A
    retcode    = spsMatSpSolve(S, A, SinvA);
    
#ifdef TRANS
    double res = *asinv;
    
    for (DSDP_INT i = 0; i < n - n % 8; i+=8) {
        res += SinvA[idx]; idx += inc;
        res += SinvA[idx]; idx += inc;
        res += SinvA[idx]; idx += inc;
        res += SinvA[idx]; idx += inc;
        res += SinvA[idx]; idx += inc;
        res += SinvA[idx]; idx += inc;
        res += SinvA[idx]; idx += inc;
        res += SinvA[idx]; idx += inc;
    }
    
    for (DSDP_INT i = n - n % 8; i < n; ++i) {
        res += SinvA[idx]; idx += inc;
    }
    *asinv = res;
    MKL_Dimatcopy('C', 'T', n, n, 1.0, SinvA, n, n);
#else
    double tmp = 0.0;
    // Transpose
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
#endif
    
    idx = 0;
    // Solve
    retcode = pardisoSolve(S, n, SinvA, fSinvASinv);
    
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

extern DSDP_INT spsSinvRkSinvSolve( spsMat *S, rkMat *A, rkMat *SinvASinv, double *asinv ) {
    
    DSDP_INT retcode = DSDP_RETCODE_OK;
    assert( S->dim == A->dim );
    assert( SinvASinv->dim == S->dim );
    assert( A->isdata && !SinvASinv->isdata );
    
    double res = 0.0, tmp;
    DSDP_INT rank = A->rank;
    SinvASinv->rank = rank;
    r1Mat *r1data = NULL;

    for (DSDP_INT i = 0; i < rank; ++i) {
        r1data = A->data[i];
        retcode = spsSinvR1SinvSolve(S, r1data, SinvASinv->data[i], &tmp);
        res += tmp;
    }
    
    *asinv = res;
    
    return retcode;
}

extern DSDP_INT spsSinvR1SinvSolve( spsMat *S, r1Mat *A, r1Mat *SinvASinv, double *asinv ) {
    
    // Routine for setting up the Schur matrix
    DSDP_INT retcode = DSDP_RETCODE_OK;
    DSDP_INT n = S->dim;
    
    double *xSinvASinv = SinvASinv->x, *xA = A->x;
    
    assert( A->sign );
    SinvASinv->sign = A->sign;
    retcode = pardisoSolve(S, 1, xA, xSinvASinv);
    SinvASinv->nnz = n;
    
    double res = 0.0;
    DSDP_INT i;
    
    if (n > 64) {
        
        for (i = 0; i < n - 7; ++i) {
            res += xA[i] * xSinvASinv[i]; i++;
            res += xA[i] * xSinvASinv[i]; i++;
            res += xA[i] * xSinvASinv[i]; i++;
            res += xA[i] * xSinvASinv[i]; i++;
            res += xA[i] * xSinvASinv[i]; i++;
            res += xA[i] * xSinvASinv[i]; i++;
            res += xA[i] * xSinvASinv[i]; i++;
            res += xA[i] * xSinvASinv[i];
        }
        
        if (i < n - 3) {
            res += xA[i] * xSinvASinv[i]; i++;
            res += xA[i] * xSinvASinv[i]; i++;
            res += xA[i] * xSinvASinv[i]; i++;
            res += xA[i] * xSinvASinv[i]; i++;
        }
        
        if (i < n - 1) {
            res += xA[i] * xSinvASinv[i]; i++;
            res += xA[i] * xSinvASinv[i]; i++;
        }
        
        if (i < n) {
            res += xA[i] * xSinvASinv[i]; i++;
        }
        
    } else {
        for (DSDP_INT i = 0; i < n; ++i) {
            res += xA[i] * xSinvASinv[i];
        }
    }
    
    *asinv = res * A->sign;    
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
    
    if (info != SPARSE_STATUS_SUCCESS) {
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

extern DSDP_INT spsMatStoreFactor( spsMat *sMat, rkMat *factor ) {
    
    // Save factorized data
    DSDP_INT retcode = DSDP_RETCODE_OK;
    sMat->factor = factor;
    return retcode;
}

extern rkMat* spsMatGetFactor( spsMat *sMat ) {
    return sMat->factor;
}

extern DSDP_INT spsMatGetRank( spsMat *sMat, DSDP_INT *rank ) {
    
    if (sMat->factor) {
        *rank = sMat->factor->rank;
    } else {
        *rank = sMat->dim;
    }
    
    return DSDP_RETCODE_OK;
}

extern DSDP_INT spsMatFillLower( spsMat *sMat, double *lowFullMat ) {
    
    DSDP_INT retcode = DSDP_RETCODE_OK;
    assert( sMat->dim > 0 );
    DSDP_INT n = sMat->dim, ni = 0;
    DSDP_INT *Ap = sMat->p, *Ai = sMat->i;
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
    assert( sMat->dim > 0 );
    
    DSDP_INT n = sMat->dim, *Ap = sMat->p, *Ai = sMat->i;
    double *Ax = sMat->x;
    
    for (DSDP_INT i = 0; i < n; ++i) {
        for (DSDP_INT j = Ap[i]; j < Ap[i + 1]; ++j) {
            fulldMat[i * n + Ai[j]] = Ax[j];
            fulldMat[Ai[j] * n + i] = Ax[j];
        }
    }
    return retcode;
}

extern DSDP_INT spsMatReset( spsMat *sMat ) {
    
    // Reset a sparse matrix to be 0
    DSDP_INT retcode = DSDP_RETCODE_OK;
    assert( sMat->dim > 0 && !sMat->factor ); // Never reset an (eigen) factorized sparse matrix
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

extern void spsMatLinvView( spsMat *S ) {
    // Lanczos debugging routine. Print P^-1 L^-1 to the screen
    
    assert( S->isFactorized );
    DSDP_INT n = S->dim;
    
    double *eye = (double *) calloc(n * n, sizeof(double));
    double *Linv = (double *) calloc(n * n, sizeof(double));
    
    for (DSDP_INT i = 0; i < n; ++i) {
        eye[n * i + i] = 1.0;
    }
    
    pardisoForwardSolve(S, n, eye, Linv, FALSE);
    
    printf("Matrix view: \n");
    for (DSDP_INT i = 0; i < n; ++i) {
        for (DSDP_INT j = 0; j < n; ++j) {
            printf("%20.12e, ", Linv[i * n + j]);
        }
        printf("\n");
    }
    
    DSDP_FREE(eye);
    DSDP_FREE(Linv);
}
