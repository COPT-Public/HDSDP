#include <stdio.h>
#include <math.h>
#include <string.h>
#include "speigs.h"
#include "spinfo.h"

static void speigs_is_diag( spint *p, spint *i, spint n, spint *is_diag ) {
    /* Find out if the matrix is diagonal */
    for ( spint j = 0; j < n; ++j ) {
        if ( (p[j + 1] - p[j] > 0) && i[j] != j ) { *is_diag = FALSE; return; }
    }
    
    *is_diag = TRUE; return;
}

static void speigs_is_rankone( spint *p, spint *i, double *x, spint n,
                               spint *is_rankone, double *work, double tol ) {
    /* Find out if the matrix is rank one*/
    *is_rankone = TRUE;
    spint j, k, sgn, c = 0, nnz = 0; double d, err = 0.0;
    
    for ( j = 0; j < n; ++j ) {
        c = j; if ( p[j + 1] - p[j] > 0 ) { break; }
    }
    
    if (i[0] != c ) { *is_rankone = FALSE; return; }
    d = x[0]; sgn = (d > 0) ? 1 : 0;
    d = ( sgn ) ? sqrt(d) : -sqrt(-d);
    
    if ( d < tol ) { *is_rankone = FALSE; return; }
    
    memset(work, 0, sizeof(double) * n);
    for ( j = p[c]; j < p[c + 1]; ++j ) {
        work[i[j]] = x[j] / d; nnz += (x[j] != 0.0);
    }
    
    if ( 2 * p[n + 1] != nnz * (nnz + 1) ) { *is_rankone = FALSE; return; }
    
    for ( j = c + 1; j < n; ++j ) {
        if ( p[j] > p[c + 1] && p[j] < p[n + 1] ) {
            if ( i[p[j]] < c ) { *is_rankone = FALSE; return; }
        }
    }
    
    if ( sgn ) {
        for ( j = 0; j < n; ++j ) {
            for ( k = p[j]; k < p[j + 1]; ++k ) {
                err += fabs(x[k] + work[j] * work[i[k]]);
            }
            if ( err > tol ) { *is_rankone = FALSE; return; }
        }
    } else {
        for ( j = 0; j < n; ++j ) {
            for ( k = p[j]; k < p[j + 1]; ++k ) {
                err += fabs(x[k] - work[j] * work[i[k]]);
            }
            if ( err > tol ) { *is_rankone = FALSE; return; }
        }
    }
    
    return;
}

static void speigs_compute_submat( spint *p, spint *i, spint n, spint *sn,
                                   spint *nnzs, spint *perm, spint *iperm ) {
    /* Find out the submatrix embedded */
    spint s = 0, j;
    for ( j = 0; j < n; ++j ) {
        if ( nnzs[j] ) {
            perm[j] = s; iperm[s] = j; ++s;
        }
    }
    *sn = s; return;
}

static spint speigs_factorize_zero( spint *p, spint *i, double *x, spint n, spint *aiwork, double *awork, spint *sn,
                                    spint *iwork, spint *liwork, double *work, spint *lwork,
                                    double *evals, double *evecs, spint *rank ) {
    /* Factorize a zero matrix. Nothing is needed. */
    spint retcode = SP_EIGS_OK;
    *rank = 0; return retcode;
}

static spint speigs_factorize_diag( spint *p, spint *i, double *x, spint n, spint *aiwork, double *awork, spint *sn,
                                    spint *iwork, spint *liwork, double *work, spint *lwork,
                                    double *evals, double *evecs, spint *rank ) {
    /* Factorize a diagonal matrix. */
    spint retcode = SP_EIGS_OK, *nnzs = aiwork, j, k;
    double *evec = evecs;
    
    for ( j = k = 0; j < n; ++j ) {
        if ( nnzs[j] ) {
            evals[k] = x[k]; evec[j] = 1.0;
            ++k; evec += n;
        }
    }
    
    *rank = k; return retcode;
}

static spint speigs_factorize_two( spint *p, spint *i, double *x, spint n, spint *aiwork, double *awork, spint *sn,
                                   spint *iwork, spint *liwork, double *work, spint *lwork,
                                   double *evals, double *evecs, spint *rank ) {
    /* Factorize a two-two matrix. Possible that diagonal elements exist */
    spint retcode = SP_EIGS_OK, *nnzs = aiwork, j, k;
    double *evec = evecs;
    
    for ( j = k = 0; j < n; ++j ) {
        if ( nnzs[j] ) {
            if ( i[p[j]] == j ) {
                evals[k] = x[k]; evecs[j] = 1.0;
                ++k; evec += n;
            } else {
                evals[i[p[j]]] = x[k]; evecs[j] = ROOT;
                ++k; evecs += n;
                evals[j] = x[k]; evecs[i[p[j]]] = -ROOT;
                ++k; evecs += n;
            }
        }
    }
    
    *rank = k; return retcode;
}

static spint speigs_factorize_rankone( spint *p, spint *i, double *x, spint n, spint *aiwork, double *awork, spint *sn,
                                       spint *iwork, spint *liwork, double *work, spint *lwork,
                                       double *evals, double *evecs, spint *rank ) {
    /* Factorize a rankone matrix. Reuse information from awork */
    spint retcode = SP_EIGS_OK;
    evals[0] = ( x[0] > 0 ) ? 1.0 : -1.0;
    memcpy(evecs, awork, sizeof(double) * n);
    *rank = 1; return retcode;
}

static spint speigs_factorize_dense( double *a, double *evals, double *evecs, spint *n, spint *liwork, spint *iwork,
                                    spint *lwork, double *work ) {
    /* Wrapper for Lapack dense factorization routine */
    spint retcode = SP_EIGS_OK;
    
    
    return retcode;
}

static spint speigs_factorize_sparse( spint *p, spint *i, double *x, spint n, spint *aiwork, double *awork, spint *sn,
                                      spint *iwork, spint *liwork, double *work, spint *lwork,
                                      double *evals, double *evecs, spint *rank ) {
    /* Factorize a sparse matrix. Use information of the submatrix */
    spint retcode = SP_EIGS_OK, *perm = aiwork + n, *iperm = aiwork + 2 * n, j, k, s = *sn;
    
    for ( j = 0; j < p[n + 1]; ++j ) {
        for ( k = p[j]; k < p[j + 1]; ++k ) {
            work[s * perm[i[k]] + perm[  j ]] = x[k];
            work[s * perm[  j ] + perm[i[k]]] = x[k];
        }
    }
    
    return retcode;
}

static spint speigs_factorize_general( spint *p, spint *i, double *x, spint n, spint *aiwork, double *awork, spint *sn,
                                       spint *iwork, spint *liwork, double *work, spint *lwork,
                                       double *evals, double *evecs, spint *rank ) {
    spint retcode = SP_EIGS_OK;
    return retcode;
}

extern spint speigs_analyze( spint *Ap, spint *Ai, double *Ax, spint *dim,
                             spint *iwork, spint *liwork, double *work, spint *lwork,
                             spint *type, spint *sn, double tol, double gthresh ) {
    /* Detect the structure of the sparse matrix to accelerate decomposition */
    spint retcode = SP_EIGS_OK;
    
    if ( !liwork || !lwork ) {
        sperr("Invalid working array length. \n");
        retcode = SP_EIGS_ERR; return retcode;
    }
    
    if ( !dim ) {
        sperr("Invalid dimension. \n");
        retcode = SP_EIGS_ERR; return retcode;
    }
    
    if ( !Ap || !Ai || !Ax || !Ax || !iwork || !work ) {
        /* Estimate size */
        *liwork = 3 * (*dim); *lwork = *dim;
        return retcode;
    }
    
    if ( *liwork < 3 * (*dim) || *lwork < *dim ) {
        sperr("Insufficient space for analysis phase. \n");
        retcode = SP_EIGS_ERR; return retcode;
    }
    
    if ( !type ) {
        sperr("Invalid type. \n");
    }
    
    /* Begin analysis */
    spint *p = Ap, *i = Ai, n = *dim, j;
    spint *nnzs = iwork + 0, *perm = iwork + n, *iperm = iwork + 2 * n;
    spint is_diag = FALSE, is_two = FALSE, is_rankone = FALSE;
    double *x = Ax;
    
    /* Case 1: Zero matrix */
    if ( p[n + 1] == 0 ) {
        *sn = 0; *type = MATRIX_TYPE_ZERO;
        return retcode;
    }
    
    /* Collect column statistics */
    for ( j = 0; j < n; ++j ) {
        nnzs[j] = p[j + 1] - p[j];
        if ( nnzs[j] >= 2 ) {
            is_diag = FALSE; is_two = FALSE;
        }
    }
    
    /* Case 2: Diagonal matrix */
    if ( is_diag ) {
        speigs_is_diag(p, i, n, &is_diag);
        if ( is_diag ) {
            *sn = n; *type = MATRIX_TYPE_DIAG;
            return retcode;
        }
    }
    
    /* Case 3: Two-two */
    if ( is_two ) {
        *sn = n; *type = MATRIX_TYPE_TWOTWO;
        return retcode;
    }
    
    /* Case 4: Rank-one */
    speigs_is_rankone(p, i, x, n, &is_rankone, work, tol);
    
    if ( is_rankone ) {
        *sn = n; *type = MATRIX_TYPE_RANKONE;
        return retcode;
    }
    
    /* Case 5: Sparse submatrix */
    speigs_compute_submat(p, i, n, sn, nnzs, perm, iperm);
    if ( *sn < gthresh * n ) {
        *type = MATRIX_TYPE_SPARSE;
        return retcode;
    }

    /* Case 6: Nothing found */
    *sn = n; *type = MATRIX_TYPE_GENERAL;
    return retcode;
}

extern spint speigs_factorize( spint *Ap, spint *Ai, double *Ax, spint *dim, spint *aiwork, double *awork,
                               spint *type, spint *sn, spint *iwork, spint *liwork, double *work, spint *lwork,
                               double *evals, double *evecs, spint *rank, double tol ) {
    /* Factorize the matrix based on the structure */
    spint retcode = SP_EIGS_OK;
    
    if ( !type ) {
        sperr("Invalid type. \n");
        retcode = SP_EIGS_ERR; return retcode;
    }
    
    if ( !sn ) {
        sperr("Invalid submatrix size. \n");
        retcode = SP_EIGS_ERR; return retcode;
    }
    
    if ( !dim ) {
        sperr("Invalid dimension. \n");
        retcode = SP_EIGS_ERR; return retcode;
    }
    
    if ( !liwork || !lwork ) {
        sperr("Invalid working array length. \n");
        retcode = SP_EIGS_ERR; return retcode;
    }
    
    
    if ( !Ap || !Ai || !aiwork || !awork || !iwork || !lwork || !evals || !evecs ) {
        *liwork = (*sn) * LAPACK_IWORK;
        *lwork = (*sn) * (*sn) + (*sn) * LAPACK_LWORK;
        return retcode;
    }
    
    if ( *liwork <= (*sn) * LAPACK_IWORK ||
         *lwork <= (*sn) * (*sn) + (*sn) * LAPACK_LWORK ) {
        sperr("Insufficient space for factorization phase. \n");
        retcode = SP_EIGS_ERR; return retcode;
    }
    
    spint *p = Ap, *i = Ai, n = *dim;
    
    double *x = Ax;
    memset(evals, 0, sizeof(double) * n);
    memset(evecs, 0, sizeof(double) * n * n);
    
    /* Begin factorize */
    switch (*type) {
        case MATRIX_TYPE_ZERO:
            retcode = speigs_factorize_zero(p, i, x, n, aiwork, awork, sn, iwork,
                                            liwork, work, lwork, evals, evecs, rank);
            break;
        case MATRIX_TYPE_DIAG:
            retcode = speigs_factorize_diag(p, i, x, n, aiwork, awork, sn, iwork,
                                            liwork, work, lwork, evals, evecs, rank);
            break;
        
        case MATRIX_TYPE_TWOTWO:
            retcode = speigs_factorize_two(p, i, x, n, aiwork, awork, sn, iwork,
                                           liwork, work, lwork, evals, evecs, rank);
            break;
            
        case MATRIX_TYPE_RANKONE:
            retcode = speigs_factorize_rankone(p, i, x, n, aiwork, awork, sn, iwork,
                                               liwork, work, lwork, evals, evecs, rank);
            break;
            
        case MATRIX_TYPE_SPARSE:
            retcode = speigs_factorize_sparse(p, i, x, n, aiwork, awork, sn, iwork,
                                              liwork, work, lwork, evals, evecs, rank);
            break;
            
        case MATRIX_TYPE_GENERAL:
            retcode = speigs_factorize_sparse(p, i, x, n, aiwork, awork, sn, iwork,
                                              liwork, work, lwork, evals, evecs, rank);
            break;
        default:
            sperr("Unknown matrix type. \n");
            retcode = SP_EIGS_ERR;
    }
    
    return retcode;
}
