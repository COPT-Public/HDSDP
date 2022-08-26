/** @file speigs.c
 *  @brief The implementation of sparse eigen decomposition routine for HDSDP
 *
 * A set of routines that factorize very sparse matrices that typically arise from semi-definite programming problems.
 * The routines detect the special structures of the matrix and accelerate the factorization procedure.
 *
 *  @author Wenzhi Gao, Shanghai University of Finance and Economics
 *  @date Aug, 24th, 2022
 *
 */

#include <stdio.h>
#include <math.h>
#include <string.h>
#include "speigs.h"
#include "spinfo.h"

/* Define Lapack-related constants */
static char   jobz = 'V';
static char   range = 'A';
static char   uplolow = 'L';
static double abstol = 0.0;


/** @brief Compute lwork and iwork
 *  @param[in] n Dimension of the matrix
 *  @param[in] sn Dimension of the submatrix
 *  @param[in] type Type of the matrix
 *  @param[out] liwork Estimated integer workspace length
 *  @param[out] lwork Estimated double workspace length
 *
 *  The function computes the space requirement for the subroutines that will be invoked in the factorization phase
 */
static void speig_get_factorize_space( spint *n, spint *sn, spint *type, spint *liwork, spint *lwork ) {
    
    switch ( *type ) {
        case MATRIX_TYPE_ZERO:    *liwork = 0; *lwork = 0; break;
        case MATRIX_TYPE_DIAG:    *liwork = 0; *lwork = 0; break;
        case MATRIX_TYPE_TWOTWO:  *liwork = 0; *lwork = 0; break;
        case MATRIX_TYPE_RANKONE: *liwork = 0; *lwork = 0; break;
        case MATRIX_TYPE_SPARSE:
            *liwork = (*sn) * LAPACK_IWORK + (*sn) * 2; // Lapack iwork + iworkup
            *lwork = 2 * (*sn) * (*sn) + (*sn) * LAPACK_LWORK + (*sn); break;
        case MATRIX_TYPE_GENERAL:
            *liwork = (*n) * LAPACK_IWORK + (*n) * 2;
            *lwork = 2 * (*n) * (*n) + (*n) * LAPACK_LWORK + (*n); break;
        default: sperr("Unknown matrix type. \n"); break;
    }
    return;
}

/** @brief Check if a matrix is diagonal
 *  @param[in] p  CSC format column pointer
 *  @param[in] i  CSC format row index
 *  @param[in] n  Dimension of the matrix
 *  @param[out] is_diag Is the matrix diagonal?
 */
static void speigs_is_diag( spint *p, spint *i, spint n, spint *is_diag ) {
    
    for ( spint j = 0; j < n; ++j ) {
        if ( (p[j + 1] - p[j] > 0) && i[p[j]] != j ) {
            *is_diag = FALSE; return;
        }
    }
    
    *is_diag = TRUE; return;
}

/** @brief Find out if a matrix is rank-one. \f$ A = \alpha  a  a^T \f$
 *  @param[in] p  CSC format column pointer
 *  @param[in] i  CSC format row index
 *  @param[in] x  CSC format matrix nonzero entries
 *  @param[in] n  Dimension of the matrix
 *  @param[out] is_rankone Is the matrix rank-one?
 *  @param[out] work Working array for rank-one detection
 *  @param[in] tol Tolerance for rank-one classification \f$||A - a a^T||_F \leq tol\f$
 *
 *  On exit, the array "work" would be filled by the rank-one factor a if A is rank-one
 */
static void speigs_is_rankone( spint *p, spint *i, double *x, spint n,
                               spint *is_rankone, double *work, double tol ) {
    /*
     A sparse rank-one matrix must look like
     
                0  0  0  0
                0  0  0  0
                0  0 [*] *
                0  0  *  *
     */
    
    spint j, k, sgn, c = 0, nnz = 0; double d, err = 0.0;
    
    for ( j = 0; j < n; ++j ) {
        c = j; if ( p[j + 1] - p[j] > 0 ) { break; }
    }
    
    if ( i[0] != c ) { *is_rankone = FALSE; return; }
    d = x[0]; sgn = (d > 0) ? 1 : 0;
    d = ( sgn ) ? sqrt(d) : -sqrt(-d);
    
    /* The leading nonzero element is actually close to 0 */
    if ( fabs(d) < tol ) { *is_rankone = FALSE; return; }
    memset(work, 0, sizeof(double) * n);
    
    for ( j = p[c]; j < p[c + 1]; ++j ) {
        work[i[j]] = x[j] / d; nnz += (x[j] != 0.0);
    }
    
    if ( 2 * p[n] != nnz * (nnz + 1) ) {
        *is_rankone = FALSE; return;
    }
    
    for ( j = c + 1; j < n; ++j ) {
        if ( work[j] && i[p[j]] != j ) {
            *is_rankone = FALSE; return;
        }
    }
    
    if ( sgn ) {
        if ( p[n] >= 5 * n ) {
            for ( j = 0; j < n; ++j ) {
                for ( k = p[j]; k < p[j + 1]; ++k ) {
                    err += fabs(x[k] - work[j] * work[i[k]]);
                }
                if ( err > tol ) { *is_rankone = FALSE; return; }
            }
        } else {
            for ( j = 0; j < n; ++j ) {
                for ( k = p[j]; k < p[j + 1]; ++k ) {
                    err += fabs(x[k] - work[j] * work[i[k]]);
                }
            }
        }
    } else {
        if ( p[n] >= 5 * n ) {
            for ( j = 0; j < n; ++j ) {
                for ( k = p[j]; k < p[j + 1]; ++k ) {
                    err += fabs(x[k] + work[j] * work[i[k]]);
                }
                if ( err > tol ) { *is_rankone = FALSE; return; }
            }
        } else {
            for ( j = 0; j < n; ++j ) {
                for ( k = p[j]; k < p[j + 1]; ++k ) {
                    err += fabs(x[k] + work[j] * work[i[k]]);
                }
            }
        }
    }
    
    if ( err > tol ) { *is_rankone = FALSE; return; }
    *is_rankone = TRUE; return;
}

/** @brief Compute the dense submatrix of a large sparse matrix
 *  @param[in] p  CSC format column pointer
 *  @param[in] i  CSC format row index
 *  @param[in] n  Dimension of the matrix
 *  @param[out] sn Dimension of the submatrix
 *  @param[in] nnzs Number of nonzeros in each column
 *  @param[out] perm Permutation that gathers nonzero elements
 *  @param[out] iperm Inverse permutation
 *
 * On exit, "perm" and "iperm" will be filled by the permutation and its inverse respectively
 */
static void speigs_compute_submat( spint *p, spint *i, spint n, spint *sn,
                                   spint *nnzs, spint *perm, spint *iperm ) {
    
    spint s = 0, j;
    for ( j = 0; j < n; ++j ) {
        if ( nnzs[j] ) {
            perm[j] = s; iperm[s] = j; ++s;
        }
    }
    *sn = s; return;
}

/** @brief Compute the eigen factorization of an all-zero matrix
 *  @param[in] p  CSC format column pointer
 *  @param[in] i  CSC format row index
 *  @param[in] x  CSC format matrix nonzero entries
 *  @param[in] n  Dimension of the matrix
 *  @param[in] aiwork Integer working array from the analysis phase
 *  @param[in] awork Double working array from the analysis phase
 *  @param[in] sn Dimension of the submatrix
 *  @param[in] iwork Integer working array for the factorization phase
 *  @param[in] liwork Length of "iwork"
 *  @param[in] work Double working array for the factorization phase
 *  @param[in] lwork Length of "work"
 *  @param[out] evals Eigen-values after factorization
 *  @param[out] evecs Eigen-vectors after factorization
 *  @param[out] rank Rank of the factorized matrix
 *  @param[in] tol Tolerance to tell if an eigen-value is 0
 *  @return retcode Status of the factorization
 *
 * On exit, "evals" and "evecs" will be overwritten by the eigen-decomposition of the matrix. "rank" is the rank of the matrix
 * Since the matrix is all-zero, no operation is needed.
 */
static spint speigs_factorize_zero( spint *p, spint *i, double *x, spint n, spint *aiwork,
                                    double *awork, spint *sn, spint *iwork, spint *liwork,
                                    double *work, spint *lwork, double *evals, double *evecs,
                                    spint *rank, double tol ) {
    
    spint retcode = SP_EIGS_OK;
    *rank = 0; return retcode;
}

/** @brief Compute the eigen factorization of a diagonal matrix
 *  @param[in] p  CSC format column pointer
 *  @param[in] i  CSC format row index
 *  @param[in] x  CSC format matrix nonzero entries
 *  @param[in] n  Dimension of the matrix
 *  @param[in] aiwork Integer working array from the analysis phase
 *  @param[in] awork Double working array from the analysis phase
 *  @param[in] sn Dimension of the submatrix
 *  @param[in] iwork Integer working array for the factorization phase
 *  @param[in] liwork Length of "iwork"
 *  @param[in] work Double working array for the factorization phase
 *  @param[in] lwork Length of "work"
 *  @param[out] evals Eigen-values after factorization
 *  @param[out] evecs Eigen-vectors after factorization
 *  @param[out] rank Rank of the factorized matrix
 *  @param[in] tol Tolerance to tell if an eigen-value is 0
 *  @return retcode Status of the factorization
 *
 * On exit, "evals" and "evecs" will be overwritten by the eigen-decomposition of the matrix. "rank" is the rank of the matrix
 * Since the matrix is diagonal, all the eigen-vectors are unit vectors and eigen-values are determined by the elements in "x"
 */
static spint speigs_factorize_diag( spint *p, spint *i, double *x, spint n, spint *aiwork,
                                    double *awork, spint *sn, spint *iwork, spint *liwork,
                                    double *work, spint *lwork, double *evals, double *evecs,
                                    spint *rank, double tol ) {
    
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

/** @brief Compute the eigen factorization of a two-two matrix
 *  @param[in] p  CSC format column pointer
 *  @param[in] i  CSC format row index
 *  @param[in] x  CSC format matrix nonzero entries
 *  @param[in] n  Dimension of the matrix
 *  @param[in] aiwork Integer working array from the analysis phase
 *  @param[in] awork Double working array from the analysis phase
 *  @param[in] sn Dimension of the submatrix
 *  @param[in] iwork Integer working array for the factorization phase
 *  @param[in] liwork Length of "iwork"
 *  @param[in] work Double working array for the factorization phase
 *  @param[in] lwork Length of "work"
 *  @param[out] evals Eigen-values after factorization
 *  @param[out] evecs Eigen-vectors after factorization
 *  @param[out] rank Rank of the factorized matrix
 *  @param[in] tol Tolerance to tell if an eigen-value is 0
 *  @return retcode Status of the factorization
 *
 * On exit, "evals" and "evecs" will be overwritten by the eigen-decomposition of the matrix. "rank" is the rank of the matrix
 * Since the matrix composes of 2 by 2 submatrices, Givens' rotation is employed to factorize the matrixf
 */
static spint speigs_factorize_two( spint *p, spint *i, double *x, spint n, spint *aiwork,
                                   double *awork, spint *sn, spint *iwork, spint *liwork,
                                   double *work, spint *lwork, double *evals, double *evecs,
                                   spint *rank, double tol ) {
    
    spint retcode = SP_EIGS_OK, *nnzs = aiwork, j, k;
    double *v = evecs, *e = evals;
    
    for ( j = k = 0; j < n; ++j ) {
        if ( nnzs[j] ) {
            if ( i[p[j]] == j ) {
                e[k] = x[k]; v[j] = 1.0;
                ++k; v += n;
            } else {
                v[j] =  ROOT; v[i[p[j]]] = ROOT;
                e[k] = x[k]; ++k; v += n;
                v[j] = -ROOT; v[i[p[j]]] = ROOT;
                e[k] = x[k]; ++k; v += n;
            }
        }
    }
    
    if ( k > n ) {
        sperr("Invalid rank for two-two matrix. \n");
        retcode = SP_EIGS_ERR; return retcode;
    }
    
    *rank = k; return retcode;
}

/** @brief Compute the eigen factorization of a rank-one matrix
 *  @param[in] p  CSC format column pointer
 *  @param[in] i  CSC format row index
 *  @param[in] x  CSC format matrix nonzero entries
 *  @param[in] n  Dimension of the matrix
 *  @param[in] aiwork Integer working array from the analysis phase
 *  @param[in] awork Double working array from the analysis phase
 *  @param[in] sn Dimension of the submatrix
 *  @param[in] iwork Integer working array for the factorization phase
 *  @param[in] liwork Length of "iwork"
 *  @param[in] work Double working array for the factorization phase
 *  @param[in] lwork Length of "work"
 *  @param[out] evals Eigen-values after factorization
 *  @param[out] evecs Eigen-vectors after factorization
 *  @param[out] rank Rank of the factorized matrix
 *  @param[in] tol Tolerance to tell if an eigen-value is 0
 *  @return retcode Status of the factorization
 *
 * On exit, "evals" and "evecs" will be overwritten by the eigen-decomposition of the matrix. "rank" is the rank of the matrix
 * Since the matrix is rank-one, the "awork" array from the analysis phase contains the eigen-decomposition
 */
static spint speigs_factorize_rankone( spint *p, spint *i, double *x, spint n, spint *aiwork,
                                       double *awork, spint *sn, spint *iwork, spint *liwork,
                                       double *work, spint *lwork, double *evals, double *evecs,
                                       spint *rank, double tol ) {
    
    spint retcode = SP_EIGS_OK;
    evals[0] = ( x[0] > 0 ) ? 1.0 : -1.0;
    memcpy(evecs, awork, sizeof(double) * n);
    *rank = 1; return retcode;
}

/** @brief Compute the eigen factorization of a general full matrix
 *  @param[in] a Dense array that contains the matrix to factorize
 *  @param[out] evals Eigen-values after factorization
 *  @param[out] evecs Eigen-vectors after factorization
 *  @param[in] n  Dimension of the dense matrix
 *  @param[in] liwork Length of the integer working array for Lapack
 *  @param[in] iwork Integer working array for Lapack
 *  @param[in] lwork Length of double working array for Lapack
 *  @param[in] work Double working array for Lapack
 *  @param[in] isuppz Auxiliary placeholder for Lapack parameter
 *  @return retcode Status of the factorization
 *
 * On exit, "evals" and "evecs" will be overwritten by the eigen-decomposition of the matrix. "rank" is the rank of the matrix
 * The routine is a wrapper of the Lapack dsyevr function
 */
static spint speigs_factorize_dense( double *a, double *evals, double *evecs, spint *n,
                                     spint *liwork, spint *iwork, spint *lwork, double *work,
                                     spint *isuppz ) {
    
    spint retcode = SP_EIGS_OK, m, info;
    dsyevr(&jobz, &range, &uplolow, n, a, n, NULL, NULL, NULL, NULL, &abstol, &m,
           evals, evecs, n, isuppz, work, lwork, iwork, liwork, &info);
    if ( info ) {
        retcode = SP_EIGS_ERR;
        sperr("Eigen-decomposition failed \n");
    }
    
    return retcode;
}

/** @brief Compute the eigen factorization of a sparse matrix admitting an easier submatrix representation
 *  @param[in] p  CSC format column pointer
 *  @param[in] i  CSC format row index
 *  @param[in] x  CSC format matrix nonzero entries
 *  @param[in] n  Dimension of the matrix
 *  @param[in] aiwork Integer working array from the analysis phase
 *  @param[in] awork Double working array from the analysis phase
 *  @param[in] sn Dimension of the submatrix
 *  @param[in] iwork Integer working array for the factorization phase
 *  @param[in] liwork Length of "iwork"
 *  @param[in] work Double working array for the factorization phase
 *  @param[in] lwork Length of "work"
 *  @param[out] evals Eigen-values after factorization
 *  @param[out] evecs Eigen-vectors after factorization
 *  @param[out] rank Rank of the factorized matrix
 *  @param[in] tol Tolerance to tell if an eigen-value is 0
 *  @return retcode Status of the factorization
 *
 * On exit, "evals" and "evecs" will be overwritten by the eigen-decomposition of the matrix. "rank" is the rank of the matrix
 * The routine uses the permutation and inverse permutation information collected in the analysis phase to formulate the submatrix,
 * factorizes the submatrix and finally recovers the decomposition using the inverse permutation
 */
static spint speigs_factorize_sparse( spint *p, spint *i, double *x, spint n, spint *aiwork,
                                      double *awork, spint *sn, spint *iwork, spint *liwork,
                                      double *work, spint *lwork, double *evals, double *evecs,
                                      spint *rank, double tol ) {
    
    spint retcode = SP_EIGS_OK, *perm = aiwork + n, *iperm = aiwork + 2 * n, j, k, q, s = *sn;
    
    for ( j = 0; j < n; ++j ) {
        for ( k = p[j]; k < p[j + 1]; ++k ) {
            work[s * perm[j] + perm[i[k]]] = x[k]; // work[s * perm[i[k]] + perm[j]] = x[k];
        }
    }
    
    /* a: submatrix; work1: eigen value; work2: eigen vector; work3: Lapack work */
    double *a = work, *work1 = work + s * s, *work2 = work1 + s, \
           *work3 = work2 + s * s, *sev, *ev;
    spint lwork2 = s * LAPACK_LWORK, liwork2 = s * LAPACK_IWORK;
    spint *isuppz = iwork + liwork2;
    
    retcode = speigs_factorize_dense(a, work1, work2, sn, &liwork2,
                                     iwork, &lwork2, work3, isuppz);
    if ( retcode != SP_EIGS_OK ) { return retcode; }
    
    /* Successfully exits */
    for ( k = j = 0; j < s; ++j ) {
        if ( fabs(work1[j]) > tol ) {
            evals[k] = work1[j];
            sev = work2 + s * j; ev = evecs + n * k;
            for ( q = 0; q < s; ++q ) {
                ev[iperm[q]] = sev[q];
            }
            ++k;
        }
    }
    
    *rank = k; return retcode;
}

/** @brief Compute the eigen factorization of a general dense matrix
 *  @param[in] p  CSC format column pointer
 *  @param[in] i  CSC format row index
 *  @param[in] x  CSC format matrix nonzero entries
 *  @param[in] n  Dimension of the matrix
 *  @param[in] aiwork Integer working array from the analysis phase
 *  @param[in] awork Double working array from the analysis phase
 *  @param[in] sn Dimension of the submatrix
 *  @param[in] iwork Integer working array for the factorization phase
 *  @param[in] liwork Length of "iwork"
 *  @param[in] work Double working array for the factorization phase
 *  @param[in] lwork Length of "work"
 *  @param[out] evals Eigen-values after factorization
 *  @param[out] evecs Eigen-vectors after factorization
 *  @param[out] rank Rank of the factorized matrix
 *  @param[in] tol Tolerance to tell if an eigen-value is 0
 *  @return retcode Status of the factorization
 *
 * On exit, "evals" and "evecs" will be overwritten by the eigen-decomposition of the matrix. "rank" is the rank of the matrix
 * The routine converts the sparse matrix into a dense array and calls Lapack directly. Slow in general
 */
static spint speigs_factorize_general( spint *p, spint *i, double *x, spint n, spint *aiwork,
                                       double *awork, spint *sn, spint *iwork, spint *liwork,
                                       double *work, spint *lwork, double *evals, double *evecs,
                                       spint *rank, double tol ) {
    
    spint retcode = SP_EIGS_OK, s = *sn, j, k;
    
    if ( n != *sn ) {
        sperr("Invalid submatrix size for general matrix. \n");
        retcode = SP_EIGS_ERR; return retcode;
    }
    
    for ( j = 0; j < n; ++j ) {
        for ( k = p[j]; k < p[j + 1]; ++k ) {
            work[n * j + i[k]] = x[k];
        }
    }
    
    /* a: submatrix; work1: eigen value; work2: eigen vector; work3: Lapack work */
    double *a = work, *work1 = work + s * s, *work2 = work1 + s, \
           *work3 = work2 + s * s, *sev, *ev;
    spint lwork2 = s * LAPACK_LWORK, liwork2 = s * LAPACK_IWORK;
    spint *isuppz = iwork + liwork2;
    retcode = speigs_factorize_dense(a, work1, work2, &n, &liwork2,
                                     iwork, &lwork2, work3, isuppz);
    if ( retcode != SP_EIGS_OK ) { return retcode; }
    
    for ( k = j = 0; j < n; ++j ) {
        if ( fabs(work1[j]) > tol ) {
            evals[k] = work1[j];
            sev = work2 + n * j; ev = evecs + n * k;
            memcpy(ev, sev, sizeof(double) * n); ++k;
        }
    }
    
    *rank = k; return retcode;
}

/**
 *  @brief Perform the analysis phase of sparse eigen-value factorization
 *  @param[in] Ap  CSC format column pointer
 *  @param[in] Ai  CSC format row index
 *  @param[in] Ax  CSC format matrix nonzero entries
 *  @param[in] dim Dimension of the matrix
 *  @param[out] iwork Integer working array for the analysis phase
 *  @param[in] liwork Length of "iwork" or the expected length of integer working array
 *  @param[out] work Double working array for the analysis phase
 *  @param[in] lwork Length of "lwork" or the expected length of double working array
 *  @param[out] type Type of the matrix
 *  @param[out] sn Size of the submatrix
 *  @param[in] tol Tolerance to classify if a matrix is rank-one by \f$||A - a a^T||_F \leq tol\f$
 *  @param[in] gthresh Threshold of (submatrix size / dim) classifying a matrix as general or sparse
 *  @return retcode Status of the analysis phase
 *
 * Perform the analysis phase of the sparse eigen-value factorization.
 *
 * If all the necessary memories are allocated, on exit, "work" and "iwork" are filled by the intermediate
 * information which can be used in the factorization phase; "type" is filled by one of the five types;
 * "sn" is filled by size of the submatrix.
 *
 * If "dim" is supplied and the rest of the working array is incomplete, "lwork" and "work" will be respectively
 * filled by the expected length of the double and integer working arrays
 *
 */
extern spint speigs_analyze( spint *Ap, spint *Ai, double *Ax, spint *dim,
                             spint *iwork, spint *liwork, double *work, spint *lwork,
                             spint *type, spint *sn, double tol, double gthresh ) {
    
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
    spint is_diag = TRUE, is_two = TRUE, is_rankone = FALSE;
    double *x = Ax;
    
    /* Case 1: Zero matrix */
    if ( p[n] == 0 ) {
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

/** @brief The jump table for eigen routines
 *
 * Currently contains six implementations of eigen routines
 */
static spint (*speig_routines[6]) ( spint*, spint*, double*, spint, spint*,
                                    double*, spint*, spint*, spint*, double*,
                                    spint*, double*, double*, spint*, double ) =
{
    &speigs_factorize_zero,
    &speigs_factorize_sparse,
    &speigs_factorize_general,
    &speigs_factorize_rankone,
    &speigs_factorize_diag,
    &speigs_factorize_two
};

/**
 *  @brief Perform the analysis phase of sparse eigen-value factorization
 *  @param[in] Ap  CSC format column pointer
 *  @param[in] Ai  CSC format row index
 *  @param[in] Ax  CSC format matrix nonzero entries
 *  @param[in] dim Dimension of the matrix
 *  @param[in] aiwork Integer working array from the analysis phase
 *  @param[in] awork Double working array from the analysis phase
 *  @param[in] type "type" from the analysis phase
 *  @param[in] sn "sn" from the analysis phase
 *  @param[in] iwork Integer working array for the factorization phase
 *  @param[in] liwork Length of "iwork" or the expected length of integer working array
 *  @param[in] work Double working array for the factorization phase
 *  @param[in] lwork Length of "lwork" or the expected length of double working array
 *  @param[out] evals Eigen-values after factorization
 *  @param[out] evecs Eigen-vectors after factorization
 *  @param[out] rank Rank of the factorized matrix
 *  @param[in] tol Tolerance to tell if an eigen-value is 0
 *  @return retcode Status of the factorization phase
 *
 * Perform the analysis phase of the sparse eigen-value factorization.
 *
 * If all the necessary memories are allocated, on exit, "work" and "iwork" are filled by the intermediate
 * information which can be used in the factorization phase; "type" is filled by one of the five types;
 * "sn" is filled by size of the submatrix.
 *
 * If "dim" is supplied and the rest of the working array is incomplete, "lwork" and "work" will be respectively
 * filled by the expected length of the double and integer working arrays
 *
 */
extern spint speigs_factorize( spint *Ap, spint *Ai, double *Ax, spint *dim, spint *aiwork,
                               double *awork, spint *type, spint *sn, spint *iwork, spint *liwork,
                               double *work, spint *lwork, double *evals, double *evecs,
                               spint *rank, double tol ) {

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
        speig_get_factorize_space(dim, sn, type, liwork, lwork);
        return retcode;
    }
    
    spint liwork2, lwork2;
    speig_get_factorize_space(dim, sn, type, &liwork2, &lwork2);
    
    if ( *liwork < liwork2 || *lwork < lwork2) {
        sperr("Insufficient space for factorization phase. \n");
        retcode = SP_EIGS_ERR; return retcode;
    }
    
    spint *p = Ap, *i = Ai, n = *dim; double *x = Ax;
    memset(evals, 0, sizeof(double) * n);
    memset(evecs, 0, sizeof(double) * n * n);
    
    /* Begin factorize */
    if ( *type >= 6 && *type <= -1 ) {
        sperr("Unknown matrix type. \n");
        retcode = SP_EIGS_ERR;
    } else {
        retcode = speig_routines[*type] (p, i, x, n, aiwork, awork, sn, iwork,
                                         liwork, work, lwork, evals, evecs, rank, tol);
    }

    return retcode;
}
