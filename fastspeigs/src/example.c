/** @file example.c
 *  @brief The example for SPEIGS package
 *
 * A little example that demonstrates how to use SPEIGS
 *
 *  @author Wenzhi Gao, Shanghai University of Finance and Economics
 *  @date Aug, 25th, 2022
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "example.h"
#include "speigs.h"

#ifdef sperr
#undef sperr
#define sperr(x) printf(x); return SP_EIGS_ERR;
#endif

static void print_mtype( spint mtype ) {
    
    switch ( mtype ) {
        case MATRIX_TYPE_ZERO:    printf("Matrix type: Zero. \n"); break;
        case MATRIX_TYPE_DIAG:    printf("Matrix type: Diag. \n"); break;
        case MATRIX_TYPE_TWOTWO:  printf("Matrix type: Two.  \n"); break;
        case MATRIX_TYPE_RANKONE: printf("Matrix type: Rank1. \n"); break;
        case MATRIX_TYPE_SPARSE:  printf("Matrix type: Sparse. \n"); break;
        case MATRIX_TYPE_GENERAL: printf("Matrix type: General. \n"); break;
        default: printf("Unknown matrix type. \n"); break;
    }
    return;
}

static void print_evd( spint n, double *e, double *V ) {
    
    spint i, j; printf("Eigen-values: \n");
    for ( i = 0; i < n; ++i ) {
        printf("e["id"] = %5.3e \n", i, e[i]);
    }
    
    printf("Eigen-vectors: \n");
    for ( i = 0; i < n; ++i ) {
        for ( j = 0; j < n; ++j ) {
            printf("%5.3e ", V[j * n + i]);
        }
        printf("\n");
    }
}

/** @brief Test an example of matrix
 *
 * Sample usage of SPEIGS routine
 */
spint test_matrix( spint n, spint *Ap, spint *Ai, double *Ax ) {
    
    spint retcode = SP_EIGS_OK, lwork, liwork, *aiwork, *fiwork, mtype, sn, rank;
    double gthresh = 0.8, tol = 1e-10, *awork, *fwork;
    
    
    printf("Analysis phase starts. "
                   "Getting memory requirement. \n");
    
    retcode = speigs_analyze(NULL, NULL, NULL, &n, NULL, &liwork,
                                 NULL, &lwork, NULL, NULL, 0.0, 0.0);
    
    if ( retcode != SP_EIGS_OK ) {
        sperr("Analysis phase failed. \n");
    }
    
    printf("Memory needed for analysis: "
              "integer array length "id", double array length "id". \n",
              liwork, lwork);
    
    aiwork = (spint *) calloc(liwork, sizeof(spint));
    awork = (double *) calloc(lwork, sizeof(double));
    
    if ( !aiwork || !awork ) {
        sperr("Memory allocation failed in analysis. \n");
    }
    
    printf("SPEIGS analysis starts. \n");
    
    retcode = speigs_analyze(Ap, Ai, Ax, &n, aiwork, &liwork, awork,
                                 &lwork, &mtype, &sn, tol, gthresh);
    if ( retcode != SP_EIGS_OK ) {
        sperr("Analysis phase failed. \n");
    }
    
    printf("Analysis phase completes. \n");
    print_mtype(mtype);
    printf("Submatrix size: "id" \n", sn);
    
    printf("Factorization phase starts. "
              "Getting memory requirement. \n");
    
    retcode = speigs_factorize(NULL, NULL, NULL, &n, aiwork, awork,
                               &mtype, &sn, NULL, &liwork, NULL, &lwork,
                               NULL, NULL, NULL, 0.0);
    
    if ( retcode != SP_EIGS_OK ) {
        sperr("Factorization phase failed. \n");
    }
    
    printf("Memory needed for factorization: "
                      "integer array length "id", double array length "id". \n",
                      liwork, lwork);
    
    fiwork = (spint *) calloc(liwork, sizeof(spint));
    fwork = (double *) calloc(lwork, sizeof(double));
    
    if ( !fiwork || !fwork ) {
        sperr("Memory allocation failed in factorization. \n");
    }
    
    printf("SPEIGS factorization starts. \n");
    
    retcode = speigs_factorize(Ap, Ai, Ax, &n, aiwork, awork,
                                   &mtype, &sn, fiwork, &liwork,
                                   fwork, &lwork, e, V, &rank, tol);
    
    if ( retcode != SP_EIGS_OK ) {
        sperr("Factorization phase failed. \n");
    }
        
    printf("Factorization phase completes. "
              "rank(A) = "id" by tol = %5.2e \n", rank, tol);
    printf("SPEIGS completes. Cleaning up. \n");
    
    free(aiwork); aiwork = NULL; free(awork); awork = NULL;
    free(fiwork); fiwork = NULL; free(fwork); fwork = NULL;
    print_evd(n, e, V);
        
    return retcode;
}


/** @brief Main function
 *
 *  Test the SPEIG routines
 */
int main() {
    
    int retcode = SP_EIGS_OK;
    printf("--------------- Test begins ---------------\n");
    
    retcode = test_matrix(n, Ap0, Ai0, Ax0);
    
    if ( retcode != SP_EIGS_OK ) {
        printf("Error on line %d. \n", __LINE__);
        return retcode;
    }
    
    printf("----------------- Passed -----------------\n");
    
    retcode = test_matrix(n, Ap1, Ai1, Ax1);
    
    if ( retcode != SP_EIGS_OK ) {
        printf("Error on line %d. \n", __LINE__);
        return retcode;
    }
    
    printf("----------------- Passed -----------------\n");
    
    retcode = test_matrix(n, Ap2, Ai2, Ax2);
    
    if ( retcode != SP_EIGS_OK ) {
        printf("Error on line %d. \n", __LINE__);
        return retcode;
    }
    
    printf("----------------- Passed -----------------\n");
    
    retcode = test_matrix(n, Ap3, Ai3, Ax3);
    
    if ( retcode != SP_EIGS_OK ) {
        printf("Error on line %d. \n", __LINE__);
        return retcode;
    }
    
    printf("----------------- Passed -----------------\n");
    
    retcode = test_matrix(n, Ap4, Ai4, Ax4);
    
    if ( retcode != SP_EIGS_OK ) {
        printf("Error on line %d. \n", __LINE__);
        return retcode;
    }
    
    printf("----------------- Passed -----------------\n");
    
    retcode = test_matrix(n, Ap5, Ai5, Ax5);
    
    if ( retcode != SP_EIGS_OK ) {
        printf("Error on line %d. \n", __LINE__);
        return retcode;
    }
    
    printf("---------------- Test ends ----------------\n");
    return SP_EIGS_OK;
}
