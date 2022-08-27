/** @file mex_speigs.c
 *  @brief Mexfile entry function for speigs
 *
 * @author Wenzhi Gao, Shanghai University of Finance and Economics
 * @date Aug, 24th, 2022
 *
 */

#include <stdio.h>
#include "speigs.h"

/** @brief Print type of matrix
 *  @param[in] mtype Type of the matrix
 */
static void print_mtype( spint mtype ) {
    
    switch ( mtype ) {
        case MATRIX_TYPE_ZERO:    mexPrintf("Matrix type: Zero. \n"); break;
        case MATRIX_TYPE_DIAG:    mexPrintf("Matrix type: Diag. \n"); break;
        case MATRIX_TYPE_TWOTWO:  mexPrintf("Matrix type: Two.  \n"); break;
        case MATRIX_TYPE_RANKONE: mexPrintf("Matrix type: Rank1. \n"); break;
        case MATRIX_TYPE_SPARSE:  mexPrintf("Matrix type: Sparse. \n"); break;
        case MATRIX_TYPE_GENERAL: mexPrintf("Matrix type: General. \n"); break;
        default: sperr("Unknown matrix type. \n"); break;
    }
    return;
}

/* [V, e] = mex_speigs(A, opts); */
#define V    plhs[0]
#define e    plhs[1]
#define A    prhs[0]
#define opts prhs[1]

/** @brief Matlab entry function
 *  @param[in] nlhs Number of left-hand-side parameters
 *  @param[out] plhs Pointers for left-hand-side parameters
 *  @param[in] nrhs Number of right-hand-side parameters
 *  @param[out] prhs Pointers for left-hand-side parameters
 *
 *  Matab entry for [V, e] = mex_speigs(A, opts);
 *  V is a n by r array that gives r eigenvectors and e is all the nonzero eigen-values.
 *  opts.gthresh specifies when submatrix permutation is used
 *  opts.tol specifies the criterion to decide if an eigen-value is 0
 *  opts.quiet hides logs during factorization
 */
extern void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {
    
    spint retcode = SP_EIGS_OK, m, n, lwork, liwork, *aiwork = NULL, *fiwork = NULL, mtype, sn, rank;
    double gthresh = 0.8, tol = 1e-10, quiet = 1.0, *awork = NULL, *fwork = NULL, *Vp, *ep;
    spint *Ap = NULL, *Ai = NULL; double *Ax = NULL;
    mxArray *param = NULL;
    
    /* ------------------- Parameter Extraction ------------------- */
    if ( nrhs <= 0 || nrhs >= 3 ) {
        sperr("Invalid number of entries. \n");
    }
    
    if ( nlhs <= 1 ) {
        sperr("Invalid number of outputs. \n");
    }
    
    if ( nrhs == 2 ) {
        param = mxGetField(opts, 0, "gthresh");
        if ( param ) { gthresh = (double) (*mxGetPr(param)); }
        param = mxGetField(opts, 0, "tol");
        if ( param ) { tol = (double) (*mxGetPr(param)); }
        param = mxGetField(opts, 0, "quiet");
        if ( param ) { quiet = (double) (*mxGetPr(param)); }
    }
    
    if ( gthresh <= 0.0 || tol <= 0.0 ) { sperr("Invalid thresholds. \n"); }
    if ( !quiet ) {
        mexPrintf("\n");
        mexPrintf("SPEIGS Parameter is invoked with"
                  " gthresh = %3.3f and tol = %5.2e \n", gthresh, tol);
    }
    
    if ( !mxIsSparse(A) ) { sperr("A must be sparse. \n"); }
    
    m = mxGetM(A); n = mxGetN(A);
    if ( m != n ) { sperr("A is not a square matrix. \n"); }
    
    /* Extract matrix data */
    Ap = mxGetJc(A); Ai = mxGetIr(A); Ax = mxGetPr(A);
    
    if ( !Ap || !Ai || !Ax ) {
        sperr("Failed to extract A. \n");
    }
    
    if ( !quiet ) {
        mexPrintf("Dimension of problem %lld, Nnz %lld. \n", m, Ap[n]);
    }
    
    /* ------------------- Analysis Phase ------------------- */
    if ( !quiet ) {
        mexPrintf("Analysis phase starts. "
                  "Getting memory requirement. \n");
    }
    retcode = speigs_analyze(NULL, NULL, NULL, &n, NULL, &liwork,
                             NULL, &lwork, NULL, NULL, 0.0, 0.0);
    if ( retcode != SP_EIGS_OK ) {
        sperr("Analysis phase failed. \n");
    }
    
    if ( !quiet ) {
        mexPrintf("Memory needed for analysis: "
                  "integer array length %lld, double array length %lld. \n",
                  liwork, lwork);
    }
    
    aiwork = (spint *) mxCalloc(liwork, sizeof(spint));
    awork = (double *) mxCalloc(lwork, sizeof(double));
    
    if ( !aiwork || !awork ) {
        sperr("Memory allocation failed in analysis. \n");
    }
    
    if ( !quiet ) {
        mexPrintf("SPEIGS analysis starts. \n");
    }
    retcode = speigs_analyze(Ap, Ai, Ax, &n, aiwork, &liwork, awork,
                             &lwork, &mtype, &sn, tol, gthresh);
    if ( retcode != SP_EIGS_OK ) {
        sperr("Analysis phase failed. \n");
    }
    
    if ( !quiet ) {
        mexPrintf("Analysis phase completes. \n");
        print_mtype(mtype);
        mexPrintf("Submatrix size: %lld \n", sn);
    }
    
    /* ------------------- Factorization Phase ------------------- */
    if ( !quiet ) {
        mexPrintf("Factorization phase starts. "
                  "Getting memory requirement. \n");
    }
    
    retcode = speigs_factorize(NULL, NULL, NULL, &n, aiwork, awork,
                               &mtype, &sn, NULL, &liwork, NULL, &lwork,
                               NULL, NULL, NULL, 0.0);
    if ( retcode != SP_EIGS_OK ) {
        sperr("Factorization phase failed. \n");
    }
    
    if ( !quiet ) {
        mexPrintf("Memory needed for factorization: "
                  "integer array length %lld, double array length %lld. \n",
                  liwork, lwork);
    }
    
    if ( liwork > 0 || fiwork > 0 ) {
        fiwork = (spint *) mxCalloc(liwork, sizeof(spint));
        fwork = (double *) mxCalloc(lwork, sizeof(double));
        if ( !fiwork || !fwork ) {
            sperr("Memory allocation failed in factorization. \n");
        }
    }
    
    if ( !quiet ) {
        mexPrintf("SPEIGS factorization starts. \n");
    }
    
    V = mxCreateDoubleMatrix(n, n, mxREAL);
    e = mxCreateDoubleMatrix(n, 1, mxREAL);
    Vp = mxGetPr(V); ep = mxGetPr(e);
    
    retcode = speigs_factorize(Ap, Ai, Ax, &n, aiwork, awork,
                               &mtype, &sn, fiwork, &liwork,
                               fwork, &lwork, ep, Vp, &rank, tol);
    
    if ( retcode != SP_EIGS_OK ) {
        sperr("Factorization phase failed. \n");
    }
    
    if ( !quiet ) {
        mexPrintf("Factorization phase completes. "
                  "rank(A) = %lld by tol = %5.2e \n", (long long int) rank, tol);
        mexPrintf("SPEIGS completes. Cleaning up. \n\n");
    }
    
    /* ------------------- Clean up ------------------- */
    mxFree(aiwork); mxFree(awork);
    mxFree(fiwork); mxFree(fwork);
    
    return;
}
