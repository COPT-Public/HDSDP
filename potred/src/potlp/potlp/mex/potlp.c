/** @file potlp.c
 *  @brief Mexfile entry function for potential reduction LP
 *
 * @author Wenzhi Gao, Shanghai University of Finance and Economics
 * @date Oct, 26th, 2022
 *
 */

#include "lp_solver.h"

/* [x, y, s] = potlp(A, b, c, params); */
#define x      plhs[0]
#define y      plhs[1]
#define s      plhs[2]
#define A      prhs[0]
#define b      prhs[1]
#define c      prhs[2]
#define params prhs[3]

/** @brief Matlab entry function
 *  @param[in] nlhs Number of left-hand-side parameters
 *  @param[out] plhs Pointers for left-hand-side parameters
 *  @param[in] nrhs Number of right-hand-side parameters
 *  @param[out] prhs Pointers for left-hand-side parameters
 *
 *  Matab entry for [x, y, s, time] = potlp(A, b, c, params);
 */
extern void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {
    
    pot_int retcode = RETCODE_OK;
    
    mwSize nCol = 0;
    mwSize nRow = 0;
    mwSize *colMatBeg = NULL;
    mwSize *colMatIdx = NULL;
    double  *colMatElem = NULL;
    double  *colObj = NULL;
    double  *rowRhs = NULL;
    mxArray *param = NULL;
    
    if ( nrhs <= 0 || nrhs >= 5 ) {
        mexErrMsgTxt("Invalid number of entries. \n");
    }
    
    if ( nlhs <= 2 ) {
        mexErrMsgTxt("Invalid number of outputs. \n");
    }
    
    /* Get A */
    if ( !mxIsSparse(A) ) {
        mexErrMsgTxt("A must be sparse. \n");
    }
    
    nRow = mxGetM(A);
    nCol = mxGetN(A);
    colMatBeg = mxGetJc(A);
    colMatIdx = mxGetIr(A);
    colMatElem = mxGetPr(A);
    
    if ( !colMatBeg || !colMatIdx || !colMatElem ) {
        mexErrMsgTxt("Failed to extract A. \n");
    }
    
    /* Get b */
    if (!mxIsDouble(b)) {
        mexErrMsgTxt("b must be a double array \n");
    }
    rowRhs = mxGetPr(b);
    
    /* Get c */
    if (!mxIsDouble(c)) {
        mexErrMsgTxt("c must be a double array \n");
    }
    colObj = mxGetPr(c);
    
    /* Get Parameters */
    int maxIter = 1000;
    int maxRuizIter = 50;
    int coefScal = 0;
    
    double relFeasTol = 1e-04;
    double relOptTol = 1e-04;
    double maxTime = 600.0;
    double compFocus = 10.0;
    
    if ( nrhs == 4 ) {
        
        param = mxGetField(params, 0, "maxIter");
        if ( param ) {
            maxIter = (int) (*mxGetPr(param));
        }
        
        param = mxGetField(params, 0, "maxRuizIter");
        if ( param ) {
            maxRuizIter = (int) (*mxGetPr(param));
        }
        
        param = mxGetField(params, 0, "coefScal");
        if ( param ) {
            coefScal = (int) (*mxGetPr(param));
        }
        
        param = mxGetField(params, 0, "relFeasTol");
        if ( param ) {
            relFeasTol = (double) (*mxGetPr(param));
        }
        
        param = mxGetField(params, 0, "relOptTol");
        if ( param ) {
            relOptTol = (double) (*mxGetPr(param));
        }
        
        param = mxGetField(params, 0, "maxTime");
        if ( param ) {
            maxTime = (double) (*mxGetPr(param));
        }
        
        param = mxGetField(params, 0, "compFocus");
        if ( param ) {
            compFocus = (double) (*mxGetPr(param));
        }
    }
    
    /* Ready to solve */
    potlp_solver *potlp = NULL;
    retcode = LPSolverCreate(&potlp);
    
    if ( retcode != RETCODE_OK ) {
        retcode = RETCODE_FAILED;
        printf("Error on line %d \n", __LINE__);
        goto exit_cleanup;
    }
    
    retcode = LPSolverInit(potlp, nCol, nRow);
    if ( retcode != RETCODE_OK ) {
        retcode = RETCODE_FAILED;
        printf("Error on line %d \n", __LINE__);
        goto exit_cleanup;
    }
    
    /* Type conversion */
    pot_int *Ap = NULL;
    pot_int *Ai = NULL;
    Ap = mxCalloc(nCol + 1, sizeof(pot_int));
    Ai = mxCalloc(colMatBeg[nCol], sizeof(pot_int));
    
    if ( !Ap || !Ai ) {
        retcode = RETCODE_FAILED;
        printf("Error on line %d \n", __LINE__);
        goto exit_cleanup;
    }
    
    for ( int i = 0; i < nCol; ++i ) {
        Ap[i + 1] = colMatBeg[i + 1];
    }
    
    for ( int i = 0; i < Ap[nCol]; ++i ) {
        Ai[i] = colMatIdx[i];
    }
    
    retcode = LPSolverSetData(potlp, Ap, Ai, colMatElem, colObj, rowRhs);
    
    mxFree(Ap); Ap = NULL;
    mxFree(Ai); Ai = NULL;
    
    if ( retcode != RETCODE_OK ) {
        retcode = RETCODE_FAILED;
        printf("Error on line %d \n", __LINE__);
        goto exit_cleanup;
    }
    
    /* Set parameters */
    potlp->intParams[INT_PARAM_MAXITER] = maxIter;
    potlp->intParams[INT_PARAM_MAXRUIZITER] = maxRuizIter;
    potlp->intParams[INT_PARAM_COEFSCALE] = coefScal;
    
    potlp->dblParams[DBL_PARAM_RELFEASTOL] = relFeasTol;
    potlp->dblParams[DBL_PARAM_RELOPTTOL] = relOptTol;
    potlp->dblParams[DBL_PARAM_TIMELIMIT] = maxTime;
    potlp->dblParams[DBL_PARAM_COMPFOCUS] = compFocus;
    
    LPSolverParamsPrint(potlp);

    retcode = LPSolverOptimize(potlp);
    if ( retcode != RETCODE_OK ) {
        retcode = RETCODE_FAILED;
        printf("Error on line %d \n", __LINE__);
        goto exit_cleanup;
    }
    
    /* Finished with solve and write solution */
    x = mxCreateDoubleMatrix(nCol, 1, mxREAL);
    y = mxCreateDoubleMatrix(nRow, 1, mxREAL);
    s = mxCreateDoubleMatrix(nCol, 1, mxREAL);
    
    double *xp = mxGetPr(x);
    double *yp = mxGetPr(y);
    double *sp = mxGetPr(s);
    
    if ( !xp || !yp || !sp ) {
        mexErrMsgTxt("Failed to allocate memory for LP solution.\n");
    }
    
    LPSolverGetSolution(potlp, xp, yp, sp);
    
exit_cleanup:
    
    if ( Ap ) {
        mxFree(Ap);
    }
    
    if ( Ai ) {
        mxFree(Ai);
    }
    
    LPSolverDestroy(&potlp);
    return;
}
