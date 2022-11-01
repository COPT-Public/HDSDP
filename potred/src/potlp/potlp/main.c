#include <stdio.h>

#include "another_lp_solver.h"
#include "data.h"

int main(int argc, const char * argv[]) {
    
    int retcode = RETCODE_OK;
    potlp_solver *potlp = NULL;
    retcode = POT_FNAME(LPSolverCreate)(&potlp);
    
    if ( retcode != RETCODE_OK ) {
        retcode = RETCODE_FAILED;
        goto exit_cleanup;
    }
    
    retcode = POT_FNAME(LPSolverInit)(potlp, nCol, nRow);
    if ( retcode != RETCODE_OK ) {
        retcode = RETCODE_FAILED;
        goto exit_cleanup;
    }
    
    retcode = POT_FNAME(LPSolverSetData)(potlp, Ap, Ai, Ax, obj, rhs);
    if ( retcode != RETCODE_OK ) {
        retcode = RETCODE_FAILED;
        goto exit_cleanup;
    }
    
    potlp->intParams[INT_PARAM_MAXITER] = 1000000;
    potlp->intParams[INT_PARAM_MAXRUIZITER] = 1000;
    potlp->intParams[INT_PARAM_CURVATURE] = 1;
    potlp->intParams[INT_PARAM_COEFSCALE] = 0;
    potlp->dblParams[DBL_PARAM_COMPFOCUS] = 10.0;
    potlp->dblParams[DBL_PARAM_RELOPTTOL] = 1e-10;
    potlp->dblParams[DBL_PARAM_RELFEASTOL] = 1e-10;
    
    POT_FNAME(LPSolverParamsPrint)(potlp);
    retcode = POT_FNAME(LPSolverOptimize)(potlp);
    if ( retcode != RETCODE_OK ) {
        retcode = RETCODE_FAILED;
        goto exit_cleanup;
    }
    
exit_cleanup:
    POT_FNAME(LPSolverDestroy)(&potlp);
    return retcode;
}

/*
 [ 0   A  0  0 -b
   A'  0 -I  0  c
   b' -c' 0 -1  0 ]
 */
