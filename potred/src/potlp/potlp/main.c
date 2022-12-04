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
    
    potlp->intParams[INT_PARAM_MAXITER] = 50000;
    potlp->intParams[INT_PARAM_MAXRUIZITER] = 10;
    potlp->intParams[INT_PARAM_MAXPCITER] = 0;
    potlp->intParams[INT_PARAM_CURVATURE] = 100;
    potlp->intParams[INT_PARAM_CURVINTERVAL] = 0;
    potlp->intParams[INT_PARAM_L2SCALE] = 0;
    potlp->intParams[INT_PARAM_COEFSCALE] = 1;
    potlp->intParams[INT_PARAM_RSCALFREQ] = 100000000;
    potlp->intParams[INT_PARAM_SCALSIMPLEX] = 1;
    potlp->dblParams[DBL_PARAM_RELOPTTOL] = 1e-08;
    potlp->dblParams[DBL_PARAM_RELFEASTOL] = 1e-08;
    
    potlp->dblParams[18] = 1.5;
    potlp->dblParams[19] = 100;
    
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
