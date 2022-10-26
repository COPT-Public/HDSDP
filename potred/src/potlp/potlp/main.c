#include <stdio.h>

#include "lp_solver.h"
#include "data_new.h"

int main(int argc, const char * argv[]) {
    
    int retcode = RETCODE_OK;
    potlp_solver *potlp = NULL;
    retcode = LPSolverCreate(&potlp);
    
    if ( retcode != RETCODE_OK ) {
        retcode = RETCODE_FAILED;
        goto exit_cleanup;
    }
    
    retcode = LPSolverInit(potlp, nCol, nRow);
    if ( retcode != RETCODE_OK ) {
        retcode = RETCODE_FAILED;
        goto exit_cleanup;
    }
    
    retcode = LPSolverSetData(potlp, Ap, Ai, Ax, obj, rhs);
    if ( retcode != RETCODE_OK ) {
        retcode = RETCODE_FAILED;
        goto exit_cleanup;
    }
    
    potlp->intParams[INT_PARAM_MAXITER] = 1000000;
    potlp->intParams[INT_PARAM_MAXRUIZITER] = 1000;
    
    retcode = LPSolverOptimize(potlp);
    if ( retcode != RETCODE_OK ) {
        retcode = RETCODE_FAILED;
        goto exit_cleanup;
    }
    
exit_cleanup:
    LPSolverDestroy(&potlp);
    return retcode;
}

/*
 [ 0   A  0  0 -b
   A'  0 -I  0  c
   b' -c' 0 -1  0 ]
 */
