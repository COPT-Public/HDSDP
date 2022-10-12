#include <stdio.h>

#include "lp_solver.h"
#include "data_new.h"

int main(int argc, const char * argv[]) {
    
    int retcode = RETCODE_OK;
    potlp_solver *potlp = NULL;
    retcode = LPSolverCreate(&potlp);
    
    if ( retcode != RETCODE_OK ) {
        retcode = RETCODE_FAILED;
    }
    
    retcode = LPSolverInit(potlp, nCol, nRow);
    retcode = LPSolverSetData(potlp, Ap, Ai, Ax, obj, rhs);
    potlp->intParams[INT_PARAM_MAXITER] = 2000000;
    potlp->intParams[INT_PARAM_MAXRUIZITER] = 100;
    retcode = LPSolverOptimize(potlp);
    
exit_cleanup:
    LPSolverDestroy(&potlp);
    return retcode;
}

/*
 [ 0   A  0  0 -b
   A'  0 -I  0  c
   b' -c' 0 -1  0 ]
 */
