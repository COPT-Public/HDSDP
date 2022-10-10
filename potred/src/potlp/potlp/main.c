#include <stdio.h>

#include "lp_solver.h"
#include "data.h"

int main(int argc, const char * argv[]) {
    
    int retcode = RETCODE_OK;
    potlp_solver *potlp = NULL;
    retcode = LPSolverCreate(&potlp);
    
    if ( retcode != RETCODE_OK ) {
        retcode = RETCODE_FAILED;
    }
    
    retcode = LPSolverInit(potlp, nCol, nRow);
    retcode = LPSolverSetData(potlp, Ap, Ai, Ax, obj, rhs);
    potlp->potIterator->intParams[INT_PARAM_MAXITER] = 5000;
    retcode = LPSolverOptimize(potlp);
    
exit_cleanup:
    LPSolverDestroy(&potlp);
    return retcode;
}
