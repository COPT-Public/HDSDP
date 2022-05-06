#ifndef heurpool_h
#define heurpool_h

#include "dsdpinitializer.h"

extern void adjustSolverParams( HSDSolver *dsdpSolver, double largeblock ) {
    
    // hamming problems
    if (dsdpSolver->m >= 50 * largeblock && dsdpSolver->m >= 15000) {
        dsdpSolver->pObjVal = 1e+05; dsdpSolver->ybound = 1e+07;
        DSDPSetDblParam(dsdpSolver, DBL_PARAM_INIT_BETA, 1.0);
        DSDPSetDblParam(dsdpSolver, DBL_PARAM_RHO, 5.0);
        DSDPSetDblParam(dsdpSolver, DBL_PARAM_RHON, 5.0);
    }
    
    if (dsdpSolver->m >= 100 * largeblock && dsdpSolver->m >= 50000) {
        dsdpSolver->pObjVal = 1e+05; dsdpSolver->ybound = 1e+05;
        DSDPSetDblParam(dsdpSolver, DBL_PARAM_INIT_BETA, 10.0);
        DSDPSetDblParam(dsdpSolver, DBL_PARAM_RHO, 5.0);
        DSDPSetDblParam(dsdpSolver, DBL_PARAM_RHON, 5.0);
    }
}

#endif /* heurpool_h */
