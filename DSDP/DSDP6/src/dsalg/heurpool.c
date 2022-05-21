#include "heurpool.h"
// Implement the Heuristics

extern void DSDP_HEURS( adjustSolverParams )
( HSDSolver *dsdpSolver, double largeblock ) {
    
    double ds;   DSDPGetStats(&dsdpSolver->dsdpStats, STAT_NUM_DENSE_MAT, &ds);
    double sps;  DSDPGetStats(&dsdpSolver->dsdpStats, STAT_NUM_SPARSE_MAT, &sps);
    double zero; DSDPGetStats(&dsdpSolver->dsdpStats, STAT_NUM_ZERO_MAT, &zero);
    double r1;   DSDPGetStats(&dsdpSolver->dsdpStats, STAT_NUM_RANKONE_MAT, &r1);
    
    // hamming problems
    if (dsdpSolver->m >= 50 * largeblock && dsdpSolver->m >= 15000) {
        dsdpSolver->pObjVal = 1e+05; dsdpSolver->ybound = 1e+07;
        DSDPSetDblParam(dsdpSolver, DBL_PARAM_INIT_BETA, 1.0);
        DSDPSetDblParam(dsdpSolver, DBL_PARAM_RHO, 5.0);
        DSDPSetDblParam(dsdpSolver, DBL_PARAM_RHON, 5.0);
    }
    
    if (dsdpSolver->m >= 100 * largeblock && dsdpSolver->m >= 50000) {
        dsdpSolver->pObjVal = 1e+10; dsdpSolver->ybound = 1e+04;
        DSDPSetDblParam(dsdpSolver, DBL_PARAM_INIT_BETA, 1000.0);
    }
    
    // prob_x_x_x problems
    if (ds > 0.7 * dsdpSolver->nBlock * dsdpSolver->m) {
        DSDPSetDblParam(dsdpSolver, DBL_PARAM_INIT_BETA, 1.0);
        DSDPSetIntParam(dsdpSolver, INT_PARAM_ACORRECTOR, 6);
    }
    
    // theta
    if (zero == 0.0 && r1 == 1.0 && ds == 0.0 && sps >= 2000) {
        DSDPSetDblParam(dsdpSolver, DBL_PARAM_INIT_BETA, 10.0);
    }
    
}

extern void DSDP_HEURS( adjustCScaler )
( HSDSolver *dsdpSolver ) {
    double tol; DSDPGetDblParam(dsdpSolver, DBL_PARAM_ABS_OPTTOL, &tol);
    if (dsdpSolver->cScaler > 1e+10) { dsdpSolver->cScaler = 1e+08;
    } else if (dsdpSolver->cScaler > 1e+05) {
        // Just use the coefficient as scaler
    } else if (dsdpSolver->cScaler < 1e-15) {
        dsdpSolver->cScaler = 1.0;
        DSDPStatUpdate(&dsdpSolver->dsdpStats, STAT_PFEAS_PROBLEM, TRUE);
    } else if (dsdpSolver->cScaler < 1e-08) { dsdpSolver->cScaler = 1e-06;
    } else { dsdpSolver->cScaler = 1.0; }
}

extern void DSDP_HEURS( adjPRelaxPenalty )
( HSDSolver *dsdpSolver ) {
    double tmp = vec_infnorm(dsdpSolver->y);
    if (tmp > 1e+07) { tmp = 1.5 * tmp; }
    else if (tmp > 1e+05) { tmp = MIN(tmp * 10, 1e+07); }
    else if (tmp > 1e+03) { tmp = MIN(tmp * 100, 1e+06); }
    else if (tmp > 10) { tmp = MIN(tmp * 200, 1e+05); }
    else { tmp = MIN(tmp * 1000, 1e+04); }
    dsdpSolver->ybound = MIN(dsdpSolver->ybound, tmp);
}

extern void DSDP_HEURS( adjDualPerturb )
( HSDSolver *dsdpSolver ) {
    double tmp, tmp2;
    DSDPGetStats(&dsdpSolver->dsdpStats, STAT_PHASE_A_ITER, &tmp);
    DSDPGetStats(&dsdpSolver->dsdpStats, STAT_ONE_NORM_C, &tmp2);
    if (tmp >= 80) {
        dsdpSolver->dperturb = 1e-07 * dsdpSolver->cScaler;
    } else if (tmp >= 50) {
        dsdpSolver->dperturb = 1e-08 * dsdpSolver->cScaler;
    } else if (tmp >= 20) {
        dsdpSolver->dperturb = 1e-10 * dsdpSolver->cScaler;
    } else {
        dsdpSolver->dperturb = 0.0;
    }
    
    dsdpSolver->dperturb = MIN(dsdpSolver->dperturb, 1e-05 );
    
}
