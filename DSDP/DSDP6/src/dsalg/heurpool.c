#include "heurpool.h"
// Implement the Heuristics

extern void DSDP_HEURS( adjustSolverParams )
( HSDSolver *dsdpSolver, double largeblock ) {
    
    double ds;   DSDPGetStats(&dsdpSolver->dsdpStats, STAT_NUM_DENSE_MAT, &ds);
    double sps;  DSDPGetStats(&dsdpSolver->dsdpStats, STAT_NUM_SPARSE_MAT, &sps);
    double zero; DSDPGetStats(&dsdpSolver->dsdpStats, STAT_NUM_ZERO_MAT, &zero);
    double r1;   DSDPGetStats(&dsdpSolver->dsdpStats, STAT_NUM_RONE_MAT, &r1);
    double nopint, nodint, impX, impyub, impylb, tmp;
    DSDPGetStats(&dsdpSolver->dsdpStats, STAT_NO_PINTERIOR, &nopint);
    DSDPGetStats(&dsdpSolver->dsdpStats, STAT_NO_DINTERIOR, &nodint);
    DSDPGetStats(&dsdpSolver->dsdpStats, STAT_IMP_BOUNDX, &impX);
    DSDPGetStats(&dsdpSolver->dsdpStats, STAT_IMP_UBOUNDY, &impyub);
    DSDPGetStats(&dsdpSolver->dsdpStats, STAT_IMP_LBOUNDY, &impylb);
    
    DSDP_INT hit = FALSE;
    
    // Case 1: Large block
    if (dsdpSolver->nBlock > 100) {
        DSDPSetIntParam(dsdpSolver, INT_PARAM_ACORRECTOR, 6);
        DSDPSetIntParam(dsdpSolver, INT_PARAM_BCORRECTOR, 0);
        DSDPSetDblParam(dsdpSolver, DBL_PARAM_INIT_BETA, 1.0);
        dsdpSolver->pObjVal = 1e+10; hit = TRUE;
    }
    
    // Case 2: Extremely dense
    if (ds > 0.7 * dsdpSolver->nBlock * dsdpSolver->m && !hit) {
        DSDPSetDblParam(dsdpSolver, DBL_PARAM_INIT_BETA, 1.0);
        DSDPSetIntParam(dsdpSolver, INT_PARAM_ACORRECTOR, 4);
        DSDPSetIntParam(dsdpSolver, INT_PARAM_GOLDSEARCH, FALSE);
        dsdpSolver->ybound = 1e+04;
    }
    
    // Case 3 Implicit trace bound
    if (impX) {
        dsdpSolver->pObjVal = 1e+05;
        // DSDPSetIntParam(dsdpSolver, INT_PARAM_GOLDSEARCH, FALSE);
        DSDPSetDblParam(dsdpSolver, DBL_PARAM_INIT_BETA, 1000.0);
        DSDPSetDblParam(dsdpSolver, DBL_PARAM_BOUND_X, impX);
        DSDPSetIntParam(dsdpSolver, INT_PARAM_GOLDSEARCH, FALSE);
        DSDPSetDblParam(dsdpSolver, DBL_PARAM_RHO, 5.0);
        DSDPSetDblParam(dsdpSolver, DBL_PARAM_RHON, 5.0);
    }
    
    // Case 4 No primal interior
    if (nopint) {
        dsdpSolver->ybound = 1e+04;
        DSDPSetDblParam(dsdpSolver, DBL_PARAM_INIT_BETA, 1000.0);
        if (impX) {
            if (ds > 0) {
                dsdpSolver->pObjVal = 1e+10;
                DSDPSetDblParam(dsdpSolver, DBL_PARAM_RHO, 3.0);
                DSDPSetDblParam(dsdpSolver, DBL_PARAM_RHON, 3.0);
            }
        }
    }
    
    // Case 5 Implicit dual bound
    if (impyub && impylb) {
        tmp = fabs(impyub); dsdpSolver->ybound = MAX(impyub, tmp);
        tmp = fabs(impylb); dsdpSolver->ybound = MAX(impylb, tmp);
        dsdpSolver->ybound *= 1.1;
        dsdpSolver->pObjVal = 1e+05;
        if (dsdpSolver->ybound < 1e+03) {
            DSDPSetDblParam(dsdpSolver, DBL_PARAM_INIT_BETA, 1.0);
            DSDPSetIntParam(dsdpSolver, INT_PARAM_BCORRECTOR, 4);
        }
    } else if (impyub || impylb) {
        // One sided bound
        if (!nopint) {
            DSDPSetIntParam(dsdpSolver, INT_PARAM_GOLDSEARCH, FALSE);
            DSDPSetDblParam(dsdpSolver, DBL_PARAM_RHO, 3.0);
            DSDPSetDblParam(dsdpSolver, DBL_PARAM_RHON, 3.0);
            DSDPSetDblParam(dsdpSolver, DBL_PARAM_INIT_BETA, 1e+06);
            DSDPSetIntParam(dsdpSolver, INT_PARAM_ACORRECTOR, 15);
        } else {
            // sensor sparse inc
            DSDPSetDblParam(dsdpSolver, DBL_PARAM_INIT_BETA, 1e+05);
            dsdpSolver->ybound = 1e+07;
            if (zero == 0.0 && sps == 1.0) {
                DSDPSetDblParam(dsdpSolver, DBL_PARAM_INIT_BETA, 1.0);
                DSDPSetIntParam(dsdpSolver, INT_PARAM_ACORRECTOR, 6);
                DSDPSetIntParam(dsdpSolver, INT_PARAM_BCORRECTOR, 0);
                dsdpSolver->ybound = 1e+04;
            }
        }
    }
    
    if (nodint) {
        if (dsdpSolver->lpDim > dsdpSolver->m) {
            dsdpSolver->ybound = 1.0;
            DSDPSetDblParam(dsdpSolver, DBL_PARAM_INIT_BETA, 1.0);
        } else {
            dsdpSolver->ybound = 10.0;
            DSDPSetIntParam(dsdpSolver, INT_PARAM_GOLDSEARCH, FALSE);
        }
        DSDPSetDblParam(dsdpSolver, DBL_PARAM_INIT_BETA, 1.0);
    }
    
    DSDPSetDblParam(dsdpSolver, DBL_PARAM_PRLX_PENTALTY, dsdpSolver->ybound);
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
        dsdpSolver->dperturb += 1e-07 * dsdpSolver->cScaler;
    } else if (tmp >= 50) {
        dsdpSolver->dperturb += 1e-08 * dsdpSolver->cScaler;
    } else if (tmp >= 20) {
        dsdpSolver->dperturb += 1e-10 * dsdpSolver->cScaler;
    } else {
        dsdpSolver->dperturb += 0.0;
    }
    
    dsdpSolver->dperturb = MIN(dsdpSolver->dperturb, 1e-05 );
    
}
