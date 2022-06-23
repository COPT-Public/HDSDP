#include "heurpool.h"
// Implement the Heuristics

extern void DSDP_HEURS( adjustSolverParams )
( HSDSolver *dsdpSolver, double largeblock ) {
    
    double ds;   DSDPGetStats(&dsdpSolver->dsdpStats, STAT_NUM_DENSE_MAT, &ds);
    double sps;  DSDPGetStats(&dsdpSolver->dsdpStats, STAT_NUM_SPARSE_MAT, &sps);
    double zero; DSDPGetStats(&dsdpSolver->dsdpStats, STAT_NUM_ZERO_MAT, &zero);
    double r1;   DSDPGetStats(&dsdpSolver->dsdpStats, STAT_NUM_RANKONE_MAT, &r1);
    double nlp =  dsdpSolver->lpDim; DSDP_INT inttmp;
    
    DSDP_INT hit = FALSE;
    
    // Case 1: large block
    if (dsdpSolver->nBlock > 100) {
        DSDPSetIntParam(dsdpSolver, INT_PARAM_ACORRECTOR, 6);
        DSDPSetIntParam(dsdpSolver, INT_PARAM_BCORRECTOR, 0);
        DSDPSetDblParam(dsdpSolver, DBL_PARAM_INIT_BETA, 1.0);
        dsdpSolver->pObjVal = 1e+10; hit = TRUE;
    }
    
    // Case 2: Extremely dense
    if (ds > 0.7 * dsdpSolver->nBlock * dsdpSolver->m && !hit) {
        printf("| Hit prob \n"); hit = TRUE;
        DSDPSetDblParam(dsdpSolver, DBL_PARAM_INIT_BETA, 1.0);
        DSDPSetIntParam(dsdpSolver, INT_PARAM_ACORRECTOR, 4);
        DSDPSetIntParam(dsdpSolver, INT_PARAM_GOLDSEARCH, FALSE);
        dsdpSolver->ybound = 1e+04;
    }
    
    // hamming problems: Case: tr(X) = alpha
    if (zero == 0.0 && dsdpSolver->m >= 50 * largeblock && dsdpSolver->m >= 15000 && !hit) {
        printf("| Hit hamming \n"); hit = TRUE;
        dsdpSolver->pObjVal = 1e+05; dsdpSolver->ybound = 1e+07;
        DSDPSetDblParam(dsdpSolver, DBL_PARAM_INIT_BETA, 1.0);
        DSDPSetDblParam(dsdpSolver, DBL_PARAM_RHO, 3.0);
        DSDPSetDblParam(dsdpSolver, DBL_PARAM_RHON, 3.0);
    }
    
    // theta: Case: tr(X) = alpha
    if (zero == 0.0 && r1 == 1.0 && ds == 0.0 && sps >= 2000 &&
        dsdpSolver->m >= dsdpSolver->n * 25 && !hit) {
        printf("| Hit theta \n"); hit = TRUE;
        DSDPSetDblParam(dsdpSolver, DBL_PARAM_INIT_BETA, 10.0);
        DSDPSetIntParam(dsdpSolver, INT_PARAM_GOLDSEARCH, TRUE);
        dsdpSolver->pObjVal = 1e+02;
    }
    
    // Gxx: Case: tr(X) = alpha
    if (zero == 0.0 && r1 == dsdpSolver->m && sps == 1.0 && !hit &&
        (dsdpSolver->n == dsdpSolver->m - 1 || dsdpSolver->n == dsdpSolver->m) &&
        dsdpSolver->m <= 10000) {
        getIntParam(dsdpSolver->param, INT_PARAM_GOLDSEARCH, &inttmp);
        printf("| Hit Gxx \n"); hit = TRUE;
        DSDPSetIntParam(dsdpSolver, INT_PARAM_GOLDSEARCH, FALSE);
        DSDPSetDblParam(dsdpSolver, DBL_PARAM_INIT_BETA, 1e+03);
        DSDPSetDblParam(dsdpSolver, DBL_PARAM_BOUND_X, 1e+04);
        DSDPSetDblParam(dsdpSolver, DBL_PARAM_RHO, 5.0);
        DSDPSetDblParam(dsdpSolver, DBL_PARAM_RHON, 5.0);
        
        if (inttmp) {
            dsdpSolver->pObjVal = 1e+05;
            dsdpSolver->ybound = 1e+05;
        } else {
            dsdpSolver->pObjVal = 1e+08;
            dsdpSolver->ybound = 1e+04;
        }
    }
    
    // Vibra and buck: one sided y bound
    if (nlp == dsdpSolver->m && dsdpSolver->n != dsdpSolver->m) {
        printf("| Hit ViBu \n"); hit = TRUE;
        DSDPSetDblParam(dsdpSolver, DBL_PARAM_RHO, 3.0);
        DSDPSetDblParam(dsdpSolver, DBL_PARAM_RHON, 3.0);
        DSDPSetDblParam(dsdpSolver, DBL_PARAM_INIT_BETA, 1e+06);
        DSDPSetIntParam(dsdpSolver, INT_PARAM_ACORRECTOR, 12);
    }
    
    // inc_xxx Case: tr(E * X) = 0.0
    if (zero == 0.0 && sps == 1.0 && ds == 0.0 &&
        (fabs(nlp / dsdpSolver->n - 4.0) < 0.5) && !hit) {
        printf("| Hit Inc \n"); hit = TRUE;
        DSDPSetDblParam(dsdpSolver, DBL_PARAM_INIT_BETA, 1.0);
        dsdpSolver->ybound = 1e+04;
    }
    
    // swissroll Case: tr(E * X) = 0.0
    if ((fabs((double) dsdpSolver->m / dsdpSolver->n - 4) < 0.5)
        && zero == 0.0 && ds ==0.0 && r1 == dsdpSolver->m && !hit) {
        printf("| Hit Swiss \n"); hit = TRUE;
        DSDPSetDblParam(dsdpSolver, DBL_PARAM_INIT_BETA, 1.0);
        dsdpSolver->ybound = 1e+04;
    }
    
    // body Case: tr(E * X) = 0.0
    if (zero == 0.0 && sps == 0.0 && ds == 1.0 && r1 == dsdpSolver->m && !hit) {
        printf("| Hit body \n"); hit = TRUE;
        DSDPSetDblParam(dsdpSolver, DBL_PARAM_INIT_BETA, 1.0);
        dsdpSolver->ybound = 1e+04;
    }
    
    // shmup: two-sided implied bound on y
    if (nlp >= dsdpSolver->m && ds == 0.0 && r1 == 0.0 && !hit) {
        printf("| Hit shmup \n"); hit = TRUE;
        DSDPSetDblParam(dsdpSolver, DBL_PARAM_INIT_BETA, 1.0);
        DSDPSetIntParam(dsdpSolver, INT_PARAM_BCORRECTOR, 4);
        dsdpSolver->pObjVal = 1e+05; dsdpSolver->ybound = 100;
    }
    
    // Large LP cone. rabmo and reimer
    if (nlp > 30 * dsdpSolver->n && nlp > dsdpSolver->m && !hit) {
        dsdpSolver->ybound = 1.0;
        printf("| Hit large LPcone \n"); hit = TRUE;
        DSDPSetDblParam(dsdpSolver, DBL_PARAM_INIT_BETA, 1.0);
    }
    
    // 1024: Case: tr(X) = alpha
    getIntParam(dsdpSolver->param, INT_PARAM_GOLDSEARCH, &inttmp);
    if (!inttmp && sps == dsdpSolver->m && r1 == 1.0 &&
        zero == 0.0 && ds == 0.0 && dsdpSolver->m >= dsdpSolver->n * 8 &&
        dsdpSolver->nBlock == 1.0 && !hit) {
        printf("| Hit 1024 \n");
        DSDPSetDblParam(dsdpSolver, DBL_PARAM_INIT_BETA, 1.0);
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
