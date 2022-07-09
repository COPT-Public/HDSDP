#include "dsdpinitializer.h"
#include "dsdpcorrector.h"
#include "dsdputils.h"
#include "dsdplog.h"
#include "heurpool.h"
/* Implement the initialization procedure of DSDP */
static char etype[] = "DSDP Initialization";

static DSDP_INT inity( HSDSolver *dsdpSolver ) {
    
    DSDP_INT retcode = DSDP_RETCODE_OK;
    vec_reset(dsdpSolver->y);
    return retcode;
}

static DSDP_INT initkappatau( HSDSolver *dsdpSolver ) {
    
    DSDP_INT retcode = DSDP_RETCODE_OK;
    dsdpSolver->tau = 1.0;
#ifdef KAPPATAU
    retcode = DSDPGetDblParam(dsdpSolver, DBL_PARAM_INIT_KAPPA, &dsdpSolver->kappa);
#else
    dsdpSolver->kappa = dsdpSolver->mu / dsdpSolver->tau;
#endif
    return retcode;
}

static DSDP_INT initmu( HSDSolver *dsdpSolver ) {
    DSDP_INT retcode = DSDP_RETCODE_OK;
    DSDPGetDblParam(dsdpSolver, DBL_PARAM_INIT_MU, &dsdpSolver->mu);
    return retcode;
}

static DSDP_INT initresi( HSDSolver *dsdpSolver ) {
    // Initialize Ry
    
    /*
     if norm(Ctau, 'fro') == 0
         Ry = - speye(n) * initbeta;
     else
         Ry = - speye(n) * max(norm(Ctau, 'fro'), 10.0) * initbeta;
     end % End if
    */
    
    DSDP_INT retcode = DSDP_RETCODE_OK;
    double beta, Cnrm = 0.0, nrm  = 0.0, tau = dsdpSolver->tau;
    DSDPGetDblParam(dsdpSolver, DBL_PARAM_INIT_BETA, &beta);
    Cnrm = DSDPConic( COPS_GET_C_FNORM )(dsdpSolver);
    dsdpSolver->Ry = (Cnrm == 0) ? (-1.0) : - MAX(Cnrm, 100.0) * tau * beta;
    DSDPGetStats(&dsdpSolver->dsdpStats, STAT_ONE_NORM_B, &nrm);
    dsdpSolver->pinfeas = nrm + 1.0;
    dsdpSolver->iterProgress[ITER_INITIALIZE] = TRUE;
    
    return retcode;
}

static void initparams( HSDSolver *dsdpSolver ) {
    // TODO: Adjust solver parameter based on problem features
    
    double dblparam, stats;
    dsdpSolver->pObjVal = 1e+10;
    
    // Adjust bound on the y variable. Set it to the initial value
    DSDPGetDblParam(dsdpSolver, DBL_PARAM_PRLX_PENTALTY, &dblparam);
    
    // If solving a feasibility problem
    DSDPGetStats(&dsdpSolver->dsdpStats, STAT_PFEAS_PROBLEM, &stats);
    dblparam = (stats) ?  MIN(0.001, dblparam) : dblparam;
    dsdpSolver->pObjVal = (stats) ? 0.0 : dsdpSolver->pObjVal;
    
    // If not using primal relaxation in Phase A
    dsdpSolver->nall = dsdpSolver->n + dsdpSolver->lpDim + dsdpSolver->m * 2;
    dsdpSolver->ybound = dblparam;
    
    // Corrector
    DSDP_INT m = dsdpSolver->m, nusercorr, ncorrA = 0;
    double largeblock;
    DSDPGetIntParam(dsdpSolver, INT_PARAM_BCORRECTOR, &nusercorr);
    DSDPGetStats(&dsdpSolver->dsdpStats, STAT_LARGEST_BLOCK, &largeblock);
    
    ncorrA = (m - 2) / largeblock;
    if (dsdpSolver->nall < 100 && ncorrA == 0) {
        ncorrA = 1;
    }
    if (ncorrA >= 1) { ncorrA ++; }
    ncorrA *= ncorrA;
    if (m < 2000 && ncorrA > 10) { ncorrA = 10; }
    nusercorr = ncorrA;
    
    if (TRUE) {
        if (largeblock >= 5 * m) {
            nusercorr = MIN(nusercorr, 0); ncorrA = 2;
        } else if (largeblock >= m) {
            nusercorr = MIN(nusercorr, 2); ncorrA = 4;
        } else {
            ncorrA = 6;
        }
        
        if (m > 20 * largeblock) {
            nusercorr = MAX(nusercorr, 12); ncorrA = 12;
        } else if (m > 5 * largeblock) {
            nusercorr = MAX(nusercorr, 10); ncorrA = 10;
        } else if (m > 2 * largeblock) {
            nusercorr = MAX(nusercorr, 8); ncorrA = 8;
        }
        nusercorr = MIN(nusercorr, 12);
    }
    ncorrA = MAX(ncorrA, 2);
    
    DSDPSetIntParam(dsdpSolver, INT_PARAM_ACORRECTOR, ncorrA);
    DSDPSetIntParam(dsdpSolver, INT_PARAM_BCORRECTOR, nusercorr);
    
    // Golden line search
    DSDPGetStats(&dsdpSolver->dsdpStats, STAT_NNZ_OBJ, &stats);
    if (stats <= 0.01 * dsdpSolver->n * dsdpSolver->n) {
        DSDPSetIntParam(dsdpSolver, INT_PARAM_GOLDSEARCH, TRUE);
    }
    
    // Some other heuristics
    DSDP_HEURS( adjustSolverParams )(dsdpSolver, largeblock);
    
    DSDPGetStats(&dsdpSolver->dsdpStats, STAT_PFEAS_PROBLEM, &stats);
    if (stats) {
        // No search heuristic for feasibility problem
        DSDPSetIntParam(dsdpSolver, INT_PARAM_GOLDSEARCH, FALSE);
    }
    DSDPParamPrint(dsdpSolver->param);
    
    DSDPGetIntParam(dsdpSolver, INT_PARAM_ACORRECTOR, &ncorrA);
    DSDPGetIntParam(dsdpSolver, INT_PARAM_BCORRECTOR, &nusercorr);
    printf("| Corrector A: %d  Corrector B: %d \n", ncorrA, nusercorr);
    
//    dsdpSolver->y->x[0] = -1e+06;
//    DSDPSetDblParam(dsdpSolver, DBL_PARAM_INIT_BETA, 0.001);
//    DSDPSetIntParam(dsdpSolver, INT_PARAM_GOLDSEARCH, FALSE);
    showBeautifulDashlines();
        
}

extern DSDP_INT dsdpInitializeA( HSDSolver *dsdpSolver ) {
    
    // Initialize iteration for DSDP solver
    DSDP_INT retcode = DSDP_RETCODE_OK;
    
    inity(dsdpSolver);  initkappatau(dsdpSolver);
    initmu(dsdpSolver); initparams(dsdpSolver);
    initresi(dsdpSolver);
    
    dsdpSolver->Pnrm = DSDP_INFINITY;
    printf("| DSDP is initialized with Ry = %3.3e * I %52s\n", dsdpSolver->Ry, "");
    return retcode;
}

extern DSDP_INT dsdpInitializeB( HSDSolver *dsdpSolver ) {
    
    // Initialize the mu parameter and the initial primal objective bound
    DSDP_INT retcode = DSDP_RETCODE_OK;
    DSDP_INT BmaxIter; double rho, initpObj, pfeas;
    DSDPGetDblParam(dsdpSolver, DBL_PARAM_RHON, &rho);
    DSDPGetDblParam(dsdpSolver, DBL_PARAM_INIT_POBJ, &initpObj);
    DSDPGetIntParam(dsdpSolver, INT_PARAM_BMAXITER, &BmaxIter);
    DSDPStatUpdate(&dsdpSolver->dsdpStats, STAT_NUM_SMALL_ITER, 0.0);
    DSDPGetStats(&dsdpSolver->dsdpStats, STAT_PFEAS_PROBLEM, &pfeas);
    
    // y = y / tau;
    vec_rscale(dsdpSolver->y, dsdpSolver->tau);
    getBslack(dsdpSolver, dsdpSolver->y, DUALVAR);
    vec_dot(dsdpSolver->dObj, dsdpSolver->y, &dsdpSolver->dObjVal);
    dsdpSolver->dPotential = DSDP_INFINITY;

    // S = S / tau;
    DSDPConic( COPS_DO_VAR_SCALE )(dsdpSolver, dsdpSolver->tau);
    
    // Primal relaxation
    if (!pfeas) DSDP_HEURS( adjPRelaxPenalty ) (dsdpSolver);
    
    printf("| Primal relaxation penalty is set to %10.3e \n", dsdpSolver->ybound);
    
    // Dual perturb
    DSDP_HEURS( adjDualPerturb ) ( dsdpSolver );
    dsdpSolver->dperturb = MAX(-dsdpSolver->Ry, dsdpSolver->dperturb);
    printf("| Perturbing dual iterations by %10.3e \n", dsdpSolver->dperturb);
    
    // Barrier
    switch (dsdpSolver->solStatus) {
        case DSDP_PD_FEASIBLE:
            // mu = min((pObj - dObj) / rho, muPrimal)
            dsdpSolver->mu = MIN((dsdpSolver->pObjVal - dsdpSolver->dObjVal) / \
                                 (rho * dsdpSolver->nall), dsdpSolver->mu);
            printf("| DSDP Phase B starts. Restarting dual-scaling %51s \n", "");
            break;
        case DSDP_PUNKNOWN_DFEAS:
            if (pfeas) {
                dsdpSolver->pObjVal = 0.0;
            }
            dsdpSolver->pObjVal = MAX(initpObj, dsdpSolver->dObjVal + 0.1 * initpObj);
            dsdpSolver->mu = (dsdpSolver->pObjVal - dsdpSolver->dObjVal) / rho;
            printf("| DSDP Phase B starts. "
                   "Trying to certificate primal feasibility in %d iterations %16s \n", BmaxIter, "");
            break;
        default:
            error(etype, "Invalid starting status for Phase B. \n");
            break;
    }
    
    printf("| Heuristic start: mu: %10.3e pObj: %10.3e  dObj: %10.3e %29s\n",
           dsdpSolver->mu, dsdpSolver->pObjVal, dsdpSolver->dObjVal, "");
    
    dsdpSolver->tau = 1.0;
    return retcode;
}
