#include "dsdpinitializer.h"
#include "dsdpcorrector.h"
#include "dsdputils.h"
/* Implement the initialization procedure of DSDP
 After the routine Ry
 */
static char etype[] = "DSDP Initialization";

static DSDP_INT isConstant( HSDSolver *dsdpSolver, DSDP_INT *isConstant ) {
    // Check whether C is a constant
    DSDP_INT retcode = DSDP_RETCODE_OK;
    *isConstant = FALSE;
    
    DSDP_INT m = dsdpSolver->m, dim = 0;
    DSDP_INT isCons = TRUE, isRank1 = FALSE;
    r1Mat *C = NULL; double num = 0.0;
    
    for (DSDP_INT i = 0; i < dsdpSolver->nBlock; ++i) {
        isCons = TRUE;
        if (dsdpSolver->sdpData[i]->types[m] == MAT_TYPE_RANKK) {
            rkMatisRank1(dsdpSolver->sdpData[i]->sdpData[m], &isRank1);
            if (isRank1) {
                C = (r1Mat *) dsdpSolver->sdpData[i]->sdpData[m];
                dim = C->dim; num = C->x[0];
                
                for (DSDP_INT j = 1; j < dim; ++j) {
                    if (fabs(num - C->x[j]) > 1e-03) {
                        isCons = FALSE; break;
                    }
                }
            } else {
                isCons = FALSE;
            }
        } else {
            isCons = FALSE;
        }
        
        if (isCons) {
            *isConstant = TRUE; break;
        }
    }
    
    return retcode;
}

static DSDP_INT inity( HSDSolver *dsdpSolver ) {
    
    DSDP_INT retcode = DSDP_RETCODE_OK;
    vec_reset(dsdpSolver->y);
    return retcode;
}

static DSDP_INT initkappatau( HSDSolver *dsdpSolver ) {
    
    DSDP_INT retcode = DSDP_RETCODE_OK;
    retcode = DSDPGetDblParam(dsdpSolver, DBL_PARAM_INIT_TAU, &dsdpSolver->tau);
#ifdef KAPPATAU
    retcode = DSDPGetDblParam(dsdpSolver, DBL_PARAM_INIT_KAPPA, &dsdpSolver->kappa);
#else
    dsdpSolver->kappa = dsdpSolver->mu / dsdpSolver->tau;
#endif
    return retcode;
}

static DSDP_INT initmu( HSDSolver *dsdpSolver ) {
    
    DSDP_INT retcode = DSDP_RETCODE_OK;
    retcode = DSDPGetDblParam(dsdpSolver, DBL_PARAM_INIT_MU, &dsdpSolver->mu);
    return retcode;
}

static DSDP_INT initresi( HSDSolver *dsdpSolver ) {
    // Initialize Ry
    
    /*
     if norm(Ctau, 'fro') == 0
         Ry = - speye(n) * initbeta;
     else
         Ry = - speye(n) * (norm(Ctau, 'fro')) * initbeta;
     end % End if
    */
    
    DSDP_INT retcode = DSDP_RETCODE_OK;
    
    double beta, Cnrm = 0.0, nrm  = 0.0, tau = dsdpSolver->tau;
    DSDPGetDblParam(dsdpSolver, DBL_PARAM_INIT_BETA, &beta);
    DSDP_INT m  = dsdpSolver->m, isCons = FALSE;
    retcode = isConstant(dsdpSolver, &isCons);
    
    if (isCons) {
        beta *= 3.0;
    }
    
    for (DSDP_INT i = 0; i < dsdpSolver->nBlock; ++i) {
        // Matrix of index m in each block is C
        retcode = getMatFnorm(dsdpSolver, i, m, &nrm);
        checkCode;
        Cnrm = nrm * nrm + Cnrm;
    }
        
    if (fabs(Cnrm) < 1e-08) {
        dsdpSolver->Ry = - beta;
    } else {
        dsdpSolver->Ry = - sqrt(Cnrm) * tau * beta;
    }
    
    dsdpSolver->Ry = MAX(dsdpSolver->Ry, -1e+08);
    dsdpSolver->Ry = MIN(dsdpSolver->Ry, -Cnrm);
    
    for (DSDP_INT i = 0; i < dsdpSolver->nBlock; ++i) {
        // S = C * tau - Ry
        retcode = addMattoS(dsdpSolver, i, m, tau);
        retcode = spsMatAdddiag(dsdpSolver->S[i],
                                - dsdpSolver->Ry,
                                dsdpSolver->symS[i]);
        checkCode;
    }
    
    // Primal infeasibility
    DSDPGetStats(&dsdpSolver->dsdpStats, STAT_ONE_NORM_B, &nrm);
    dsdpSolver->pinfeas = nrm + 1.0;
    
    // Iteration Monitor
    dsdpSolver->iterProgress[ITER_INITIALIZE] = TRUE;
    
    return retcode;
}

static void initparams( HSDSolver *dsdpSolver ) {
    // TODO: Adjust solver parameter based on problem features
    
    double dblparam, stats; DSDP_INT intparam;
    
    // Adjust bound on the y variable. Set it to the initial value
    DSDPGetDblParam(dsdpSolver, DBL_PARAM_PRLX_PENTALTY, &dblparam);
    
    // If solving a feasibility problem
    DSDPGetStats(&dsdpSolver->dsdpStats, STAT_PFEAS_PROBLEM, &stats);
    dblparam = (stats) ?  MIN(0.001, dblparam) : dblparam;
    dsdpSolver->pObjVal = (stats) ? 0.0 : dsdpSolver->pObjVal;
    
    // If not using primal relaxation in Phase A
    DSDPGetIntParam(dsdpSolver, INT_PARAM_PRELAX, &intparam);
    if (!intparam) {
        dblparam = DSDP_INFINITY;
    }
    
    dsdpSolver->ybound = dblparam;
    
    // Conditions for using conservative HSD start and bound
    
    
    // Conditions for using aggressive dual infeasibility elimination scheme and fixed bound
    
}

static void adjpObj( HSDSolver *dsdpSolver ) {
    // Adjust primal objecive
    double pObj = dsdpSolver->pObjVal;
    if (pObj > 1e+10) {
        pObj *= 1.1;
    } else if (pObj < -1e+05) {
        pObj += 1e+03;
    } else if (pObj < 0) {
        pObj *= 0.5;
    }
    
    pObj = MAX(1e+08, pObj);
    pObj = MIN(1e+15, pObj);
    dsdpSolver->pObjVal = pObj;
}

extern DSDP_INT dsdpInitializeA( HSDSolver *dsdpSolver ) {
    
    // Initialize iteration for DSDP solver
    // TODO: Add special case if C is all-constant
    
    DSDP_INT retcode = DSDP_RETCODE_OK;
    
    retcode = inity(dsdpSolver); checkCode;
    retcode = initkappatau(dsdpSolver); checkCode;
    retcode = initmu(dsdpSolver); checkCode;
    retcode = initresi(dsdpSolver); checkCode;
    initparams(dsdpSolver); checkCode;
    // dsdpSolver->mu *= - dsdpSolver->Ry * sqrt(dsdpSolver->n);
    dsdpSolver->Pnrm = DSDP_INFINITY;
    printf("| DSDP is initialized with Ry = %3.3e * I %52s|\n", dsdpSolver->Ry, "");
    
    return retcode;
}

extern DSDP_INT dsdpInitializeB( HSDSolver *dsdpSolver ) {
    
    // Initialize the mu parameter and the initial primal objective bound
    DSDP_INT retcode = DSDP_RETCODE_OK;
    
    assert( dsdpSolver->eventMonitor[EVENT_IN_PHASE_B] );
    
    DSDP_INT BmaxIter; double rho, initpObj, tmp, tmp2;
    retcode = DSDPGetDblParam(dsdpSolver, DBL_PARAM_RHON, &rho);
    retcode = DSDPGetDblParam(dsdpSolver, DBL_PARAM_INIT_POBJ, &initpObj);
    retcode = DSDPGetIntParam(dsdpSolver, INT_PARAM_BMAXITER, &BmaxIter);

    // Reset the number of iterations taking small stepsize
    DSDPStatUpdate(&dsdpSolver->dsdpStats, STAT_NUM_SMALL_ITER, 0.0);
    
    // y = y / tau;
    vec_rscale(dsdpSolver->y, dsdpSolver->tau);
    vec_dot(dsdpSolver->dObj, dsdpSolver->y, &dsdpSolver->dObjVal);
    dsdpSolver->dPotential = DSDP_INFINITY;
    
    // S = S / tau;
    for (DSDP_INT i = 0; i < dsdpSolver->nBlock; ++i) {
        retcode = spsMatRscale(dsdpSolver->S[i], dsdpSolver->tau);
    }
    
    // Re-initialize bound of the dual feasible solution
    tmp = vec_infnorm(dsdpSolver->y);
    
    if (tmp > 1e+07) {
        dsdpSolver->ybound = 1.5 * tmp;
    } else if (tmp > 1e+05) {
        tmp = MIN(tmp * 10, 1e+07);
    } else if (tmp > 1e+03) {
        tmp = MIN(tmp * 100, 1e+06);
    } else if (tmp > 10) {
        tmp = MIN(tmp * 200, 1e+05);
    } else {
        tmp = MIN(tmp * 1000, 1e+04);
    }
    
    dsdpSolver->ybound = tmp;
    
    printf("| Primal relaxation penalty is set to %10.3e \n", dsdpSolver->ybound);
    
    DSDPGetStats(&dsdpSolver->dsdpStats, STAT_PHASE_A_ITER, &tmp);
    DSDPGetStats(&dsdpSolver->dsdpStats, STAT_ONE_NORM_C, &tmp2);
    if (tmp >= 80) {
        dsdpSolver->dperturb = 5e-06 * tmp2;
    } else if (tmp >= 50) {
        dsdpSolver->dperturb = 1e-06 * tmp2;
    } else if (tmp >= 30) {
        dsdpSolver->dperturb = 1e-07 * tmp2;
    } else {
        dsdpSolver->dperturb = 0.0;
    }
    
    dsdpSolver->dperturb = 0.0;
    
    printf("| Perturbing dual iterations by %10.3e \n", dsdpSolver->dperturb);
    switch (dsdpSolver->solStatus) {
            
        case DSDP_PD_FEASIBLE:
            // mu = min((pObj - dObj) / rho, muPrimal)
            adjpObj(dsdpSolver);
            dsdpSolver->mu = MIN((dsdpSolver->pObjVal - dsdpSolver->dObjVal) / rho, dsdpSolver->mu);
            printf("| DSDP Phase B starts. Restarting dual-scaling %51s |\n", "");
            break;
        case DSDP_PUNKNOWN_DFEAS:
            // pObj  = max(dsdpParam{5}, dObj + dsdpParam{5} / 10);
            dsdpSolver->pObjVal = MAX(initpObj, dsdpSolver->dObjVal + 0.1 * initpObj);
            dsdpSolver->mu = (dsdpSolver->pObjVal - dsdpSolver->dObjVal) / rho;
            printf("| DSDP Phase B starts. "
                   "Trying to certificate primal feasibility in %d iterations %16s |\n", BmaxIter, "");
            break;
        default:
            error(etype, "Invalid starting status for Phase B. \n");
            break;
    }
    
    DSDPGetStats(&dsdpSolver->dsdpStats, STAT_PFEAS_PROBLEM, &tmp);
    if (tmp) {
        dsdpSolver->pObjVal = 0.0;
    }
    
    printf("| Heuristic start: mu: %10.3e pObj: %10.3e  dObj: %10.3e %29s |\n",
           dsdpSolver->mu, dsdpSolver->pObjVal, dsdpSolver->dObjVal, "");
    
    dsdpSolver->tau = 1.0;
    
    return retcode;
}
