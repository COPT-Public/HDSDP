#include "dsdpinitializer.h"
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
    retcode = DSDPGetDblParam(dsdpSolver, DBL_PARAM_INIT_KAPPA, &dsdpSolver->kappa);
    retcode = DSDPGetDblParam(dsdpSolver, DBL_PARAM_INIT_TAU, &dsdpSolver->tau);
    return retcode;
}

static DSDP_INT initmu( HSDSolver *dsdpSolver ) {
    
    DSDP_INT retcode = DSDP_RETCODE_OK;
    retcode = DSDPGetDblParam(dsdpSolver, DBL_PARAM_INIT_MU, &dsdpSolver->mu);
    return retcode;
}

static DSDP_INT initresi( HSDSolver *dsdpSolver ) {
    // Initialize Ry TODO: and ry
    
    /*
     if norm(Ctau, 'fro') == 0
         Ry = - speye(n) * initbeta;
     else
         Ry = - speye(n) * (norm(Ctau, 'fro')) * initbeta;
     end % End if
    */
    
    DSDP_INT retcode = DSDP_RETCODE_OK;
    
    double beta;
    retcode = DSDPGetDblParam(dsdpSolver, DBL_PARAM_INIT_BETA, &beta);
    double Cnrm = 0.0, nrm  = 0.0, tau = dsdpSolver->tau;
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
    
    for (DSDP_INT i = 0; i < dsdpSolver->nBlock; ++i) {
        // S = C * tau - Ry
        retcode = addMattoS(dsdpSolver, i, m, tau);
        retcode = spsMatAdddiag(dsdpSolver->S[i],
                                - dsdpSolver->Ry,
                                dsdpSolver->symS[i]);
        checkCode;
    }
    
    // Iteration Monitor
    dsdpSolver->iterProgress[ITER_INITIALIZE] = TRUE;
    
    return retcode;
}

extern DSDP_INT dsdpInitializeA( HSDSolver *dsdpSolver ) {
    
    // Initialize iteration for DSDP solver
    // TODO: Add special case if C is all-constant
    
    DSDP_INT retcode = DSDP_RETCODE_OK;
    
    retcode = inity(dsdpSolver); checkCode;
    retcode = initkappatau(dsdpSolver); checkCode;
    retcode = initmu(dsdpSolver); checkCode;
    retcode = initresi(dsdpSolver); checkCode;
    
    // dsdpSolver->mu *= - dsdpSolver->Ry * sqrt(dsdpSolver->n);
    dsdpSolver->Pnrm = DSDP_INFINITY;
    printf("| DSDP is initialized with Ry = %3.3e * I %52s|\n", dsdpSolver->Ry, "");
    
    return retcode;
}

extern DSDP_INT dsdpInitializeB( HSDSolver *dsdpSolver ) {
    
    // Initialize the mu parameter and the initial primal objective bound
    DSDP_INT retcode = DSDP_RETCODE_OK;
    
    assert( dsdpSolver->eventMonitor[EVENT_IN_PHASE_B] );
    
    DSDP_INT BmaxIter; double rho, initpObj;
    retcode = DSDPGetDblParam(dsdpSolver, DBL_PARAM_RHO, &rho);
    retcode = DSDPGetDblParam(dsdpSolver, DBL_PARAM_INIT_POBJ, &initpObj);
    retcode = DSDPGetIntParam(dsdpSolver, INT_PARAM_BMAXITER, &BmaxIter);

    // Reset the number of iterations taking small stepsize
    DSDPStatUpdate(&dsdpSolver->dsdpStats, STAT_NUM_SMALL_ITER, 0.0);
    dsdpSolver->dObjVal = dsdpSolver->dObjVal / dsdpSolver->tau;
    dsdpSolver->dPotential = DSDP_INFINITY;
    
    // y = y / tau;
    vec_rscale(dsdpSolver->y, dsdpSolver->tau);
    
    // S = S / tau;
    for (DSDP_INT i = 0; i < dsdpSolver->nBlock; ++i) {
        retcode = spsMatRscale(dsdpSolver->S[i], dsdpSolver->tau);
    }
    
    switch (dsdpSolver->solStatus) {
            
        case DSDP_PD_FEASIBLE:
            // mu = min((pObj - dObj) / rho, muPrimal)
            dsdpSolver->mu = MIN((dsdpSolver->pObjVal - dsdpSolver->dObjVal) / rho, dsdpSolver->mu);
            printf("| DSDP Phase B starts. Restarting dual-scaling %51s |\n", "");
            break;
        case DSDP_PUNKNOWN_DFEAS:
            //pObj  = max(dsdpParam{5}, dObj + dsdpParam{5} / 10);
            dsdpSolver->pObjVal = MAX(initpObj, dsdpSolver->dObjVal + 0.1 * initpObj);
            dsdpSolver->mu = (dsdpSolver->pObjVal - dsdpSolver->dObjVal) / rho;
            printf("| DSDP Phase B starts. "
                   "Trying to certificate primal feasibility in %d iterations %16s |\n", BmaxIter, "");
            break;
        default:
            error(etype, "Invalid starting status for Phase B. \n");
            break;
    }
    
    printf("| Heuristic start: mu: %10.3e pObj: %10.3e  dObj: %10.3e %29s |\n",
           dsdpSolver->mu, dsdpSolver->pObjVal, dsdpSolver->dObjVal, "");
    
    return retcode;
}
