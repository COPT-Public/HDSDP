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
    
    DSDP_INT m = dsdpSolver->m;
    DSDP_INT dim = 0;
    r1Mat *C = NULL;
    DSDP_INT isCons = TRUE;
    double num = 0.0;
    
    for (DSDP_INT i = 0; i < dsdpSolver->nBlock; ++i) {
        isCons = TRUE;
        if (dsdpSolver->sdpData[i]->types[m] == MAT_TYPE_RANK1) {
            C = (r1Mat *) dsdpSolver->sdpData[i]->sdpData[m];
            dim = C->dim;
            num = C->x[0];
            for (DSDP_INT j = 1; j < dim; ++j) {
                if (fabs(num - C->x[j]) > 1e-03) {
                    isCons = FALSE;
                    break;
                }
            }
        } else {
            isCons = FALSE;
        }
        
        if (isCons) {
            *isConstant = TRUE;
            break;
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
    dsdpSolver->kappa = 1.0;
    dsdpSolver->tau   = 1.0;
    return retcode;
}

static DSDP_INT initmu( HSDSolver *dsdpSolver ) {
    
    DSDP_INT retcode = DSDP_RETCODE_OK;
    dsdpSolver->mu = dsdpSolver->param->initMu;
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
    
    double beta = dsdpSolver->param->initBeta;
    double Cnrm = 0.0;
    double nrm  = 0.0;
    double tau  = dsdpSolver->tau;
    DSDP_INT m  = dsdpSolver->m;
    
    DSDP_INT isCons = FALSE;
    retcode = isConstant(dsdpSolver, &isCons);
    
    if (isCons) {
        beta *= 10.0;
    }
    
    for (DSDP_INT i = 0; i < dsdpSolver->nBlock; ++i) {
        // Matrix of index m in each block is C
        retcode = getMatnrm(dsdpSolver, i, m, &nrm);
        checkCode;
        Cnrm += nrm * nrm;
    }
    
    Cnrm = sqrt(Cnrm);
    
    if (fabs(Cnrm) < 1e-08) {
        dsdpSolver->Ry = - beta;
    } else {
        dsdpSolver->Ry = - Cnrm * tau * beta;
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
    
    printf("| DSDP is initialized with Ry = %3.3e * I \n", dsdpSolver->Ry);
    
    return retcode;
}

extern DSDP_INT dsdpInitializeB( HSDSolver *dsdpSolver ) {
    
    // Initialize the mu parameter and the initial primal objective bound
    DSDP_INT retcode = DSDP_RETCODE_OK;
    
    assert( dsdpSolver->eventMonitor[EVENT_IN_PHASE_B] );
    
    // Reset the number of iterations taking small stepsize
    dsdpSolver->smallIter = 0;
    dsdpSolver->dObjVal = dsdpSolver->dObjVal / dsdpSolver->tau;
    
    // y = y / tau;
    vec_rscale(dsdpSolver->y, dsdpSolver->tau);
    
    // S = S / tau;
    for (DSDP_INT i = 0; i < dsdpSolver->nBlock; ++i) {
        retcode = spsMatRscale(dsdpSolver->S[i], dsdpSolver->tau);
    }
    
    switch (dsdpSolver->solStatus) {
            
        case DSDP_PD_FEASIBLE:
            // mu = min((pObj - dObj) / rho, muPrimal)
            dsdpSolver->mu = MIN((dsdpSolver->pObjVal - dsdpSolver->dObjVal) / \
                                  dsdpSolver->param->rho, dsdpSolver->mu);
            printf("| DSDP Phase B starts. Restarting dual-scaling \n");
            break;
        case DSDP_PUNKNOWN_DFEAS:
            //pObj  = max(dsdpParam{5}, dObj + dsdpParam{5} / 10);
            dsdpSolver->pObjVal = MAX(dsdpSolver->param->initpObj,
                                      dsdpSolver->dObjVal + 0.1 * dsdpSolver->param->initpObj);
            dsdpSolver->mu = (dsdpSolver->pObjVal - dsdpSolver->dObjVal) / dsdpSolver->param->rho;
            printf("| DSDP Phase B starts. "
                   "Trying to certificate primal feasibility in %d iterations \n",
                   dsdpSolver->param->BmaxIter);
            break;
        default:
            error(etype, "Invalid starting status for Phase B. \n");
            break;
    }
    
    printf("| Heuristic start: mu: %10.3e pObj: %10.3e  dObj: %10.3e \n",
           dsdpSolver->mu, dsdpSolver->pObjVal, dsdpSolver->dObjVal);
    
    return retcode;
}
