#include "dsdpinitializer.h"
#include "dsdputils.h"
/* Implement the initialization procedure of DSDP
 After the routine Ry
 */
static char etype[] = "DSDP Initialization";

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

extern DSDP_INT dsdpInitialize( HSDSolver *dsdpSolver ) {
    
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
