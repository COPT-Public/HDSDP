#include "dsdpinitializer.h"
#include "dsdputils.h"
/* Implement the initialization procedure of DSDP
 After the routine Ry
 */
static char etype[] = "DSDP Initialization";

static DSDP_INT inity( HSDSolver *dsdpSolver ) {
    
    DSDP_INT retcode = DSDP_RETCODE_OK;
    memset(dsdpSolver->y, 0, sizeof(double) * dsdpSolver->m);
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
    DSDP_INT retcode = DSDP_RETCODE_OK;
    
    double beta = dsdpSolver->param->initBeta;
    double Cnrm = 0.0;
    double nrm = 0.0;
    double tau = dsdpSolver->tau;
    DSDP_INT m  = dsdpSolver->m;
    
    for (DSDP_INT i = 0; i < dsdpSolver->nBlock; ++i) {
        retcode = getMatnrm(dsdpSolver, i, m, &nrm);
        Cnrm += nrm * nrm;
    }
    
    Cnrm = sqrt(Cnrm);
    
    if (fabs(Cnrm) < 1e-08) {
        dsdpSolver->Ry = - beta;
    } else {
        dsdpSolver->Ry = - Cnrm * tau * beta;
    }
    
    for (DSDP_INT i = 0; i < dsdpSolver->nBlock; ++i) {
        retcode = addMattoS(dsdpSolver, i, m, tau);
        retcode = spsMatAdddiag(dsdpSolver->S[i], - dsdpSolver->Ry);
    }
    
    return retcode;
}

extern DSDP_INT dsdpInitialize( HSDSolver *dsdpSolver ) {
    
    // Initialize iteration for DSDP solver
    DSDP_INT retcode = DSDP_RETCODE_OK;
    
    retcode = inity(dsdpSolver);
    retcode = initkappatau(dsdpSolver);
    retcode = initmu(dsdpSolver);
    retcode = initresi(dsdpSolver);
    
    printf("* DSDP is initialized. \n");
    
    return retcode;
}
