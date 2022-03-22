#include <stdio.h>
#include <stdlib.h>
#include "test.h"
#include "sparsemat.h"
#include "dsdplanczos.h"

DSDP_INT test_lanczos(void) {
    
    DSDP_INT retcode = DSDP_RETCODE_OK;
    
    DSDP_INT p[] = {0, 2, 3};
    DSDP_INT i[] = {0, 1, 1};
    DSDP_INT p2[] = {0, 0, 0};
    DSDP_INT i2[] = {0, 0, 0};
    double x[] = {94.816148385236914, 0.00036219860919956828, 94.825133003301544};
    double x2[] = {4.950904159595997397039741981344e-09,
                   1.949168277611792912239890851467e-04,
                  -7.147300985018342485248499240669e-05};
    double x3[3] = {0.0};
    
    spsMat S, dS, spaux;
    DSDPLanczos lczSolver;
    
    dsdpLanczosInit(&lczSolver);
    dsdpLanczosAlloc(&lczSolver, 2);
    
    S.p = p; dS.p = p; spaux.p = p2;
    S.i = i; dS.i = i; spaux.i = i2;
    S.x = x; dS.x = x2; spaux.x = x3;
    S.dim = 2; dS.dim = 2; spaux.dim = 2;
    S.nnz = 3; dS.nnz = 3; spaux.nnz = 3;
    S.isFactorized = TRUE;
    
    
    
    spsMatSymbolic(&S);
    spsMatFactorize(&S);
    spsMatInvView(&S);
    spsMatLinvView(&S);
    
    double alpha = 0.0;
    dsdpGetAlpha(&lczSolver, &S, &dS, &spaux, &alpha);
    dsdpLanczosFree(&lczSolver);
    
    return retcode;
}
