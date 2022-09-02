#include <stdio.h>
#include "rotsolve.h"

#define check_code(code) if (retcode != DSDP_RETCODE_OK) \
                            { printf("Optimization failed on line %d", __LINE__); \
                              return retcode;}

void printDual( DSDP_INT n, double *y ) {
    
    for (DSDP_INT i = 0; i < n; ++i) {
        printf("y[%d] = %+10.6e, \n", i + 1, y[i]);
    }
    
    return;
}

DSDP_INT Asp[37] = {
    0, 1, 2, 4, 5, 6, 7, 9, 11, 12, 14, 16, 18,
    19, 20, 21, 22, 24, 26, 27, 29, 32, 34,
    36, 38, 39, 41, 43, 45, 47, 49, 51, 52,
    53, 54, 55, 55
};

DSDP_INT Asi[55] = {
    54, 53, 51, 52, 50, 49, 48, 44, 47, 43,
    46, 42, 39, 45, 38, 41, 37, 40, 36, 35,
    34, 33, 26, 32, 25, 31, 24, 18, 30, 17,
    23, 29, 16, 22, 15, 28, 14, 21, 13, 9,
    27, 8, 20, 7, 19, 6, 12, 5, 11, 4, 10,
    3, 2, 1, 0
};

double Asx[55] = {
    -1.00, -1.00, -1.00, -1.00, -1.00, -1.00, -1.00, -1.00,
    -1.00, -1.00, -1.00, -1.00, -1.00, -1.00, -1.00, -1.00,
    -1.00, -1.00, -1.00, -1.00, -1.00, -1.00, -1.00, -1.00,
    -1.00, -1.00, -1.00, -1.00, -1.00, -1.00, -1.00, -1.00,
    -1.00, -1.00, -1.00, -1.00, -1.00, -1.00, -1.00, -1.00,
    -1.00, -1.00, -1.00, -1.00, -1.00, -1.00, -1.00, -1.00,
    -1.00, -1.00, -1.00, -1.00, -1.00, -1.00, -1.00
};

DSDP_INT Alp[36] = {
    0,
    1, 1,
    2, 2,
    3, 3, 3, 3, 3,
    4, 4,
    5, 5, 5,
    6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6,
    7, 7,
    8, 8, 8,
    9, 9, 9, 9,
    10
};

DSDP_INT Ali[MATRIX_DIM] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
double   Alx[MATRIX_DIM] = {-1.0, -2.0, -1.0, -2.0, -2.0, -1.0, -2.0, -2.0, -2.0, -1.0};
double lpObj[1] = {-1.0};
double btmp[NUM_CONSTR] = {0.0};
double sol[NUM_CONSTR]  = {0.0};

DSDP_INT solve_rot( double *b ) {
    
    DSDP_INT retcode = DSDP_RETCODE_OK;
    HSDSolver *dsdp = NULL;
    
    retcode = DSDPCreate(&dsdp, NULL);
    check_code(retcode);
    
    retcode = DSDPSetDim(dsdp, MATRIX_DIM, 1, NUM_CONSTR, NUM_LPVAR, NULL);
    check_code(retcode);
    
    retcode = DSDPSetSDPConeData(dsdp, 0, MATRIX_DIM, Asp, Asi, Asx);
    check_code(retcode);
    
    retcode = DSDPSetLPData(dsdp, NUM_LPVAR, Alp, Ali, Alx, lpObj);
    check_code(retcode);
    
    retcode = DSDPSetObj(dsdp, b);
    check_code(retcode);
    
    DSDPSetDblParam(dsdp, DBL_PARAM_REL_OPTTOL, 1e-10);
    DSDPSetDblParam(dsdp, DBL_PARAM_ABS_OPTTOL, 1e-10);
    retcode = DSDPOptimize(dsdp);
    check_code(retcode);
    
    retcode = DSDPGetDual(dsdp, sol, NULL);
    printDual(NUM_CONSTR, sol);
    check_code(retcode);
    
    return retcode;
}
