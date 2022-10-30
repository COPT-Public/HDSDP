#include <sys/time.h>
#include "pot_utils.h"

/* TODO: Add compatibility for Windows platform */
static double my_clock( void ) {
    
    struct timeval t;
    gettimeofday(&t, NULL);
    return (1.0e-6 * t.tv_usec + t.tv_sec);
}

extern double potUtilGetTimeStamp( void ) {
    
    return my_clock();
}

extern void potUtilGetDefaultParams( double dblParams[NUM_DBL_PARAM], int intParams[NUM_INT_PARAM] ) {
    
    /* Absolute optimality tolerance */
    dblParams[DBL_PARAM_RELOPTTOL] = 1e-04;
    /* Absolute feasibility tolerance */
    dblParams[DBL_PARAM_RELFEASTOL] = 1e-04;
    /* Time limit  */
    dblParams[DBL_PARAM_TIMELIMIT] = 600.0;
    /* Objective focus */
    dblParams[DBL_PARAM_COMPFOCUS] = 10.0;
    /* Maximum iteration */
    intParams[INT_PARAM_MAXITER] = 10000;
    /* Maximum maximum Ruiz iteration */
    intParams[INT_PARAM_MAXRUIZITER] = 100;
    /* Switch of coefficient scaling */
    intParams[INT_PARAM_COEFSCALE] = 0;
    /* Switch of curvature usage */
    intParams[INT_PARAM_CURVATURE] = 1;
    
    return;
}

extern void potUtilPrintParams( double dblParams[NUM_DBL_PARAM], int intParams[NUM_INT_PARAM] ) {
    
    printf("Parameter summary \n");
    printf("MaxIter     is set to %d \n", intParams[INT_PARAM_MAXITER]);
    printf("RuizMaxIter is set to %d \n", intParams[INT_PARAM_MAXRUIZITER]);
    printf("CoeffScal   is set to %d \n", intParams[INT_PARAM_COEFSCALE]);
    printf("Curvature   is set to %d \n", intParams[INT_PARAM_CURVATURE]);
    printf("RelFeasTol  is set to %3.3e \n", dblParams[DBL_PARAM_RELFEASTOL]);
    printf("RelOptTol   is set to %3.3e \n", dblParams[DBL_PARAM_RELOPTTOL]);
    printf("TimeLimit   is set to %.0fs \n", dblParams[DBL_PARAM_TIMELIMIT]);
    printf("compFocus   is set to %3.3e \n", dblParams[DBL_PARAM_COMPFOCUS]);
    
    return;
}

