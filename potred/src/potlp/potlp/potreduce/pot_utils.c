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
    /* Maximum iteration */
    intParams[INT_PARAM_MAXITER] = 10000;
    /* Maximum maximum Ruiz iteration */
    intParams[INT_PARAM_MAXRUIZITER] = 100;
    /* Switch of coefficient scaling */
    intParams[INT_PARAM_COEFSCALE] = 0;
    
    return;
}

