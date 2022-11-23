#include "pot_utils.h"
#include "pot_param.h"

#include <math.h>
#include <sys/time.h>

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
    /* Maximum Ruiz iteration */
    intParams[INT_PARAM_MAXRUIZITER] = 10;
    /* Maximum PC iteration */
    intParams[INT_PARAM_MAXPCITER] = 10;
    /* Switch of coefficient scaling */
    intParams[INT_PARAM_COEFSCALE] = 0;
    /* Number of curvature usage */
    intParams[INT_PARAM_CURVATURE] = 20;
    /* Frequency of switching the weight of the residuals */
    intParams[INT_PARAM_RSCALFREQ] = -1;
    /* Whether to scale simplex to start from all-one */
    intParams[INT_PARAM_SCALSIMPLEX] = 1;
    
    /* Internal parameters */
    dblParams[DBL_IPARAM_RESTARTRATE] = 20.0;
    dblParams[DBL_IPARAM_RESTARTMAX] = 1e+04;
    
    return;
}

extern void potUtilPrintParams( double dblParams[NUM_DBL_PARAM], int intParams[NUM_INT_PARAM] ) {
    
    printf("External parameter summary \n");
    printf("MaxIter     is set to %d \n", intParams[INT_PARAM_MAXITER]);
    printf("RuizMaxIter is set to %d \n", intParams[INT_PARAM_MAXRUIZITER]);
    printf("PCMaxIter   is set to %d \n", intParams[INT_PARAM_MAXPCITER]);
    printf("CoeffScal   is set to %d \n", intParams[INT_PARAM_COEFSCALE]);
    printf("Curvature   is set to %d \n", intParams[INT_PARAM_CURVATURE]);
    printf("CInterVal   is set to %d \n", intParams[INT_PARAM_CURVINTERVAL]);
    printf("RScalFreq   is set to %d \n", intParams[INT_PARAM_RSCALFREQ]);
    printf("ScalSpx     is set to %d \n", intParams[INT_PARAM_SCALSIMPLEX]);
    printf("RelFeasTol  is set to %3.3e \n", dblParams[DBL_PARAM_RELFEASTOL]);
    printf("RelOptTol   is set to %3.3e \n", dblParams[DBL_PARAM_RELOPTTOL]);
    printf("TimeLimit   is set to %.0fs \n", dblParams[DBL_PARAM_TIMELIMIT]);
    
    return;
}

extern void potUtilPrintIParams( double dblParams[NUM_DBL_PARAM], int intParams[NUM_INT_PARAM] ) {
    
    printf("\nInternal parameter summary \n");
    printf("RestartRate is set to %3.3e \n", dblParams[DBL_IPARAM_RESTARTRATE]);
    printf("RestartMax  is set to %3.3e \n", dblParams[DBL_IPARAM_RESTARTMAX]);
    
    return;
}


/* Debugging */
extern void potUtilPrintDblContent( int n, double *d ) {
    
    for ( int i = 0; i < n; ++i ) {
        printf("%5.3e, ", d[i]);
    }
    printf("\n");
    return;
}

extern void potUtilPrintIntContent( int n, int *d ) {
    
    for ( int i = 0; i < n; ++i ) {
        printf("%5d, ", d[i]);
    }
    printf("\n");
    return;
}

extern void potUtilPrintDblSum( int n, double *d ) {
    
    double ds = 0.0;
    
    for ( int i = 0; i < n; ++i ) {
        ds += d[i];
    }
    
    printf("Sum = %10.6e \n", ds);
    return;
}

extern int potUtilVerifyNeighbour( int n, double *d, double r ) {
    /* Check if log(d_i) >= r */
    int inNeighbour = 1;
    for ( int i = 0; i < n; ++i ) {
        if ( log(d[i]) < r ) {
            inNeighbour = 0;
            break;
        }
    }
    
    return inNeighbour;
}
