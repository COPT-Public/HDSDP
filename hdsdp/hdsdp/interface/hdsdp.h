/** @file   hdsdp.h
 *  @brief  Define the utilities and macros for HDSDP solver
 *  @author Wenzhi Gao
 */
#ifndef hdsdp_h
#define hdsdp_h

/* Macro for MATLAB */
#ifdef MATLAB_MEX_FILE
#include "mex.h"
#define calloc mxCalloc
#define printf mexPrintf
#define free mxFree
#else
#endif

typedef enum {
    
    HDSDP_RETCODE_OK,
    HDSDP_RETCODE_FAILED,
    HDSDP_RETCODE_MEMORY,
    
} hdsdp_retcode;

typedef enum {
    
    HDSDP_UNKNOWN,
    HDSDP_OPTIMAL,
    HDSDP_MAXITER,
    HDSDP_INFEAS_OR_UNBOUNDED,
    HDSDP_TIMELIMIT,
    HDSDP_USER_INTERRUPT,
    HDSDP_INTERNAL_ERROR,
    HDSDP_NUMERICAL
    
} hdsdp_status;

typedef struct hdsdp_solver_internal hdsdp;

// Integer Parameters
#define INT_PARAM_MAXITER       0

// Double Parameters
#define DBL_PARAM_RELFEASTOL    0
#define DBL_PARAM_RELOPTTOL     1
#define DBL_PARAM_TIMELIMIT     2

// Version information
#define VERSION_MAJOR           2
#define VERSION_MINOR           0
#define VERSION_TECHNICAL       0

// Build date
#define BUILD_DATE_YEAR         2022
#define BUILD_DATE_MONTH        11
#define BUILD_DATE_DAY          24


#endif /* hdsdp_h */
