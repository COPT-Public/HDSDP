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

//#define HDSDP_LANCZOS_DEBUG
//#define HDSDP_CONIC_DEBUG
#define HDSDP_KKT_DEBUG
#define HDSDP_LINSYS_DEBUG
#define HDSDP_CONJGRAD_DEBUG
#define KKT_ACCURACY (1e-08)
#define HDSDP_SPARSE_CONE_THRESHOLD  (0.3)
#define HDSDP_SPARSE_SCHUR_THRESHOLD (0.3)

#ifndef MEMWATCH
#define MEMWATCH
#ifdef HEADERPATH
#include "external/memwatch.h"
#else
#include "memwatch.h"
#endif
#endif

typedef enum {
    
    HDSDP_RETCODE_OK,
    HDSDP_RETCODE_FAILED,
    HDSDP_RETCODE_MEMORY,
    
} hdsdp_retcode;

typedef enum {
    
    HDSDP_UNKNOWN,
    HDSDP_DUAL_FEASIBLE,
    HDSDP_DUAL_OPTIMAL,
    HDSDP_PRIMAL_DUAL_OPTIMAL,
    HDSDP_MAXITER,
    HDSDP_SUSPECT_INFEAS_OR_UNBOUNDED,
    HDSDP_INFEAS_OR_UNBOUNDED,
    HDSDP_TIMELIMIT,
    HDSDP_USER_INTERRUPT,
    HDSDP_INTERNAL_ERROR,
    HDSDP_NUMERICAL
    
} hdsdp_status;

typedef struct hdsdp_solver_internal hdsdp;

#define HDSDP_INFINITY          1e+30

// Integer Parameters
#define INT_PARAM_MAXITER       0

// Double Parameters
#define DBL_PARAM_RELFEASTOL    0
#define DBL_PARAM_RELOPTTOL     1
#define DBL_PARAM_ABSFEASTOL    2
#define DBL_PARAM_ABSOPTTOL     3
#define DBL_PARAM_TIMELIMIT     4
#define DBL_PARAM_POTRHOVAL     5
#define DBL_PARAM_HSDGAMMA      6
#define DBL_PARAM_DUALBND       7
#define DBL_PARAM_BARMUSTART    8
#define DBL_PARAM_DUALSTART     9

// Version information
#define VERSION_MAJOR           1
#define VERSION_MINOR           0
#define VERSION_TECHNICAL       0

// Build date
#define BUILD_DATE_YEAR         2023
#define BUILD_DATE_MONTH        07
#define BUILD_DATE_DAY          21

#ifdef __cplusplus
extern "C" {
#endif

extern hdsdp_retcode HDSDPCreate( hdsdp **pHSolver );
extern hdsdp_retcode HDSDPInit( hdsdp *HSolver, int nRows, int nCones );
extern hdsdp_retcode HDSDPSetCone( hdsdp *HSolver, int iCone, void *userCone );
extern void HDSDPSetDualObjective( hdsdp *HSolver, double *dObj );
extern hdsdp_retcode HDSDPOptimize( hdsdp *HSolver, int dOptOnly );
extern hdsdp_retcode HDSDPGetRowDual( hdsdp *HSolver, double *pObjVal, double *dObjVal, double *dualVal );
extern hdsdp_retcode HDSDPGetConeValues( hdsdp *HSolver, int iCone, double *conePrimal, double *coneDual );
extern void HDSDPCheckSolution( hdsdp *HSolver, double diErrors[6] );
extern void HDSDPClear( hdsdp *HSolver );
extern void HDSDPDestroy( hdsdp **pHSolver );

#ifdef __cplusplus
}
#endif

#endif /* hdsdp_h */
