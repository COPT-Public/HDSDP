/** @file hdsdp.c
 *  @brief A homogeneous dual-scaling interior point solver
 */

#ifdef HEADERPATH
#include "interface/def_hdsdp.h"
#include "interface/hdsdp.h"
#include "interface/hdsdp_conic.h"
#include "interface/hdsdp_utils.h"
#include "interface/hdsdp_user_data.h"
#include "interface/hdsdp_schur.h"
#include "interface/hdsdp_algo.h"
#else
#include "def_hdsdp.h"
#include "hdsdp.h"
#include "hdsdp_conic.h"
#include "hdsdp_utils.h"
#include "hdsdp_user_data.h"
#include "hdsdp_schur.h"
#include "hdsdp_algo.h"
#endif

static void HDSDPIGetStatistics( hdsdp *HSolver ) {
    /* Get features for Heurstics */
    
    /* Get sum of conic dimensions. Used to convert the identity matrix to Frobenius norm */
    int sumConeDims = 0;
    for ( int iCone = 0; iCone < HSolver->nCones; ++iCone ) {
        sumConeDims += HConeGetDim(HSolver->HCones[iCone]);
    }
    
    set_int_feature(HSolver, HDSDP_INTFEATURE_N_SUMCONEDIMS, sumConeDims);
    set_int_feature(HSolver, HDSDP_INTFEATURE_N_ROWS, HSolver->nRows);
    set_int_feature(HSolver, HDSDP_INTFEATURE_N_CONES, HSolver->nCones);
    
    return;
}

static void HDSDPIPrintStatistics( hdsdp *HSolver ) {
    
    hdsdp_printf("\nStatistics: \n");
    
    print_int_feature(HSolver, HDSDP_INTFEATURE_N_ROWS, "Number of rows");
    print_int_feature(HSolver, HDSDP_INTFEATURE_N_CONES, "Number of cones");
    print_int_feature(HSolver, HDSDP_INTFEATURE_N_SUMCONEDIMS, "Sum of cone dimensions");
    
    return;
}

extern void HDSDPIGetDefaultParams( hdsdp *HSolver ) {
    
    set_int_param(HSolver, INT_PARAM_MAXITER, 100);
    
    set_dbl_param(HSolver, DBL_PARAM_ABSOPTTOL, 1e-08);
    set_dbl_param(HSolver, DBL_PARAM_ABSFEASTOL, 1e-08);
    set_dbl_param(HSolver, DBL_PARAM_RELOPTTOL, 1e-08);
    set_dbl_param(HSolver, DBL_PARAM_RELFEASTOL, 1e-08);
    
    set_dbl_param(HSolver, DBL_PARAM_TIMELIMIT, 3600.0);
    set_dbl_param(HSolver, DBL_PARAM_POTRHOVAL, 5.0);
    set_dbl_param(HSolver, DBL_PARAM_HSDGAMMA, 0.7);
    
    return;
}

extern void HDSDPIPrintParams( hdsdp *HSolver ) {
    
    printf("\nParameters\n");
    print_int_param(HSolver, INT_PARAM_MAXITER, "Maximum iteration");
    
    print_dbl_param(HSolver, DBL_PARAM_ABSOPTTOL, "Abs optimality");
    print_dbl_param(HSolver, DBL_PARAM_RELOPTTOL, "Rel optimality");
    print_dbl_param(HSolver, DBL_PARAM_ABSFEASTOL, "Abs feasibility");
    print_dbl_param(HSolver, DBL_PARAM_RELFEASTOL, "Rel feasibility");
    print_dbl_param(HSolver, DBL_PARAM_TIMELIMIT, "Time limit");
    print_dbl_param(HSolver, DBL_PARAM_POTRHOVAL, "Potential param");
    print_dbl_param(HSolver, DBL_PARAM_HSDGAMMA, "Mu reduction rate");
    
    printf("\n");
    
    return;
}

/* Define HDSDP Solver interface */
extern hdsdp_retcode HDSDPCreate( hdsdp **pHSolver ) {
    
    hdsdp_retcode retcode = HDSDP_RETCODE_OK;
    
    if ( !pHSolver ) {
        retcode = HDSDP_RETCODE_FAILED;
        return retcode;
    }
    
    hdsdp *HSolver = NULL;
    HDSDP_INIT(HSolver, hdsdp, 1);
    HDSDP_MEMCHECK(HSolver);
    
    *pHSolver = HSolver;
    
    HUtilStartCtrlCCheck();
    
exit_cleanup:
    return retcode;
}

extern hdsdp_retcode HDSDPInit( hdsdp *HSolver, int nRows, int nCones ) {
    /* Initialize the dual potential reduction solver */
    hdsdp_retcode retcode = HDSDP_RETCODE_OK;
    
    HSolver->nRows = nRows;
    HSolver->nCones = nCones;
    
    HDSDP_INIT(HSolver->rowRHS, double, nRows);
    HDSDP_MEMCHECK(HSolver->rowRHS);
        
    HDSDP_INIT(HSolver->HCones, hdsdp_cone *, nCones);
    HDSDP_MEMCHECK(HSolver->HCones);
    
    for ( int iCone = 0; iCone < nCones; ++iCone ) {
        HDSDP_CALL(HConeCreate(&HSolver->HCones[iCone]));
    }
    
    HDSDP_CALL(HKKTCreate(&HSolver->HKKT));
    
    /* Allocate vectors */
    HDSDP_INIT(HSolver->dRowDual, double, nRows);
    HDSDP_MEMCHECK(HSolver->dRowDual);
    
    HDSDP_INIT(HSolver->dRowDualStep, double, nRows);
    HDSDP_MEMCHECK(HSolver->dRowDualStep);
    
    HDSDP_INIT(HSolver->dMinvRowRHS, double, nRows);
    HDSDP_MEMCHECK(HSolver->dMinvRowRHS);
    
    HDSDP_INIT(HSolver->dMinvASinv, double, nRows);
    HDSDP_MEMCHECK(HSolver->dMinvASinv);
    
    HDSDP_INIT(HSolver->dMinvASinvRdSinv, double, nRows);
    HDSDP_MEMCHECK(HSolver->dMinvASinvRdSinv);
    
    HDSDP_INIT(HSolver->dMinvASinvCSinv, double, nRows);
    HDSDP_MEMCHECK(HSolver->dMinvASinvCSinv);
    
    HDSDP_INIT(HSolver->dHAuxiVec1, double, nRows);
    HDSDP_MEMCHECK(HSolver->dHAuxiVec1);
    
    HDSDP_INIT(HSolver->dHAuxiVec2, double, nRows);
    HDSDP_MEMCHECK(HSolver->dHAuxiVec2);
    
    HSolver->dBarrierMu = 1e+10;
    HSolver->comp = HDSDP_INFINITY;
    HSolver->compRel = HDSDP_INFINITY;
    
    HSolver->HStatus = HDSDP_UNKNOWN;
    
exit_cleanup:
    return retcode;
}

extern hdsdp_retcode HDSDPSetCone( hdsdp *HSolver, int iCone, void *userCone ) {
    
    hdsdp_retcode retcode = HDSDP_RETCODE_OK;
    
    HDSDP_CALL(HConeSetData(HSolver->HCones[iCone], userCone));
    
exit_cleanup:
    return retcode;
}

extern void HDSDPSetDualStart( hdsdp *HSolver, double *dStart ) {
    
    if ( dStart ) {
        HDSDP_MEMCPY(HSolver->dRowDual, dStart, double, HSolver->nRows);
    }
    
    return;
}

extern hdsdp_retcode HDSDPOptimize( hdsdp *HSolver, int dOptOnly ) {
    
    hdsdp_retcode retcode = HDSDP_RETCODE_OK;
    
    /* Start optimization */
    hdsdp_printf("\nHDSDP: software for semi-definite programming \n");
    printf("Wenzhi Gao, Dongdong Ge, Yinyu Ye, 2023\n");
    printf("---------------------------------------------\n");
    HSolver->dTimeBegin = HUtilGetTimeStamp();
    
    HDSDPIGetDefaultParams(HSolver);
    
    /* Process conic data */
    hdsdp_printf("Calling pre-solver \n");
    printf("  Processing the cones \n");
    for ( int iCone = 0; iCone < HSolver->nCones; ++iCone ) {
        HDSDP_CALL(HConeProcData(HSolver->HCones[iCone]));
        HDSDP_CALL(HConePresolveData(HSolver->HCones[iCone]));
    }
    
    printf("  Starting KKT analysis \n");
    /* Prepare KKT solver */
    HDSDP_CALL(HKKTInit(HSolver->HKKT, HSolver->nRows, HSolver->nCones, HSolver->HCones));
    
    /* Get statistics */
    printf("  Collecting statistcs \n");
    HDSDPIGetStatistics(HSolver);
    
    /* End pre-solver */
    hdsdp_printf("Pre-solver ends. Elapsed time: %3.1f seconds \n", HUtilGetTimeStamp() - HSolver->dTimeBegin);
    
    HDSDPIPrintStatistics(HSolver);
    HDSDPIPrintParams(HSolver);
    
exit_cleanup:
    return retcode;
}

extern hdsdp_retcode HDSDPGetRowDual( hdsdp *HSolver, double *pObjVal, double *dObjVal, double *dualVal ) {
    
    hdsdp_retcode retcode = HDSDP_RETCODE_OK;
    
    if ( pObjVal ) {
        *pObjVal = HSolver->pObjVal;
    }
    
    if ( dObjVal ) {
        *dObjVal = HSolver->dObjVal;
    }
    
    if ( dualVal ) {
        HDSDP_MEMCPY(dualVal, HSolver->dRowDual, double, HSolver->nRows);
    }
    
exit_cleanup:
    return retcode;
}

extern hdsdp_retcode HDSDPGetConeValues( hdsdp *HSolver, int iCone, double *conePrimal, double *coneDual ) {
    
    hdsdp_retcode retcode = HDSDP_RETCODE_OK;
    
    
exit_cleanup:
    return retcode;
}

#define DIMACS_ERROR_1 0
#define DIMACS_ERROR_2 1
#define DIMACS_ERROR_3 2
#define DIMACS_ERROR_4 3
#define DIMACS_ERROR_5 4
#define DIMACS_ERROR_6 5
extern void HDSDPCheckSolution( hdsdp *HSolver, double dErrs[6] ) {
    /* Compute and report 6 DIMACS error */
    
    dErrs[DIMACS_ERROR_1] = 1.0;
    dErrs[DIMACS_ERROR_2] = 1.0;
    dErrs[DIMACS_ERROR_3] = 1.0;
    dErrs[DIMACS_ERROR_4] = 1.0;
    dErrs[DIMACS_ERROR_5] = 1.0;
    dErrs[DIMACS_ERROR_6] = 1.0;
    
    hdsdp_printf("DIMACS errors: %5.2e %5.2e %5.2e %5.2e %5.2e %5.2e \n",
                 dErrs[DIMACS_ERROR_1], dErrs[DIMACS_ERROR_2], dErrs[DIMACS_ERROR_3],
                 dErrs[DIMACS_ERROR_4], dErrs[DIMACS_ERROR_5], dErrs[DIMACS_ERROR_6]);
    
    double dMaxDimacsErr = 0.0;
    
    for ( int iElem = 0; iElem < 6; ++iElem ) {
        dMaxDimacsErr = HDSDP_MAX(dMaxDimacsErr, dErrs[iElem]);
    }
    
    hdsdp_printf("Maximum error: %5.2e.", dMaxDimacsErr);
    
    return;
}

extern void HDSDPClear( hdsdp *HSolver ) {
    
    if ( !HSolver ) {
        return;
    }
    
    HDSDP_FREE(HSolver->rowRHS);
    
    for ( int iCone = 0; iCone < HSolver->nCones; ++iCone ) {
        HConeDestroy(&HSolver->HCones[iCone]);
    }
    
    HDSDP_FREE(HSolver->HCones);
    HKKTDestroy(&HSolver->HKKT);
    
    HDSDP_FREE(HSolver->dRowDual);
    HDSDP_FREE(HSolver->dRowDualStep);
    HDSDP_FREE(HSolver->dMinvRowRHS);
    HDSDP_FREE(HSolver->dMinvASinv);
    HDSDP_FREE(HSolver->dMinvASinvRdSinv);
    HDSDP_FREE(HSolver->dMinvASinvCSinv);
    HDSDP_FREE(HSolver->dHAuxiVec1);
    HDSDP_FREE(HSolver->dHAuxiVec2);
    
    HDSDP_ZERO(HSolver, hdsdp, 1);
    
    return;
}

extern void HDSDPDestroy( hdsdp **pHSolver ) {
    
    if ( !pHSolver ) {
        return;
    }
    
    HDSDPClear(*pHSolver);
    HDSDP_FREE(*pHSolver);
    
    hdsdp_printf("HDSDP ends. Exiting ... \n");
    
    return;
}
