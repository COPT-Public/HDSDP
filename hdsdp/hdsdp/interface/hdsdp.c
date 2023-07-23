/** @file hdsdp.c
 *  @brief A homogeneous dual-scaling interior point solver
 */

#ifdef HEADERPATH
#include "interface/hdsdp.h"
#include "interface/hdsdp_conic.h"
#include "interface/hdsdp_utils.h"
#include "interface/hdsdp_user_data.h"
#include "interface/hdsdp_schur.h"
#else
#include "hdsdp.h"
#include "hdsdp_conic.h"
#include "hdsdp_utils.h"
#include "hdsdp_user_data.h"
#include "hdsdp_schur.h"
#endif

struct hdsdp_solver_internal {
    
    char *coneModelName[100];
    
    /* Logging */
    int HPhaseA;
    
    /* User data */
    int nRows;
    double *rowRHS;
    
    /* Cones */
    int nCones;
    hdsdp_cone **HCones;
    
    /* KKT solver */
    hdsdp_kkt *HKKT;
    
    /* Step and vector */
    double *dRowDual;
    double *dRowDualStep;
    double dBarHsdTau;
    double dBarHsdTauStep;
    double *dMinvRowRHS;
    double *dMinvASinv;
    double *dMinvRdSinv;
    double *MinvASinvCSinv;
    double *dHAuxiVec1;
    double *dHAuxiVec2;
    
    /* Monitor */
    int nIterCount;
    double dBarrierMu;
    double dProxNorm;
    double dPotentialVal;
    double dBarrierVal;
    double dPStep;
    double dDStep;
    double dResidual;
    
    /* Convergence criterion */
    double dPotentialRho;
    double pObjVal;
    double dObjVal;
    double pInfeas;
    double dInfeas;
    double compl;
    double pInfeasRel;
    double dInfeasRel;
    double complRel;
    
    double dTimeBegin;
    
    hdsdp_status HStatus;
    
    /* Parameters */
    int HIntParams[20];
    int HIntFeatures[10];
    double HDblParams[20];
};

static void HDSDPIGetFeatures( hdsdp *HSolver ) {
    /* Get features for Heurstics */
    
    return;
}

static void HDSDPIPrintStatistics( hdsdp *HSolver ) {
    
    hdsdp_printf("Instance statistics: \n");
    
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
    
    HDSDP_INIT(HSolver->dMinvRdSinv, double, nRows);
    HDSDP_MEMCHECK(HSolver->dMinvRdSinv);
    
    HDSDP_INIT(HSolver->MinvASinvCSinv, double, nRows);
    HDSDP_MEMCHECK(HSolver->MinvASinvCSinv);
    
    HDSDP_INIT(HSolver->dHAuxiVec1, double, nRows);
    HDSDP_MEMCHECK(HSolver->dHAuxiVec1);
    
    HDSDP_INIT(HSolver->dHAuxiVec2, double, nRows);
    HDSDP_MEMCHECK(HSolver->dHAuxiVec2);
    
    HSolver->dBarrierMu = 1e+10;
    HSolver->compl = HDSDP_INFINITY;
    HSolver->complRel = HDSDP_INFINITY;
    
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
    hdsdp_printf("HDSDP: software for semi-definite programming. \n");
    HSolver->dTimeBegin = HUtilGetTimeStamp();
    
    HUtilGetDefaultParams(HSolver->HIntParams, HSolver->HDblParams);
    
    /* Process conic data */
    hdsdp_printf("Call conic data pre-solver. \n");
    for ( int iCone = 0; iCone < HSolver->nCones; ++iCone ) {
        HDSDP_CALL(HConeProcData(HSolver->HCones[iCone]));
        HDSDP_CALL(HConePresolveData(HSolver->HCones[iCone]));
    }
    hdsdp_printf("Conic pre-solver ends. Elapsed time: %f seconds \n", HUtilGetTimeStamp() - HSolver->dTimeBegin);
    
    /* Prepare KKT solver */
    hdsdp_printf("Call KKT pre-solver. \n");
    HDSDP_CALL(HKKTInit(HSolver->HKKT, HSolver->nRows, HSolver->nCones, HSolver->HCones));
    hdsdp_printf("KKT pre-solver ends. Elapsed time: %f seconds \n", HUtilGetTimeStamp() - HSolver->dTimeBegin);
    
    /* Get features */
    hdsdp_printf("Call feature pre-solver. \n");
    HDSDPIGetFeatures(HSolver);
    hdsdp_printf("Feature pre-solver ends. Elapsed time: %f seconds \n", HUtilGetTimeStamp() - HSolver->dTimeBegin);
    
    HUtilPrintParams(HSolver->HIntParams, HSolver->HDblParams);
    
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
    HDSDP_FREE(HSolver->dMinvRdSinv);
    HDSDP_FREE(HSolver->MinvASinvCSinv);
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
