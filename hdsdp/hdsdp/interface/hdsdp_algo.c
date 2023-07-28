#ifdef HEADERPATH
#include "interface/def_hdsdp.h"
#include "interface/hdsdp.h"
#include "interface/hdsdp_utils.h"
#include "interface/hdsdp_user_data.h"
#include "interface/hdsdp_file_io.h"
#include "interface/hdsdp_conic.h"
#include "interface/hdsdp_schur.h"
#include "interface/hdsdp_algo.h"
#else
#include "def_hdsdp.h"
#include "hdsdp.h"
#include "hdsdp_utils.h"
#include "hdsdp_user_data.h"
#include "hdsdp_file_io.h"
#include "hdsdp_conic.h"
#include "hdsdp_schur.h"
#include "hdsdp_algo.h"
#endif

#include <math.h>

/* Implement the major algorithms for HDSDP.
 
   To maximize the practical efficiency of HDSDP, we implement
   six major algorithms and one infeasibility (ray) detection module.
 
   The six algorithms are respectively
 
   1. Infeasible-start dual interior point method
   2. Dual embedding
   3. Dual potential reduction
   4. Primal-dual interior point method for linear programs
   5. Infeasible dual corrector step
   6. Feasible dual corrector step
 
   and the detection routine is invoked when
 
   1. the dual objective exceeds some large threshold value
   2. the final solution check shows large primal residual
 
   By default we start the algorithm with infeasible-start dual interior point method,
   which does not involve the potential expensive operation trace( C * S^-1 * C * S^-1 ).
   If we find dual infeasibility stagnates, self-dual embedding is invoked to detect dual infeasibility.
   If either infeasible start or embedding method finds a dual feasible solution, based on whether a primal
   solution is needed, the feasible solution is
 
   1. normalized and handed over to dual potential reduction for primal-dual optimality.
   2. kept in the algorithm that gets it for just dual optimality
 
 */

#ifdef HDSDP_ALGO_DEBUG
#define algo_debug printf
#else
#define algo_debug(...)
#endif

#define HDSDP_ALGO_DUAL_INFEAS     (0)
#define HDSDP_ALGO_DUAL_HSD        (1)
#define HDSDP_ALGO_DUAL_POTENTIAL  (2)
static void HDSDP_SetStart( hdsdp *HSolver, int SDPMethod, int dOnly ) {
    
    /* Set staring point for the dual iterations. */
    HDSDP_ZERO(HSolver->dRowDual, double, HSolver->nRows);
    HSolver->dBarHsdTau = 1.0;
    
    double objFroNorm = get_dbl_feature(HSolver, DBL_FEATURE_OBJFRONORM);
    objFroNorm = HDSDP_MAX(objFroNorm, 100.0);
    
    if ( SDPMethod == HDSDP_ALGO_DUAL_HSD ) {
        
        if ( dOnly ) {
            HSolver->dBarrierMu = 1e+08;
            HSolver->dResidual = - objFroNorm * get_dbl_param(HSolver, DBL_PARAM_DUALSTART);
        } else {
            HSolver->dBarrierMu = 1e+08;
            HSolver->dResidual = - objFroNorm * 1e+01;
        }
        
    } else if ( SDPMethod == HDSDP_ALGO_DUAL_INFEAS ) {
        
        /* Initialization for dual infeasible algorithm */
        HSolver->dResidual = - objFroNorm * get_dbl_param(HSolver, DBL_PARAM_DUALSTART) * get_dbl_feature(HSolver, DBL_FEATURE_OBJSCALING);
        HSolver->pInfeas = 1.0 + get_dbl_feature(HSolver, DBL_FEATURE_RHSFRONORM);
        HSolver->pObjInternal = get_dbl_param(HSolver, DBL_PARAM_POBJSTART);
        HSolver->dBarrierMu = HSolver->pObjInternal - HSolver->dObjInternal - \
                                HSolver->dResidual * get_dbl_param(HSolver, DBL_PARAM_TRXESTIMATE);
        HSolver->dBarrierMu = HSolver->dBarrierMu / HSolver->dAllConeDims;
        
    } else {
        hdsdp_printf("Invalid starting strategy. \n");
    }
    
    hdsdp_printf("Initialize with dual residual %3.1e\n", - HSolver->dResidual);
    
    for ( int iCone = 0; iCone < HSolver->nCones; ++iCone ) {
        HConeSetStart(HSolver->HCones[iCone], HSolver->dResidual);
    }
    
    return;
}

static void HDSDP_ResetStart( hdsdp *HSolver ) {
    
    HDSDP_ZERO(HSolver->dRowDual, double, HSolver->nRows);
    HSolver->dBarHsdTau = 1.0;
    
    double objFroNorm = get_dbl_feature(HSolver, DBL_FEATURE_OBJFRONORM);
    HSolver->dResidual = - HDSDP_MAX(objFroNorm, 1e+02);
    HSolver->dResidual = - HSolver->dResidual * 1e+10;
    HSolver->dResidual = HDSDP_MAX(HSolver->dResidual, -1e+15);
    
    hdsdp_printf("Reset with dual residual %3.1e\n", - HSolver->dResidual);
    for ( int iCone = 0; iCone < HSolver->nCones; ++iCone ) {
        HConeSetStart(HSolver->HCones[iCone], HSolver->dResidual);
    }
    
    return;
}

static void HDSDP_PrintHeader( hdsdp *HSolver, int SDPMethod ) {
    
    switch (SDPMethod) {
        case HDSDP_ALGO_DUAL_HSD:
            hdsdp_printf("HDSDP starts. Using self-dual method \n\n");
            hdsdp_printf("    %5s  %12s  %12s  %8s  %8s  %5s  %6s   %5s \n",
                         "nIter", "pObj", "dObj", "dInf", "Mu", "Step", "Tau", "T [H]");
            break;
        case HDSDP_ALGO_DUAL_INFEAS:
            hdsdp_printf("HDSDP starts. Using infeasible dual method \n\n");
            hdsdp_printf("    %5s  %12s  %12s  %8s  %8s  %5s  %6s   %5s \n",
                         "nIter", "pObj", "dObj", "dInf", "Mu", "Step", "|P|", "T [D]");
            break;
        case HDSDP_ALGO_DUAL_POTENTIAL:
            hdsdp_printf("HDSP re-starts. Using feasible dual method \n\n");
            hdsdp_printf("    %5s  %12s  %12s  %8s  %8s  %5s  %6s   %5s\n",
                         "nIter", "pObj", "dObj", "pInf", "Mu", "Step", "|P|", "T [P]");
            break;
        default:
            assert( 0 );
            break;
    }
    
    return;
}

static void HDSDP_PrintLog( hdsdp *HSolver, int SDPMethod ) {
    
    double dRhsScal = get_dbl_feature(HSolver, DBL_FEATURE_RHSSCALING);
    double dObjScal = get_dbl_feature(HSolver, DBL_FEATURE_OBJSCALING);
    double elapsedTime = HUtilGetTimeStamp() - HSolver->dTimeBegin;
    
    /* Get infeasibilities */
    double nSumCones = (double) get_int_feature(HSolver, INT_FEATURE_N_SUMCONEDIMS);
    double pdObjScal = 1.0 / (dRhsScal * dObjScal * HSolver->dBarHsdTau);
    HSolver->dInfeas =  sqrt(nSumCones) * fabs(HSolver->dResidual) / (dRhsScal * HSolver->dBarHsdTau);
    
    HSolver->dObjInternal = 0.0;
    for ( int iRow = 0; iRow < HSolver->nRows; ++iRow ) {
        HSolver->dObjInternal += HSolver->rowRHS[iRow] * HSolver->dRowDual[iRow];
    }
    
    HSolver->dObjVal = HSolver->dObjInternal * pdObjScal;
    HSolver->pObjVal = HSolver->pObjInternal * pdObjScal;
    
    switch (SDPMethod) {
        case HDSDP_ALGO_DUAL_HSD:
            hdsdp_printf("    %5d  %+12.5e  %+12.5e  %8.2e  %8.2e  %5.2f  %5.1e  %4.1f\n", HSolver->nIterCount + 1,
                   HDSDP_INFINITY, HSolver->dObjVal, HSolver->dInfeas,
                   HSolver->dBarrierMu, HSolver->dDStep, HSolver->dBarHsdTau, elapsedTime);
            break;
        case HDSDP_ALGO_DUAL_INFEAS:
            hdsdp_printf("    %5d  %+12.5e  %+12.5e  %8.2e  %8.2e  %5.2f  %5.1e  %4.1f \n", HSolver->nIterCount + 1,
                   HSolver->pObjVal * pdObjScal, HSolver->dObjVal * pdObjScal, HSolver->dInfeas, HSolver->dBarrierMu,
                   HSolver->dDStep, HSolver->dProxNorm, elapsedTime);
            break;
        case HDSDP_ALGO_DUAL_POTENTIAL:
            hdsdp_printf("    %5d  %+12.5e  %+12.5e  %8.2e  %8.2e  %5.2f  %5.1e  %4.1f \n", HSolver->nIterCount + 1,
                   HSolver->pObjVal * pdObjScal, HSolver->dObjVal * pdObjScal, HSolver->pInfeas, HSolver->dBarrierMu,
                   HSolver->dDStep, HSolver->dProxNorm, elapsedTime);
            break;
        default:
            assert( 0 );
            break;
    }
    
    return;
}

static int HDSDP_CheckIsInterior( hdsdp *HSolver, double barHsdTau, double *rowDual ) {
    
    int isInterior = 0;
    
    for ( int iCone = 0; iCone < HSolver->nCones; ++iCone ) {
        HConeCheckIsInterior(HSolver->HCones[iCone], 1.0, rowDual, &isInterior);
        if ( !isInterior ) {
            return 0;
        }
    }
    
    if ( HSolver->whichMethod != HDSDP_ALGO_DUAL_HSD ) {
        HConeCheckIsInterior(HSolver->HBndCone, 1.0, rowDual, &isInterior);
    }
    
    if ( !isInterior ) {
        return 0;
    }
    
    return 1;
}

static hdsdp_retcode HDSDP_GetLogBarrier( hdsdp *HSolver, double barHsdTau, double *rowDual, int whichBuffer, double *logdet ) {
    
    hdsdp_retcode retcode = HDSDP_RETCODE_OK;
    
    double dLogDet = 0.0;
    double dBarrierVal = 0.0;
    
    for ( int iCone = 0; iCone < HSolver->nCones; ++iCone ) {
        HDSDP_CALL(HConeGetLogBarrier(HSolver->HCones[iCone], barHsdTau, rowDual, whichBuffer, &dLogDet));
        dBarrierVal -= dLogDet;
    }
    
    if ( HSolver->whichMethod != HDSDP_ALGO_DUAL_HSD ) {
        HDSDP_CALL(HConeGetLogBarrier(HSolver->HBndCone, barHsdTau, rowDual, whichBuffer, &dLogDet));
        dBarrierVal -= dLogDet;
    }
    
    *logdet = dBarrierVal;
    
exit_cleanup:
    return retcode;
}

static void HDSDP_SetResidual( hdsdp *HSolver, double dResi ) {
        
    for ( int iCone = 0; iCone < HSolver->nCones; ++iCone ) {
        HConeReduceResi(HSolver->HCones[iCone], dResi);
    }
    
    return;
}

/* Macros in MATLAB style */
#define b    HSolver->rowRHS
#define y    HSolver->dRowDual
#define dy   HSolver->dRowDualStep
#define d1   HSolver->dMinvRowRHS
#define d2   HSolver->dMinvASinv
#define d3   HSolver->dMinvASinvRdSinv
#define d4   HSolver->dMinvASinvCSinv
#define dd1  HSolver->dHAuxiVec1
#define bTy  HSolver->dObjInternal
#define mu   HSolver->dBarrierMu
#define tau  HSolver->dBarHsdTau
#define dtau HSolver->dBarHsdTauStep
static void HDSDP_HSD_BuildStep( hdsdp *HSolver ) {
    
    /* Prepare intermediate results */
    double dHsdTraceCSinvCSinv = 0.0;
    double dHsdTraceCSinv = 0.0;
    double dHsdTraceCSinvRdCSinv = 0.0;
    
    HSolver->dBarHsdTauStep = 0.0;
    
    /* First get dual objective */
    double dOldObjVal = HSolver->dObjInternal;
    bTy = 0.0;
    
    for ( int iRow = 0; iRow < HSolver->nRows; ++iRow ) {
        bTy += b[iRow] * y[iRow];
    }
    
    HSolver->dObjImprove = bTy - dOldObjVal;
    
    /* Then compute dd1 = b - mu * ASinvCSinv */
    double *dASinvCSinvVec = HSolver->HKKT->dASinvCSinvVec;
    
    for ( int iRow = 0; iRow < HSolver->nRows; ++iRow ) {
        dd1[iRow] = b[iRow] - mu * dASinvCSinvVec[iRow];
    }
    
    /* Compute dTau */
    HKKTExport(HSolver->HKKT, NULL, NULL, NULL,
               &dHsdTraceCSinvCSinv, &dHsdTraceCSinv, &dHsdTraceCSinvRdCSinv, NULL);
    
    double dTauStepEnumerator = - bTy + mu / tau + mu * (dHsdTraceCSinv - dHsdTraceCSinvRdCSinv);
    double dTauStepDenominator = mu * dHsdTraceCSinvCSinv + mu / (tau * tau);
    double tauOverMu = tau / mu;
    
    for ( int iRow = 0; iRow < HSolver->nRows; ++iRow ) {
        dTauStepEnumerator -= dd1[iRow] * (d1[iRow] * tauOverMu - d2[iRow] + d3[iRow]);
        dTauStepDenominator += dd1[iRow] * (d1[iRow] / mu + d4[iRow]);
    }
    
    if ( fabs(dTauStepDenominator) < 1e-12 ) {
        dtau = 0.0;
    } else {
        dtau = dTauStepEnumerator / dTauStepDenominator;
    }
    
    /* Get dy assembled */
    for ( int iRow = 0; iRow < HSolver->nRows; ++iRow ) {
        dy[iRow] = d1[iRow] * (tau + dtau) / mu + d4[iRow] * dtau - d2[iRow] + d3[iRow];
    }
    
    return;
}

static hdsdp_retcode HDSDP_HSD_RatioTest( hdsdp *HSolver, double *dMaxDist ) {
    
    hdsdp_retcode retcode = HDSDP_RETCODE_OK;
    
    double dMaxStep = HDSDP_INFINITY;
    double dStep = 0.0;
    
    /* Get step of tau  */
    dStep = HSolver->dBarHsdTau / HSolver->dBarHsdTauStep;
    if ( dStep < 0.0 ) {
        dMaxStep = HDSDP_MIN(dMaxStep, -dStep);
    }
    
    /* Aggregate step sizes from different cones */
    for ( int iCone = 0; iCone < HSolver->nCones; ++iCone ) {
        /* We eleminate dual infeasibility with full Newton's step */
        HDSDP_CALL(HConeRatioTest(HSolver->HCones[iCone], HSolver->dBarHsdTauStep,
                                  HSolver->dRowDualStep, 1.0, BUFFER_DUALVAR, &dStep));
        
        dMaxStep = HDSDP_MIN(dMaxStep, dStep);
    }
    
    if ( HSolver->whichMethod != HDSDP_ALGO_DUAL_HSD ) {
        HDSDP_CALL(HConeRatioTest(HSolver->HBndCone, HSolver->dBarHsdTauStep,
                                  HSolver->dRowDualStep, 1.0, BUFFER_DUALVAR, &dStep));
    }
    
    if ( dMaxStep < 1e-02 ) {
        HSolver->nSmallStep += 1;
        if ( HSolver->nSmallStep > 2 ) {
            hdsdp_printf("HDSDP stagates at the cone boundary. \n");
            retcode = HDSDP_RETCODE_FAILED;
            goto exit_cleanup;
        }
    }
    
exit_cleanup:
    *dMaxDist = dMaxStep;
    return retcode;
}

static int HDSDP_PhaseA_ProxMeasure( hdsdp *HSolver ) {
    
    int isPFeasible = 0;
    double dTraceSinv = 0.0;
    double pObjNew = HSolver->dObjInternal;
    
    HKKTExport(HSolver->HKKT, NULL, NULL, NULL, NULL, NULL, NULL, &dTraceSinv);
    for ( int iRow = 0; iRow < HSolver->nRows; ++iRow ) {
        HSolver->dHAuxiVec1[iRow] = \
            HSolver->dMinvRowRHS[iRow] / HSolver->dBarrierMu - HSolver->dMinvASinv[iRow];
        HSolver->dHAuxiVec2[iRow] = \
            HSolver->rowRHS[iRow] / HSolver->dBarrierMu - HSolver->HKKT->dASinvVec[iRow];
    }
    
    HSolver->dProxNorm = 0.0;
    for ( int iRow = 0; iRow < HSolver->nRows; ++iRow ) {
        HSolver->dProxNorm += HSolver->dHAuxiVec1[iRow] * HSolver->dHAuxiVec2[iRow];
    }
    
    if ( HSolver->dProxNorm < 0.0 ) {
        algo_debug("Proximity norm is negative \n");
        HSolver->dProxNorm = 1.0;
        isPFeasible = 0;
        return isPFeasible;
    } else {
        HSolver->dProxNorm = sqrt(HSolver->dProxNorm);
    }
    
    for ( int iRow = 0; iRow < HSolver->nRows; ++iRow ) {
        HSolver->dHAuxiVec2[iRow] = HSolver->dHAuxiVec1[iRow] - HSolver->dRowDual[iRow];
    }
    
    /* Check feasibility B <- dResiCoef * I + dACoefScal * A' * y + C * dCCoef */
    HConeCheckIsInteriorExpert(HSolver->HBndCone, 1.0, 1.0, HSolver->dHAuxiVec2,
                               -HSolver->dResidual, BUFFER_DUALCHECK, &isPFeasible);
    
    for ( int iCone = 0; iCone < HSolver->nCones && isPFeasible; ++iCone ) {
        HConeCheckIsInteriorExpert(HSolver->HCones[iCone], 1.0, 1.0, HSolver->dHAuxiVec2,
                                   -HSolver->dResidual, BUFFER_DUALCHECK, &isPFeasible);
        if ( !isPFeasible ) {
            break;
        }
    }
    
    if ( isPFeasible ) {
        
        double dRelGap = 0.0;
        for ( int iRow = 0; iRow < HSolver->nRows; ++iRow ) {
            dRelGap += HSolver->dHAuxiVec1[iRow] * \
                        (HSolver->HKKT->dASinvRdSinvVec[iRow] + HSolver->HKKT->dASinvVec[iRow]);
        }
        
        dRelGap += HSolver->dAllConeDims;
        dRelGap += dTraceSinv * HSolver->dResidual;
        pObjNew += dRelGap * HSolver->dBarrierMu;
        
        if ( dRelGap < 0 ) {
            algo_debug("Dropped invalid solution %10.6e \n", dRelGap);
            isPFeasible = 0;
        } else {
            algo_debug("Found new primal bound %10.6e \n", pObjNew);
            HSolver->pObjInternal = pObjNew;
            /* Record solution */
            HSolver->dInaccBarrierMaker = HSolver->dBarrierMu;
            HDSDP_MEMCPY(HSolver->dInaccRowDualMaker, HSolver->dRowDual, double, HSolver->nRows);
            HDSDP_MEMCPY(HSolver->dInaccRowDualStepMaker, HSolver->dHAuxiVec2, double, HSolver->nRows);
        }
    }
    
    return isPFeasible;
}

static hdsdp_retcode HDSDP_PhaseA_AdaptiveResi( hdsdp *HSolver, double *dResiRate ) {
    
    hdsdp_retcode retcode = HDSDP_RETCODE_OK;
    
    /* First do a ratio test to get alpha_c such that s + alpha_c * A^T * dy_c >= 0 */
    double dAlphaC = 0.0;
    double dAlphaInf = 0.0;
    double dAdaRatio = 0.0;
    double dMaxStep = HDSDP_INFINITY;
    double dStep = 0.0;
    
    for ( int iRow = 0; iRow < HSolver->nRows; ++iRow ) {
        HSolver->dHAuxiVec1[iRow] = -HSolver->dMinvASinv[iRow];
    }
    
    for ( int iCone = 0; iCone < HSolver->nCones; ++iCone ) {
        HDSDP_CALL(HConeRatioTest(HSolver->HCones[iCone], 0.0, HSolver->dHAuxiVec1, 0.0, BUFFER_DUALVAR, &dStep));
        dMaxStep = HDSDP_MIN(dMaxStep, dStep);
    }
    
//    HDSDP_CALL(HConeRatioTest(HSolver->HBndCone, 0.0, HSolver->dHAuxiVec1, 0.0, BUFFER_DUALVAR, &dStep));
//    dMaxStep = HDSDP_MIN(dMaxStep, dStep);
    
    dAlphaC = HDSDP_MIN(0.98 * dMaxStep, 1.0) ;
    
    /* Perform line search to guarantee validity */
    int isInterior = 0;
    while ( !isInterior && dAlphaC > 1e-02 * dMaxStep ) {
        
        for ( int iCone = 0; iCone < HSolver->nCones; ++iCone ) {
            HDSDP_CALL(HConeAddStepToBufferAndCheck(HSolver->HCones[iCone], dAlphaC, BUFFER_DUALCHECK, &isInterior));
            if ( !isInterior ) {
                break;
            }
        }
        
        HDSDP_CALL(HConeAddStepToBufferAndCheck(HSolver->HBndCone, dAlphaC, BUFFER_DUALCHECK, &isInterior));
        
        if ( !isInterior ) {
            dAlphaC = dAlphaC * 0.8;
        }
    }
    
    /* Now we have alpha_c in hand. Do the next ratio test for
        s' + alpha * (rd - A^T dy_r) */
    dMaxStep = HDSDP_INFINITY;
    for ( int iCone = 0; iCone < HSolver->nCones; ++iCone ) {
        HDSDP_CALL(HConeRatioTest(HSolver->HCones[iCone], 0.0, HSolver->dMinvASinvRdSinv, 1.0, BUFFER_DUALCHECK, &dStep));
        dMaxStep = HDSDP_MIN(dMaxStep, dStep);
    }
    
    HDSDP_CALL(HConeRatioTest(HSolver->HBndCone, 0.0, HSolver->dMinvASinvRdSinv, 1.0, BUFFER_DUALCHECK, &dStep));
    dMaxStep = HDSDP_MIN(dMaxStep, dStep);
    
    dAlphaInf = dMaxStep;
    dAdaRatio = dAlphaInf / dAlphaC;
    dAdaRatio = HDSDP_MIN(0.98 * dAdaRatio, 1.0);
    
    /* Finally adjust the adaptive ratio based on proximity norm */
    if ( HSolver->dProxNorm < 1.0 ) {
        dAdaRatio = HDSDP_MAX(0.9, dAdaRatio);
    } else if ( HSolver->dProxNorm < 10.0 ) {
        dAdaRatio = HDSDP_MAX(0.3, dAdaRatio);
    } else {
        dAdaRatio = HDSDP_MAX(0.1, dAdaRatio);
    }
    
    *dResiRate = dAdaRatio;
 
exit_cleanup:
    return retcode;
}

static void HDSDP_Infeasible_BuildStep( hdsdp *HSolver, double dAdaRatio ) {
    
    /* Assemble dual step dy = 1/mu * dy_b - dy_c + gamma * dy_r */
    for ( int iRow = 0; iRow < HSolver->nRows; ++iRow ) {
        HSolver->dRowDualStep[iRow] = \
        HSolver->dMinvRowRHS[iRow] / HSolver->dBarrierMu - \
        HSolver->dMinvASinv[iRow] + dAdaRatio * HSolver->dMinvASinvRdSinv[iRow];
    }
    
    return;
}

static hdsdp_retcode HDSDP_Infeasible_RatioTest( hdsdp *HSolver, double dAdaRatio, double *dMaxDist ) {
    
    hdsdp_retcode retcode = HDSDP_RETCODE_OK;
    
    double dMaxStep = HDSDP_INFINITY;
    double dStep = 0.0;
    
    /* Aggregate step sizes from different cones */
    for ( int iCone = 0; iCone < HSolver->nCones; ++iCone ) {
        /* We eleminate dual infeasibility with full Newton's step */
        HDSDP_CALL(HConeRatioTest(HSolver->HCones[iCone], 0.0,
                                  HSolver->dRowDualStep, dAdaRatio, BUFFER_DUALVAR, &dStep));
        dMaxStep = HDSDP_MIN(dMaxStep, dStep);
    }
    
    HDSDP_CALL(HConeRatioTest(HSolver->HBndCone, 0.0, HSolver->dRowDualStep,
                              dAdaRatio, BUFFER_DUALVAR, &dStep));
    dMaxStep = HDSDP_MIN(dMaxStep, dStep);
    
    if ( dMaxStep < 1e-02 ) {
        HSolver->nSmallStep += 1;
        if ( HSolver->nSmallStep > 2 ) {
            retcode = HDSDP_RETCODE_FAILED;
            goto exit_cleanup;
        }
    }
    
exit_cleanup:
    *dMaxDist = dMaxStep;
    return retcode;
}

static hdsdp_retcode HDSDP_Infeasible_Corrector( hdsdp *HSolver, int lastStep ) {
    
    hdsdp_retcode retcode = HDSDP_RETCODE_OK;
    
    /* Implement the centrality corrector for HDSDP
       The corrector step resembles the main infeasible start dual method, but the Schur complement
       is reused. Also a line-search is applied to reduce the barrier function,
     
       If corrector is called in the last step, the algorithm will perform aggressive update to
       reduce the barrier function, giving stronger centrality guarantees.
     */
    
    (void) lastStep;
    int nMaxCorr = get_int_param(HSolver, INT_PARAM_CORRECTORA);
    
    algo_debug("Starting infeasible dual corrector \n");
    
    int isInterior = 0;
    
    isInterior = HDSDP_CheckIsInterior(HSolver, 1.0, HSolver->dRowDual);
    
    if ( !isInterior ) {
        algo_debug("Incumbent dual solution is infeasible");
        retcode = HDSDP_RETCODE_FAILED;
        goto exit_cleanup;
    }
    
    /* Compute value of the barrier function */
    double dBarrierVal = 0.0;
    double dNewBarrierVal = 0.0;
    double dMaxStep = 0.0;
    double dStep = 0.0;
    double dAlphaC = 0.0;
    double dAlphaInf = 0.0;
    double dAdaRatio = 0.0;
    double dAdaRatioMax = 0.8;
    
    
    HDSDP_CALL(HDSDP_GetLogBarrier(HSolver, 0.0, NULL, BUFFER_DUALVAR, &dBarrierVal));
    
    for ( int nCorr = 0; nCorr < nMaxCorr; ++nCorr ) {
        
        if ( !HSolver->dResidual ) {
            break;
        }
        
        /* Build the Schur complement elements */
        HDSDP_CALL(HKKTBuildUp(HSolver->HKKT, KKT_TYPE_CORRECTOR));
        HKKTBuildUpExtraCone(HSolver->HKKT, HSolver->HBndCone, KKT_TYPE_CORRECTOR);
        
        /* Solve for the directions */
        HKKTExport(HSolver->HKKT, HSolver->dMinvASinv, HSolver->dMinvASinvRdSinv,
                   NULL, NULL, NULL, NULL, NULL);
        
        HDSDP_CALL(HKKTSolve(HSolver->HKKT, HSolver->dMinvASinv, NULL));
        
        if ( dAdaRatioMax ) {
            HDSDP_CALL(HKKTSolve(HSolver->HKKT, HSolver->dMinvASinvRdSinv, NULL));
        }
        
        /* Build direction */
        for ( int iRow = 0; iRow < HSolver->nRows; ++iRow ) {
            HSolver->dRowDualStep[iRow] = - HSolver->dMinvASinv[iRow];
        }
        
        /* Ratio test */
        HDSDP_Infeasible_RatioTest(HSolver, 0.0, &dMaxStep);
        dMaxStep = HDSDP_MIN(0.8 * dMaxStep, 1.0);
        
        /* Guarantee feasibility */
        while ( 1 ) {
            /* Update dual */
            for ( int iRow = 0; iRow < HSolver->nRows; ++iRow ) {
                HSolver->dHAuxiVec1[iRow] = HSolver->dRowDual[iRow] + dMaxStep * HSolver->dRowDualStep[iRow];
            }
            
            isInterior = HDSDP_CheckIsInterior(HSolver, 1.0, HSolver->dHAuxiVec1);
            
            if ( !isInterior ) {
                dMaxStep *= 0.5;
                break;
            }
            
            if ( isInterior || dMaxStep < 5e-03 ) {
                break;
            }
        }
        
        if ( dMaxStep < 5e-03 ) {
            break;
        }
        
        /* Compute the barrier function */
        HDSDP_CALL(HDSDP_GetLogBarrier(HSolver, 0.0, NULL, BUFFER_DUALVAR, &dNewBarrierVal));
        
        if ( dNewBarrierVal > dBarrierVal ) {
            dMaxStep *= 0.5;
            for ( int iRow = 0; iRow < HSolver->nRows; ++iRow ) {
                HSolver->dHAuxiVec1[iRow] = HSolver->dRowDual[iRow] + dMaxStep * HSolver->dRowDualStep[iRow];
            }
            isInterior = HDSDP_CheckIsInterior(HSolver, 1.0, HSolver->dHAuxiVec1);
            assert( isInterior );
        }
        
        dAlphaC = dMaxStep;
        
        /* Reduce infeasibility */
        dMaxStep = HDSDP_INFINITY;
        for ( int iCone = 0; iCone < HSolver->nCones; ++iCone ) {
            HDSDP_CALL(HConeRatioTest(HSolver->HCones[iCone], 0.0, HSolver->dMinvASinvRdSinv, 1.0, BUFFER_DUALVAR, &dStep));
            dMaxStep = HDSDP_MIN(dMaxStep, dStep);
        }
        HDSDP_CALL(HConeRatioTest(HSolver->HBndCone, 0.0, HSolver->dMinvASinvRdSinv, 1.0, BUFFER_DUALVAR, &dStep));
        dMaxStep = HDSDP_MIN(dMaxStep, dStep);
        
        
        dAlphaInf = dMaxStep;
        dAdaRatio = dAlphaInf / dAlphaC;
        dAdaRatio = HDSDP_MIN(1.0, dAdaRatioMax * dAdaRatio);
        
        /* Guarantee feasibility after further reduction */
        double dResi = HSolver->dResidual;
        for ( int iTry = 0; dAdaRatio > 1e-03; ++iTry ) {
            HDSDP_SetResidual(HSolver, dResi * (1 - dAlphaC * dAdaRatio));
            for ( int iRow = 0; iRow < HSolver->nRows; ++iRow ) {
                HSolver->dHAuxiVec1[iRow] = HSolver->dRowDual[iRow] + dAlphaC * \
                                            (dAdaRatio * HSolver->dMinvASinvRdSinv[iRow] - HSolver->dMinvASinv[iRow]);
            }
            
            isInterior = HDSDP_CheckIsInterior(HSolver, 1.0, HSolver->dHAuxiVec1);
            
            if ( isInterior ) {
                break;
            }
            
            dAdaRatio *= 0.8;
        }
        
        /* Adjust parameter heuristically */
        if ( dAlphaC * dAdaRatio < 5e-04 ) {
            dAdaRatioMax = 0.0;
        } else if ( dAlphaC * dAdaRatio < 0.1 ) {
            dAdaRatioMax = dAdaRatioMax * 0.9;
        }
        
        if ( dAlphaC * dAdaRatio > 0.8 ) {
            HSolver->dBarrierMu *= 0.8;
            dAdaRatioMax = HDSDP_MIN(dAdaRatioMax * 2.0, 0.9);
        } else if ( dAlphaC * dAdaRatio > 0.3 ) {
            HSolver->dBarrierMu *= 0.95;
            dAdaRatioMax = HDSDP_MIN(dAdaRatioMax * 2.0, 0.8);
        }
        
        algo_debug("Infeasible corrector step %d. AdaRatio %f\n", nCorr + 1, dAdaRatio);
        
        /* Save progress */
        if ( dAdaRatioMax ) {
            HSolver->dResidual = dResi * (1 - dAlphaC * dAdaRatio);
            HDSDP_MEMCPY(HSolver->dRowDual, HSolver->dHAuxiVec1, double, HSolver->nRows);
            dBarrierVal = -HDSDP_INFINITY;
        } else {
            isInterior = HDSDP_CheckIsInterior(HSolver, 1.0, HSolver->dRowDual);
            assert( isInterior );
        }
        
        if ( dAdaRatioMax == 0.0 ) {
            break;
        }
    }
    
exit_cleanup:
    return retcode;
}

static hdsdp_retcode HDSDP_PhaseA_BarInfeasSolve( hdsdp *HSolver, int dOnly ) {
    
    hdsdp_retcode retcode = HDSDP_RETCODE_OK;
    /*
     Implement the infeasible-start dual potential reduction method.
     An adaptive stepsize strategy is applied to find eliminate dual infeasibility efficiently.
     Also the output from the KKT solver will be reused as if in the corrector method of DSDP.
     
     Primal solution, if available. Will be recorded by this phase.
     If it is hard to eliminate dual infeasibility, the embedding will be invoked to provide a certificate.
        
     The Phase A of HDSDP implements an infeasible-start dual scaling interior point method.
     Its main idea is similar to the primal-dual method and tries to take Newton's step towards the
     perturbed KKT system
   
                    A * (x + dx) = b
                  -A^T * dy - ds = -gamma * rd
      mu * s^-1 * ds * s^-1 + dx = mu * s^-1 - x,
     
     where gamma is what we call a "residual reduction" rate.
   
     Taking gamma -> 0.0 puts less emphasis on dual infeasibility and more on centrality,
     while gamma -> 0.0 recovers the full Newton's step
   
     The goal of Phase A is to
   
     1. eliminate dual infeasibility as fast as possible
     2. maintain centrality for the Phase B to start
     3. get a rough primal solution value bound
     
     To fulfill the above two purposes, HDSDP implements a heuristic strategy to adaptively determine the
     residual reduction rate. We briefly describe it here.
    
     Given a dual infeasible solution (y, S) such that rd = -A^T * y - s + c,
     the dual update (with step length alpha) can be written as
   
     s + alpha * ds = s + alpha * A^T * dy_c + alpha * gamma * (rd - A^T * dy_r)
                    
     where M * dy_c = Asinv and M * dy_r = Asinvrdsinv
   
     Phase A would first take gamma = 0.0 and pretend a corrector step is taken. Then we get
     alpha_c such that s + alpha_c * A^T * dy_c >= 0, where alpha_c can either take the maximum step length
     towards conic boundary, or we can further do line-search to reduce the barrier function
     log det (s + alpha_c * A^T * dy_c).
   
     After getting alpha_c and s' = [s + alpha_c * A^T * dy_c], we compute the maximum step length alpha_inf such that
     
       [s + alpha * A^T * dy_c] + alpha_inf * (rd - A^T * dy_r)
     = s' + alpha_inf * (rd - A^T * dy_r),
     
     and the correspondence alpha_inf = alpha_c * gamma would provide gamma = alpha_inf / alpha_c.
   
     Finally, we plug in gamma in the Newton's system and compute the maximum step length, with a
     0.95 shrinkage as our final step.
     
     */
    
    (void) dOnly;
    
    HSolver->whichMethod = HDSDP_ALGO_DUAL_INFEAS;
    
    int nMaxIter = get_int_param(HSolver, INT_PARAM_MAXITER);
    double dAbsfeasTol = get_dbl_param(HSolver, DBL_PARAM_ABSFEASTOL);
    double dRelfeasTol = get_dbl_param(HSolver, DBL_PARAM_RELFEASTOL);
    double dTimeLimit = get_dbl_param(HSolver, DBL_PARAM_TIMELIMIT);
    double nSumCones = (double) get_int_feature(HSolver, INT_FEATURE_N_SUMCONEDIMS);
    double dObjScal = get_dbl_feature(HSolver, DBL_FEATURE_OBJSCALING);
    double dObjOneNorm = get_dbl_feature(HSolver, DBL_FEATURE_OBJONENORM);
    double dFeasTol = HDSDP_MIN(dAbsfeasTol, dRelfeasTol * (1 + dObjOneNorm));
    double dAdaRatio = 0.0;
    dFeasTol = dFeasTol * dObjScal / sqrt(nSumCones);
    
    /* We are starting from scratch */
    HDSDP_SetStart(HSolver, HDSDP_ALGO_DUAL_INFEAS, 0);
    
    int isInteriorPoint = 0;
    int pObjFound = 0;
    
    /* Initial check */
    for ( int iCone = 0; iCone < HSolver->nCones; ++iCone ) {
        HDSDP_CALL(HConeCheckIsInterior(HSolver->HCones[iCone], HSolver->dBarHsdTau,
                                        HSolver->dRowDual, &isInteriorPoint));
        if ( !isInteriorPoint ) {
            break;
        }
        
        HConeCheckIsInterior(HSolver->HBndCone, HSolver->dBarHsdTau, HSolver->dRowDual, &isInteriorPoint);
    }
    
    if ( !isInteriorPoint ) {
        hdsdp_printf("Initial point is not in the cone. Adding slack value.\n");
        HSolver->dResidual *= 100.0;
        HDSDP_ResetStart(HSolver);
    }
        
    HDSDP_PrintHeader(HSolver, HDSDP_ALGO_DUAL_INFEAS);
    
    /* Enter the main iteration */
    while ( 1 ) {
        
        /* Build up the Schur complement */
        HDSDP_CALL(HKKTBuildUp(HSolver->HKKT, KKT_TYPE_INFEASIBLE));
        /* Schur complement now contains an extra bound cone */
        HDSDP_CALL(HKKTBuildUpExtraCone(HSolver->HKKT, HSolver->HBndCone, KKT_TYPE_INFEASIBLE));
        /* Export the information needed */
        HKKTExport(HSolver->HKKT, HSolver->dMinvASinv, HSolver->dMinvASinvRdSinv,
                   NULL, NULL, NULL, NULL, NULL);
        
        /* Factorize the KKT system */
        HDSDP_CALL(HKKTFactorize(HSolver->HKKT));
        
        /* Solve the KKT system */
        HDSDP_CALL(HKKTSolve(HSolver->HKKT, HSolver->rowRHS, HSolver->dMinvRowRHS));
        HDSDP_CALL(HKKTSolve(HSolver->HKKT, HSolver->dMinvASinv, NULL));
        HDSDP_CALL(HKKTSolve(HSolver->HKKT, HSolver->dMinvASinvRdSinv, NULL));
        
        /* Compute proximity norm and verify primal feasibility */
        pObjFound += HDSDP_PhaseA_ProxMeasure(HSolver);
        
        /* Update barrier parameter */
        double dBarrierTarget = HSolver->pObjInternal - HSolver->dObjInternal;
        dBarrierTarget = dBarrierTarget - HSolver->dResidual * get_dbl_param(HSolver, DBL_PARAM_TRXESTIMATE);
        dBarrierTarget = dBarrierTarget / (5.0 * HSolver->dAllConeDims);
        
        if ( HSolver->dProxNorm < 1.0 ) {
            HSolver->dBarrierMu = HSolver->dBarrierMu * 0.005;
        } else if ( HSolver->dProxNorm < 5.0 ) {
            HSolver->dBarrierMu = HSolver->dBarrierMu * 0.01;
            HSolver->dBarrierMu = HDSDP_MAX(HSolver->dBarrierMu, dBarrierTarget);
        } else if ( HSolver->dProxNorm < 10.0 ) {
            HSolver->dBarrierMu = HSolver->dBarrierMu * 0.1;
            HSolver->dBarrierMu = HDSDP_MAX(HSolver->dBarrierMu, dBarrierTarget);
        } else {
            HSolver->dBarrierMu = HSolver->dBarrierMu * 0.95;
            HSolver->dBarrierMu = HDSDP_MAX(HSolver->dBarrierMu, dBarrierTarget);
        }
        
        /* Compute adaptive ration gamma */
        HDSDP_CALL(HDSDP_PhaseA_AdaptiveResi(HSolver, &dAdaRatio));
        algo_debug("Adaptive ratio: %f\n", dAdaRatio);
        
        /* Assemble direction */
        HDSDP_Infeasible_BuildStep(HSolver, dAdaRatio);
        
        /* Final ratio test */
        HDSDP_Infeasible_RatioTest(HSolver, dAdaRatio, &HSolver->dDStep);
        HSolver->dDStep = HDSDP_MIN(0.95 * HSolver->dDStep, 1.0);
        algo_debug("Stepsize: %f\n", HSolver->dDStep);
        
        /* Take step */
        for ( int iRow = 0; iRow < HSolver->nRows; ++iRow ) {
            HSolver->dRowDual[iRow] += HSolver->dDStep * HSolver->dRowDualStep[iRow];
        }
        HSolver->dResidual = HSolver->dResidual * (1.0 - dAdaRatio * HSolver->dDStep);
        HDSDP_SetResidual(HSolver, HSolver->dResidual);
        
        /* Restart if there is no valid primal bound */
        if ( HSolver->nIterCount == 3 && !pObjFound ) {
            hdsdp_printf("Increasing dual infeasibility \n");
            HDSDP_ResetStart(HSolver);
            for ( int iCone = 0; iCone < HSolver->nCones; ++iCone ) {
                HDSDP_CALL(HConeCheckIsInterior(HSolver->HCones[iCone], HSolver->dBarHsdTau,
                                                HSolver->dRowDual, &isInteriorPoint));
                if ( !isInteriorPoint ) {
                    retcode = HDSDP_RETCODE_FAILED;
                    goto exit_cleanup;
                }
            }
        }
        
        /* Print log */
        HDSDP_CALL(HDSDP_Infeasible_Corrector(HSolver, 0));
        HDSDP_PrintLog(HSolver, HDSDP_ALGO_DUAL_INFEAS);
        
        /* Convergence check */
        if ( fabs(HSolver->dResidual) < dFeasTol ) {
            HSolver->HStatus = HDSDP_DUAL_FEASIBLE;
            break;
        }
        
        if ( HUtilCheckCtrlC() ) {
            HSolver->HStatus = HDSDP_USER_INTERRUPT;
            break;
        }
        
        if ( HSolver->dDStep < 1e-03 ) {
            HSolver->HStatus = HDSDP_SUSPECT_INFEAS_OR_UNBOUNDED;
        }
        
        if ( HUtilGetTimeStamp() - HSolver->dTimeBegin >= dTimeLimit ) {
            HSolver->HStatus = HDSDP_TIMELIMIT;
            break;
        }
        
        HSolver->nIterCount += 1;
        if ( HSolver->nIterCount >= nMaxIter ) {
            HSolver->HStatus = HDSDP_MAXITER;
            break;
        }
    }
    
exit_cleanup:
    return retcode;
}

static hdsdp_retcode HDSDP_PhaseA_BarHsdSolve( hdsdp *HSolver, int dOnly ) {
    
    hdsdp_retcode retcode = HDSDP_RETCODE_OK;
    
    /* Implement the self-dual embedding method. In the current form, HDSDP does not try to extract
       primal solution from it and it serves as a certificate tool of dual (in)feasibility.
     
       It can be either invoked after the infeasible method or independently invoked. HDSDP will use solver
       solution status to tell which is the case.
     */
    
    HSolver->whichMethod = HDSDP_ALGO_DUAL_HSD;
    
    int nMaxIter = get_int_param(HSolver, INT_PARAM_MAXITER);
    double dBarrierGamma = get_dbl_param(HSolver, DBL_PARAM_HSDGAMMA);
    double dAbsfeasTol = get_dbl_param(HSolver, DBL_PARAM_ABSFEASTOL);
    double dAbsoptTol = get_dbl_param(HSolver, DBL_PARAM_ABSOPTTOL);
    double dRelfeasTol = get_dbl_param(HSolver, DBL_PARAM_RELFEASTOL);
    double dReloptTol = get_dbl_param(HSolver, DBL_PARAM_RELOPTTOL);
    double dTimeLimit = get_dbl_param(HSolver, DBL_PARAM_TIMELIMIT);
    double dObjOneNorm = get_dbl_feature(HSolver, DBL_FEATURE_OBJONENORM);
    double dObjScal = get_dbl_feature(HSolver, DBL_FEATURE_OBJSCALING);
    double nSumCones = (double) get_int_feature(HSolver, INT_FEATURE_N_SUMCONEDIMS);
    
    /* Here we need to transform between two convergence criteria
       Recall that for absolute measure we have
       
            sqrt(sumCones) * |Rd| < absFeasTol * dObjScal
        =>  |Rd| < absFeasTol * dObjScal / sqrt(sumCones)
            Mu < absOptTol
     
       and for relative measure we have
                
            sqrt(sumCones) * |Rd| < relFeasTol * (1 + |C|) * dObjScal
        =>  |Rd| < relFeasTol * (1 + |C|) * dObjScal / sqrt(sumCones)
            Mu < relOptTol * (1 + |dObj|)
     */
    
    /* If we only need a dual feasible solution */
    if ( !dOnly ) {
        dAbsoptTol = 1e+20;
        dReloptTol = 1e+20;
    }
    
    double dFeasTol = HDSDP_MIN(dAbsfeasTol, dRelfeasTol * (1.0 + dObjOneNorm));
    dFeasTol = dFeasTol * dObjScal / sqrt(nSumCones);
    dAbsoptTol = dAbsoptTol * 1e-04;
    dReloptTol = dAbsoptTol * 1e-04;
    
    /* If we are starting from scratch */
    if ( HSolver->HStatus == HDSDP_UNKNOWN ) {
        HDSDP_SetStart(HSolver, HDSDP_ALGO_DUAL_HSD, dOnly);
    }
    
    HDSDP_PrintHeader(HSolver, HDSDP_ALGO_DUAL_HSD);
    
    /* Enter the main iteration */
    while ( 1 ) {

        /* Assert initial iteration is in the cone */
        int isInteriorPoint = 0;
        for ( int iCone = 0; iCone < HSolver->nCones; ++iCone ) {
            HDSDP_CALL(HConeCheckIsInterior(HSolver->HCones[iCone], HSolver->dBarHsdTau,
                                            HSolver->dRowDual, &isInteriorPoint));
            if ( !isInteriorPoint ) {
                break;
            }
        }
        
        if ( !isInteriorPoint ) {
            if ( HSolver->nIterCount == 0 ) {
                hdsdp_printf("Initial point is not in the cone. Adding slack value.\n");
                HSolver->dResidual *= 100.0;
                HDSDP_ResetStart(HSolver);
                HSolver->nIterCount += 1;
                continue;
            } else {
                hdsdp_printf("Iteration %d is not in the cone. \n", HSolver->nIterCount);
                retcode = HDSDP_RETCODE_FAILED;
                goto exit_cleanup;
                break;
            }
        }
        
        /* Afterwards, we set up the Schur complement system */
        HDSDP_CALL(HKKTBuildUp(HSolver->HKKT, KKT_TYPE_HOMOGENEOUS));
        HKKTRegularize(HSolver->HKKT, 1e-05);
        
        /* Then we export information needed */
        HKKTExport(HSolver->HKKT, HSolver->dMinvASinv, HSolver->dMinvASinvRdSinv,
                   HSolver->dMinvASinvCSinv, NULL, NULL, NULL, NULL);
        
        /* Factorize the KKT system */
        HDSDP_CALL(HKKTFactorize(HSolver->HKKT));
        
        /* Solve the KKT system */
        HDSDP_CALL(HKKTSolve(HSolver->HKKT, HSolver->rowRHS, HSolver->dMinvRowRHS));
        HDSDP_CALL(HKKTSolve(HSolver->HKKT, HSolver->dMinvASinv, NULL));
        HDSDP_CALL(HKKTSolve(HSolver->HKKT, HSolver->dMinvASinvRdSinv, NULL));
        HDSDP_CALL(HKKTSolve(HSolver->HKKT, HSolver->dMinvASinvCSinv, NULL));
        
        /* Assemble the iterations */
        HDSDP_HSD_BuildStep(HSolver);
        
        /* Do ratio test */
        HDSDP_CALL(HDSDP_HSD_RatioTest(HSolver, &HSolver->dDStep));
        
        /* Choose step size */
        if ( HSolver->dDStep > 1.0 ) {
            HSolver->dDStep = HDSDP_MIN(0.7 * HSolver->dDStep, 1.0);
        } else if ( HSolver->dDStep > 0.5 ) {
            HSolver->dDStep = HDSDP_MIN(0.5 * HSolver->dDStep, 1.0);
        } else if ( HSolver->dDStep > 0.2 ) {
            HSolver->dDStep = HDSDP_MIN(0.3 * HSolver->dDStep, 1.0);
        } else {
            HSolver->dDStep = HDSDP_MIN(0.2 * HSolver->dDStep, 1.0);
        }
        
        /* Print log */
        HDSDP_PrintLog(HSolver, HDSDP_ALGO_DUAL_HSD);
        
        /* Take step */
        HSolver->dBarHsdTau += HSolver->dDStep * HSolver->dBarHsdTauStep;
        for ( int iRow = 0; iRow < HSolver->nRows; ++iRow ) {
            HSolver->dRowDual[iRow] += HSolver->dDStep * HSolver->dRowDualStep[iRow];
        }
        HSolver->dResidual = HSolver->dResidual * (1.0 - HSolver->dDStep);
        HDSDP_SetResidual(HSolver, HSolver->dResidual);
        
        /* Reduce barrier parameter */
        double dBarrierMuTarget = 0.0;
        if ( HSolver->dBarrierMu > 1e-12 ) {
            if ( HSolver->dDStep > 0.8 && HSolver->dBarHsdTau > 1.0 ) {
                dBarrierMuTarget = 0.1 * HSolver->dBarrierMu;
                dBarrierMuTarget = HDSDP_MAX(dBarrierMuTarget, -0.1 * HSolver->dResidual / HSolver->dBarHsdTau);
                HSolver->dBarrierMu = HDSDP_MIN(HSolver->dBarrierMu, dBarrierMuTarget);
            } else {
                dBarrierMuTarget = dBarrierGamma * HSolver->dBarrierMu;
                dBarrierMuTarget = HDSDP_MAX(dBarrierMuTarget, -0.1 * HSolver->dResidual / HSolver->dBarHsdTau);
                HSolver->dBarrierMu = HDSDP_MIN(HSolver->dBarrierMu, dBarrierMuTarget);
            }
        } else {
            dBarrierMuTarget = 0.8 * HSolver->dBarrierMu;
            HSolver->dBarrierMu = HDSDP_MIN(HSolver->dBarrierMu, dBarrierMuTarget);
        }
        
        /* Convergence check */
        if ( fabs(HSolver->dResidual) < dFeasTol * HSolver->dBarHsdTau && HSolver->dBarrierMu < dAbsoptTol &&
            (HSolver->dBarrierMu < dReloptTol * (1 + 2.0 * fabs(HSolver->dObjVal)) &&
            fabs(HSolver->dObjImprove) < 1e-05 * (fabs(HSolver->dObjInternal) + 1.0))) {
            
            /* If we do not need a dual solution */
            if ( dOnly ) {
                HSolver->HStatus = HDSDP_DUAL_OPTIMAL;
            } else {
                HSolver->HStatus = HDSDP_DUAL_FEASIBLE;
            }
            
            break;
        }
        
        /* Infeasibility check */
        if ( HSolver->dBarHsdTau <= 1e-10 ) {
            HSolver->HStatus = HDSDP_SUSPECT_INFEAS_OR_UNBOUNDED;
            break;
        }
        
        if ( HUtilCheckCtrlC() ) {
            HSolver->HStatus = HDSDP_USER_INTERRUPT;
            break;
        }
        
        if ( HUtilGetTimeStamp() - HSolver->dTimeBegin >= dTimeLimit ) {
            HSolver->HStatus = HDSDP_TIMELIMIT;
            break;
        }
        
        HSolver->nIterCount += 1;
        if ( HSolver->nIterCount >= nMaxIter ) {
            HSolver->HStatus = HDSDP_MAXITER;
            break;
        }
    }
    
exit_cleanup:
    
    if ( retcode != HDSDP_RETCODE_OK ) {
        HSolver->HStatus = HDSDP_NUMERICAL;
    }
    
    return retcode;
}

static hdsdp_retcode HDSDP_PhaseB_BarDualPotentialSolve( hdsdp *HSolver ) {
    
    hdsdp_retcode retcode = HDSDP_RETCODE_OK;
    /* Implement the same potential reduction algorithm as in DSDP 5.8 */
    
    HSolver->whichMethod = HDSDP_ALGO_DUAL_POTENTIAL;
    
    
exit_cleanup:
    return retcode;
}

static hdsdp_retcode HDSDP_PhaseC_BarPrimalInfeasCheck( hdsdp *HSolver ) {
    
    hdsdp_retcode retcode = HDSDP_RETCODE_OK;
    /* Check existence of improving ray */
    
    
    
exit_cleanup:
    return retcode;
}

static hdsdp_retcode HDSDP_Conic_Solve( hdsdp *HSolver, int dOnly ) {
    
    hdsdp_retcode retcode = HDSDP_RETCODE_OK;
    
    HDSDP_CALL(HDSDP_PhaseA_BarInfeasSolve(HSolver, dOnly));
    
    if ( HSolver->HStatus == HDSDP_SUSPECT_INFEAS_OR_UNBOUNDED ) {
        HDSDP_CALL(HDSDP_PhaseA_BarHsdSolve(HSolver, dOnly));
    }
    
    if ( HSolver->HStatus == HDSDP_DUAL_FEASIBLE ) {
        HDSDP_CALL(HDSDP_PhaseB_BarDualPotentialSolve(HSolver));
    }
    
exit_cleanup:
    return retcode;
}

static hdsdp_retcode HDSDP_PureLP_Solve( hdsdp *HSolver ) {
    
    hdsdp_retcode retcode = HDSDP_RETCODE_OK;
    
    
    
    
exit_cleanup:
    return retcode;
}

extern hdsdp_retcode HDSDPSolve( hdsdp *HSolver, int dOnly ) {
    
    hdsdp_retcode retcode = HDSDP_RETCODE_OK;
    
    HDSDP_CALL(HDSDP_Conic_Solve(HSolver, dOnly));
    
exit_cleanup:
    return retcode;
}
