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
 
   The four algorithms are respectively
 
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

#define HDSDP_ALGO_DUAL_INFEAS     (0)
#define HDSDP_ALGO_DUAL_HSD        (1)
#define HDSDP_ALGO_DUAL_POTENTIAL  (2)
static void HDSDP_SetStart( hdsdp *HSolver ) {
    
    
    return;
}

static void HDSDP_ResetStart( hdsdp *HSolver ) {
    
    return;
}

static void HDSDP_PrintHeader( hdsdp *HSolver, int SDPMethod ) {
    
    switch (SDPMethod) {
        case HDSDP_ALGO_DUAL_HSD:
            printf("%5s  %8s  %8s  %6s  %6s  %6s %10s \n", "nIter", "pObj", "dObj", "dInf", "Tau", "Step", "Time [H]");
            break;
        case HDSDP_ALGO_DUAL_INFEAS:
            printf("%5s  %8s  %8s  %6s  %6s  %10s \n", "nIter", "pObj", "dObj", "dInf", "Step", "Time [D]");
            break;
        case HDSDP_ALGO_DUAL_POTENTIAL:
            printf("%5s  %8s  %8s  %6s  %6s  %10s \n", "nIter", "pObj", "dObj", "pInf", "Step", "Time [P]");
            break;
        default:
            assert( 0 );
            break;
    }
    
    return;
}

static void HDSDP_PrintLog( hdsdp *HSolver, int SDPMethod ) {
    
    double elapsedTime = HUtilGetTimeStamp() - HSolver->dTimeBegin;
    
    switch (SDPMethod) {
        case HDSDP_ALGO_DUAL_HSD:
            printf("%5d  %8.3e  %8.3e  %6.2e  %6.2e  %6.2e  %10.1f\n", HSolver->nIterCount,
                   HDSDP_INFINITY, HSolver->dObjVal, HSolver->dInfeas,
                   HSolver->dBarHsdTau, HSolver->dDStep, elapsedTime);
            break;
        case HDSDP_ALGO_DUAL_INFEAS:
            printf("%5d  %8.3e  %8.3e  %6.2e  %6.2e  %10.1f \n", HSolver->nIterCount,
                   HSolver->pObjVal, HSolver->dObjVal, HSolver->dInfeas, HSolver->dDStep, elapsedTime);
            break;
        case HDSDP_ALGO_DUAL_POTENTIAL:
            printf("%5d  %8.3e  %8.3e  %6.2e  %6.2e  %10.1f \n", HSolver->nIterCount,
                   HSolver->pObjVal, HSolver->dObjVal, HSolver->pInfeas, HSolver->dDStep, elapsedTime);
            break;
        default:
            assert( 0 );
            break;
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
#define bTy  HSolver->dObjVal
#define mu   HSolver->dBarrierMu
#define tau  HSolver->dBarHsdTau
#define dtau HSolver->dBarHsdTauStep
static void HDSDP_BuildHsdStep( hdsdp *HSolver ) {
    
    /* Prepare intermediate results */
    double dHsdTraceCSinvCSinv = 0.0;
    double dHsdTraceCSinv = 0.0;
    double dHsdTraceCSinvRdCSinv = 0.0;
    
    HSolver->dBarHsdTauStep = 0.0;
    
    /* First get dual objective */
    bTy = 0.0;
    
    for ( int iRow = 0; iRow < HSolver->nRows; ++iRow ) {
        bTy += b[iRow] * y[iRow];
    }
    
    /* Then compute dd1 = b - mu * ASinvCSinv */
    double *dASinvCSinvVec = HSolver->HKKT->dASinvCSinvVec;
    
    for ( int iRow = 0; iRow < HSolver->nRows; ++iRow ) {
        dd1[iRow] = b[iRow] - mu * dASinvCSinvVec[iRow];
    }
    
    /* Compute dTau */
    HKKTExport(HSolver->HKKT, NULL, NULL, NULL,
               &dHsdTraceCSinvCSinv, &dHsdTraceCSinv, &dHsdTraceCSinvRdCSinv);
    
    double dTauStepEnumerator = - bTy + mu / tau + mu * (dHsdTraceCSinv - dHsdTraceCSinvCSinv);
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

static hdsdp_retcode HDSDP_RatioTest( hdsdp *HSolver, int SDPMethod, double *dMaxDist ) {
    
    hdsdp_retcode retcode = HDSDP_RETCODE_OK;
    
    double dMaxStep = HDSDP_INFINITY;
    double dStep = 0.0;
    
    /* Get step of tau  */
    if ( SDPMethod == HDSDP_ALGO_DUAL_HSD && HSolver->dBarHsdTauStep ) {
        dStep = HSolver->dBarHsdTau / HSolver->dBarHsdTauStep;
        if ( dStep < 0.0 ) {
            dMaxStep = HDSDP_MIN(dMaxStep, dStep);
        }
    }
    
    if ( dMaxStep < 1e-03 ) {
        goto exit_cleanup;
    }
    
    /* Aggregate step sizes from different cones */
    for ( int iCone = 0; iCone < HSolver->nCones; ++iCone ) {
        /* We eleminate dual infeasibility with full Newton's step */
        HDSDP_CALL(HConeRatioTest(HSolver->HCones[iCone], HSolver->dBarHsdTauStep,
                                  HSolver->dRowDualStep, 1.0, &dStep));
        dMaxStep = HDSDP_MIN(dMaxStep, dStep);
    }
    
exit_cleanup:
    *dMaxDist = dMaxStep;
    return retcode;
}

static hdsdp_retcode HDSDP_PhaseC_BarPrimalInfeasCheck( hdsdp *HSolver ) {
    
    hdsdp_retcode retcode = HDSDP_RETCODE_OK;
    
    
    
    
exit_cleanup:
    return retcode;
}

static hdsdp_retcode HDSDP_PhaseA_BarInfeasSolve( hdsdp *HSolver, int dOnly ) {
    
    hdsdp_retcode retcode = HDSDP_RETCODE_OK;
    
    
    
exit_cleanup:
    return retcode;
}

static hdsdp_retcode HDSDP_PhaseA_BarHsdSolve( hdsdp *HSolver, int dOnly ) {
    
    hdsdp_retcode retcode = HDSDP_RETCODE_OK;
    
    (void) dOnly;
    /* Implement the self-dual embedding method. In the current form, HDSDP does not try to extract
       primal solution from it and it serves as a certificate tool of dual (in)feasibility.
     
       It can be either invoked after the infeasible method or independently invoked. HDSDP will use solver
       solution status to tell which is the case.
     */
    
    int nIterLeft = get_int_param(HSolver, INT_PARAM_MAXITER) - HSolver->nIterCount;
    double dBarrierGamma = get_dbl_param(HSolver, DBL_PARAM_HSDGAMMA);
    double dAbsfeasTol = get_dbl_param(HSolver, DBL_PARAM_ABSFEASTOL);
    double dAbsoptTol = get_dbl_param(HSolver, DBL_PARAM_ABSOPTTOL);
    double dRelfeasTol = get_dbl_param(HSolver, DBL_PARAM_RELFEASTOL);
    double dReloptTol = get_dbl_param(HSolver, DBL_PARAM_RELOPTTOL);
    
    /* If we are starting from scratch */
    if ( HSolver->HStatus == HDSDP_UNKNOWN ) {
        HDSDP_SetStart(HSolver);
    }
    
    HDSDP_PrintHeader(HSolver, HDSDP_ALGO_DUAL_HSD);
    /* Enter the main iteration */
    for ( int iBarHsdIter = 0; iBarHsdIter < nIterLeft; ++iBarHsdIter ) {

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
            if ( iBarHsdIter == 0 ) {
                printf("Initial point is not in the cone. Adding slack value.\n");
                HSolver->dResidual *= 100.0;
                iBarHsdIter -= 1;
                HDSDP_ResetStart(HSolver);
                continue;
            } else {
                printf("Iteration %d is not in the cone. \n", iBarHsdIter);
                retcode = HDSDP_RETCODE_FAILED;
                goto exit_cleanup;
                break;
            }
        }
        
        /* Afterwards, we set up the Schur complement system */
        HDSDP_CALL(HKKTBuildUp(HSolver->HKKT, KKT_TYPE_HOMOGENEOUS));
        
        /* Then we export information needed */
        HKKTExport(HSolver->HKKT, HSolver->dMinvASinv, HSolver->dMinvASinvRdSinv, HSolver->dMinvASinvCSinv, NULL, NULL, NULL);
        
        /* Factorize the KKT system */
        HDSDP_CALL(HKKTFactorize(HSolver->HKKT));
        
        /* Solve the KKT system */
        HDSDP_CALL(HKKTSolve(HSolver->HKKT, HSolver->rowRHS, HSolver->dMinvRowRHS));
        HDSDP_CALL(HKKTSolve(HSolver->HKKT, HSolver->dMinvASinv, NULL));
        HDSDP_CALL(HKKTSolve(HSolver->HKKT, HSolver->dMinvASinvRdSinv, NULL));
        HDSDP_CALL(HKKTSolve(HSolver->HKKT, HSolver->dMinvASinvCSinv, NULL));
        
        /* Assemble the iterations */
        HDSDP_BuildHsdStep(HSolver);
        
        /* Do ratio test */
        HDSDP_CALL(HDSDP_RatioTest(HSolver, HDSDP_ALGO_DUAL_HSD, &HSolver->dDStep));
        
        /* Choose step */
        HSolver->dDStep = HDSDP_MIN(0.9 * HSolver->dDStep, 1.0);
        
        /* Take step */
        HSolver->dBarHsdTau = HSolver->dDStep * HSolver->dBarHsdTauStep;
        for ( int iRow = 0; iRow < HSolver->nRows; ++iRow ) {
            HSolver->dRowDual[iRow] += HSolver->dDStep * HSolver->dRowDualStep[iRow];
        }
        
        HSolver->dResidual = HSolver->dResidual * (1.0 - HSolver->dDStep);
        
        for ( int iCone = 0; iCone < HSolver->nCones; ++iCone ) {
            HConeReduceResi(HSolver->HCones[iCone], 1.0 - HSolver->dDStep);
        }
        
        /* Reduce barrier parameter */
        double dBarrierMuTarget =  dBarrierGamma * HSolver->dBarrierMu;
        
        if ( HSolver->dResidual > dAbsoptTol )
        dBarrierGamma = HDSDP_MAX(dBarrierMuTarget, fabs(HSolver->dResidual / HSolver->dBarHsdTau));
        HSolver->dBarrierMu = HDSDP_MIN(HSolver->dBarrierMu, dBarrierGamma);
        
        /* Print log */
        HDSDP_PrintLog(HSolver, HDSDP_ALGO_DUAL_HSD);
    }
    
exit_cleanup:
    return retcode;
}

static hdsdp_retcode HDSDP_PhaseB_BarDualPotentialSolve( hdsdp *HSolver ) {
    
    hdsdp_retcode retcode = HDSDP_RETCODE_OK;
    
    
    
    
exit_cleanup:
    return retcode;
}

static hdsdp_retcode HDSDP_Conic_Solve( hdsdp *HSolver ) {
    
    hdsdp_retcode retcode = HDSDP_RETCODE_OK;
    
    
    
    
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
    
    HDSDP_CALL(HDSDP_PhaseA_BarHsdSolve(HSolver, 0));
    
exit_cleanup:
    return retcode;
}
