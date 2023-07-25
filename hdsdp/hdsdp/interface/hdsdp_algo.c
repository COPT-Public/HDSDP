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

#define HDSDP_ALGO_DUAL_INFEAS     (0)
#define HDSDP_ALGO_DUAL_HSD        (1)
#define HDSDP_ALGO_DUAL_POTENTIAL  (2)
static void HDSDP_SetStart( hdsdp *HSolver, int SDPMethod, int dOnly ) {
    
    /* Set staring point for the dual iterations. */
    HDSDP_ZERO(HSolver->dRowDual, double, HSolver->nRows);
    HSolver->dBarHsdTau = 1.0;
    
    if ( SDPMethod == HDSDP_ALGO_DUAL_HSD ) {
        if ( dOnly ) {
            HSolver->dBarrierMu = 1.0;
        } else {
            HSolver->dBarrierMu = 1e+08;
        }
    } else {
        HSolver->dBarrierMu = get_dbl_param(HSolver, DBL_PARAM_BARMUSTART);
    }
    
    double objFroNorm = get_dbl_feature(HSolver, DBL_FEATURE_OBJFRONORM);
    
    if ( dOnly ) {
        HSolver->dResidual = - objFroNorm * get_dbl_param(HSolver, DBL_PARAM_DUALSTART);
    } else {
        HSolver->dResidual = - objFroNorm * 1e+01;
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
    HSolver->dResidual = - objFroNorm * 1e+05;
    
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
            hdsdp_printf("    %5s  %12s  %12s  %8s  %8s  %5s  %5s \n",
                         "nIter", "pObj", "dObj", "dInf", "Mu", "Step", "T [H]");
            break;
        case HDSDP_ALGO_DUAL_INFEAS:
            hdsdp_printf("HDSDP starts. Using infeasible dual method \n\n");
            hdsdp_printf("    %5s  %12s  %12s  %8s  %8s  %5s  %5s \n",
                         "nIter", "pObj", "dObj", "dInf", "Mu", "Step", "T [D]");
            break;
        case HDSDP_ALGO_DUAL_POTENTIAL:
            hdsdp_printf("HDSDP starts. Using feasible dual method \n\n");
            hdsdp_printf("    %5s  %12s  %12s  %8s  %8s  %5s  %5s \n",
                         "nIter", "pObj", "dObj", "pInf", "Mu", "Step", "T [P]");
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
    double nSumCones = (double) get_dbl_feature(HSolver, INT_FEATURE_N_SUMCONEDIMS);
    double pdObjScal = dRhsScal * dObjScal / HSolver->dBarHsdTau;
    HSolver->dInfeas = dRhsScal * sqrt(nSumCones) * fabs(HSolver->dResidual) / HSolver->dBarHsdTau;
    
    HSolver->dObjVal = HSolver->dObjInternal * pdObjScal;
    HSolver->pObjVal = HSolver->pObjInternal * pdObjScal;
    
    switch (SDPMethod) {
        case HDSDP_ALGO_DUAL_HSD:
            hdsdp_printf("    %5d  %+12.5e  %+12.5e  %8.2e  %8.2e  %5.2f  %4.1f  %8.2e\n", HSolver->nIterCount + 1,
                   HDSDP_INFINITY, HSolver->dObjVal, HSolver->dInfeas,
                   HSolver->dBarrierMu, HSolver->dDStep, elapsedTime, HSolver->dBarHsdTau);
            break;
        case HDSDP_ALGO_DUAL_INFEAS:
            hdsdp_printf("    %5d  %+12.5e  %+12.5e  %8.2e  %8.2e  %5.2f  %4.1f \n", HSolver->nIterCount + 1,
                   HSolver->pObjVal * pdObjScal, HSolver->dObjVal * pdObjScal, HSolver->dInfeas, HSolver->dBarrierMu,
                   HSolver->dDStep, elapsedTime);
            break;
        case HDSDP_ALGO_DUAL_POTENTIAL:
            hdsdp_printf("    %5d  %+12.5e  %+12.5e  %8.2e  %8.2e  %5.2f  %4.1f \n", HSolver->nIterCount + 1,
                   HSolver->pObjVal * pdObjScal, HSolver->dObjVal * pdObjScal, HSolver->pInfeas, HSolver->dBarrierMu,
                   HSolver->dDStep, elapsedTime);
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
#define bTy  HSolver->dObjInternal
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

static hdsdp_retcode HDSDP_RatioTest( hdsdp *HSolver, int SDPMethod, double *dMaxDist ) {
    
    hdsdp_retcode retcode = HDSDP_RETCODE_OK;
    
    double dMaxStep = HDSDP_INFINITY;
    double dStep = 0.0;
    
    /* Get step of tau  */
    if ( SDPMethod == HDSDP_ALGO_DUAL_HSD && HSolver->dBarHsdTauStep ) {
        dStep = HSolver->dBarHsdTau / HSolver->dBarHsdTauStep;
        if ( dStep < 0.0 ) {
            dMaxStep = HDSDP_MIN(dMaxStep, -dStep);
        }
    }
    
    if ( dMaxStep < 1e-03 ) {
        retcode = HDSDP_RETCODE_FAILED;
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

static hdsdp_retcode HDSDP_PhaseA_BarInfeasSolve( hdsdp *HSolver, int dOnly ) {
    
    hdsdp_retcode retcode = HDSDP_RETCODE_OK;
    /* Implement the infeasible-start dual potential reduction method.
       An adaptive stepsize strategy is applied to find eliminate dual infeasibility efficiently.
       Also the output from the KKT solver will be reused as if in the corrector method of DSDP.
     
       Primal solution, if available. Will be recorded by this phase.
       If it is hard to eliminate dual infeasibility, the embedding will be invoked to provide a certificate
     */
    
    
    
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
    
    int nMaxIter = get_int_param(HSolver, INT_PARAM_MAXITER);
    double dBarrierGamma = get_dbl_param(HSolver, DBL_PARAM_HSDGAMMA);
    double dAbsfeasTol = get_dbl_param(HSolver, DBL_PARAM_ABSFEASTOL);
    double dAbsoptTol = get_dbl_param(HSolver, DBL_PARAM_ABSOPTTOL);
    double dRelfeasTol = get_dbl_param(HSolver, DBL_PARAM_RELFEASTOL);
    double dReloptTol = get_dbl_param(HSolver, DBL_PARAM_RELOPTTOL);
    double dTimeLimit = get_dbl_param(HSolver, DBL_PARAM_TIMELIMIT);
    double dObjOneNorm = get_dbl_feature(HSolver, DBL_FEATURE_OBJONENORM);
    double dObjScal = get_dbl_feature(HSolver, DBL_FEATURE_OBJSCALING);
    double nSumCones = (double) get_dbl_feature(HSolver, INT_FEATURE_N_SUMCONEDIMS);
    
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
    
    /* If we only solve the dual */
    if ( !dOnly ) {
        dAbsoptTol = 1e+20;
        dReloptTol = 1e+20;
    }
    
    double dFeasTol = HDSDP_MIN(dAbsfeasTol, dRelfeasTol * (1 + dObjOneNorm));
    dFeasTol = dFeasTol * sqrt(nSumCones) * dObjScal;
    
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
        HKKTExport(HSolver->HKKT, HSolver->dMinvASinv, HSolver->dMinvASinvRdSinv, HSolver->dMinvASinvCSinv,
                   NULL, NULL, NULL);
        
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
        for ( int iCone = 0; iCone < HSolver->nCones; ++iCone ) {
            HConeReduceResi(HSolver->HCones[iCone], 1.0 - HSolver->dDStep);
        }
        
        /* Reduce barrier parameter */
        double dBarrierMuTarget = 0.0;
        if ( HSolver->dBarrierMu > 1e-12 ) {
            if ( HSolver->dDStep > 0.6 && HSolver->dBarHsdTau > 1.0 ) {
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
            HSolver->dBarrierMu < dReloptTol * (1 + 2.0 * fabs(HSolver->dObjVal)) ) {
            
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
    
    
    
    
exit_cleanup:
    return retcode;
}

static hdsdp_retcode HDSDP_PhaseC_BarPrimalInfeasCheck( hdsdp *HSolver ) {
    
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
    
    HDSDP_CALL(HDSDP_PhaseA_BarHsdSolve(HSolver, dOnly));
    
exit_cleanup:
    return retcode;
}
