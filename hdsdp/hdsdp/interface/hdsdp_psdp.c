#include "hdsdp_psdp.h"

#ifdef HEADERPATH
#include "interface/hdsdp_psdp.h"
#include "interface/hdsdp_utils.h"
#include "interface/hdsdp_conic.h"
#include "interface/def_hdsdp.h"
#include "linalg/dense_opts.h"
#include "linalg/vec_opts.h"
#else
#include "hdsdp_psdp.h"
#include "hdsdp_utils.h"
#include "hdsdp_conic.h"
#include "def_hdsdp.h"
#include "dense_opts.h"
#include "vec_opts.h"
#endif

#include <math.h>

static void denseMatTripMultiplypy( int nCol, double *S, double *X, double *aux, double *XSX ) {
    
    /* Routine for multiplying three dense matrices X * S * X and adding it to buffer
       Check dataMatDenseKKT3ComputeSinvASinvImpl for more details */
    
    double *XCol = NULL;
    double *SXCol = NULL;
    double *SX = aux;
    
    HDSDP_ZERO(SX, double, nCol * nCol);
    for ( int iCol = 0; iCol < nCol; ++iCol ) {
        XCol = X + nCol * iCol;
        SXCol = SX + nCol * iCol;
        fds_symv(nCol, 1.0, S, XCol, 0.0, SXCol);
    }
    
    for ( int iCol = 0; iCol < nCol; ++iCol ) {
        SXCol = SX + nCol * iCol;
        for ( int iRow = 0; iRow < iCol; ++iRow ) {
            XCol = X + nCol * iRow;
            double dDotVal = dot(&nCol, SXCol, &HIntConstantOne, XCol, &HIntConstantOne);
            FULL_ENTRY(XSX, nCol, iCol, iRow) += dDotVal;
            FULL_ENTRY(XSX, nCol, iRow, iCol) += dDotVal;
        }
        
        XCol = X + nCol * iCol;
        FULL_ENTRY(XSX, nCol, iCol, iCol) += \
        dot(&nCol, SXCol, &HIntConstantOne, XCol, &HIntConstantOne);
    }

    return;
}

static void sdpPrimalConeILanczosMultiply( void *Hpsdp, double *dLhsVec, double *dRhsVec ) {
    
    /* Ratio test subroutine */
    hdsdp_psdp *psdp = (hdsdp_psdp *) Hpsdp;
    HFpLinsysBSolve(psdp->XFactor, 1, dLhsVec, dRhsVec);
    fds_symv(psdp->nCol, -1.0, psdp->dPrimalXStep, dRhsVec, 0.0, psdp->dPrimalAuxiVec1);
    HFpLinsysFSolve(psdp->XFactor, 1, psdp->dPrimalAuxiVec1, dRhsVec);
    
    return;
}

extern hdsdp_retcode HPSDPCreate( hdsdp_psdp **pHpsdp ) {
    
    hdsdp_retcode retcode = HDSDP_RETCODE_OK;
    HDSDP_NULLCHECK(pHpsdp);
    
    hdsdp_psdp *Hpsdp = NULL;
    HDSDP_INIT(Hpsdp, hdsdp_psdp, 1);
    HDSDP_MEMCHECK(Hpsdp);
    HDSDP_ZERO(Hpsdp, hdsdp_psdp, 1);
    
    *pHpsdp = Hpsdp;
    
exit_cleanup:
    return retcode;
}

extern hdsdp_retcode HPSDPInit( hdsdp_psdp *Hpsdp, hdsdp *HSolver ) {
    /* Initialize primal solver from dual solver */
    hdsdp_retcode retcode = HDSDP_RETCODE_OK;
    
    if ( HSolver->nCones > 1 || HSolver->dInfeas > 0.0 ||
        HSolver->HCones[0]->cone != HDSDP_CONETYPE_DENSE_SDP ) {
        retcode = HDSDP_RETCODE_FAILED;
        goto exit_cleanup;
    }
    
    Hpsdp->HSolver = HSolver;
    Hpsdp->HCone = HSolver->HCones[0];
    
    /* Space in dual solver is reused in PSDP */
    Hpsdp->nRow = HSolver->nRows;
    Hpsdp->nCol = HConeGetDim(Hpsdp->HCone);
    Hpsdp->nCones = 1;
    Hpsdp->rowRHS = HSolver->rowRHS;
    Hpsdp->HKKT = HSolver->HKKT;
    Hpsdp->dRowDual = HSolver->dRowDual;
    Hpsdp->dRowDualStep = HSolver->dRowDualStep;
    Hpsdp->dPrimalAuxiVec1 = HSolver->dHAuxiVec1;
    Hpsdp->dPrimalAuxiVec2 = HSolver->dHAuxiVec2;
    
    Hpsdp->dBarrierMu = HSolver->dBarrierMu;
    
    HDSDP_CALL(HFpLinsysCreate(&Hpsdp->XFactor, Hpsdp->nCol, HDSDP_LINSYS_DENSE_DIRECT));
    HDSDP_CALL(HFpLinsysSymbolic(Hpsdp->XFactor, NULL, NULL));
    
    /* Use new Lanczos solver for primal ratio test */
    HDSDP_CALL(HLanczosCreate(&Hpsdp->Lanczos));
    HDSDP_CALL(HLanczosInit(Hpsdp->Lanczos, Hpsdp->nCol, 30));
    HLanczosSetData(Hpsdp->Lanczos, Hpsdp, sdpPrimalConeILanczosMultiply);
    
    /* Primal memory allocation */
    HDSDP_INIT(Hpsdp->dPrimalKKTRhs, double, Hpsdp->nRow);
    HDSDP_MEMCHECK(Hpsdp->dPrimalKKTRhs);
    HDSDP_INIT(Hpsdp->dPrimalX, double, Hpsdp->nCol * Hpsdp->nCol);
    HDSDP_MEMCHECK(Hpsdp->dPrimalX);
    HDSDP_INIT(Hpsdp->dPrimalScalX, double, Hpsdp->nCol * Hpsdp->nCol);
    HDSDP_MEMCHECK(Hpsdp->dPrimalScalX);
    HDSDP_INIT(Hpsdp->dDualS, double, Hpsdp->nCol * Hpsdp->nCol);
    HDSDP_MEMCHECK(Hpsdp->dDualS);
    HDSDP_INIT(Hpsdp->dPrimalXStep, double, Hpsdp->nCol * Hpsdp->nCol);
    HDSDP_MEMCHECK(Hpsdp->dPrimalXStep);
    HDSDP_INIT(Hpsdp->dPrimalMatBuffer, double, Hpsdp->nCol * Hpsdp->nCol);
    HDSDP_MEMCHECK(Hpsdp->dPrimalMatBuffer);
    
    /* Register primal matrix in KKT solver */
    HDSDPGetConeValues(Hpsdp->HSolver, 0, Hpsdp->dPrimalX, NULL, Hpsdp->dPrimalMatBuffer);
    HDSDP_MEMCPY(Hpsdp->dPrimalScalX, Hpsdp->dPrimalX, double, Hpsdp->nCol * Hpsdp->nCol);
    Hpsdp->HKKT->dPrimalX = &Hpsdp->dPrimalX;
    Hpsdp->HKKT->iKKTCone = 0;
    
    hdsdp_printf("Primal solver is initialized. \n");
    
exit_cleanup:
    return retcode;
}


extern hdsdp_retcode HPSDPOptimize( hdsdp_psdp *Hpsdp ) {
    
    hdsdp_retcode retcode = HDSDP_RETCODE_OK;
    
    
    double *dPrimalX = Hpsdp->dPrimalX;
    double *dPrimalScalX = Hpsdp->dPrimalScalX;
    double *dDualS = Hpsdp->dDualS;
    double *dPrimalXStep = Hpsdp->dPrimalXStep;
    double *dPrimalMatBuffer = Hpsdp->dPrimalMatBuffer;
    double *dPrimalKKTRhs = Hpsdp->dPrimalKKTRhs;
    
    double *dRowDual = Hpsdp->dRowDual;
    double *dRowDualStep = Hpsdp->dRowDualStep;
    double *dPrimalInfeas = Hpsdp->dPrimalAuxiVec1;
    double dBarrierMu = Hpsdp->dBarrierMu;
    
    hdsdp_cone *cone = Hpsdp->HCone;
    hdsdp *HSolver  = Hpsdp->HSolver;
    
    double dAbsfeasTol = get_dbl_param(HSolver, DBL_PARAM_ABSFEASTOL);
    double dAbsoptTol = get_dbl_param(HSolver, DBL_PARAM_ABSOPTTOL);
    double dRelfeasTol = get_dbl_param(HSolver, DBL_PARAM_RELFEASTOL);
    double dReloptTol = get_dbl_param(HSolver, DBL_PARAM_RELOPTTOL);
    double dTimeLimit = get_dbl_param(HSolver, DBL_PARAM_TIMELIMIT);
    double dObjOneNorm = get_dbl_feature(HSolver, DBL_FEATURE_OBJONENORM);
    double dObjScal = get_dbl_feature(HSolver, DBL_FEATURE_OBJSCALING);
    double dRhsScal = get_dbl_feature(HSolver, DBL_FEATURE_RHSSCALING);
    double pdScal = dObjScal * dRhsScal;
    
    /* Hack into data structure */
    hdsdp_cone_sdp_dense *coneData = (hdsdp_cone_sdp_dense *) Hpsdp->HCone->coneData;
    double *dDualSStep = coneData->dualStep;
    
    int nIter = 0;
    
    /* Build KKT system */
    HDSDP_CALL(HKKTBuildUp(Hpsdp->HKKT, KKT_TYPE_PRIMAL));
    hdsdp_printf("KKT Build up: %f s", HUtilGetTimeStamp() - HSolver->dTimeBegin);
    HDSDP_CALL(HKKTFactorize(Hpsdp->HKKT));
    hdsdp_printf("KKT Factorize: %f s", HUtilGetTimeStamp() - HSolver->dTimeBegin);
    
    /* Enter primal main loop */
    for ( nIter = 0; nIter < 5; ++nIter ) {
        
        /* Compute primal infeasibility rp = A * X - b */
        HDSDP_ZERO(dPrimalInfeas, double, Hpsdp->nRow);
        HConeComputeATimesX(cone, dPrimalX, dPrimalInfeas);
        
        double dPrimalInfeasNorm = 0.0;
        
        for ( int iRow = 0; iRow < Hpsdp->nRow; ++iRow ) {
            double dTmp = Hpsdp->rowRHS[iRow] - dPrimalInfeas[iRow];
            dPrimalInfeasNorm += dTmp * dTmp;
        }
        dPrimalInfeasNorm = sqrt(dPrimalInfeasNorm);
        
        /* Get current dual solution from cone. Dual solution is always feasible */
        HDSDPGetConeValues(Hpsdp->HSolver, 0, NULL, Hpsdp->dDualS, NULL);
                
        /* Prepare primal KKT RHS */
        /* Get buffer matrix X * S * X */
        HDSDP_ZERO(dPrimalMatBuffer, double, Hpsdp->nCol * Hpsdp->nCol);
        denseMatTripMultiplypy(Hpsdp->nCol, dDualS, dPrimalX, dPrimalXStep, dPrimalMatBuffer);
        HDSDP_ZERO(dPrimalKKTRhs, double, Hpsdp->nRow);
        HConeComputeATimesX(cone, dPrimalMatBuffer, dPrimalKKTRhs);
        
        for ( int iRow = 0; iRow < Hpsdp->nRow; ++iRow ) {
            dPrimalKKTRhs[iRow] -= dBarrierMu * dPrimalInfeas[iRow];
        }
        
        /* Solve the KKT system */
        HDSDP_CALL(HKKTSolve(Hpsdp->HKKT, dPrimalKKTRhs, dRowDualStep));
        
        /* Get Primal direction */
        /* Compute X * S * X + XScal * dS * XScal */
        /* Dual ratio test: necessary for getting dS from buffer */
        double dMaxDualStep = 0.0;
        HDSDP_CALL(HConeRatioTest(Hpsdp->HCone, 0.0, dRowDualStep, 1.0, BUFFER_DUALVAR, &dMaxDualStep));
        HUtilMatSymmetrize(Hpsdp->nCol, dDualSStep);
        denseMatTripMultiplypy(Hpsdp->nCol, dDualSStep, dPrimalScalX, dPrimalXStep, dPrimalMatBuffer);
        
        /* Compute X - 1 / mu * (X * S * X + XScal * dS * XScal) */
        HDSDP_MEMCPY(dPrimalXStep, dPrimalX, double, Hpsdp->nCol * Hpsdp->nCol);
        int nElem = Hpsdp->nCol * Hpsdp->nCol;
        double oneOverMu = - 1.0 / dBarrierMu;
        axpy(&nElem, &oneOverMu, dPrimalMatBuffer, &HIntConstantOne, dPrimalXStep, &HIntConstantOne);
        
        /* Factorize the primal matrix */
        if ( nIter == 0 ) {
            HDSDP_CALL(HFpLinsysNumeric(Hpsdp->XFactor, NULL, NULL, dPrimalX));
        }
        
        /* Ratio test */
        double dMaxPrimalStep = 0.0;
        HDSDP_CALL(HLanczosSolve(Hpsdp->Lanczos, NULL, &dMaxPrimalStep));
        
        dMaxPrimalStep = HDSDP_MIN(0.95 * dMaxPrimalStep, 1.0);
        dMaxDualStep = HDSDP_MIN(0.95 * dMaxDualStep, 1.0);
        
        /* Take step */
        /* Primal step */
        axpy(&nElem, &dMaxPrimalStep, dPrimalXStep, &HIntConstantOne, dPrimalX, &HIntConstantOne);
        
        /* Dual step */
        axpy(&Hpsdp->nRow, &dMaxDualStep, dRowDualStep, &HIntConstantOne, dRowDual, &HIntConstantOne);
        HConeUpdate(cone, 1.0, dRowDual);
        
        int isInterior = 0;
        HDSDP_CALL(HConeCheckIsInterior(cone, 1.0, dRowDual, &isInterior));
        
        if ( !isInterior ) {
            retcode = HDSDP_RETCODE_FAILED;
            goto exit_cleanup;
        }
        
        HFpLinsysPsdCheck(Hpsdp->XFactor, NULL, NULL, dPrimalX, &isInterior);
        
        if ( !isInterior ) {
            retcode = HDSDP_RETCODE_FAILED;
            goto exit_cleanup;
        }
        
        /* Update barrier parameter */
        double dDualObj = dot(&Hpsdp->nRow, Hpsdp->rowRHS, &HIntConstantOne, dRowDual, &HIntConstantOne);
        double dPrimalObj = HConeComputeTraceCX(cone, dPrimalX);
        double dCompl = dot(&nElem, dPrimalX, &HIntConstantOne, dDualS, &HIntConstantOne);
        dCompl = dCompl / Hpsdp->nCol;
        
        dBarrierMu = HDSDP_MIN(dBarrierMu, dCompl / (5.0 * Hpsdp->nCol));
        
        /* Synchronize data to HSolver */
        HSolver->pObjVal = dPrimalObj;
        HSolver->dObjVal = dDualObj;
        HSolver->pInfeas = dPrimalInfeasNorm;
        HSolver->dBarrierMu = dBarrierMu;
        HSolver->dDStep = dMaxDualStep;
        
        double elapsedTime = HUtilGetTimeStamp() - HSolver->dTimeBegin;
        Hpsdp->HSolver->nIterCount += 1;
        
        hdsdp_printf("    %5d  %+15.8e  %+15.8e  %8.2e  %8.2e  %5.2f %5.2f     %4.1f \n", Hpsdp->HSolver->nIterCount + 1,
               HSolver->pObjVal, HSolver->dObjVal, HSolver->pInfeas, HSolver->dBarrierMu,
               dMaxPrimalStep, HSolver->dDStep, elapsedTime);
        
        if ( HSolver->comp < (fabs(HSolver->pObjVal) + fabs(HSolver->dObjVal) + 1.0) * dReloptTol &&
            HSolver->comp < dAbsoptTol / pdScal ) {
            break;
        }
    }
 
exit_cleanup:
    return retcode;
}

extern void HPSDPClear( hdsdp_psdp *Hpsdp ) {
    
    if ( !Hpsdp ) {
        return;
    }
    
    HDSDP_FREE(Hpsdp->dPrimalKKTRhs);
    HDSDP_FREE(Hpsdp->dPrimalX);
    HDSDP_FREE(Hpsdp->dPrimalScalX);
    HDSDP_FREE(Hpsdp->dDualS);
    HDSDP_FREE(Hpsdp->dPrimalXStep);
    HDSDP_FREE(Hpsdp->dPrimalMatBuffer);
    
    HLanczosDestroy(&Hpsdp->Lanczos);
    HFpLinsysDestroy(&Hpsdp->XFactor);
    
    return;
}

extern void HPSDPDestroy( hdsdp_psdp **pHpsdp ) {
    
    if ( !pHpsdp ) {
        return;
    }
    
    HPSDPClear(*pHpsdp);
    HDSDP_FREE(*pHpsdp);
    
    return;
}
