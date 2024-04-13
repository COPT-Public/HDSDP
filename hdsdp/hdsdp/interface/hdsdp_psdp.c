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

static void sdpPrimalConeILanczosMultiply( void *Hpsdp, double *dLhsVec, double *dRhsVec ) {
    
    /* Ratio test subroutine */
    hdsdp_psdp *psdp = (hdsdp_psdp *) Hpsdp;
    int iLanczos = psdp->iLanczos;
    HFpLinsysBSolve(psdp->XFactors[iLanczos], 1, dLhsVec, dRhsVec);
    fds_symv(HConeGetDim(psdp->HCones[iLanczos]), -1.0, psdp->dPrimalXStep[iLanczos], dRhsVec, 0.0, psdp->dPrimalAuxiVec1);
    HFpLinsysFSolve(psdp->XFactors[iLanczos], 1, psdp->dPrimalAuxiVec1, dRhsVec);
    
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
    
    /* We need initial dual feasible solution */
    int nLpCones = get_int_feature(HSolver, INT_FEATURE_N_LPCONES);
    
    if ( HSolver->dInfeas > 0.0 || nLpCones > 0 ) {
        retcode = HDSDP_RETCODE_FAILED;
        goto exit_cleanup;
    }
    
    Hpsdp->HSolver = HSolver;
    Hpsdp->HCones = HSolver->HCones;
    
    /* Space in dual solver is reused in PSDP */
    Hpsdp->nRow = HSolver->nRows;
    Hpsdp->nCones = HSolver->nCones;
    Hpsdp->rowRHS = HSolver->rowRHS;
    Hpsdp->HKKT = HSolver->HKKT;
    Hpsdp->dRowDual = HSolver->dRowDual;
    Hpsdp->dRowDualStep = HSolver->dRowDualStep;
    Hpsdp->dPrimalAuxiVec1 = HSolver->dHAuxiVec1;
    Hpsdp->dPrimalAuxiVec2 = HSolver->dHAuxiVec2;
    
    Hpsdp->dBarrierMu = HSolver->dBarrierMu;
    
    HDSDP_INIT(Hpsdp->XFactors, hdsdp_linsys *, Hpsdp->nCones);
    HDSDP_MEMCHECK(Hpsdp->XFactors);
    
    /* Use new Lanczos solver for primal ratio test */
    HDSDP_INIT(Hpsdp->Lanczos, hdsdp_lanczos *, Hpsdp->nCones);
    HDSDP_MEMCHECK(Hpsdp->Lanczos);
    
    for ( int iCone = 0; iCone < Hpsdp->nCones; ++iCone ) {
        int nCol = HConeGetDim(Hpsdp->HCones[iCone]);
        HDSDP_CALL(HFpLinsysCreate(&Hpsdp->XFactors[iCone], nCol, HDSDP_LINSYS_DENSE_DIRECT));
        HDSDP_CALL(HFpLinsysSymbolic(Hpsdp->XFactors[iCone], NULL, NULL));
        HDSDP_CALL(HLanczosCreate(&Hpsdp->Lanczos[iCone]));
        HDSDP_CALL(HLanczosInit(Hpsdp->Lanczos[iCone], nCol, 30));
        HLanczosSetData(Hpsdp->Lanczos[iCone], Hpsdp, sdpPrimalConeILanczosMultiply);
    }
    
    /* Primal memory allocation */
    HDSDP_INIT(Hpsdp->dPrimalKKTRhs, double, Hpsdp->nRow);
    HDSDP_MEMCHECK(Hpsdp->dPrimalKKTRhs);
    
    HDSDP_INIT(Hpsdp->dPrimalX, double *, Hpsdp->nCones);
    HDSDP_MEMCHECK(Hpsdp->dPrimalX);
    
    HDSDP_INIT(Hpsdp->dPrimalScalX, double *, Hpsdp->nCones);
    HDSDP_MEMCHECK(Hpsdp->dPrimalScalX);
    
    HDSDP_INIT(Hpsdp->dPrimalXStep, double *, Hpsdp->nCones);
    HDSDP_MEMCHECK(Hpsdp->dPrimalXStep);
    
    HDSDP_INIT(Hpsdp->dPrimalMatBuffer, double *, Hpsdp->nCones);
    HDSDP_MEMCHECK(Hpsdp->dPrimalMatBuffer);
    
    for ( int iCone = 0; iCone < Hpsdp->nCones; ++iCone ) {
        int nCol = HConeGetDim(Hpsdp->HCones[iCone]);
        HDSDP_INIT(Hpsdp->dPrimalX[iCone], double, nCol * nCol);
        HDSDP_MEMCHECK(Hpsdp->dPrimalX[iCone]);
        HDSDP_INIT(Hpsdp->dPrimalScalX[iCone], double, nCol * nCol);
        HDSDP_MEMCHECK(Hpsdp->dPrimalScalX[iCone]);
        HDSDP_INIT(Hpsdp->dPrimalXStep[iCone], double, nCol * nCol);
        HDSDP_MEMCHECK(Hpsdp->dPrimalXStep[iCone]);
        HDSDP_INIT(Hpsdp->dPrimalMatBuffer[iCone], double, nCol * nCol);
        HDSDP_MEMCHECK(Hpsdp->dPrimalMatBuffer[iCone]);
        HDSDPGetConeValues(Hpsdp->HSolver, iCone, Hpsdp->dPrimalX[iCone], NULL, Hpsdp->dPrimalMatBuffer[iCone]);
        HDSDP_MEMCPY(Hpsdp->dPrimalScalX[iCone], Hpsdp->dPrimalX[iCone], double, nCol * nCol);
    }
    
    int isInterior = 0;
    for ( int iCone = 0; iCone < Hpsdp->nCones; ++iCone ) {
        HFpLinsysPsdCheck(Hpsdp->XFactors[iCone], NULL, NULL, Hpsdp->dPrimalX[iCone], &isInterior);
    }
    
    if ( !isInterior ) {
        retcode = HDSDP_RETCODE_FAILED;
        goto exit_cleanup;
    }
    
    /* Register primal matrix in KKT solver */
    HKKTRegisterPSDP(Hpsdp->HKKT, Hpsdp->dPrimalScalX);
    hdsdp_printf("Primal solver is initialized. \n");
    
exit_cleanup:
    return retcode;
}


extern hdsdp_retcode HPSDPOptimize( hdsdp_psdp *Hpsdp ) {
    
    hdsdp_retcode retcode = HDSDP_RETCODE_OK;
    
    double *dPrimalKKTRhs = Hpsdp->dPrimalKKTRhs;
    
    double **dPrimalX = Hpsdp->dPrimalX;
    double **dPrimalScalX = Hpsdp->dPrimalScalX;
    double **dPrimalXStep = Hpsdp->dPrimalXStep;
    double **dPrimalMatBuffer = Hpsdp->dPrimalMatBuffer;
    
    double *dRowDual = Hpsdp->dRowDual;
    double *dRowDualStep = Hpsdp->dRowDualStep;
    double *dPrimalInfeas = Hpsdp->dPrimalAuxiVec1;
    double dBarrierMu = Hpsdp->dBarrierMu;
    
    hdsdp_cone **cones = Hpsdp->HCones;
    hdsdp *HSolver  = Hpsdp->HSolver;
    
    double dAbsoptTol = get_dbl_param(HSolver, DBL_PARAM_ABSOPTTOL);
    double dRelfeasTol = get_dbl_param(HSolver, DBL_PARAM_RELFEASTOL);
    double dReloptTol = get_dbl_param(HSolver, DBL_PARAM_RELOPTTOL);
    double dObjScal = get_dbl_feature(HSolver, DBL_FEATURE_OBJSCALING);
    double dRhsScal = get_dbl_feature(HSolver, DBL_FEATURE_RHSSCALING);
    double dSumDims = HSolver->dAllConeDims - HSolver->nRows * 2;
    double pdScal = 1.0 / (dObjScal * dRhsScal);
    
    int nIter = 0;
    int nMaxIter = 20;
    double dPrimalInfeasNorm = 0.0;
    
    /* Build KKT system */
    HDSDP_CALL(HKKTBuildUp(Hpsdp->HKKT, KKT_TYPE_PRIMAL));
    // hdsdp_printf("KKT Build up: %f s \n", HUtilGetTimeStamp() - HSolver->dTimeBegin);
    HDSDP_CALL(HKKTFactorize(Hpsdp->HKKT));
    // hdsdp_printf("KKT Factorize: %f s \n", HUtilGetTimeStamp() - HSolver->dTimeBegin);
    
    /* Enter primal main loop */
    for ( nIter = 0; nIter < nMaxIter; ++nIter ) {
        
        /* Compute primal infeasibility rp = A * X - b */
        HDSDP_ZERO(dPrimalInfeas, double, Hpsdp->nRow);
        for ( int iCone = 0; iCone < Hpsdp->nCones; ++iCone ) {
            HConeComputeATimesX(cones[iCone], dPrimalX[iCone], dPrimalInfeas);
        }
        dPrimalInfeasNorm = 0.0;
        for ( int iRow = 0; iRow < Hpsdp->nRow; ++iRow ) {
            double dTmp = Hpsdp->rowRHS[iRow] - dPrimalInfeas[iRow];
            dPrimalInfeasNorm += dTmp * dTmp;
        }
        dPrimalInfeasNorm = sqrt(dPrimalInfeasNorm);
                
        /* Prepare primal KKT RHS */
        /* Get buffer matrix X * S * X */
        int iDualMat = 1;
        for ( int iCone = 0; iCone < Hpsdp->nCones; ++iCone) {
            int nCol = HConeGetDim(cones[iCone]);
            HDSDP_ZERO(dPrimalMatBuffer[iCone], double, nCol * nCol);
            HConeBuildPrimalXSXDirection(cones[iCone], Hpsdp->HKKT, dPrimalX[iCone], dPrimalMatBuffer[iCone], iDualMat);
        }
        
        HDSDP_ZERO(dPrimalKKTRhs, double, Hpsdp->nRow);
        for ( int iCone = 0; iCone < Hpsdp->nCones; ++iCone) {
            HConeComputeATimesX(cones[iCone], dPrimalMatBuffer[iCone], dPrimalKKTRhs);
        }
        
        for ( int iRow = 0; iRow < Hpsdp->nRow; ++iRow ) {
            dPrimalKKTRhs[iRow] -= dBarrierMu * dPrimalInfeas[iRow];
        }
        
        /* Solve the KKT system */
        HDSDP_CALL(HKKTSolve(Hpsdp->HKKT, dPrimalKKTRhs, dRowDualStep));
        
        /* Get Primal direction */
        /* Compute X * S * X + XScal * dS * XScal */
        /* Dual ratio test: necessary for getting dS from buffer */
        double dMaxDualStep = HDSDP_INFINITY;
        double dDualStep = 0.0;
        double oneOverMu = - 1.0 / dBarrierMu;
        
        iDualMat = 0;
        for ( int iCone = 0; iCone < Hpsdp->nCones; ++iCone ) {
            int nCol = HConeGetDim(cones[iCone]);
            HDSDP_CALL(HConeRatioTest(cones[iCone], 0.0, dRowDualStep, 1.0, BUFFER_DUALVAR, &dDualStep));
            dMaxDualStep = HDSDP_MIN(dMaxDualStep, dDualStep);
            HConeBuildPrimalXSXDirection(cones[iCone], Hpsdp->HKKT, dPrimalScalX[iCone], dPrimalMatBuffer[iCone], iDualMat);
            /* Compute X - 1 / mu * (X * S * X + XScal * dS * XScal) */
            HDSDP_MEMCPY(dPrimalXStep[iCone], dPrimalX[iCone], double, nCol * nCol);
            int nElem = nCol * nCol;
            axpy(&nElem, &oneOverMu, dPrimalMatBuffer[iCone], &HIntConstantOne, dPrimalXStep[iCone], &HIntConstantOne);
        }
        
        /* Ratio test */
        double dMaxPrimalStep = HDSDP_INFINITY;
        double dPrimalStep = 0.0;
        for ( int iCone = 0; iCone < Hpsdp->nCones; ++iCone ) {
            Hpsdp->iLanczos = iCone;
            HDSDP_CALL(HLanczosSolve(Hpsdp->Lanczos[iCone], NULL, &dPrimalStep));
            dMaxPrimalStep = HDSDP_MIN(dMaxPrimalStep, dPrimalStep);
        }
        
        dMaxPrimalStep = HDSDP_MIN(0.7 * dMaxPrimalStep, 1.0);
        dMaxDualStep = HDSDP_MIN(0.7 * dMaxDualStep, 1.0);
        
        /* Take step */
        /* Dual step */
        axpy(&Hpsdp->nRow, &dMaxDualStep, dRowDualStep, &HIntConstantOne, dRowDual, &HIntConstantOne);
        
        /* Primal step */
        for ( int iCone = 0; iCone < Hpsdp->nCones; ++iCone ) {
            int nCol = HConeGetDim(cones[iCone]);
            int nElem = nCol * nCol;
            axpy(&nElem, &dMaxPrimalStep, dPrimalXStep[iCone], &HIntConstantOne, dPrimalX[iCone], &HIntConstantOne);
            HConeUpdate(cones[iCone], 1.0, dRowDual);
        }
        
        int isInterior = 0;
        
        for ( int iCone = 0; iCone < Hpsdp->nCones; ++iCone ) {
            HDSDP_CALL(HConeCheckIsInterior(cones[iCone], 1.0, dRowDual, &isInterior));
            if ( !isInterior ) {
                retcode = HDSDP_RETCODE_FAILED;
                goto exit_cleanup;
            }
        }
        
        for ( int iCone = 0; iCone < Hpsdp->nCones; ++iCone ) {
            HFpLinsysPsdCheck(Hpsdp->XFactors[iCone], NULL, NULL, dPrimalX[iCone], &isInterior);
            if ( !isInterior ) {
                retcode = HDSDP_RETCODE_FAILED;
                goto exit_cleanup;
            }
        }
        
        /* Update barrier parameter */
        double dDualObj = dot(&Hpsdp->nRow, Hpsdp->rowRHS, &HIntConstantOne, dRowDual, &HIntConstantOne);
        double dPrimalObj = 0.0;
        double dCompl = 0.0;
        for ( int iCone = 0; iCone < Hpsdp->nCones; ++iCone ) {
            dPrimalObj += HConeComputeTraceCX(cones[iCone], dPrimalX[iCone]);
            dCompl += HConeComputeXDotS(cones[iCone], dPrimalX[iCone]);
        }
        
        dCompl = dCompl / dSumDims;
        dBarrierMu = HDSDP_MIN(dBarrierMu, dCompl / (1.25 * dSumDims));
        
        /* Synchronize data to HSolver */
        HSolver->pObjInternal = dPrimalObj;
        HSolver->dObjInternal = dDualObj;
        HSolver->dObjVal = HSolver->dObjInternal * pdScal;
        HSolver->pObjVal = HSolver->pObjInternal * pdScal;
        HSolver->pInfeas = dPrimalInfeasNorm / (1 + get_dbl_feature(HSolver, DBL_FEATURE_RHSONENORM));
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
        
        if ( dMaxPrimalStep < 1e-02 && dMaxDualStep < 1e-02 ) {
            break;
        }
        
        if ( dCompl > HSolver->comp ) {
            break;
        }
        
        HSolver->comp = dCompl;
    }
    
    double dPDGap = HSolver->pObjVal - HSolver->dObjVal;
    dPDGap = dPDGap / (fabs(HSolver->pObjVal) + fabs(HSolver->dObjVal) + 1.0);
    double dPrimalInf = dPrimalInfeasNorm / (1.0 + get_dbl_feature(HSolver, DBL_FEATURE_RHSONENORM));
    double dCompl = HSolver->comp;
    dCompl = dCompl / (fabs(HSolver->pObjVal) + fabs(HSolver->dObjVal) + 1.0);
    hdsdp_printf("Primal solver ends. Primal. inf: %10.5e  Gap: %10.5e  Compl: %10.5e \n", dPrimalInf, dPDGap, dCompl);
    
exit_cleanup:
    
    return retcode;
}

extern void HPSDPClear( hdsdp_psdp *Hpsdp ) {
    
    if ( !Hpsdp ) {
        return;
    }
    
    HDSDP_FREE(Hpsdp->dPrimalKKTRhs);
    
    for ( int iCone = 0; iCone < Hpsdp->nCones; ++iCone ) {
        HLanczosDestroy(&Hpsdp->Lanczos[iCone]);
        HFpLinsysDestroy(&Hpsdp->XFactors[iCone]);
        HDSDP_FREE(Hpsdp->dPrimalX[iCone]);
        HDSDP_FREE(Hpsdp->dPrimalXStep[iCone]);
        HDSDP_FREE(Hpsdp->dPrimalScalX[iCone]);
        HDSDP_FREE(Hpsdp->dPrimalMatBuffer[iCone]);
    }
    
    HDSDP_FREE(Hpsdp->dPrimalX);
    HDSDP_FREE(Hpsdp->dPrimalScalX);
    HDSDP_FREE(Hpsdp->dPrimalXStep);
    HDSDP_FREE(Hpsdp->dPrimalMatBuffer);
    HDSDP_FREE(Hpsdp->Lanczos);
    HDSDP_FREE(Hpsdp->XFactors);
    
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
