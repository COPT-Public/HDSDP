/** @file hdsdp\_conic.h
 *  @brief Implement HDSDP general conic interface
 */

#ifdef HEADERPATH
#include "interface/hdsdp_utils.h"
#include "interface/def_hdsdp_user_data.h"
#include "interface/hdsdp_user_data.h"
#include "interface/def_hdsdp_conic.h"
#include "interface/hdsdp_conic_sdp.h"
#include "interface/hdsdp_conic_bound.h"
#else
#include "hdsdp_utils.h"
#include "def_hdsdp_user_data.h"
#include "hdsdp_user_data.h"
#include "def_hdsdp_conic.h"
#include "hdsdp_conic_sdp.h"
#include "hdsdp_conic_bound.h"
#endif

extern hdsdp_retcode HConeCreate( hdsdp_cone **pHCone ) {
    
    hdsdp_retcode retcode = HDSDP_RETCODE_OK;
    HDSDP_NULLCHECK(pHCone);
    hdsdp_cone *HCone = NULL;
    HDSDP_INIT(HCone, hdsdp_cone, 1);
    HDSDP_MEMCHECK(HCone);
    HDSDP_ZERO(HCone, hdsdp_cone, 1);
    *pHCone = HCone;
    
exit_cleanup:
    
    return retcode;
}

extern hdsdp_retcode HConeSetData( hdsdp_cone *HCone, user_data *usrData ) {
    
    hdsdp_retcode retcode = HDSDP_RETCODE_OK;
    
    HCone->usrData = usrData;
    HCone->cone = HUserDataChooseCone(usrData);
    
    switch ( HCone->cone ) {
        case HDSDP_CONETYPE_SCALAR_BOUND:
            set_func_pointer(HCone->coneCreate, sBoundConeCreateImpl);
            set_func_pointer(HCone->coneProcData, sBoundConeProcDataImpl);
            set_func_pointer(HCone->conePresolveData, sBoundConePresolveDummyImpl);
            set_func_pointer(HCone->coneDestroyData, sBoundConeDestroyImpl);
            set_func_pointer(HCone->coneSetStart, sBoundConeSetStartDummyImpl);
            set_func_pointer(HCone->coneUpdate, sBoundConeUpdateImpl);
            set_func_pointer(HCone->coneRatioTest, sBoundConeRatioTestImpl);
            set_func_pointer(HCone->coneGetSymNnz, sBoundConeGetSymNnzImpl);
            set_func_pointer(HCone->coneAddSymNz, sBoundConeAddSymNnzImpl);
            set_func_pointer(HCone->coneGetKKTMap, sBoundConeGetSymMappingImpl);
            set_func_pointer(HCone->coneGetDim, sBoundConeGetDimImpl);
            set_func_pointer(HCone->coneGetObjNorm, sBoundConeGetObjNormDummyImpl);
            set_func_pointer(HCone->coneGetCoeffNorm, sBoundConeGetCoeffNormDummyImpl);
            set_func_pointer(HCone->coneScal, sBoundConeScalDummyImpl);
            set_func_pointer(HCone->coneBuildSchur, sBoundConeGetKKT);
            set_func_pointer(HCone->coneBuildSchurFixed, sBoundConeGetKKTByFixedStrategy);
            set_func_pointer(HCone->coneGetBarrier, sBoundConeGetBarrier);
            set_func_pointer(HCone->coneAxpyBufferAndCheck, sBoundConeAddStepToBufferAndCheck);
            set_func_pointer(HCone->coneInteriorCheck, sBoundConeInteriorCheck);
            set_func_pointer(HCone->coneInteriorCheckExpert, sBoundConeInteriorCheckExpert);
            set_func_pointer(HCone->coneReduceResi, sBoundConeReduceResidual);
            set_func_pointer(HCone->coneSetPerturb, sBoundConeSetPerturb);
            set_func_pointer(HCone->conePFeasCheck, NULL);
            set_func_pointer(HCone->conePRecover, sBoundConeGetPrimal);
            set_func_pointer(HCone->coneView, sBoundConeViewImpl);
            break;
        case HDSDP_CONETYPE_LP:
            set_func_pointer(HCone->coneCreate, NULL);
            set_func_pointer(HCone->coneProcData, NULL);
            set_func_pointer(HCone->conePresolveData, NULL);
            set_func_pointer(HCone->coneDestroyData, NULL);
            set_func_pointer(HCone->coneSetStart, NULL);
            set_func_pointer(HCone->coneUpdate, NULL);
            set_func_pointer(HCone->coneRatioTest, NULL);
            set_func_pointer(HCone->coneGetSymNnz, NULL);
            set_func_pointer(HCone->coneAddSymNz, NULL);
            set_func_pointer(HCone->coneGetKKTMap, NULL);
            set_func_pointer(HCone->coneGetObjNorm, NULL);
            set_func_pointer(HCone->coneGetDim, NULL);
            set_func_pointer(HCone->coneGetObjNorm, NULL);
            set_func_pointer(HCone->coneGetCoeffNorm, NULL);
            set_func_pointer(HCone->coneScal, NULL);
            set_func_pointer(HCone->coneBuildSchur, NULL);
            set_func_pointer(HCone->coneBuildSchurFixed, NULL);
            set_func_pointer(HCone->coneGetBarrier, NULL);
            set_func_pointer(HCone->coneAxpyBufferAndCheck, NULL);
            set_func_pointer(HCone->coneInteriorCheck, NULL);
            set_func_pointer(HCone->coneInteriorCheckExpert, NULL);
            set_func_pointer(HCone->coneReduceResi, NULL);
            set_func_pointer(HCone->coneSetPerturb, sBoundConeSetPerturb);
            set_func_pointer(HCone->conePFeasCheck, NULL);
            set_func_pointer(HCone->conePRecover, NULL);
            set_func_pointer(HCone->coneView, NULL);
            break;
        case HDSDP_CONETYPE_DENSE_SDP:
            set_func_pointer(HCone->coneCreate, sdpDenseConeCreateImpl);
            set_func_pointer(HCone->coneProcData, sdpDenseConeProcDataImpl);
            set_func_pointer(HCone->conePresolveData, sdpDenseConePresolveImpl);
            set_func_pointer(HCone->coneDestroyData, sdpDenseConeDestroyImpl);
            set_func_pointer(HCone->coneSetStart, sdpDenseConeSetStartImpl);
            set_func_pointer(HCone->coneUpdate, sdpDenseConeUpdateImpl);
            set_func_pointer(HCone->coneRatioTest, sdpDenseConeRatioTestImpl);
            set_func_pointer(HCone->coneGetSymNnz, sdpDenseConeGetSymNnzImpl);
            set_func_pointer(HCone->coneAddSymNz, sdpDenseConeAddSymNnzImpl);
            set_func_pointer(HCone->coneGetKKTMap, sdpDenseConeGetSymMapping);
            set_func_pointer(HCone->coneGetDim, sdpDenseConeGetDim);
            set_func_pointer(HCone->coneGetObjNorm, sdpDenseConeGetObjNorm);
            set_func_pointer(HCone->coneGetCoeffNorm, sdpDenseConeGetCoeffNorm);
            set_func_pointer(HCone->coneScal, sdpDenseConeScal);
            set_func_pointer(HCone->coneBuildSchur, sdpDenseConeGetKKT);
            set_func_pointer(HCone->coneBuildSchurFixed, sdpDenseConeGetKKTByFixedStrategy);
            set_func_pointer(HCone->coneGetBarrier, sdpDenseConeGetBarrier);
            set_func_pointer(HCone->coneAxpyBufferAndCheck, sdpDenseConeAddStepToBufferAndCheck);
            set_func_pointer(HCone->coneInteriorCheck, sdpDenseConeInteriorCheck);
            set_func_pointer(HCone->coneInteriorCheckExpert, sdpDenseConeInteriorCheckExpert);
            set_func_pointer(HCone->coneReduceResi, sdpDenseConeReduceResidual);
            set_func_pointer(HCone->coneSetPerturb, sdpDenseConeSetPerturb);
            set_func_pointer(HCone->conePFeasCheck, NULL);
            set_func_pointer(HCone->conePRecover, NULL);
            set_func_pointer(HCone->coneView, sdpDenseConeViewImpl);
            break;
        case HDSDP_CONETYPE_SPARSE_SDP:
            set_func_pointer(HCone->coneCreate, sdpSparseConeCreateImpl);
            set_func_pointer(HCone->coneProcData, sdpSparseConeProcDataImpl);
            set_func_pointer(HCone->conePresolveData, sdpSparseConePresolveImpl);
            set_func_pointer(HCone->coneDestroyData, sdpSparseConeDestroyImpl);
            set_func_pointer(HCone->coneSetStart, sdpSparseConeSetStartImpl);
            set_func_pointer(HCone->coneUpdate, sdpSparseConeUpdateImpl);
            set_func_pointer(HCone->coneRatioTest, sdpSparseConeRatioTestImpl);
            set_func_pointer(HCone->coneGetDim, sdpSparseConeGetDim);
            set_func_pointer(HCone->coneGetCoeffNorm, sdpSparseConeGetCoeffNorm);
            set_func_pointer(HCone->coneGetSymNnz, sdpSparseConeGetSymNnzImpl);
            set_func_pointer(HCone->coneAddSymNz, sdpSparseConeAddSymNnzImpl);
            set_func_pointer(HCone->coneGetKKTMap, sdpSparseConeGetSymMapping);
            set_func_pointer(HCone->coneGetObjNorm, sdpSparseConeGetObjNorm);
            set_func_pointer(HCone->coneScal, sdpSparseConeScal);
            set_func_pointer(HCone->coneBuildSchur, sdpSparseConeGetKKT);
            set_func_pointer(HCone->coneBuildSchurFixed, sdpSparseConeGetKKTByFixedStrategy);
            set_func_pointer(HCone->coneGetBarrier, sdpSparseConeGetBarrier);
            set_func_pointer(HCone->coneAxpyBufferAndCheck, sdpSparseConeAddStepToBufferAndCheck);
            set_func_pointer(HCone->coneInteriorCheck, sdpSparseConeInteriorCheck);
            set_func_pointer(HCone->coneInteriorCheckExpert, sdpSparseConeInteriorCheckExpert);
            set_func_pointer(HCone->coneReduceResi, sdpSparseConeReduceResidual);
            set_func_pointer(HCone->coneSetPerturb, sdpSparseConeSetPerturb);
            set_func_pointer(HCone->conePFeasCheck, NULL);
            set_func_pointer(HCone->conePRecover, NULL);
            set_func_pointer(HCone->coneView, sdpSparseConeViewImpl);
            break;
        case HDSDP_CONETYPE_SOCP:
            retcode = HDSDP_RETCODE_FAILED;
            goto exit_cleanup;
        default:
            retcode = HDSDP_RETCODE_FAILED;
            goto exit_cleanup;
    }
    
exit_cleanup:
    
    return retcode;
}

extern void HConeClear( hdsdp_cone *HCone ) {
    
    if ( !HCone ) {
        return;
    }
    
    HCone->coneDestroyData(&HCone->coneData);
    HDSDP_ZERO(HCone, hdsdp_cone, 1);
    
    return;
}

extern void HConeDestroy( hdsdp_cone **pHCone ) {
    
    if ( !pHCone ) {
        return;
    }
    
    HConeClear(*pHCone);
    HDSDP_FREE(*pHCone);
    
    return;;
}

extern void HConeView( hdsdp_cone *HCone ) {
     
    HCone->coneView(HCone->coneData);
    
    return;
}

extern hdsdp_retcode HConeProcData( hdsdp_cone *HCone ) {
    
    hdsdp_retcode retcode = HDSDP_RETCODE_OK;
    user_data *usrData = (user_data *) HCone->usrData;
    
    HDSDP_CALL(HCone->coneCreate(&HCone->coneData));
    HDSDP_CALL(HCone->coneProcData(HCone->coneData, usrData->nConicRow, usrData->nConicCol,
                                   usrData->coneMatBeg, usrData->coneMatIdx, usrData->coneMatElem));
    
exit_cleanup:
    
    return retcode;
}

extern hdsdp_retcode HConePresolveData( hdsdp_cone *HCone ) {
    
    hdsdp_retcode retcode = HDSDP_RETCODE_OK;
    HDSDP_CALL(HCone->conePresolveData(HCone->coneData));
    
exit_cleanup:
    
    return retcode;
}

/* Conic algorithm interface */
extern void HConeSetStart( hdsdp_cone *HCone, double dConeStartVal ) {
    
    HCone->coneSetStart(HCone->coneData, dConeStartVal);
    return;
}

extern void HConeUpdate( hdsdp_cone *HCone, double barHsdTau, double *rowDual ) {
    
    HCone->coneUpdate(HCone->coneData, barHsdTau, rowDual);
    return;
}

extern hdsdp_retcode HConeRatioTest( hdsdp_cone *HCone, double barHsdTauStep, double *rowDualStep, double dAdaRatio, int whichBuffer, double *maxStep ) {
    
    return HCone->coneRatioTest(HCone->coneData, barHsdTauStep, rowDualStep, dAdaRatio, whichBuffer, maxStep);
}

/* Schur complement and algorithm iterates */
extern int64_t HConeGetSymNnz( hdsdp_cone *HCone ) {
    
    return HCone->coneGetSymNnz(HCone->coneData);
}

extern void HConeAddSymNz( hdsdp_cone *HCone, int iCol, int *schurMatCol ) {
    
    HCone->coneAddSymNz(HCone->coneData, iCol, schurMatCol);
    return;
}

extern void HConeGetSymMapping( hdsdp_cone *HCone, int iCol, int *schurMatCol ) {
    
    HCone->coneGetKKTMap(HCone->coneData, iCol, schurMatCol);
    return;
}

extern int HConeGetDim( hdsdp_cone *HCone ) {
    
    return HCone->coneGetDim(HCone->coneData);
}

extern double HConeGetCoeffNorm( hdsdp_cone *HCone, int whichNorm ) {
    
    return HCone->coneGetCoeffNorm(HCone->coneData, whichNorm);
}

extern double HConeGetObjNorm( hdsdp_cone *HCone, int whichNorm ) {
    
    return HCone->coneGetObjNorm(HCone->coneData, whichNorm);
}

extern hdsdp_retcode HConeBuildSchurComplement( hdsdp_cone *HCone, void *schurMat, int typeKKT ) {
    /* There are three types of KKT system we set up during the algorithm iterations,
       which we respectively denote by
     
     1. KKT_TYPE_INFEASIBLE
     KKT system in the infeasible-start self-dual model.
     Objective C relevant terms are not set up.
     
     2. KKT_TYPE_INFEAS_CORRECTOR
     KKT system in the infeasible start corrector
     The Schur complement matrix M is not set up and only dASinvRdSinvVec, dASinv are setup
     to generate corrector step
     
     3. KKT_TYPE_HOMOGENEOUS
     KKT system in the homogeneous self-dual model
     All the information involing Rd, A, and C will be set up
     (Most time-consuming and will only be invoked when the problem is considered dual infeasible)
     */
    
    return HCone->coneBuildSchur(HCone->coneData, schurMat, typeKKT);
}

extern hdsdp_retcode HConeBuildSchurComplementFixed( hdsdp_cone *HCone, void *schurMat, int typeKKT, int kktStrategy ) {
    
    return HCone->coneBuildSchurFixed(HCone->coneData, schurMat, typeKKT, kktStrategy);
}

/* Barrier, projection and recovery */
extern hdsdp_retcode HConeGetLogBarrier( hdsdp_cone *HCone, double barHsdTau, double *rowDual, int whichBuffer, double *logdet ) {
    
    hdsdp_retcode retcode = HDSDP_RETCODE_OK;
    HDSDP_CALL(HCone->coneGetBarrier(HCone->coneData, barHsdTau, rowDual, whichBuffer, logdet));
    
exit_cleanup:
    return retcode;
}

extern hdsdp_retcode HConeAddStepToBufferAndCheck( hdsdp_cone *HCone, double dStep, int whichBuffer, int *isInterior ) {
    
    hdsdp_retcode retcode = HDSDP_RETCODE_OK;
    HDSDP_CALL(HCone->coneAxpyBufferAndCheck(HCone->coneData, dStep, whichBuffer, isInterior));
    
exit_cleanup:
    return retcode;
}

extern hdsdp_retcode HConeCheckIsInterior( hdsdp_cone *HCone, double barHsdTau, double *rowDual, int *isInterior ) {
    
    hdsdp_retcode retcode = HDSDP_RETCODE_OK;
    HDSDP_CALL(HCone->coneInteriorCheck(HCone->coneData, barHsdTau, rowDual, isInterior));
    
exit_cleanup:
    return retcode;
}

extern hdsdp_retcode HConeCheckIsInteriorExpert( hdsdp_cone *HCone, double dCCoef, double dACoefScal, double *dACoef, double dEyeCoef, int whichBuffer, int *isInterior ) {
    
    hdsdp_retcode retcode = HDSDP_RETCODE_OK;
    HDSDP_CALL(HCone->coneInteriorCheckExpert(HCone->coneData, dCCoef, dACoefScal, dACoef, dEyeCoef, whichBuffer, isInterior));
    
exit_cleanup:
    return retcode;
}

extern void HConeReduceResi( hdsdp_cone *HCone, double resiReduction ) {
    
    HCone->coneReduceResi(HCone->coneData, resiReduction);
    return;
}

extern void HConeSetPerturb( hdsdp_cone *HCone, double dPerturb ) {
    
    HCone->coneSetPerturb(HCone->coneData, dPerturb);
    return;
}

extern int HConePFeasSolFound( hdsdp_cone *HCone, double barHsdTauStep, double *rowDualStep ) {
    
    return HCone->conePFeasCheck(HCone->coneData, barHsdTauStep, rowDualStep);
}

extern void HConeGetPrimal( hdsdp_cone *HCone, double dBarrierMu, double *dRowDual, double *dRowDualStep, double *dConePrimal, double *dConePrimal2 ) {
    
    HCone->conePRecover(HCone->coneData, dBarrierMu, dRowDual, dRowDualStep, dConePrimal, dConePrimal2);
    return;
}

extern void HConeScalByConstant( hdsdp_cone *HCone, double dScal ) {
    
    HCone->coneScal(HCone->coneData, dScal);
    return;
}
