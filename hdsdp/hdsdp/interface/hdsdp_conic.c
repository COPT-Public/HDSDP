/** @file hdsdp\_conic.h
 *  @brief Implement HDSDP general conic interface
 */

#ifdef HEADERPATH
#include "interface/hdsdp_utils.h"
#include "interface/def_hdsdp_user_data.h"
#include "interface/hdsdp_user_data.h"
#include "interface/def_hdsdp_conic.h"
#include "interface/hdsdp_conic_sdp.h"
#else
#include "hdsdp_utils.h"
#include "def_hdsdp_user_data.h"
#include "hdsdp_user_data.h"
#include "def_hdsdp_conic.h"
#include "hdsdp_conic_sdp.h"
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
        case HDSDP_CONETYPE_BOUND:
            HCone->coneCreate = NULL;
            HCone->coneProcData = NULL;
            HCone->conePresolveData = NULL;
            HCone->coneDestroyData = NULL;
            HCone->coneSetStart = NULL;
            HCone->coneUpdate = NULL;
            HCone->coneRatioTest = NULL;
            HCone->coneGetSymNnz = NULL;
            HCone->coneAddSymNz = NULL;
            HCone->coneBuildSchur = NULL;
            HCone->coneGetBarrier = NULL;
            HCone->conePFeasCheck = NULL;
            HCone->conePRecover = NULL;
            HCone->coneScal = NULL;
            break;
        case HDSDP_CONETYPE_LP:
            HCone->coneCreate = NULL;
            HCone->coneProcData = NULL;
            HCone->conePresolveData = NULL;
            HCone->coneDestroyData = NULL;
            HCone->coneSetStart = NULL;
            HCone->coneUpdate = NULL;
            HCone->coneRatioTest = NULL;
            HCone->coneGetSymNnz = NULL;
            HCone->coneAddSymNz = NULL;
            HCone->coneBuildSchur = NULL;
            HCone->coneGetBarrier = NULL;
            HCone->conePFeasCheck = NULL;
            HCone->conePRecover = NULL;
            HCone->coneScal = NULL;
            break;
        case HDSDP_CONETYPE_DENSE_SDP:
            HCone->coneCreate = sdpDenseConeCreateImpl;
            HCone->coneProcData = sdpDenseConeProcDataImpl;
            HCone->conePresolveData = sdpDenseConePresolveImpl;
            HCone->coneDestroyData = sdpDenseConeDestroyImpl;
            HCone->coneSetStart = sdpDenseConeSetStartImpl;
            HCone->coneUpdate = sdpDenseConeUpdateImpl;
            HCone->coneRatioTest = sdpDenseConeRatioTestImpl;
            HCone->coneGetSymNnz = sdpDenseConeGetSymNnzImpl;
            HCone->coneAddSymNz = sdpDenseConeAddSymNnzImpl;
            HCone->coneGetObjNorm = sdpDenseConeGetObjNorm;
            HCone->coneBuildSchur = NULL;
            HCone->coneGetBarrier = NULL;
            HCone->conePFeasCheck = NULL;
            HCone->conePRecover = NULL;
            HCone->coneScal = NULL;
            HCone->coneView = sdpDenseConeViewImpl;
            break;
        case HDSDP_CONETYPE_SPARSE_SDP:
            HCone->coneCreate = sdpSparseConeCreateImpl;
            HCone->coneProcData = sdpSparseConeProcDataImpl;
            HCone->conePresolveData = sdpSparseConePresolveImpl;
            HCone->coneDestroyData = sdpSparseConeDestroyImpl;
            HCone->coneSetStart = sdpSparseConeSetStartImpl;
            HCone->coneUpdate = sdpSparseConeUpdateImpl;
            HCone->coneRatioTest = sdpSparseConeRatioTestImpl;
            HCone->coneGetSymNnz = sdpSparseConeGetSymNnzImpl;
            HCone->coneAddSymNz = sdpSparseConeAddSymNnzImpl;
            HCone->coneBuildSchur = sdpSparseConeGetObjNorm;
            HCone->coneGetBarrier = NULL;
            HCone->conePFeasCheck = NULL;
            HCone->conePRecover = NULL;
            HCone->coneScal = NULL;
            HCone->coneView = sdpSparseConeViewImpl;
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

extern double HConeRatioTest( hdsdp_cone *HCone, double barHsdTauStep, double *rowDualStep ) {
    
    return HCone->coneRatioTest(HCone->coneData, barHsdTauStep, rowDualStep);
}

/* Schur complement and algorithm iterates */
extern int64_t HConeGetSymNnz( hdsdp_cone *HCone ) {
    
    return HCone->coneGetSymNnz(HCone->coneData);
}

extern void HConeAddSymNz( hdsdp_cone *HCone, int *schurMatCol ) {
    
    HCone->coneAddSymNz(HCone->coneData, schurMatCol);
    return;
}

extern void HConeBuildSchurComplement( hdsdp_cone *HCone, int iCol, void *schurMat ) {
    
    HCone->coneBuildSchur(HCone->coneData, iCol, schurMat);
    return;
}

/* Barrier, projection and recovery */
extern void HConeGetLogBarrier( hdsdp_cone *HCone, double barHsdTau, double *rowDual ) {
    
    HCone->coneGetBarrier(HCone->coneData, barHsdTau, rowDual);
    return;
}

extern int HConePFeasSolFound( hdsdp_cone *HCone, double barHsdTauStep, double *rowDualStep ) {
    
    return HCone->conePFeasCheck(HCone->coneData, barHsdTauStep, rowDualStep);
}

extern void HConePVarRecover( hdsdp_cone *HCone, double *pVarArr ) {
    
    HCone->conePRecover(HCone->coneData, pVarArr);
    return;
}

extern void HConeScalByConstant( hdsdp_cone *HCone, double dScal ) {
    
    HCone->coneScal(HCone->coneData, dScal);
    return;
}
