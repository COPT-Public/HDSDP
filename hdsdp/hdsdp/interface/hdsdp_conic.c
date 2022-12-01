/** @file hdsdp\_conic.h
 *  @brief Implement HDSDP general conic interface
 */

#include "hdsdp_utils.h"
#include "def_hdsdp_conic.h"

static void HConeISetDenseSDPMethods( hdsdp_cone *HCone ) {
    
    HCone->coneInitData = NULL;
    HCone->coneSetData = NULL;
    HCone->coneProcData = NULL;
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
    
    return;
}

static void HConeISetSparseSDPMethods( hdsdp_cone *HCone ) {
    
    HCone->coneInitData = NULL;
    HCone->coneSetData = NULL;
    HCone->coneProcData = NULL;
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
    
    return;
}

static void HConeISetLPMethods( hdsdp_cone *HCone ) {
    
    HCone->coneInitData = NULL;
    HCone->coneSetData = NULL;
    HCone->coneProcData = NULL;
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
    
    return;
}

extern hdsdp_retcode HConeCreate( hdsdp_cone **pHCone ) {
    
    hdsdp_retcode retcode = HDSDP_RETCODE_OK;
    
    if ( !pHCone ) {
        retcode = HDSDP_RETCODE_FAILED;
        goto exit_cleanup;
    }
    
    hdsdp_cone *HCone = NULL;
    HDSDP_INIT(HCone, hdsdp_cone, 1);
    
    if ( !HCone ) {
        retcode = HDSDP_RETCODE_MEMORY;
        goto exit_cleanup;
    }
    
    HDSDP_ZERO(HCone, hdsdp_cone, 1);
    *pHCone = HCone;
    
exit_cleanup:
    
    return retcode;
}

extern hdsdp_retcode HConeInit( hdsdp_cone *HCone, cone_type cone ) {
    
    hdsdp_retcode retcode = HDSDP_RETCODE_OK;
    
    if ( !HCone ) {
        retcode = HDSDP_RETCODE_FAILED;
        goto exit_cleanup;
    }
    
    switch ( cone ) {
        case HDSDP_CONETYPE_LP:
            HConeISetLPMethods(HCone);
            break;
        case HDSDP_CONETYPE_DENSE_SDP:
            HConeISetDenseSDPMethods(HCone);
            break;
        case HDSDP_CONETYPE_SPARSE_SDP:
            HConeISetSparseSDPMethods(HCone);
            break;
        case HDSDP_CONETYPE_SOCP:
            retcode = HDSDP_RETCODE_FAILED;
            goto exit_cleanup;
        default:
            retcode = HDSDP_RETCODE_FAILED;
            goto exit_cleanup;
    }
    
    HCone->cone = cone;
    
exit_cleanup:
    
    return retcode;
}

/* Conic data interface */
extern hdsdp_retcode HConeInitData( hdsdp_cone *HCone ) {
    
    return HCone->coneInitData(&HCone->coneData);
}


extern hdsdp_retcode HConeSetData( hdsdp_cone *HCone, void *coneData ) {
    
    return HCone->coneSetData(HCone->coneData, coneData);
}

extern hdsdp_retcode HConeProcData( hdsdp_cone *HCone ) {
    
    return HCone->coneProcData(HCone->coneData);
}

extern void HConeDestroyData( hdsdp_cone *HCone ) {
    
    HCone->coneDestroyData(&HCone->coneData);
    return;
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

extern void HConeBuildSchurComplement( hdsdp_cone *HCone, void *schurMat ) {
    
    HCone->coneBuildSchur(HCone->coneData, schurMat);
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
