#ifndef hdsdp_conic_h
#define hdsdp_conic_h

#ifdef HEADERPATH
#include "interface/def_hdsdp_conic.h"
#include "interface/hdsdp_user_data.h"
#else
#include "def_hdsdp_conic.h"
#include "hdsdp_user_data.h"
#endif

#ifdef __cplusplus
extern "C" {
#endif

#define KKT_TYPE_INFEASIBLE  (0)
#define KKT_TYPE_CORRECTOR   (1)
#define KKT_TYPE_HOMOGENEOUS (2)

#define ABS_NORM    (1)
#define FRO_NORM    (2)

#define BUFFER_DUALVAR   (0)
#define BUFFER_DUALCHECK (1)
#define BUFFER_DUALSTEP  (2)
extern hdsdp_retcode HConeCreate( hdsdp_cone **pHCone );
extern hdsdp_retcode HConeSetData( hdsdp_cone *HCone, user_data *coneData );
extern hdsdp_retcode HConeProcData( hdsdp_cone *HCone );
extern hdsdp_retcode HConePresolveData( hdsdp_cone *HCone );
extern void HConeClear( hdsdp_cone *HCone );
extern void HConeDestroy( hdsdp_cone **pHCone );
extern void HConeView( hdsdp_cone *HCone );

extern void HConeSetStart( hdsdp_cone *HCone, double dConeStartVal );
extern void HConeUpdate( hdsdp_cone *HCone, double barHsdTau, double *rowDual );
extern hdsdp_retcode HConeRatioTest( hdsdp_cone *HCone, double barHsdTauStep, double *rowDualStep, double dAdaRatio, int whichBuffer, double *maxStep );
extern int64_t HConeGetSymNnz( hdsdp_cone *HCone );
extern void HConeAddSymNz( hdsdp_cone *HCone, int iCol, int *schurMatCol );
extern void HConeGetSymMapping( hdsdp_cone *HCone, int iCol, int *schurMatCol );
extern int HConeGetDim( hdsdp_cone *HCone );
extern double HConeGetCoeffNorm( hdsdp_cone *HCone, int whichNorm );
extern double HConeGetObjNorm( hdsdp_cone *HCone, int whichNorm );
extern hdsdp_retcode HConeBuildSchurComplement( hdsdp_cone *HCone, void *schurMat, int typeKKT );
extern hdsdp_retcode HConeBuildSchurComplementFixed( hdsdp_cone *HCone, void *schurMat, int typeKKT, int kktStrategy );
extern hdsdp_retcode HConeGetLogBarrier( hdsdp_cone *HCone, double barHsdTau, double *rowDual, int whichBuffer, double *logdet );
extern hdsdp_retcode HConeAddStepToBufferAndCheck( hdsdp_cone *HCone, double dStep, int whichBuffer, int *isInterior );
extern hdsdp_retcode HConeCheckIsInterior( hdsdp_cone *HCone, double barHsdTau, double *rowDual, int *isInterior );
extern hdsdp_retcode HConeCheckIsInteriorExpert( hdsdp_cone *HCone, double dCCoef, double dACoefScal, double *dACoef, double dEyeCoef, int whichBuffer, int *isInterior );
extern void HConeReduceResi( hdsdp_cone *HCone, double resiReduction );
extern void HConeSetPerturb( hdsdp_cone *HCone, double dPerturb );
extern int HConePFeasSolFound( hdsdp_cone *HCone, double barHsdTauStep, double *rowDualStep );
extern void HConeGetPrimal( hdsdp_cone *HCone, double dBarrierMu, double *dRowDual, double *dRowDualStep, double *dConePrimal, double *dConePrimal2 );
extern void HConeScalByConstant( hdsdp_cone *HCone, double dScal );

#ifdef __cplusplus
}
#endif

#endif /* hdsdp_conic_h */
