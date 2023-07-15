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

extern hdsdp_retcode HConeCreate( hdsdp_cone **pHCone );
extern hdsdp_retcode HConeSetData( hdsdp_cone *HCone, user_data *coneData );
extern hdsdp_retcode HConeProcData( hdsdp_cone *HCone );
extern hdsdp_retcode HConePresolveData( hdsdp_cone *HCone );
extern void HConeClear( hdsdp_cone *HCone );
extern void HConeDestroy( hdsdp_cone **pHCone );
extern void HConeView( hdsdp_cone *HCone );

extern void HConeSetStart( hdsdp_cone *HCone, double dConeStartVal );
extern void HConeUpdate( hdsdp_cone *HCone, double barHsdTau, double *rowDual );
extern double HConeRatioTest( hdsdp_cone *HCone, double barHsdTauStep, double *rowDualStep );
extern int64_t HConeGetSymNnz( hdsdp_cone *HCone );
extern void HConeAddSymNz( hdsdp_cone *HCone, int iCol, int *schurMatCol );
extern void HConeGetSymMapping( hdsdp_cone *HCone, int iCol, int *schurMatCol );
extern int HConeGetDim( hdsdp_cone *HCone );
extern hdsdp_retcode HConeBuildSchurComplement( hdsdp_cone *HCone, void *schurMat, int newKKT );
extern hdsdp_retcode HConeGetLogBarrier( hdsdp_cone *HCone, double barHsdTau, double *rowDual, double *logdet );
extern int HConePFeasSolFound( hdsdp_cone *HCone, double barHsdTauStep, double *rowDualStep );
extern void HConePVarRecover( hdsdp_cone *HCone, double *pVarArr );
extern void HConeScalByConstant( hdsdp_cone *HCone, double dScal );

#ifdef __cplusplus
}
#endif

#endif /* hdsdp_conic_h */
