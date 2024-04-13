#ifndef hdsdp_conic_bound_h
#define hdsdp_conic_bound_h

#ifdef HEADERPATH
#include "interface/def_hdsdp_conic.h"
#else
#include "def_hdsdp_conic.h"
#endif

#ifdef __cplusplus
extern "C" {
#endif

extern hdsdp_retcode sBoundConeCreateImpl( hdsdp_cone_bound_scalar **pCone );
extern hdsdp_retcode sBoundConeProcDataImpl( hdsdp_cone_bound_scalar *cone, int nRow, int nCol, int *coneMatBeg, int *coneMatIdx, double *coneMatElem );
extern hdsdp_retcode sBoundConePresolveDummyImpl( hdsdp_cone_bound_scalar *cone );
extern void sBoundConeSetStartDummyImpl( hdsdp_cone_bound_scalar *cone, double dummy );
extern double sBoundConeGetObjNormDummyImpl( hdsdp_cone_bound_scalar *cone, int dummy );
extern double sBoundConeGetCoeffNormDummyImpl( hdsdp_cone_bound_scalar *cone, int dummy );
extern void sBoundConeScalDummyImpl( hdsdp_cone_bound_scalar *cone );
extern void sBoundConeUpdateImpl( hdsdp_cone_bound_scalar *cone, double barHsdTau, double *rowDual );
extern hdsdp_retcode sBoundConeRatioTestImpl( hdsdp_cone_bound_scalar *cone, double barHsdTauStep, double *rowDualStep, double dAdaRatio, int whichBuffer, double *maxStep );
extern int sBoundConeGetDimImpl( hdsdp_cone_bound_scalar *cone );
extern hdsdp_retcode sBoundConeGetKKT( hdsdp_cone_bound_scalar *cone, int iCone, void *kkt, int typeKKT );
extern hdsdp_retcode sBoundConeGetKKTByFixedStrategy( hdsdp_cone_bound_scalar *cone, int iCone, void *kkt, int typeKKT, int ikktStrategy );
extern void sBoundConeBuildPrimalXSXDirection( hdsdp_cone_bound_scalar *cone, void *kkt, double *dPrimalScalMatrix, double *dPrimalXSXBuffer, int iDualMat );
extern double sBoundConeXDotS( hdsdp_cone_bound_scalar *cone, double *dConePrimal );
extern int64_t sBoundConeGetSymNnzImpl( hdsdp_cone_bound_scalar *cone );
extern void sBoundConeAddSymNnzImpl( hdsdp_cone_bound_scalar *cone, int iCol, int *schurMatCol );
extern void sBoundConeGetSymMappingImpl( hdsdp_cone_bound_scalar *cone, int dummy1, int *dummy2 );
extern hdsdp_retcode sBoundConeInteriorCheck( hdsdp_cone_bound_scalar *cone, double barHsdTau, double *rowDual, int *isInterior );
extern hdsdp_retcode sBoundConeInteriorCheckExpert( hdsdp_cone_bound_scalar *cone, double dCCoef, double dACoefScal, double *dACoef, double dEyeCoef, int whichBuffer, int *isInterior );
extern void sBoundConeReduceResidual( hdsdp_cone_bound_scalar *cone, double dummy );
extern void sBoundConeSetPerturb( hdsdp_cone_bound_scalar *cone, double dummy );
extern hdsdp_retcode sBoundConeGetBarrier( hdsdp_cone_bound_scalar *cone, double barHsdTau, double *rowDual, int whichBuffer, double *logdet );
extern hdsdp_retcode sBoundConeAddStepToBufferAndCheck( hdsdp_cone_bound_scalar *cone, double dStep, int whichBuffer, int *isInterior );
extern void sBoundConeGetPrimal( hdsdp_cone_bound_scalar *cone, double dBarrierMu, double *dRowDual, double *dRowDualStep, double *dBoundLowerPrimal, double *dBoundUpperPrimal );
extern void sBoundConeClearImpl( hdsdp_cone_bound_scalar *cone );
extern void sBoundConeDestroyImpl( hdsdp_cone_bound_scalar **pCone );
extern void sBoundConeViewImpl( hdsdp_cone_bound_scalar *cone );

#ifdef __cplusplus
}
#endif

#endif /* hdsdp_conic_bound_h */
