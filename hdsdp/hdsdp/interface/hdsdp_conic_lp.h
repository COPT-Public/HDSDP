#ifndef hdsdp_conic_lp_h
#define hdsdp_conic_lp_h

#include <stdio.h>

#ifdef HEADERPATH
#include "interface/def_hdsdp_conic.h"
#else
#include "def_hdsdp_conic.h"
#endif

#ifdef __cplusplus
extern "C" {
#endif

extern hdsdp_retcode LPConeCreateImpl( hdsdp_cone_lp **pCone );
extern hdsdp_retcode LPConeProcDataImpl( hdsdp_cone_lp *cone, int nRow, int nCol, int *coneMatBeg, int *coneMatIdx, double *coneMatElem );
extern hdsdp_retcode LPConePresolveImpl( hdsdp_cone_lp *cone );
extern void LPConeSetStartImpl( hdsdp_cone_lp *cone, double rResi );
extern double LPConeGetObjNorm( hdsdp_cone_lp *cone, int whichNorm );
extern double LPConeGetCoeffNorm( hdsdp_cone_lp *cone, int whichNorm );
extern void LPConeScal( hdsdp_cone_lp *cone, double dScal );
extern void LPConeUpdateImpl( hdsdp_cone_lp *cone, double barHsdTau, double *rowDual );
extern hdsdp_retcode LPConeRatioTestImpl( hdsdp_cone_lp *cone, double barHsdTauStep, double *rowDualStep, double dAdaRatio,
                                         int whichBuffer, double *maxStep );
extern int LPConeGetDim( hdsdp_cone_lp *cone );
extern hdsdp_retcode LPConeGetKKT( hdsdp_cone_lp *cone, int iCone, void *kkt, int typeKKT );
extern hdsdp_retcode LPConeGetKKTFixedStrategy( hdsdp_cone_lp *cone, int iCone, void *kkt, int typeKKT, int ikktStrategy );
extern void LPConeBuildPrimalXSXDirection( hdsdp_cone_lp *cone, void *kkt, double *dPrimalScalMatrix, double *dPrimalXSXBuffer, int iDualMat );
extern double LPConeXDotS( hdsdp_cone_lp *cone, double *dConePrimal );
extern int64_t LPConeGetSymNnzImpl( hdsdp_cone_lp *cone );
extern void LPConeAddSymNnzImpl( hdsdp_cone_sdp_dense *cone, int iCol, int *schurMatCol );
extern void LPConeGetSymMapping( hdsdp_cone_lp *cone, int iCol, int *schurMatCol );
extern hdsdp_retcode LPConeInteriorCheck( hdsdp_cone_lp *cone, double barHsdTau, double *rowDual, int *isInterior );
extern hdsdp_retcode LPConeInteriorCheckExpert( hdsdp_cone_lp *cone, double dCCoef, double dACoefScal, double *dACoef,
                                                double dEyeCoef, int whichBuffer, int *isInterior );
extern void LPConeReduceResidual( hdsdp_cone_lp *cone, double resiReduction );
extern void LPConeSetPerturb( hdsdp_cone_lp *cone, double dDualPerturb );
extern hdsdp_retcode LPConeGetBarrier( hdsdp_cone_lp *cone, double barHsdTau, double *rowDual, int whichBuffer, double *logdet );
extern hdsdp_retcode LPConeAddStepToBufferAndCheck( hdsdp_cone_lp *cone, double dStep, int whichBuffer, int *isInterior );
extern hdsdp_retcode LPConeGetPrimal( hdsdp_cone_lp *cone, double dBarrierMu, double *dRowDual, double *dRowDualStep, double *dConePrimal, double *dAuxiMat );
extern void LPConeGetDual( hdsdp_cone_lp *cone, double *dConeDual, double *ddummy );
extern double LPConeTraceCX( hdsdp_cone_lp *cone, double *dConePrimal );
extern void LPConeATimesX( hdsdp_cone_lp *cone, double *dPrimalX, double *dATimesX );
extern void LPConeClearImpl( hdsdp_cone_lp *cone );
extern void LPConeDestroyImpl( hdsdp_cone_lp **pCone );
extern void LPConeViewImpl( hdsdp_cone_lp *cone );
extern void LPConeGetStatsImpl( hdsdp_cone_lp *cone, double *rowRHS, int coneIntFeatures[20], double coneDblFeatures[20] );

#ifdef __cplusplus
}
#endif


#endif /* hdsdp_conic_lp_h */
