#ifndef hdsdp_conic_sdp_h
#define hdsdp_conic_sdp_h

#ifdef HEADERPATH
#include "interface/def_hdsdp_conic.h"
#else
#include "def_hdsdp_conic.h"
#endif

#ifdef __cplusplus
extern "C" {
#endif

/* Dense SDP cone */
extern hdsdp_retcode sdpDenseConeCreateImpl( hdsdp_cone_sdp_dense **pCone );
extern hdsdp_retcode sdpDenseConeProcDataImpl( hdsdp_cone_sdp_dense *cone, int nRow, int nCol,
                                               int *coneMatBeg, int *coneMatIdx, double *coneMatElem );
extern hdsdp_retcode sdpDenseConePresolveImpl( hdsdp_cone_sdp_dense *cone );
extern void sdpDenseConeSetStartImpl( hdsdp_cone_sdp_dense *cone, double rResi );
extern double sdpDenseConeGetObjNorm( hdsdp_cone_sdp_dense *cone, int whichNorm );
extern double sdpDenseConeGetCoeffNorm( hdsdp_cone_sdp_dense *cone, int whichNorm );
extern void sdpDenseConeUpdateImpl( hdsdp_cone_sdp_dense *cone, double barHsdTau, double *rowDual );
extern hdsdp_retcode sdpDenseConeRatioTestImpl( hdsdp_cone_sdp_dense *cone, double barHsdTauStep, double *rowDualStep, double dAdaRatio, int whichBuffer, double *maxStep );
extern int64_t sdpDenseConeGetSymNnzImpl( hdsdp_cone_sdp_dense *cone );
extern void sdpDenseConeAddSymNnzImpl( hdsdp_cone_sdp_dense *cone, int iCol, int *schurMatCol );
extern void sdpDenseConeGetSymMapping( hdsdp_cone_sdp_dense *cone, int iCol, int *schurMatCol );
extern int sdpDenseConeGetDim( hdsdp_cone_sdp_dense *cone );
extern double sdpDenseConeGetObjNorm( hdsdp_cone_sdp_dense *cone, int whichNorm );
extern double sdpDenseConeGetCoeffNorm( hdsdp_cone_sdp_dense *cone, int whichNorm );
extern void sdpDenseConeScal( hdsdp_cone_sdp_dense *cone, double dScal );
extern hdsdp_retcode sdpDenseConeGetKKT( hdsdp_cone_sdp_dense *cone, int iCone, void *kkt, int typeKKT );
extern hdsdp_retcode sdpDenseConeGetKKTByFixedStrategy( hdsdp_cone_sdp_dense *cone, int iCone, void *kkt, int typeKKT, int kktStrategy );
extern void sdpDenseConeBuildPrimalXSXDirection( hdsdp_cone_sdp_dense *cone, void *kkt, double *dPrimalScalMatrix, double *dPrimalXSXBuffer, int iDualMat );
extern hdsdp_retcode sdpDenseConeInteriorCheck( hdsdp_cone_sdp_dense *cone, double barHsdTau, double *rowDual, int *isInterior );
extern hdsdp_retcode sdpDenseConeInteriorCheckExpert( hdsdp_cone_sdp_dense *cone, double dCCoef, double dACoefScal, double *dACoef, double dEyeCoef, int whichBuffer, int *isInterior );
extern void sdpDenseConeReduceResidual( hdsdp_cone_sdp_dense *cone, double resiReduction );
extern void sdpDenseConeSetPerturb( hdsdp_cone_sdp_dense *cone, double dDualPerturb );
extern hdsdp_retcode sdpDenseConeGetBarrier( hdsdp_cone_sdp_dense *cone, double barHsdTau, double *rowDual, int whichBuffer, double *logdet );
extern hdsdp_retcode sdpDenseConeAddStepToBufferAndCheck( hdsdp_cone_sdp_dense *cone, double dStep, int whichBuffer, int *isInterior );
extern void sdpDenseConeGetPrimal( hdsdp_cone_sdp_dense *cone, double dBarrierMu, double *dRowDual, double *dRowDualStep, double *dConePrimal, double *dAuxiMat );
extern void sdpDenseConeGetDual( hdsdp_cone_sdp_dense *cone, double *dConeDual, double *ddummy );
extern double sdpDenseConeTraceCX( hdsdp_cone_sdp_dense *cone, double *dConePrimel );
extern void sdpDenseConeATimesX( hdsdp_cone_sdp_dense *cone, double *dPrimalX, double *dATimesX );
extern double sdpDenseConeXDotS( hdsdp_cone_sdp_dense *cone, double *dConePrimal );
extern void sdpDenseConeClearImpl( hdsdp_cone_sdp_dense *cone );
extern void sdpDenseConeDestroyImpl( hdsdp_cone_sdp_dense **pCone );
extern void sdpDenseConeFeatureDetectImpl( hdsdp_cone_sdp_dense *cone, double *rowRHS, int coneIntFeatures[20], double coneDblFeatures[20] );
extern void sdpDenseConeViewImpl( hdsdp_cone_sdp_dense *cone );

/* Sparse SDP cone */
extern hdsdp_retcode sdpSparseConeCreateImpl( hdsdp_cone_sdp_sparse **pCone );
extern hdsdp_retcode sdpSparseConeProcDataImpl( hdsdp_cone_sdp_sparse *cone, int nRow, int nCol,
                                               int *coneMatBeg, int *coneMatIdx, double *coneMatElem );
extern hdsdp_retcode sdpSparseConePresolveImpl( hdsdp_cone_sdp_sparse *cone );
extern hdsdp_retcode sdpSparseConeIGetKKT( hdsdp_cone_sdp_sparse *cone, void *kkt );
extern void sdpSparseConeSetStartImpl( hdsdp_cone_sdp_sparse *cone, double rResi );
extern double sdpSparseConeGetObjNorm( hdsdp_cone_sdp_sparse *cone, int whichNorm );
extern double sdpSparseConeGetCoeffNorm( hdsdp_cone_sdp_sparse *cone, int whichNorm );
extern void sdpSparseConeScal( hdsdp_cone_sdp_sparse *cone, double dScal );
extern void sdpSparseConeUpdateImpl( hdsdp_cone_sdp_sparse *cone, double barHsdTau, double *rowDual );
extern hdsdp_retcode sdpSparseConeRatioTestImpl( hdsdp_cone_sdp_sparse *cone, double barHsdTauStep, double *rowDualStep, double dAdaRatio, int whichBuffer, double *maxStep );
extern int sdpSparseConeGetDim( hdsdp_cone_sdp_sparse *cone );
extern hdsdp_retcode sdpSparseConeGetKKT( hdsdp_cone_sdp_sparse *cone, int iCone, void *kkt, int typeKKT );
extern hdsdp_retcode sdpSparseConeGetKKTByFixedStrategy( hdsdp_cone_sdp_sparse *cone, int iCone, void *kkt, int typeKKT, int kktStrategy );
extern void sdpSparseConeBuildPrimalXSXDirection( hdsdp_cone_sdp_sparse *cone, void *kkt, double *dPrimalScalMatrix, double *dPrimalXSXBuffer, int iDualMat );
extern double sdpSparseConeXDotS( hdsdp_cone_sdp_sparse *cone, double *dConePrimal );
extern int64_t sdpSparseConeGetSymNnzImpl( hdsdp_cone_sdp_sparse *cone );
extern void sdpSparseConeAddSymNnzImpl( hdsdp_cone_sdp_sparse *cone, int iCol, int *schurMatCol );
extern void sdpSparseConeGetSymMapping( hdsdp_cone_sdp_sparse *cone, int iCol, int *schurMatCol );
extern hdsdp_retcode sdpSparseConeInteriorCheck( hdsdp_cone_sdp_sparse *cone, double barHsdTau, double *rowDual, int *isInterior );
extern hdsdp_retcode sdpSparseConeInteriorCheckExpert( hdsdp_cone_sdp_sparse *cone, double dCCoef, double dACoefScal, double *dACoef, double dEyeCoef, int whichBuffer, int *isInterior );
extern void sdpSparseConeReduceResidual( hdsdp_cone_sdp_sparse *cone, double resiReduction );
extern void sdpSparseConeSetPerturb( hdsdp_cone_sdp_sparse *cone, double dDualPerturb );
extern hdsdp_retcode sdpSparseConeGetBarrier( hdsdp_cone_sdp_sparse *cone, double barHsdTau, double *rowDual, int whichBuffer, double *logdet );
extern hdsdp_retcode sdpSparseConeAddStepToBufferAndCheck( hdsdp_cone_sdp_sparse *cone, double dStep, int whichBuffer, int *isInterior );
extern void sdpSparseConeGetPrimal( hdsdp_cone_sdp_sparse *cone, double dBarrierMu, double *dRowDual, double *dRowDualStep, double *dConePrimal, double *dAuxiMat );
extern void sdpSparseConeGetDual( hdsdp_cone_sdp_sparse *cone, double *dConeDual, double *ddummy );
extern double sdpSparseConeTraceCX( hdsdp_cone_sdp_sparse *cone, double *dConePrimel );
extern void sdpSparseConeATimesX( hdsdp_cone_sdp_sparse *cone, double *dPrimalX, double *dATimesX );
extern double sdpSparseConeXDotS( hdsdp_cone_sdp_sparse *cone, double *dConePrimal );
extern void sdpSparseConeClearImpl( hdsdp_cone_sdp_sparse *cone );
extern void sdpSparseConeDestroyImpl( hdsdp_cone_sdp_sparse **pCone );
extern void sdpSparseConeFeatureDetectImpl( hdsdp_cone_sdp_dense *cone, double *rowRHS, int coneIntFeatures[20], double coneDblFeatures[20] );
extern void sdpSparseConeViewImpl( hdsdp_cone_sdp_sparse *cone );

#ifdef __cplusplus
}
#endif

#endif /* hdsdp_conic_sdp_h */
