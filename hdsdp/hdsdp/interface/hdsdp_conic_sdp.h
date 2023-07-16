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
extern hdsdp_retcode sdpDenseConeRatioTestImpl( hdsdp_cone_sdp_dense *cone, double barHsdTauStep, double *rowDualStep, double dAdaRatio, double *maxStep );
extern int64_t sdpDenseConeGetSymNnzImpl( hdsdp_cone_sdp_dense *cone );
extern void sdpDenseConeAddSymNnzImpl( hdsdp_cone_sdp_dense *cone, int iCol, int *schurMatCol );
extern void sdpDenseConeGetSymMapping( hdsdp_cone_sdp_dense *cone, int iCol, int *schurMatCol );
extern int sdpDenseConeGetDim( hdsdp_cone_sdp_dense *cone );
extern hdsdp_retcode sdpDenseConeGetKKT( hdsdp_cone_sdp_dense *cone, void *kkt, int newKKT );
extern hdsdp_retcode sdpDenseConeGetBarrier( hdsdp_cone_sdp_dense *cone, double barHsdTau, double *rowDual, double *logdet );
extern void sdpDenseConeClearImpl( hdsdp_cone_sdp_dense *cone );
extern void sdpDenseConeDestroyImpl( hdsdp_cone_sdp_dense **pCone );
extern void sdpDenseConeViewImpl( hdsdp_cone_sdp_dense *cone );

/* Sparse SDP cone */
extern hdsdp_retcode sdpSparseConeCreateImpl( hdsdp_cone_sdp_sparse **pCone );
extern hdsdp_retcode sdpSparseConeProcDataImpl( hdsdp_cone_sdp_sparse *cone, int nRow, int nCol,
                                               int *coneMatBeg, int *coneMatIdx, double *coneMatElem );
extern hdsdp_retcode sdpSparseConePresolveImpl( hdsdp_cone_sdp_sparse *cone );
extern hdsdp_retcode sdpSparseConeIGetKKT( hdsdp_cone_sdp_sparse *cone, void *kkt );
extern void sdpSparseConeSetStartImpl( hdsdp_cone_sdp_sparse *cone, double rResi );
extern double sdpSparseConeGetObjNorm( hdsdp_cone_sdp_sparse *cone, int whichNorm );
extern double sdpDenseConeGetCoeffNorm( hdsdp_cone_sdp_dense *cone, int whichNorm );
extern void sdpSparseConeUpdateImpl( hdsdp_cone_sdp_sparse *cone, double barHsdTau, double *rowDual );
extern hdsdp_retcode sdpSparseConeRatioTestImpl( hdsdp_cone_sdp_sparse *cone, double barHsdTauStep, double *rowDualStep, double dAdaRatio, double *maxStep );
extern int sdpSparseConeGetDim( hdsdp_cone_sdp_sparse *cone );
extern hdsdp_retcode sdpSparseConeGetKKT( hdsdp_cone_sdp_sparse *cone, void *kkt, int newKKT );
extern int64_t sdpSparseConeGetSymNnzImpl( hdsdp_cone_sdp_sparse *cone );
extern void sdpSparseConeAddSymNnzImpl( hdsdp_cone_sdp_sparse *cone, int iCol, int *schurMatCol );
extern void sdpSparseConeGetSymMapping( hdsdp_cone_sdp_sparse *cone, int iCol, int *schurMatCol );
extern hdsdp_retcode sdpSparseConeGetBarrier( hdsdp_cone_sdp_sparse *cone, double barHsdTau, double *rowDual, double *logdet );
extern void sdpSparseConeClearImpl( hdsdp_cone_sdp_sparse *cone );
extern void sdpSparseConeDestroyImpl( hdsdp_cone_sdp_sparse **pCone );
extern void sdpSparseConeViewImpl( hdsdp_cone_sdp_sparse *cone );

#ifdef __cplusplus
}
#endif

#endif /* hdsdp_conic_sdp_h */
