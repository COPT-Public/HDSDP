/** @file hdsdp\_sdpdata.h
 *  @brief Implenent the SDP Coeffcient data operations
 */
#ifndef hdsdp_sdpdata_h
#define hdsdp_sdpdata_h

#ifdef HEADERPATH
#include "interface/hdsdp.h"
#include "linalg/def_hdsdp_sdpdata.h"
#else
#include "hdsdp.h"
#include "def_hdsdp_sdpdata.h"
#endif

#ifdef __cplusplus
extern "C" {
#endif

extern hdsdp_retcode sdpDataMatCreate( sdp_coeff **psdpCoeff );
extern hdsdp_retcode sdpDataMatSetData( sdp_coeff *sdpCoeff, int nSDPCol, int dataMatNnz, int *dataMatIdx, double *dataMatElem );
extern int sdpDataMatGetRank( sdp_coeff *sdpCoeff );
extern void sdpDataMatScal( sdp_coeff *sdpCoeff, double scal );

extern double sdpDataMatNorm( sdp_coeff *sdpCoeff, int type );
extern hdsdp_retcode sdpDataMatBuildUpEigs( sdp_coeff *sdpCoeff, double *dAuxFullMatrix );
extern int sdpDataMatGetNnz( sdp_coeff *sdpCoeff );
extern void sdpDataMatDump( sdp_coeff *sdpCoeff, double *dFullMatrix );
extern void sdpDataMatGetMatNz( sdp_coeff *sdpCoeff, int *iMatSpsPattern );
extern void sdpDataMatAddToBuffer( sdp_coeff *sdpCoeff, double dElem, int *iMatSpsPattern, double *dBuffer );
extern sdp_coeff_type sdpDataMatGetType( sdp_coeff *sdpCoeff );
extern void sdpDataMatClear( sdp_coeff *sdpCoeff );
extern void sdpDataMatDestroy( sdp_coeff **psdpCoeff );
extern void sdpDataMatView( sdp_coeff *sdpCoeff );
extern int sdpDataMatIsEye( sdp_coeff *sdpCoeff, double *dEyeMultiple );
extern int sdpDataMatIsUnitCol( sdp_coeff *sdpCoeff, int *iUnitCol );

/* KKT operations */
extern void sdpDataMatKKT2SolveRankOne( sdp_coeff *sdpCoeff, hdsdp_linsys *dualFactor,
                                        double *dInvMatrix, double *dInvFactor, double *sign );
extern double sdpDataMatKKT2QuadForm( sdp_coeff *sdpCoeff, double *dQuadVector, double *dAuxiVec );
extern double sdpDataMatKKT2TraceASinv( sdp_coeff *sdpCoeff, double *dSinvAVec );

extern double sdpDataMatKKT3ComputeSinvASinv( sdp_coeff *sdpCoeff, hdsdp_linsys *dualFactor, double *dInvMatrix,
                                             double *dAuxiMat, double *dSinvASinvBuffer );
extern double sdpDataMatKKT3TraceABuffer( sdp_coeff *sdpCoeff, double *dSinvASinvBuffer, double *dAuxiMat );

extern double sdpDataMatKKT4ComputeASinv( sdp_coeff *sdpCoeff, hdsdp_linsys *dualFactor, double *dInvMatrix,
                                          double *dAuxiMat, double dResidual, double *dASinvBuffer );
extern double sdpDataMatKKT4TraceASinvBuffer( sdp_coeff *sdpCoeff, hdsdp_linsys *dualFactor, double *dInvMatrix,
                                             double *dASinvBuffer, double *dAuxiMat );

extern double sdpDataMatKKT5TraceASinvBSinv( sdp_coeff *sdpCoeff, sdp_coeff *sdpCoeff2, double *Sinv, double *dAuxiMat );
extern double sdpDataMatKKT5SinvADotSinv( sdp_coeff *sdpCoeff, hdsdp_linsys *dualFactor, double *dInvMatrix, double *dAuxiMat);

#ifdef __cplusplus
}
#endif

#endif /* hdsdp_sdpdata_h */
