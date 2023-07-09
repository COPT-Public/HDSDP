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
extern double sdpDataMatDot( sdp_coeff *sdpCoeff, double *dFullMatrix );
extern void sdpDataMatScal( sdp_coeff *sdpCoeff, double scal );

#define ABS_NORM (1)
#define FRO_NORM (2)
extern double sdpDataMatNorm( sdp_coeff *sdpCoeff, int type );
extern hdsdp_retcode sdpDataMatBuildUpEigs( sdp_coeff *sdpCoeff, double *dAuxFullMatrix );
extern int sdpDataMatGetNnz( sdp_coeff *sdpCoeff );
extern void sdpDataMatDump( sdp_coeff *sdpCoeff, double *dFullMatrix );
extern inline void sdpDataMatGetMatNz( sdp_coeff *sdpCoeff, int *iMatSpsPattern );
extern sdp_coeff_type sdpDataMatGetType( sdp_coeff *sdpCoeff );
extern void sdpDataMatClear( sdp_coeff *sdpCoeff );
extern void sdpDataMatDestroy( sdp_coeff **psdpCoeff );
extern void sdpDataMatView( sdp_coeff *sdpCoeff );

#ifdef __cplusplus
}
#endif

#endif /* hdsdp_sdpdata_h */
