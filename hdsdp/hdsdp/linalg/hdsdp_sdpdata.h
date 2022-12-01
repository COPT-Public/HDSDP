/** @file hdsdp\_sdpdata.h
 *  @brief Implenent the SDP Coeffcient data operations
 */
#ifndef hdsdp_sdpdata_h
#define hdsdp_sdpdata_h

#include "hdsdp.h"
#include "def_hdsdp_sdpdata.h"

#ifdef __cplusplus
extern "C" {
#endif

#define ABS_NORM 1
#define FRO_NORM 2

extern void dataMatScalZeroImpl( void *A, double alpha );
extern void dataMatScalSparseImpl( void *A, double alpha );
extern void dataMatScalDenseImpl( void *A, double alpha );
extern void dataMatScalRankOneSparseImpl( void *A, double alpha );
extern void dataMatScalRankOneDenseImpl( void *A, double alpha );

extern double dataMatNormZeroImpl( void *A, int type );
extern double dataMatNormSparseImpl( void *A, int type );
extern double dataMatNormDenseImpl( void *A, int type );
extern double dataMatNormRankOneSparseImpl( void *A, int type );
extern double dataMatNormRankOneDenseImpl( void *A, int type );

extern int dataMatGetNnzZeroImpl( void *A );
extern int dataMatGetNnzSparseImpl( void *A );
extern int dataMatGetNnzDenseImpl( void *A );
extern int dataMatGetNnzRankOneSparseImpl( void *A );
extern int dataMatGetNnzRankOneDenseImpl( void *A );

extern void dataMatDumpZeroImpl( void *A, double *v );
extern void dataMatDumpSparseImpl( void *A, double *v );
extern void dataMatDumpDenseImpl( void *A, double *v );
extern void dataMatDumpRankOneSparseImpl( void *A, double *v );
extern void dataMatDumpRankOneDenseImpl( void *A, double *v );

#ifdef __cplusplus
}
#endif

#endif /* hdsdp_sdpdata_h */
