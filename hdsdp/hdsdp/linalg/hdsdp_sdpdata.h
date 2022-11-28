/** @file hdsdp\_sdpdata.h
 *  @brief Implenent the SDP Coeffcient data operations
 */
#ifndef hdsdp_sdpdata_h
#define hdsdp_sdpdata_h

#include "hdsdp.h"

/** @struct eigFactor
 *  @brief The eigen decomposition structure
 */
typedef struct {
    
    int     nCol;
    int     rank;
    double *eVals;
    double *eVecs;
    
} eigFactor;

typedef struct {
    
    int        nSDPCol;
    void      *dataMat;
    eigFactor *eig;
    
    void (*dataMataApB)          ( void *, double, void * );
    double (*dataMatDot)         ( void *, double * );
    void (*dataMatScal)          ( void *, double );
    void (*dataMatNorm)          ( void *, int );
    hdsdp_retcode (*dataMatEig)  ( void *, eigFactor ** );
    int (*dataMatGetNnz)         ( void * );
    void (*dataMatDump)          ( void *, double * );
    
} sdpCoeffMat;

/* Implementations of the SDP coefficient matrix */
typedef struct {
    
    int nSDPCol;
    
} sdpZeroData;

typedef struct {
    
    int     nSDPCol;
    int     nTriMatElem;
    int    *triMatCol;
    int    *triMatRow;
    double *triMatElem;
    
} sdpSparseData;


typedef struct {
    
    int     nSDPCol;
    double *dsMatElem;
    
} sdpDenseData;


typedef struct {
    
    int     nSDPCol;
    double  spR1FactorSign; ///< Include scale, may not equal to +1.0 or -1.0
    int     nSpR1FactorElem;
    int    *spR1MatIdx;
    double *spR1MatElem;
    
} sdpRankOneSparseData;


typedef struct {
    
    int     nSDPCol;
    double  r1FactorSign; ///< Include scale, may not equal to +1.0 or -1.0
    double *r1MatFactor;
    
} sdpRankOneDenseData;

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
