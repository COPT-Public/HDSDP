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
    double  spR1FactorSign;
    int     nSpR1FactorElem;
    int    *spR1MatIdx;
    double *spR1MatElem;
    
} sdpRankOneSparseData;


typedef struct {
    
    int     nSDPCol;
    double  r1FactorSign;
    double *r1MatFactor;
    
} sdpRankOneDenseData;

extern void dataMatScalSparse( void *A, double alpha );
extern void dataMatScalDense( void *A, double alpha );
extern void dataMatScalRankOneSparse( void *A, double alpha );
extern void dataMatScalRankOneDense( void *A, double alpha );

extern double dataMatNormSparse( void *A, int type );
extern double dataMatNormDense( void *A, int type );
extern double dataMatNormRankOneSparse( void *A, int type );
extern double dataMatNormRankOneDense( void *A, int type );

extern int dataMatGetNnzSparse( void *A );
extern int dataMatGetNnzDense( void *A );
extern int dataMatGetNnzRankOneSparse( void *A );
extern int dataMatGetNnzRankOneDense( void *A );

extern void dataMatDumpSparse( void *A, double *v );
extern void dataMatDumpDense( void *A, double *v );
extern void dataMatDumpRankOneSparse( void *A, double *v );
extern void dataMatDumpRankOneDense( void *A, double *v );

#endif /* hdsdp_sdpdata_h */
