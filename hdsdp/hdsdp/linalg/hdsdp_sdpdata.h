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
    int     nSpR1FactorElem;
    int    *spR1MatIdx;
    double *spR1MatElem;
    
} sdpRankOneSparseData;


typedef struct {
    
    int     nSDPCol;
    double *r1MatFactor;
    
} sdpRankOneDenseData;


#endif /* hdsdp_sdpdata_h */
