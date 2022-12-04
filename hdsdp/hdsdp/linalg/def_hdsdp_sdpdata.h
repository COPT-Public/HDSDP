#ifndef def_hdsdp_sdpdata_h
#define def_hdsdp_sdpdata_h

#include "hdsdp.h"

/** @struct eigFactor
 *  @brief The eigen decomposition structure
 */
typedef struct {
    
    int     nCol;
    int     rank;
    double *eVals;
    double *eVecs;
    
} eig_factor;

/* Implementations of the SDP coefficient matrix */
typedef enum {
    
    SDP_COEFF_ZERO,
    SDP_COEFF_SPARSE,
    SDP_COEFF_DENSE,
    SDP_COEFF_SPR1,
    SDP_COEFF_DSR1
    
} sdp_coeff_type;

typedef struct {
    
    int        nSDPCol;
    
    sdp_coeff_type dataType;
    void      *dataMat;
    
    eig_factor *eig;
    
    void (*dataMataApB)          ( void *, double, void * );
    double (*dataMatDot)         ( void *, double * );
    void (*dataMatScal)          ( void *, double );
    void (*dataMatNorm)          ( void *, int );
    hdsdp_retcode (*dataMatEig)  ( void *, eig_factor ** );
    int (*dataMatGetNnz)         ( void * );
    void (*dataMatDump)          ( void *, double * );
    
} sdp_coeff;

typedef struct {
    
    int nSDPCol;
    
} sdp_coeff_zero;

typedef struct {
    
    int     nSDPCol;
    int     nTriMatElem;
    int    *triMatCol;
    int    *triMatRow;
    double *triMatElem;
    
} sdp_coeff_sparse;


typedef struct {
    
    int     nSDPCol;
    double *dsMatElem;
    
} sdp_coeff_dense;


typedef struct {
    
    int     nSDPCol;
    double  spR1FactorSign; ///< Include scale, may not equal to +1.0 or -1.0
    int     nSpR1FactorElem;
    int    *spR1MatIdx;
    double *spR1MatElem;
    
} sdp_coeff_spr1;


typedef struct {
    
    int     nSDPCol;
    double  r1FactorSign; ///< Include scale, may not equal to +1.0 or -1.0
    double *r1MatFactor;
    
} sdp_coeff_dsr1;


#endif /* def_hdsdp_sdpdata_h */
