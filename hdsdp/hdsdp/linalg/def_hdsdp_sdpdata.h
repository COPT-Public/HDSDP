#ifndef def_hdsdp_sdpdata_h
#define def_hdsdp_sdpdata_h

#include "interface/hdsdp.h"

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
    
    int eigRank;
    double *eigVals;
    double *eigVecs;
    
    hdsdp_retcode (*create) ( void **, int, int, int *, double * );
    double (*dot)         ( void *, double * );
    void (*scal)          ( void *, double );
    double (*norm)        ( void *, int );
    hdsdp_retcode (*eig)  ( void *, int *, double *, double **, double ** );
    int (*getNnz)         ( void * );
    void (*dump)          ( void *, double * );
    void (*destroy)       ( void ** );
    void (*view)          ( void * );
    
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
