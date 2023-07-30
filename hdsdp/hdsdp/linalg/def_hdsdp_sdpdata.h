#ifndef def_hdsdp_sdpdata_h
#define def_hdsdp_sdpdata_h

#ifdef HEADERPATH
#include "interface/hdsdp.h"
#include "linalg/hdsdp_linsolver.h"
#else
#include "hdsdp.h"
#include "hdsdp_linsolver.h"
#endif

/* Implementations of the SDP coefficient matrix
 In HDSDP, we implement five data structures for LP coefficient matrix
 
 1. Zero matrix
 2. Sparse matrix
 3. Dense matrix
 4. Sparse rank one matrix
 5. Denser rank one matrix
 
 Each type of coefficient will be asociated with a data type, its eigen-decomposition
 and several methods that make it work in the HDSDP solver.
 
 */
typedef enum {
    
    SDP_COEFF_ZERO,
    SDP_COEFF_SPARSE,
    SDP_COEFF_DENSE,
    SDP_COEFF_SPR1,
    SDP_COEFF_DSR1
    
} sdp_coeff_type;

struct sdp_coeff_s {
    
    int        nSDPCol;
    
    sdp_coeff_type dataType;
    void      *dataMat;
    
    int eigRank;
    double *eigVals;
    double *eigVecs;
    
    hdsdp_retcode (*create) ( void **, int, int, int *, double * );
    void (*scal)          ( void *, double );
    double (*norm)        ( void *, int );
    hdsdp_retcode (*eig)  ( void *, int *, double *, double **, double ** );
    int (*getnnz)         ( void * );
    void (*dump)          ( void *, double * );
    void (*getmatnz)      ( void *, int * );
    void (*add2buffer)    ( void *, double, int *, double *);
    void (*destroy)       ( void ** );
    void (*view)          ( void * );
    int  (*iseye)         ( void *, double * );
    int  (*isunitcol)     ( void *, int * );
    
    /* Schur complement operations */
    /* KKT strategy 2 */
    void   (*kkt2r1solve)   ( void *, hdsdp_linsys *, double *, double *, double * );
    double (*kkt2quadform)  ( void *, double *, double * );
    double (*kkt2r1asinv)   ( void *, double * );
    
    /* KKT strategy 3 */
    double (*kkt3sinvAsinv) ( void *, hdsdp_linsys *, double *, double *, double * );
    double (*kkt3AdotB)     ( void *, double *, double * );
    
    /* KKT strategy 4 */
    double (*kkt4Asinv)     ( void *, hdsdp_linsys *, double *, double *, double, double * );
    double (*kkt4AdotsinvB) ( void *, hdsdp_linsys *, double *, double *, double * );
    
    /* KKT strategy 5 */
    double (*kkt5AsinvBsinv) ( void *, struct sdp_coeff_s *, double *, double * );
    double (*kkt5SinvAdotSinv)  ( void *, hdsdp_linsys *, double *, double * );
};

typedef struct sdp_coeff_s sdp_coeff;

typedef struct {
    
    /* In a zero matrix. There is nothing but an integer recording matrix dimension */
    
    int nSDPCol;
    
} sdp_coeff_zero;

typedef struct {
    
    /* In a sparse matrix, we adopt the triplet format i, j, x */
    int     nSDPCol;
    int     nTriMatElem;
    int    *triMatCol;
    int    *triMatRow;
    double *triMatElem;
    
} sdp_coeff_sparse;


typedef struct {
    
    /* In a dense matrix, we store an n * (n + 1) / 2 array in packed format */
    int     nSDPCol;
    double *dsMatElem;
    
} sdp_coeff_dense;


typedef struct {
    
    /* In the sparse rank 1 structure, we store the number of
       SDP columns as well as the coefficient sign, so that
     
       A = s * a * a',
     
       where a is represented using a triplet format.
       Also a full vector is used to save the effort to expand the sparse vector
     
     */
    
    int     nSDPCol;
    double  spR1FactorSign; ///< Include scale, may not equal to +1.0 or -1.0
    int     nSpR1FactorElem;
    int    *spR1MatIdx;
    double *spR1MatElem;
    double *spR1MatFactor;
    
} sdp_coeff_spr1;


typedef struct {
    
    /* A dense rank 1 matrix is stored in a similar way to the sparse version
       There is no sparse representation of the factor
     */

    int     nSDPCol;
    double  r1FactorSign; ///< Include scale, may not equal to +1.0 or -1.0
    double *r1MatFactor;
    
} sdp_coeff_dsr1;


#endif /* def_hdsdp_sdpdata_h */
