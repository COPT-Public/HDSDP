#ifndef dsdpdata_h
#define dsdpdata_h

/* Implement problem data interface of DSDP-HSD solver */
#include <stdio.h>
#include "dsdphsd.h"
#include "sparsemat.h"
#include "densemat.h"
#include "rankonemat.h"
#include "vec.h"
#include "cs.h"

#define MAT_TYPE_UNKNOWN 0
#define MAT_TYPE_SPARSE  1
#define MAT_TYPE_DENSE   2
#define MAT_TYPE_RANK1   3
#define MAT_TYPE_ZERO    4

/*
 SDP data stucture
 
 In the current version we accept CSC, Packed and Rank1 representations of SDP data
 and CSC representation of LP linear constraint matrix.
 
 The dual objective b is shared by the data and is supplied by setObjective method
*/

typedef struct {
    
    DSDP_INT  blockId;       // Index of the block
    DSDP_INT  dimy;          // Dimension of dual variable y
    DSDP_INT  dimS;          // Dimension of dual variable S
    
    DSDP_INT  nzeroMat;      // Number of zero matrices
    DSDP_INT  nspsMat;       // Number of sparse matrices
    DSDP_INT  *spsMatIdx;    // Index of the sparse matrices
    DSDP_INT  ndenseMat;     // Number of dense matrices
    DSDP_INT  *denseMatIdx;  // Index of the dense matrices
    DSDP_INT  nr1Mat;        // Number of rank 1 matrices
    DSDP_INT  *r1MatIdx;     // Index of rank1 matrices
    
    DSDP_INT  *types;        // Types of matrices
    void      **sdpData;     // Data of different types
    double    scaler;        // Scaler for presolving
    
} sdpMat;


/*
 LP data stucture
 
 We only support CSC for linear constraint matrix Al
*/

typedef struct {
    
    DSDP_INT dimy;       // Dimension of dual variable y
    DSDP_INT dims;       // Dimension of dual variable s
    
    // Sparse LP data in CSC format
    cs       *lpdata;
    double   *xscale;    // Scaler for presolving
    
} lpMat;


#ifdef __cplusplus
extern "C" {
#endif

extern DSDP_INT lpMatInit     ( lpMat  *lpData );
extern DSDP_INT lpMatSetDim   ( lpMat  *lpData, DSDP_INT dimy, DSDP_INT dims );
extern DSDP_INT lpMatSetData  ( lpMat *lpData, DSDP_INT *Ap, DSDP_INT *Ai, double *Ax );
extern DSDP_INT lpMataATy     ( double alpha, lpMat *lpData, vec *y, double *ATy );
extern DSDP_INT lpMatFree     ( lpMat *lpData );

extern DSDP_INT sdpMatInit    ( sdpMat *sdpData );
extern DSDP_INT sdpMatAlloc   ( sdpMat *sdpData );
extern DSDP_INT sdpMatSetDim  ( sdpMat *sdpData, DSDP_INT dimy, DSDP_INT dimS, DSDP_INT blockId );
extern DSDP_INT sdpMatSetHint ( sdpMat *sdpData, DSDP_INT *hint );
extern DSDP_INT sdpMatSetData ( sdpMat *sdpData, DSDP_INT *Ap, DSDP_INT *Ai, double *Ax );
extern DSDP_INT sdpMatFree    ( sdpMat *sdpData );

#ifdef __cplusplus
}
#endif


#endif /* dsdpdata_h */
