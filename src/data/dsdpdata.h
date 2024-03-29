#ifndef dsdpdata_h
#define dsdpdata_h

/* Implement problem data interface of DSDP-HSD solver */
#include "structs.h"

#define MAT_TYPE_UNKNOWN 0
#define MAT_TYPE_SPARSE  1
#define MAT_TYPE_DENSE   2
#define MAT_TYPE_RANKK   3
#define MAT_TYPE_ZERO    4

#ifdef superDebug
#define denseThresh      (1.0)
#else
#define denseThresh      (0.6)
#endif

#define rankThresh       (0.1)

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
    DSDP_INT  nnzAmat;       // Number of nonzero A matrices
    DSDP_INT  nspsMat;       // Number of sparse matrices
    DSDP_INT  *spsMatIdx;    // Index of the sparse matrices
    DSDP_INT  ndenseMat;     // Number of dense matrices
    DSDP_INT  *denseMatIdx;  // Index of the dense matrices
    DSDP_INT  nrkMat;        // Number of rank 1 matrices
    DSDP_INT  *rkMatIdx;     // Index of rank1 matrices
    
    DSDP_INT  *types;        // Types of matrices
    void      **sdpData;     // Data of different types
    double    scaler;        // Scaler for presolving
    
    DSDP_INT  *nzIdx;        // Position of nonzero matrices
    DSDP_INT  *schurspIdx;   // Assembling order in the schur matrix if the block is sparse
    
} sdpMat;


/*
 LP data stucture
 
 We only support CSC for linear constraint matrix Al
*/

typedef struct {
    
    DSDP_INT dimy;       // Dimension of dual variable y
    DSDP_INT dims;       // Dimension of dual variable s
    
    // Sparse LP data in CSC format
    DSDP_INT *Ap;
    DSDP_INT *Ai;
    double   *Ax;
    DSDP_INT nnz;
    double   *xscale;    // Scaler for presolving
    
} lpMat;


#ifdef __cplusplus
extern "C" {
#endif

extern void lpMatInit        ( lpMat  *lpData );
extern void lpMatSetDim      ( lpMat  *lpData, DSDP_INT dimy, DSDP_INT dims );
extern DSDP_INT lpMatSetData ( lpMat *lpData, DSDP_INT *Ap, DSDP_INT *Ai, double *Ax );
extern void lpMatAx          ( lpMat *lpData, vec *x, vec *Ax );
extern double lpMatcx        ( vec *c, vec *x );
extern void lpMatFree        ( lpMat *lpData );
extern void lpMatView        ( lpMat *lpData );

extern void sdpMatInit           ( sdpMat *sdpData );
extern DSDP_INT sdpMatAlloc      ( sdpMat *sdpData );
extern void sdpMatSetDim         ( sdpMat *sdpData, DSDP_INT dimy, DSDP_INT dimS, DSDP_INT blockId );
extern DSDP_INT sdpMatSetData    ( sdpMat *sdpData, DSDP_INT *Ap, DSDP_INT *Ai, double *Ax, double *cnnz );
extern DSDP_INT sdpMatScatterNnz ( sdpMat *sdpData, DSDP_INT start, DSDP_INT col, DSDP_INT *colNnz );
extern void sdpMatSetSchurIndex  ( sdpMat *sdpData, DSDP_INT start, DSDP_INT col, DSDP_INT *csum, DSDP_INT ishift );

extern double sdpMatGetCFnorm   ( sdpMat *sdpData );
extern double sdpMatGetAOneNorm ( sdpMat *sdpData );
extern double sdpMatGetCOneNorm ( sdpMat *sdpData );
extern void   sdpMatRScaleC     ( sdpMat *sdpData, double r );
extern void   sdpMatAX          ( sdpMat *sdpData, dsMat *X, vec *AX );
extern double sdpMatCX          ( sdpMat *sdpData, dsMat *X );
extern void   sdpMatATy         ( sdpMat *sdpData, double ycoef, vec *y, double tau, spsMat *S, DSDP_INT *sumHash );

extern void sdpMatFree ( sdpMat *sdpData );

#ifdef __cplusplus
}
#endif


#endif /* dsdpdata_h */
