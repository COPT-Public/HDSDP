#ifndef dsdpdata_h
#define dsdpdata_h

/* Implement problem data interface of DSDP-HSD solver */

#include <stdio.h>
#include "dsdphsd.h"
#include "vec.h"

/*
 SDP data stucture
 
 In the current version we accept CSC, Packed and Rank1 representations of SDP data
 and CSC representation of LP linear constraint matrix.
 
 The dual objective b is shared by the data and is supplied by setObjective method
*/

typedef struct {
    
    DSDP_INT dimy;       // Dimension of dual variable y
    DSDP_INT dimS;       // Dimension of dual variable S
    
    // Sparse SDP data transpose( A ) in CSC format
    DSDP_INT  nSpMat;    // Number of sparse matrices
    DSDP_INT *ASDPspBeg; // Beginning column index
    DSDP_INT *ASDPspIdx; // Row index
    double   *ASDPspVal; // Values
    
    // Dense SDP data A
    DSDP_INT  nDsMat;    // Number of dense matrices
    double   *ASDPds;    // Array representing dense matrices
    
    // Rank 1 SDP data A
    DSDP_INT  nR1Mat;    // Number of rank-one matrices
    double   *ASDPr1;    // Rank 1 data
    
} sdpMat;


/*
 LP data stucture
 
 We only support CSC for linear constraint matrix Al
*/

typedef struct {
    
    DSDP_INT dimy;       // Dimension of dual variable y
    DSDP_INT dims;       // Dimension of dual variable s
    
    // Sparse LP data in CSC format
    DSDP_INT *ALPBeg;
    DSDP_INT *ALPIdx;
    double   *ALPVal;
    
} lpMat;

#endif /* dsdpdata_h */
