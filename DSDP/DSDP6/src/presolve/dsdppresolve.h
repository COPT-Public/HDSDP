#ifndef dsdppresolve_h
#define dsdppresolve_h
/* Implement the presolving interface for DSDP. Currently the presolving procedure contains
 
 1. Coefficient normalization for LP (by column) and SDP (by block)
 2. rank-1 structure detection for dense matrices and sparse matrices with only 1 element
 
*/
#include <stdio.h>
#include "dsdphsd.h"
#include "sparsemat.h"
#include "densemat.h"
#include "rankonemat.h"
#include "dsdpdata.h"


#ifdef __cplusplus
extern "C" {
#endif


extern DSDP_INT preSDPMatPScale( sdpMat *dataMat, vec *pScaler );
extern DSDP_INT preSDPMatDScale( sdpMat *dataMat );


#ifdef __cplusplus
}
#endif


#endif /* dsdppresolve_h */
