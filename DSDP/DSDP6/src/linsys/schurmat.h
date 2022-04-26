#ifndef schurmat_h
#define schurmat_h

/* Implement the Schur matrix related operations */
#include <stdio.h>
#include "dsdphsd.h"
#include "dsdplapack.h"
#include "vec.h"
#include "structs.h"

#define SCHUR_TYPE_UNKNOWN 0
#define SCHUR_TYPE_DENSE   1
#define SCHUR_TYPE_SPARSE  2

#ifdef __cplusplus
extern "C" {
#endif

extern void     schurMatInit       ( schurMat *sMat );
extern DSDP_INT schurMatAlloc      ( schurMat *sMat, DSDP_INT dim );
extern DSDP_INT schurMatselectType ( schurMat *sMat, DSDP_INT stype, DSDP_INT *nnzIdx );
extern void     schurMatFree       ( schurMat *sMat );
extern void     schurMatGetdiag    ( schurMat *sMat, vec *diag );
extern void     schurMatAdddiag    ( schurMat *sMat, double d );
extern void     schurMatFactorize  ( schurMat *sMat );
extern void     schurMatSolve      ( schurMat *sMat, DSDP_INT nrhs, double *B, double *aux );
extern void     schurMatMx         ( schurMat *sMat, vec *x, vec *Ax );
// extern void     schurMatPrint      ( schurMat *sMat );
extern void     schurMatReset      ( schurMat *sMat, DSDP_INT cleanM );

#ifdef __cplusplus
}
#endif



#endif /* schurmat_h */
