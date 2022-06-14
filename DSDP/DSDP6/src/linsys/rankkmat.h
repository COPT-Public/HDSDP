#ifndef rankkmat_h
#define rankkmat_h

/* Implement the rank-k matrix based on rank-one matrix */
#include <stdio.h>
#include "dsdphsd.h"
#include "vec.h"

#ifdef __cplusplus
extern "C" {
#endif

extern DSDP_INT rkMatInit               ( rkMat  *R                                );
extern DSDP_INT rkMatAllocIter          ( rkMat  *R,  DSDP_INT n                   );
extern DSDP_INT rkMatAllocAndSetData    ( rkMat  *R,  DSDP_INT n, DSDP_INT rank,
                                          double *eigvals, double *eigvecs         );
extern DSDP_INT rkMatAllocAndSelectData ( rkMat  *R,  DSDP_INT n, DSDP_INT rank,
                                          double thresh, double *eigvals,
                                          double *eigvecs                          );
extern void rkMatrkTrace            ( rkMat  *R1, rkMat *R2, double *trace     );
extern void rkMatdenseTrace         ( rkMat  *R,  dsMat  *A, double *trace     );
extern void rkMatdenseUpdate        ( dsMat  *dAMat, rkMat *rkBMat             );
extern void rkMatdiagTrace          ( rkMat  *R,  double diag, double *trace   );
extern DSDP_INT rkMatCountNnz           ( rkMat  *R                                );
extern DSDP_INT rkMatFree               ( rkMat  *R                                );
extern DSDP_INT rkMatFnorm              ( rkMat  *R,  double *fnrm                 );
extern DSDP_INT rkMatScale              ( rkMat  *R,  double a                     );
extern DSDP_INT rkMatRscale             ( rkMat  *R,  double r                     );
extern DSDP_INT rkMatisRank1            ( rkMat  *R,  DSDP_INT *isRank1            );
extern DSDP_INT rkMatGetRank            ( rkMat  *R                                );
extern r1Mat*   rkMatGetBase            ( rkMat  *R,  DSDP_INT i                   );
extern void rkMatCheckSparsity          ( rkMat *R, DSDP_INT *isdense, double thresh );
extern void rkMatGetSymbolic            ( rkMat *R, DSDP_INT *hash, DSDP_INT *nnzs );
extern DSDP_INT rkMatView               ( rkMat  *R                                );

#ifdef __cplusplus
}
#endif

#endif /* rankkmat_h */
