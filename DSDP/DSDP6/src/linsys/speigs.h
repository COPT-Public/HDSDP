#ifndef speigs_h
#define speigs_h

// Implement the sparse eigen-decomposition used in DSDP
#include "structs.h"

#define EIG_FACTOR_FEAST   1
#define EIG_FACTOR_LAPACK  2

#define LWORK   30
#define IWORK   12

#ifdef __cplusplus
extern "C" {
#endif

extern void     speigInit ( speigfac *eigfac );
extern DSDP_INT speigAlloc( speigfac *eigfac, DSDP_INT nmax );
extern DSDP_INT speigSfac ( speigfac *eigfac, spsMat *A, double *eigvals, double *eigvecs );
extern void     speigDfac ( speigfac *eigfac, dsMat  *A, double *eigvals, double *eigvecs );
extern void     speigFree ( speigfac *eigfac );

#ifdef __cplusplus
}
#endif

#endif /* speigs_h */
