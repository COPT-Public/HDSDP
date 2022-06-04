#ifndef speigs_h
#define speigs_h

// Implement the sparse eigen-decomposition used in DSDP
#include "dsdphsd.h"
#include "structs.h"

#define EIG_FACTOR_FEAST   1
#define EIG_FACTOR_LAPACK  2

#define LWORK   30
#define IWORK   12

typedef struct {
    
    DSDP_INT  nmax;    // Maximum possible dimension of matrix to factorize
    DSDP_INT  lwork;   // Length of working space
    DSDP_INT  liwork;  // Length of working space
    DSDP_INT  factorMethod;
    
    double   *dwork;   // Double working space
    double   *dworkmat;// Double working space
    double   *dworkevc;// Double working space
    double   *dworkevl;// Double working space
    
    DSDP_INT *perm;    // Permutation vector
    DSDP_INT *pinv;    // Inverse of the permutation
    DSDP_INT *colnnz;  // Column nnz counter
    DSDP_INT *iwork;   // Integer working space
    DSDP_INT *iworkup; // Another integer working space
    
} speigfac;

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
