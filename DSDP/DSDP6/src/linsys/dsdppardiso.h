#ifndef dsdppardiso_h
#define dsdppardiso_h

/* A wrapper for the pardiso linear system solver in terms of CSparse */
#include "dsdphsd.h"

#define PARDISO_OK      ( 0)
#define PARDISOINDEX    (64)      // Pardiso working array length
#define PARDISO_SYM_FAC (12)      // Symbolic analysis and factorization
#define PARDISO_SOLVE   (33)      // Solve linear system
#define PARDISO_FREE    (-1)      // Free internal data structure

// Pardiso default parameters
static DSDP_INT maxfct = 1; // Maximum number of factors
static DSDP_INT mnum   = 1; // The matrix used for the solution phase
static DSDP_INT mtype  = 2; // Real and symmetric positive definite
static DSDP_INT msglvl = 0; // Print no information
static DSDP_INT idummy = 0; // Dummy variable for taking up space

// Pardiso cholesky solver
// A good choice of iterative refinement parameter is 3
static DSDP_INT PARDISO_PARAMS_CHOLESKY[PARDISOINDEX] = {
    
    1, /* Non-default value */ 3, /* P Nested dissection */ 0, /* Reserved          */
    0, /* No CG             */ 0, /* No user permutation */ 0, /* No overwriting    */
    0, /* Refinement report */ 0, /* Two steps of ItRef  */ 0, /* Reserved          */
    100,/* NO perturb       */ 1, /* Disable scaling     */ 0, /* No transpose      */
    1, /* Disable matching  */ 0, /* Report on pivots    */ 0, /* Output            */
    0, /* Output            */ 0, /* Output              */-1, /* No report         */
    0, /* No report         */ 0, /* Output              */ 1, /* Pivoting          */
    0, /* nPosEigVals       */ 0, /* nNegEigVals         */ 0, /* Classic factorize */
    0,                         0,                           0,
    0,                         0,                           0,
    0,                         0,                           0,
    0,                         1, /* 0-based solve       */ 0,
    0,                         0,                           0,
    0,                         0,                           0,
    0,                         0,                           0,
    0,                         0,                           0,
    0,                         0,                           0,
    0,                         0,                           0,
    0,                         1, /* Getting diagonal    */ 0,
    0,                         0,                           0,
    0,                         0,                           0,
    0
};

// Pardiso CG solver
static DSDP_INT PARDISO_PARAMS_CG[PARDISOINDEX] = {
    
    1, /* Non-default value */ 3, /* P Nested dissection */ 0, /* Reserved          */
    1, /* CG                */ 0, /* No user permutation */ 0, /* No overwriting    */
    0, /* Refinement report */ 0, /* Two steps of ItRef  */ 0, /* Reserved          */
    100,/* NO perturb       */ 1, /* Disable scaling     */ 0, /* No transpose      */
    1, /* Disable matching  */ 0, /* Report on pivots    */ 0, /* Output            */
    0, /* Output            */ 0, /* Output              */-1, /* No report         */
    0, /* No report         */ 0, /* Output              */ 1, /* Pivoting          */
    0, /* nPosEigVals       */ 0, /* nNegEigVals         */ 0, /* Classic factorize */
    0,                         0,                           0,
    0,                         0,                           0,
    0,                         0,                           0,
    0,                         1, /* 0-based solve       */ 0,
    0,                         0,                           0,
    0,                         0,                           0,
    0,                         0,                           0,
    0,                         0,                           0,
    0,                         0,                           0,
    0,                         0,                           0,
    0,                         1, /* Getting diagonal    */ 0,
    0,                         0,                           0,
    0,                         0,                           0,
    0
};

extern void pardisoinit ( void *, const DSDP_INT *, DSDP_INT * );

extern void pardiso     ( void     *, DSDP_INT    *, DSDP_INT *, DSDP_INT *, DSDP_INT *, DSDP_INT *,
                          double   *, DSDP_INT    *, DSDP_INT *, DSDP_INT *, DSDP_INT *, DSDP_INT *,
                          DSDP_INT *, double      *, double   *, DSDP_INT * );

extern void pardiso_getdiag ( const void     * pt[64],
                              void           * df,
                              void           * da,
                              const DSDP_INT * mnum,
                              DSDP_INT       * error );

#endif
