#ifndef dsdppardiso_h
#define dsdppardiso_h

/* A wrapper for the pardiso linear system solver in terms of CSparse */
#include "dsdphsd.h"

#define PARDISO_OK       (  0)
#define PARDISOINDEX     ( 64)      // Pardiso working array length
#define PARDISO_SYM      ( 11)      // Pardiso symbolic analysis
#define PARDISO_FAC      ( 22)      // Pardiso numerical factorization
#define PARDISO_SYM_FAC  ( 12)      // Symbolic analysis and factorization
#define PARDISO_SOLVE    ( 33)      // Solve linear system
#define PARDISO_FORWARD  (331)      // Pardiso forward solve
#define PARDISO_BACKWARD (333)      // Pardiso backward solve
#define PARDISO_FREE     ( -1)      // Free internal data structure

// Pardiso default parameters
static DSDP_INT maxfct = 1; // Maximum number of factors
static DSDP_INT mnum   = 1; // The matrix used for the solution phase
static DSDP_INT mtype  = 2; // Real and symmetric positive definite
static DSDP_INT msglvl = 0; // Print no information
static DSDP_INT idummy = 0; // Dummy variable for taking up space

// Pardiso cholesky solver
// A good choice of iterative refinement parameter is 3

#define SYMBOLIC 0
// TODO: Disable matrix checker iparam[26] when the debugging is done
static DSDP_INT PARDISO_PARAMS_CHOLESKY[PARDISOINDEX] = {
    
    1, /* Non-default value */ SYMBOLIC, /* P Nested dissection */ 0, /* Reserved          */
    0, /* No CG             */ 0, /* No user permutation */ 0, /* No overwriting    */
    0, /* Refinement report */ 0, /* Auto ItRef step     */ 0, /* Reserved          */
    6, /* Perturb           */ 0, /* Disable scaling     */ 0, /* No transpose      */
    0, /* Disable matching  */ 0, /* Report on pivots    */ 0, /* Output            */
    0, /* Output            */ 0, /* Output              */-1, /* No report         */
    0, /* No report         */ 0, /* Output              */ 1, /* Pivoting          */
    0, /* nPosEigVals       */ 0, /* nNegEigVals         */ 0, /* Classic factorize */
    0,                         0,                           0, /* Matrix checker    */
    0,                         0,                           0,
    0,                         0,                           0,
    0,                         1, /* 0-based solve       */ 0,
    0,                         0,                           0,
    0,                         0,                           0,
    0,                         0,                           0,
    0,                         0,                           0,
    0,                         0,                           0,
    0,                         0,                           0,
    0,                         1, /* Get diagonal        */ 0,
    0,                         0,                           0,
    0,                         0,                           0,
    0
};

static DSDP_INT PARDISO_PARAMS_PARTIAL_SOLVE[PARDISOINDEX] = {
    
    1, /* Non-default value */ SYMBOLIC, /* P Nested dissection */ 0, /* Reserved          */
    0, /* No CG             */ 0, /* No user permutation */ 0, /* No overwriting    */
    0, /* Refinement report */ 1, /* No ItRef step       */ 0, /* Reserved          */
    6, /* Perturb           */ 0, /* Disable scaling     */ 0, /* No transpose      */
    0, /* Disable matching  */ 0, /* Report on pivots    */ 0, /* Output            */
    0, /* Output            */ 0, /* Output              */-1, /* No report         */
    0, /* No report         */ 0, /* Output              */ 1, /* Pivoting          */
    0, /* nPosEigVals       */ 0, /* nNegEigVals         */ 0, /* Classic factorize */
    0,                         0,                           0, /* Matrix checker    */
    0,                         0,                           0,
    1, /* Partial solve     */ 0,                           0,
    0,                         1, /* 0-based solve       */ 0,
    0,                         0,                           0,
    0,                         0,                           0,
    0,                         0,                           0,
    0,                         0,                           0,
    0,                         0,                           0,
    0,                         0,                           0,
    0,                         1, /* Get diagonal        */ 0,
    0,                         0,                           0,
    0,                         0,                           0,
    0
};

static DSDP_INT PARDISO_PARAMS_PSD_CHECK[PARDISOINDEX] = {
    
    1, /* Non-default value */ SYMBOLIC, /* P Nested dissection */ 0, /* Reserved          */
    0, /* No CG             */ 0, /* No user permutation */ 0, /* No overwriting    */
    0, /* Refinement report */ 1, /* No ItRef step       */ 0, /* Reserved          */
    100, /* Perturb           */ 0, /* Disable scaling     */ 0, /* No transpose      */
    0, /* Disable matching  */ 0, /* Report on pivots    */ 0, /* Output            */
    0, /* Output            */ 0, /* Output              */-1, /* No report         */
    0, /* No report         */ 0, /* Output              */ 1, /* Pivoting          */
    0, /* nPosEigVals       */ 0, /* nNegEigVals         */ 0, /* Classic factorize */
    0,                         0,                           0, /* Matrix checker    */
    0,                         0,                           0,
    0, /* Partial solve     */ 0,                           0,
    0,                         1, /* 0-based solve       */ 0,
    0,                         0,                           0,
    0,                         0,                           0,
    0,                         0,                           0,
    0,                         0,                           0,
    0,                         0,                           0,
    0,                         0,                           0,
    0,                         1, /* Get diagonal        */ 0,
    0,                         0,                           0,
    0,                         0,                           0,
    0
};

static DSDP_INT PARDISO_PARAMS_FORWARD_BACKWORD[PARDISOINDEX] = {
    
    1, /* Non-default value */ SYMBOLIC, /* P Nested dissection */ 0, /* Reserved          */
    0, /* No CG             */ 0, /* No user permutation */ 1, /* Overwriting       */
    0, /* Refinement report */ 0, /* Auto ItRef step     */ 0, /* Reserved          */
    8, /* Perturb           */ 0, /* Disable scaling     */ 0, /* No transpose      */
    0, /* Disable matching  */ 0, /* Report on pivots    */ 0, /* Output            */
    0, /* Output            */ 0, /* Output              */-1, /* No report         */
    0, /* No report         */ 0, /* Output              */ 1, /* Pivoting          */
    0, /* nPosEigVals       */ 0, /* nNegEigVals         */ 0, /* Classic factorize */
    0,                         0,                           0, /* Matrix checker    */
    0,                         0,                           0,
    0,                         0,                           0,
    0,                         1, /* 0-based solve       */ 0,
    0,                         0,                           0,
    0,                         0,                           0,
    0,                         0,                           0,
    0,                         0,                           0,
    0,                         0,                           0,
    0,                         0,                           0,
    0,                         0, /* Get diagonal        */ 0,
    0,                         0,                           0,
    0,                         0,                           0,
    0
};

static DSDP_INT PARDISO_PARAMS_FORWARD_BACKWORD_LANCZOS[PARDISOINDEX] = {
    
    1, /* Non-default value */ SYMBOLIC, /* P Nested dissection */ 0, /* Reserved          */
    0, /* No CG             */ 0, /* No user permutation */ 0, /* Overwriting       */
    0, /* Refinement report */ 0, /* Auto ItRef step     */ 0, /* Reserved          */
    8, /* Perturb           */ 0, /* Disable scaling     */ 0, /* No transpose      */
    0, /* Disable matching  */ 0, /* Report on pivots    */ 0, /* Output            */
    0, /* Output            */ 0, /* Output              */-1, /* No report         */
    0, /* No report         */ 0, /* Output              */ 1, /* Pivoting          */
    0, /* nPosEigVals       */ 0, /* nNegEigVals         */ 0, /* Classic factorize */
    0,                         0,                           0, /* Matrix checker    */
    0,                         0,                           0,
    0,                         0,                           0,
    0,                         1, /* 0-based solve       */ 0,
    0,                         0,                           0,
    0,                         0,                           0,
    0,                         0,                           0,
    0,                         0,                           0,
    0,                         0,                           0,
    0,                         0,                           0,
    0,                         0, /* Get diagonal        */ 0,
    0,                         0,                           0,
    0,                         0,                           0,
    0
};


#ifdef DSDP64
#define pardiso pardiso_64
#endif

#ifdef __cplusplus
extern "C" {
#endif

extern void pardisoinit ( void *, const DSDP_INT *, DSDP_INT * );

extern void pardiso     ( void     *, DSDP_INT    *, DSDP_INT *, DSDP_INT *, DSDP_INT *, DSDP_INT *,
                          double   *, DSDP_INT    *, DSDP_INT *, DSDP_INT *, DSDP_INT *, DSDP_INT *,
                          DSDP_INT *, double      *, double   *, DSDP_INT * );

extern void pardiso_getdiag ( const void     * pt[PARDISOINDEX],
                              void           * df,
                              void           * da,
                              const DSDP_INT * mnum,
                              DSDP_INT       * error );

#ifdef __cplusplus
}
#endif
#endif
