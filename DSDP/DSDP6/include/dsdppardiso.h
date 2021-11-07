#ifndef dsdppardiso_h
#define dsdppardiso_h

/*
 The code implements the linear system solver substitution for an initial version of DSDP6 solver
 Here the PARDISO solver is applied to replace the sparse solver routines in DSDP as well as the conjugate
 gradient solver
 
 July 8th, 2021
 By Gwz, Shanghai University of Finance and Economics
 */

#include "dsdphsd.h"

#define PARDISOINDEX 64

#define POSITIVEDEFINITE  2
#define INDEFINITE      (-2)

void pardisoinit ( void *, const DSDP_INT *, DSDP_INT * );

void pardiso     ( void     *, DSDP_INT    *, DSDP_INT *, DSDP_INT *, DSDP_INT *, DSDP_INT *,
                   double   *, DSDP_INT    *, DSDP_INT *, DSDP_INT *, DSDP_INT *, DSDP_INT *,
                   DSDP_INT *, double      *, double   *, DSDP_INT * );

void pardiso_getdiag( const void     * pt[64],
                      void           * df,
                      void           * da,
                      const DSDP_INT * mnum,
                      DSDP_INT       * error );

// Pardiso cholesky solver
static DSDP_INT PARDISO_PARAMS_CHOLESKY[PARDISOINDEX] = {
    
    1, /* Non-default value */ 3, /* P Nested dissection */ 0, /* Reserved          */
    0, /* No CG             */ 0, /* No user permitation */ 0, /* No overwriting    */
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
    1, /* No CG             */ 0, /* No user permitation */ 0, /* No overwriting    */
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

#endif
