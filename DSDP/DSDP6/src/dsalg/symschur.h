#ifndef symschur_h
#define symschur_h

/* DSDP Schur matrix setup
   Advanced version incorporating permutation
   Only involving setup (no perturbation or solve is available)
 */

#include "dsdphsd.h"
#include "dsdpdata.h"
#include "structs.h"

#define SCHUR_M1 1
#define SCHUR_M2 2
#define SCHUR_M3 3
#define SCHUR_M4 4
#define SCHUR_M5 5

typedef struct {
    
    DSDP_INT m;       // Dimension of the dual matrix
    DSDP_INT nblock;  // Number of blocks
    DSDP_INT Mready;  // Ready to setup M ?
    
    DSDP_INT **perms; // Permutation for Schur matrix setup
    DSDP_INT **MX;    // Use MX technique to compute the row from the permutation
    spsMat   **S;     // Dual matrices
    double   **Sinv;  // Inverse matrix
    dsMat    **B;     // Arrays to store inverse of dual matrices. Filled by dsaux
    rkMat    **rkaux; // Auxiliary rank-k matrix
    
    sdpMat   **Adata; // Data arrays
    double   *Ry;     // Dual infeasibility
    dsMat    *M;      // Schur matrix
    
    // Auxiliary
    vec      *asinv;
    vec      *asinvrysinv;
    vec      *asinvcsinv;
    double   *csinvrysinv;
    double   *csinvcsinv;
    double   *csinv;   
    
    DSDP_INT *useTwo;  // Only two strategies used ?
    
} DSDPSchur;

#define KAPPA (1.5)  // Ratio between sparse and dense memory access efficiency

#endif /* symschur_h */
