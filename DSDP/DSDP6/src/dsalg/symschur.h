#ifndef symschur_h
#define symschur_h

/* DSDP Schur matrix setup: Advanced version incorporating permutation
   Only involving setup (no perturbation or solve)
 */

#include "dsdphsd.h"
#include "dsdpdata.h"
#include "structs.h"

#define SCHUR_M1 1
#define SCHUR_M2 2
#define SCHUR_M3 3
#define SCHUR_M4 4
#define SCHUR_M5 5

#define KAPPA (1.5)  // Ratio between sparse and dense memory access efficiency

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
    schurMat *M;      // Schur matrix
    
    // Auxiliary
    DSDP_INT *phaseA;  // Which phase DSDP is in
    vec      *asinv;
    vec      *asinvrysinv;
    vec      *asinvcsinv;
    double   *csinvrysinv;
    double   *csinvcsinv;
    double   *csinv;
    double   *rysinv;
    
    double    scaler;   // Scaling factor for S
    double   *schurAux;
    DSDP_INT *useTwo;   // Only two strategies used ?
    DSDP_INT *buildhsd; // Is HSD being used ?
    DSDP_INT *gold;     // Is Golden search being used ?
    DSDP_INT *lpSet;    // Is LP cone present ?
    
} symM;


#ifdef __cplusplus
extern "C" {
#endif

extern DSDP_INT symSchurMatInit    ( symM *M );
extern void symSchurMatSetDim  ( symM *M, DSDP_INT m, DSDP_INT nblock );
extern DSDP_INT symSchurMatAlloc   ( symM *M );
extern void symSchurMatRegister( symM *M, spsMat **S, dsMat **B, sdpMat **Adata, schurMat *Msdp,
                                     vec *asinv, vec *asinvrysinv, vec *asinvcsinv, double *csinvrysinv,
                                     double *csinv, double *csinvcsinv, double *rysinv, double *Ry,
                                     rkMat **rkaux, DSDP_INT *phaseA, DSDP_INT *buildHsd,
                                     DSDP_INT *gold, DSDP_INT *lpSet );
extern DSDP_INT symSchurMatFree    ( symM *M );
extern DSDP_INT DSDPSchurReorder   ( symM *M );
extern DSDP_INT DSDPCheckSchurType ( symM *M );
extern DSDP_INT DSDPSchurSetup     ( symM *M );
extern void asinvSetup         ( symM *M );
extern void arysinvSetup       ( symM *M );

#ifdef __cplusplus
}
#endif

#endif /* symschur_h */
