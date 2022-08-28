#ifndef dsdpsolver_h
#define dsdpsolver_h

#include <string.h>
#include "structs.h"
#include "dsdpdata.h"
#include "dsdpparam.h"
#include "dsdpstats.h"
#include "symschur.h"

#define IterStep 20
#define nEvent   20

#define LIKELY(x)   __builtin_expect(!!(x), 1)
#define UNLIKELY(x) __builtin_expect(!!(x), 0)

struct hdsdp {
    
    // Model name
    char SDPModel[100];   // Name of SDP model
    
    // Problem data
    sdpMat   **sdpData;   // SDP data A after transformation.
    DSDP_INT *isSDPset;
    
    lpMat    *lpData;     // LP data A after transformation
    vec      *lpObj;        // Linear coefficient
    DSDP_INT isLPset;
    
    vec      *dObj;       // Dual objective b
    
    // Dimension data
    DSDP_INT n;           // Total number of variables (SDP and LP)
    DSDP_INT m;           // Dimension of dual
    DSDP_INT nall;        // Total dimension of cones
    DSDP_INT nBlock;      // Number of SDP blocks
    DSDP_INT lpDim;       // Dimension of LP
    
    // Iteration Monitor
    DSDP_INT iterProgress[IterStep];
    DSDP_INT eventMonitor[nEvent];
    
    // Iterator
    double  pObjVal;      // Primal objective
    double  dObjVal;      // Dual objective
    double  mu;           // Duality gap
    
    spsMat   **S;         // SDP dual iteration matrix
    DSDP_INT **symS;      // Symbolic structure of S matrices
    vec      *s;          // LP dual iteration matrix
    vec      *x;          // LP primal iteration vector
    
    vec      *asinv;      // Store the sum of trace(A * Sinv) for different A_i (excluding C),
    vec      *asinvrysinv;
    double   csinvcsinv;  // Store trace(C * Sinv * C * Sinv)
    double   csinv;       // Store csinv
    double   csinvrysinv; // Store trace(C * Sinv * Ry * Sinv)
    
    symM     *M;          // Advanced Schur matrix setup
    schurMat *Msdp;       // Schur matrix for SDP (dense or sparse)
    CGSolver *cgSolver;   // Internal CG solvers
    vec      *Mdiag;      // Diagonal elements of Schur matrix for pre-conditioning
    double   Mscaler;     // Scaler of Schur complement
    
    vec      *u;          // A_ls^(-2)c + AS^(-1)CS^(-1)
    vec      *b1;         // Auxiliary array 1
    vec      *b2;         // Auxiliary array 2
    vec      *d1;         // Embedding
    vec      *d12;        // Auxiliary for d1 
    vec      *d2;         // Affine scaling
    vec      *d3;         // Dual centering
    vec      *d4;         // Dual infeasibility
    
    vec      *y;          // y
    
    double   tau;         // tau
    double   kappa;       // kappa
    double   alpha;       // Stepsize
    
    // Proximity measure
    double   ybound;      // Bound of dual variable
    double   dperturb;    // Dual perturbation
    double   Pnrm;
    double   dPotential;  // Value of dual potential function
    
    // Residuals
    double    rysinv;     // Complementarity residual
    double    Ry;         // SDP Dual infeasibility
    double    drate;      // Rate for eliminating dual infeasibility
    
    // Step matrix
    spsMat   **dS;        // SDP step matrix
    
    lczstep **lczs;   // Lanczos solver
    dsMat    **dsaux;     // Used for Schur matrix setup
    rkMat    **rkaux;     // Use for Schur matrix setup
    spsMat   **Scker;     // Used for checking positive definiteness
    vec      *scker;      // Used as a buffer
    vec      *ds;         // LP step matrix
    vec      *dy;         // Dual step matrix
    double   dtau;        // Tau step
    double   dkappa;      // Kappa step
    
    vec      *sl;         // Dual slack y - l
    vec      *slcker;
    vec      *su;         // Dual slack u - y
    vec      *sucker;
    double   pinfeas;     // Primal infeasibility
    
    // Solver parameter
    dsdpparam *param;      // Solver parameters
    
    // Solver status
    DSDP_INT insStatus;   // Solver instance status
    DSDP_INT solStatus;   // Solver solution status
    
    // Logging
    DSDP_INT verbosity;   // Solver information
    
    // Primal variable
    vec      *pScaler;    // Primal scaling coefficient
    double    cScaler;    // C is scaled by
    vec      *ymaker;     // y
    vec      *dymaker;    // dy
    double    mumaker;    // mu
    vec      *ymaker2;
    vec      *dymaker2;
    double    mumaker2;
    
    // Solver statistics
    DSDPStats dsdpStats;  // Solver statistics
    double    startTime;  // Solver start time
    
};

typedef struct hdsdp HSDSolver;


#endif /* dsdpsolver_h */
