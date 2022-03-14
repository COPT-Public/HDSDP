#ifndef dsdpsolver_h
#define dsdpsolver_h

#include "structs.h"
#include "dsdpdata.h"
#include "dsdpcg.h"
#include "dsdpparam.h"
#include "dsdpstats.h"
#include "symschur.h"

#define IterStep 20
#define nEvent   20

typedef struct {
    
    // Problem data
    sdpMat   **sdpData;   // SDP data A after transformation.
    DSDP_INT *isSDPset;
    
    lpMat    *lpData;     // LP data A after transformation
    vec      *lpObj;        // Linear coefficient
    DSDP_INT isLPset;
    
    vec      *dObj;       // Dual objective b
    
    // Dimension data
    DSDP_INT n;           // Total number of variables
    DSDP_INT m;           // Dimension of dual
    DSDP_INT nBlock;      // Number of SDP blocks ( = 1) currently
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
    
    vec      *asinv;      // Store the sum of trace(A * Sinv) for different A_i (including C),
                          // setup as a by product
    double   csinvcsinv;  // Store trace(C * Sinv * C * Sinv)
    double   csinv;       // Store csinv
    double   csinvrysinv; // Store trace(C * Sinv * Ry * Sinv)
    
    
    DSDPSchur *M;         // Advanced Schur matrix setup
    dsMat    *Msdp;       // Schur matrix for SDP (dense)
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
    double   Pnrm;
    double   dPotential;  // Value of dual potential function
    
    // Residuals
    double   rtk;         // Complementarity residual
    double    Ry;         // SDP Dual infeasibility
    vec      *ry;         // LP dual infeasibility
    
    // Step matrix
    spsMat   **dS;        // SDP step matrix
    
    /* TODO: Remove spaux by applying Lanczos iterations */
    spsMat   **spaux;     // Used for maximum step computation
    dsMat    **dsaux;     // Used for Schur matrix setup
    rkMat    **rkaux;     // Use for Schur matrix setup
    spsMat   **Scker;     // Used for checking positive definiteness
    vec      *ds;         // LP step matrix
    vec      *dy;         // Dual step matrix
    double   dtau;        // Tau step
    double   dkappa;      // Kappa step
    
    double   pinfeas;     // Primal infeasibility
    
    // Solver parameter
    hsdParam *param;      // Solver parameters
    
    // Solver status
    DSDP_INT insStatus;   // Solver instance status
    DSDP_INT solStatus;   // Solver solution status
    
    // Logging
    DSDP_INT verbosity;   // Solver information
    
    // Primal variable
    vec      *pScaler;    // Primal scaling coefficient
    vec      *ymaker;     // y
    vec      *dymaker;    // dy
    double   mumaker;     // mu
    
    // Solver statistics
    DSDPStats dsdpStats;  // Solver statistics
    
} HSDSolver;

typedef HSDSolver Solver;

#ifdef __cplusplus
extern "C" {
#endif
// Solver interface
extern DSDP_INT DSDPCreate( Solver **dsdpSolver );

extern DSDP_INT DSDPSetDim( Solver    *dsdpSolver,
                            DSDP_INT  nVars,
                            DSDP_INT  nBlock,
                            DSDP_INT  nConstrs,
                            DSDP_INT  lpDim,
                            DSDP_INT  *nNzs );

extern DSDP_INT DSDPSetLPData( Solver    *dsdpSolver,
                               DSDP_INT  nCol,
                               DSDP_INT  *Ap,
                               DSDP_INT  *Ai,
                               double    *Ax,
                               double    *lpObj );

extern DSDP_INT DSDPSetSDPConeData( Solver    *dsdpSolver,
                                    DSDP_INT  blockid,
                                    DSDP_INT  coneSize,
                                    DSDP_INT  *typehint,
                                    DSDP_INT  *Asdpp,
                                    DSDP_INT  *Asdpi,
                                    double    *Asdpx );

extern DSDP_INT DSDPSetObj   ( HSDSolver *dsdpSolver, double *dObj );
extern DSDP_INT DSDPOptimize ( Solver *dsdpSolver );

extern DSDP_INT DSDPSetDblParam ( HSDSolver *dsdpSolver, DSDP_INT pName, double   dblVal  );
extern DSDP_INT DSDPSetIntParam ( HSDSolver *dsdpSolver, DSDP_INT pName, DSDP_INT intVal  );
extern DSDP_INT DSDPGetDblParam ( HSDSolver *dsdpSolver, DSDP_INT pName, double   *dblVal );
extern DSDP_INT DSDPGetIntParam ( HSDSolver *dsdpSolver, DSDP_INT pName, DSDP_INT *intVal );

extern DSDP_INT DSDPDestroy  ( Solver *dsdpSolver );

extern void     DSDPPrintVersion (void);

#ifdef __cplusplus
}
#endif


#endif /* dsdpsolver_h */
