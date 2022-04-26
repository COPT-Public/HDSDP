#ifndef dsdpcg_h
#define dsdpcg_h

#include "dsdphsd.h"
#include "structs.h"

typedef struct {
    
    schurMat *M;        // LHS data
    vec      *r;        // Residual
    vec      *rnew;     // Workspace array
    vec      *d;        // Workspace array
    vec      *Pinvr;    // Workspace array
    vec      *Md;       // Workspace array
    vec      *x;        // CG solution vector
    vec      *aux;      // CG auxiliary array
    
    DSDP_INT pType;     // Pre-conditioner type
    schurMat *cholPre;  // Cholesky pre-conditioner
    vec      *vecPre;   // Diagonal pre-conditioner
        
    double   tol;       // Relative tolerance of CG
    double   resinrm;   // Residual norm
    DSDP_INT dim;       // Dimension of linear system
    DSDP_INT niter;     // Number of iterations
    DSDP_INT maxiter;   // Maximum number of iterations
    DSDP_INT status;    // Solution status
    DSDP_INT reuse;     // Reuse Cholesky pre-conditioner
    DSDP_INT nused;     // # iterations current pre-conditioner is already used
    DSDP_INT nfailed;   // Number of non-successfull solves
    
} CGSolver;

#define CG_PRECOND_DIAG      101
#define CG_PRECOND_CHOL      102

#define CG_STATUS_SOLVED     103
#define CG_STATUS_MAXITER    104
#define CG_STATUS_FAILED     105
#define CG_STATUS_INDEFINITE 106
#define CG_STATUS_UNKNOWN    107

#ifdef __cplusplus
extern "C" {
#endif

extern DSDP_INT dsdpCGInit            ( CGSolver *cgSolver                                     );
extern DSDP_INT dsdpCGAlloc           ( CGSolver *cgSolver, DSDP_INT m                         );
extern DSDP_INT dsdpCGFree            ( CGSolver *cgSolver                                     );
extern DSDP_INT dsdpCGSetTol          ( CGSolver *cgSolver, double tol                         );
extern DSDP_INT dsdpCGSetMaxIter      ( CGSolver *cgSolver, DSDP_INT maxiter                   );
extern DSDP_INT dsdpCGSetM            ( CGSolver *cgSolver, schurMat *M                        );
extern DSDP_INT dsdpCGSetDPre         ( CGSolver *cgSolver, vec   *preCond                     );
extern DSDP_INT dsdpCGSetCholPre      ( CGSolver *cgSolver, schurMat *preCond                  );
extern DSDP_INT dsdpCGSetPType        ( CGSolver *cgSolver, DSDP_INT pType                     );
extern DSDP_INT dsdpCGprepareP        ( CGSolver *cgSolver                                     );
extern DSDP_INT dsdpCGGetStatus       ( CGSolver *cgSolver, DSDP_INT *status                   );
extern DSDP_INT dsdpCGGetSolStatistic ( CGSolver *cgSolver, DSDP_INT *status, double *resinorm );
extern DSDP_INT dsdpCGSetPreReuse     ( CGSolver *cgSolver, DSDP_INT reuse                     );
extern DSDP_INT dsdpCGStoreRHS        ( CGSolver *cgSolver, vec *bin                           );
extern DSDP_INT dsdpCGRestoreRHS      ( CGSolver *cgSolver, vec *bout                          );
extern void     dsdpCGStoreLdiag      ( CGSolver *cgSolver                                     );
extern void     dsdpCGRestoreLdiag    ( CGSolver *cgSolver                                     );
extern DSDP_INT dsdpCGSolve           ( CGSolver *cgSolver, vec *b, vec *x0                    );

#ifdef __cplusplus
}
#endif

#endif /* dsdpcg_h */
