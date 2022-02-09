#ifndef dsdpcg_h
#define dsdpcg_h

#include "dsdphsd.h"
#include "structs.h"

typedef struct {
    
    dsMat    *M;
    vec      *r;
    vec      *rnew;
    vec      *d;
    vec      *Pinvr;
    vec      *Md;
    vec      *x;
    
    DSDP_INT pType;
    void     *preCond;
    
    double   tol;
    double   resinrm;
    DSDP_INT dim;
    DSDP_INT niter;
    DSDP_INT maxiter;
    DSDP_INT status;
    
} CGSolver;

#define CG_PRECOND_DIAG    101
#define CG_PRECOND_CHOL    102

#define CG_STATUS_MAXITER  103
#define CG_STATUS_SOLVED   104
#define CG_STATUS_FAILED   105
#define CG_STATUS_UNKNOWN  106

#ifdef __cplusplus
extern "C" {
#endif

extern DSDP_INT dsdpCGInit            ( CGSolver *cgSolver                                     );
extern DSDP_INT dsdpCGAlloc           ( CGSolver *cgSolver, DSDP_INT m                         );
extern DSDP_INT dsdpCGFree            ( CGSolver *cgSolver                                     );
extern DSDP_INT dsdpCGSetTol          ( CGSolver *cgSolver, double tol                         );
extern DSDP_INT dsdpCGSetMaxIter      ( CGSolver *cgSolver, DSDP_INT maxiter                   );
extern DSDP_INT dsdpCGSetM            ( CGSolver *cgSolver, dsMat *M                           );
extern DSDP_INT dsdpCGSetPreCond      ( CGSolver *cgSolver, DSDP_INT pType, void *preCond      );
extern DSDP_INT dsdpGetCGSolStatistic ( CGSolver *cgSolver, DSDP_INT *status, double *resinorm );
extern DSDP_INT dsdpCGSolve           ( CGSolver *cgSolver, vec *b, vec *x0                    );

#ifdef __cplusplus
}
#endif

#endif /* dsdpcg_h */
