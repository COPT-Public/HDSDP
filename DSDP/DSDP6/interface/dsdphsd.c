#include <stdio.h>
#include <stdlib.h>
#include "dsdphsd.h"
#include "dsdpdata.h"
#include "vec.h"
#include "sparsemat.h"
#include "densemat.h"
#include "rankonemat.h"
#include "dsdpparam.h"

static char *etype = "DSDP Interface";

/* The data structure of DSDP-HSD solver */
typedef struct {
    
    // Problem data
    sdpMat   *Asdp;       // SDP data
    spsMat   **AsdpData;  // SDP data A after transformation
    spsMat   *CsdpData;   // SDP data C after transformation
    lpMat    *Alp;        // LP data
    vec      *dObj;       // Dual objective b
    
    // Dimension data
    DSDP_INT m;           // Dimension of dual
    DSDP_INT nBlock;      // Number of SDP blocks ( = 1) currently
    
    
    // Iterator
    spsMat   *S;          // SDP dual iteration matrix
    dsMat    *X;          // Primal solution matrix
    vec      *s;          // LP dual iteration matrix
    vec      *x;          // LP primal iteration matrix
    
    dsMat    *SinvASinv;  // Store Sinv * A * Sinv for different A_i (including C)
    vec      *asinv;      // Store trace(A * Sinv) for different A_i (including C),
                          // setup as a by product
    
    dsMat    *M;          // Schur matrix
    vec      *u;          // A_ls^(-2)c + AS^(-1)CS^(-1)
    vec      *d1;         // Temporary iteration array
    vec      *d2;         // Temporary iteration array
    
    vec      *yp;         // Full iteration for y
    double   taup;        // Full iteration for tau
    
    double   alpha;       // Stepsize
    
    // Step matrix
    spsMat   *dS;         // SDP step matrix
    vec      *ds;         // LP step matrix
    double   *dkappa;     // Kappa step
    
    // Solver parameter
    hsdParam *param;      // Solver parameters
    
} HSDSolver;


/* DSDP internal methods */
static DSDP_INT DSDPIInit( HSDSolver *dsdpSolver ) {
    
    // Allocate memory for the internal solver
    DSDP_INT retcode = DSDP_RETCODE_OK;
    dsdpSolver->param = &defaultParam;
    return retcode;
}

static DSDP_INT DSDPIAllocate( HSDSolver       *dsdpSolver,
                                const DSDP_INT  ndualVars,
                                const DSDP_INT  nlpConstrs,
                                const DSDP_INT  nsdpConstrsDs,
                                const DSDP_INT  nsdpConstrSps,
                                const DSDP_INT  nsdpConstrRk1 ) {
    
    // Allocate memory for the internal solver
    DSDP_INT retcode = DSDP_RETCODE_OK;
    
    // Setup problem data
    return retcode;
}

static DSDP_INT DSDPIPresolve( HSDSolver *dsdpSolver ) {
    
    // Do presolve
    DSDP_INT retcode = DSDP_RETCODE_OK;
    
    
    return retcode;
}

static DSDP_INT DSDPIPostsolve( HSDSolver *dsdpSolver ) {
    
    DSDP_INT retcode = DSDP_RETCODE_OK;
    // Post-solver
    
    return retcode;
}


extern DSDP_INT DSDPCreate( HSDSolver **dsdpSolver ) {
    
    /* Create solver */
    DSDP_INT retcode = DSDP_RETCODE_OK;
    HSDSolver *solver = NULL;
    solver = (HSDSolver *) calloc(1, sizeof(HSDSolver));
    retcode = DSDPIInit(solver);
    return retcode;
}




