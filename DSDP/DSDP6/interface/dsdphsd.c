#include <stdio.h>
#include <stdlib.h>
#include "dsdphsd.h"
#include "dsdpdata.h"
#include "vec.h"
#include "sparsemat.h"
#include "densemat.h"
#include "rankonemat.h"
#include "dsdppresolve.h"
#include "dsdpparam.h"

static char *etype = "DSDP Interface";

#define IterStep 20

/* The data structure of DSDP-HSD solver */
typedef struct {
    
    // Problem data
    sdpMat   **sdpData;   // SDP data A after transformation.
    DSDP_INT *isSDPset;
    
    lpMat    *lpData;     // LP data A after transformation
    vec      *clp;        // Linear coefficient
    DSDP_INT isLPset;
    
    vec      *dObj;       // Dual objective b
    
    // Dimension data
    DSDP_INT m;           // Dimension of dual
    DSDP_INT nBlock;      // Number of SDP blocks ( = 1) currently
    DSDP_INT lpDim;       // Dimension of LP
    
    // Iteration Monitor
    DSDP_INT iterProgress[IterStep];
    
    // Iterator
    spsMat   **S;         // SDP dual iteration matrix
    vec      *s;          // LP dual iteration matrix
    vec      *x;          // LP primal iteration matrix
    
    dsMat    **SinvASinv; // Store Sinv * A * Sinv for different A_i
    vec      **asinv;     // Store trace(A * Sinv) for different A_i (including C),
                          // setup as a by product
    
    dsMat    *M;          // Schur matrix
    vec      *u;          // A_ls^(-2)c + AS^(-1)CS^(-1)
    vec      *d1;         // Temporary iteration array
    vec      *d2;         // Temporary iteration array
    
    vec      *yp;         // Full iteration for y
    double   taup;        // Full iteration for tau
    
    double   alpha;       // Stepsize
    
    // Step matrix
    spsMat   **dS;        // SDP step matrix
    vec      *ds;         // LP step matrix
    double   dtau;        // Tau step
    double   dkappa;      // Kappa step
    
    // Solver parameter
    hsdParam *param;      // Solver parameters
    
    // Solver status
    DSDP_INT insStatus;   // Solver instance status
    DSDP_INT solStatus;   // Solver solution status
    
    // Logging
    DSDP_INT verbosity;   // Solver information
    
    // Primal variable
    vec      *pScaler;    // Primal scaling coefficient
    dsMat    **X;         // Primal solution matrix
    DSDP_INT isXComputed; // Whether X has been computed
    
    
} HSDSolver;


/* DSDP internal methods */
static DSDP_INT DSDPIInit( HSDSolver *dsdpSolver ) {
    
    // Allocate memory for the internal solver
    DSDP_INT retcode = DSDP_RETCODE_OK;
    assert( dsdpSolver->insStatus == DSDP_STATUS_UNINIT );
    
    if (dsdpSolver->insStatus != DSDP_STATUS_UNINIT) {
        error(etype, "Instance has been initialized. \n");
    }
    
    // Problem data
    dsdpSolver->sdpData = NULL;
    dsdpSolver->clp     = NULL;
    dsdpSolver->lpData  = NULL;
    dsdpSolver->dObj    = NULL;
    
    dsdpSolver->isLPset  = FALSE;
    dsdpSolver->isSDPset = NULL;
    
    // Dimension data
    dsdpSolver->m      = 0;
    dsdpSolver->nBlock = 0;
    dsdpSolver->lpDim  = 0;
    
    // IterProgress monitor
    memset(dsdpSolver->iterProgress, 0,
           sizeof(DSDP_INT) * IterStep);
    
    // Iterator
    dsdpSolver->S = NULL;
    dsdpSolver->s = NULL;
    dsdpSolver->x = NULL;
    
    dsdpSolver->SinvASinv = NULL;
    dsdpSolver->asinv     = NULL;
    
    dsdpSolver->M      = NULL;
    dsdpSolver->u      = NULL;
    dsdpSolver->d1     = NULL;
    dsdpSolver->d2     = NULL;
    
    dsdpSolver->yp     = NULL;
    dsdpSolver->taup   = 0.0;
    dsdpSolver->alpha  = 0.0;
    
    // Step matrix
    dsdpSolver->dS     = NULL;
    dsdpSolver->ds     = NULL;
    dsdpSolver->dkappa = 0.0;
    dsdpSolver->dtau   = 0.0;
    
    dsdpSolver->param     = &defaultParam;
    dsdpSolver->insStatus = DSDP_STATUS_INIT_UNSET;
    dsdpSolver->solStatus = DSDP_UNKNOWN;
    
    // Verbosity
    dsdpSolver->verbosity = 0;
    
    // Primal variable
    dsdpSolver->pScaler     = NULL;
    dsdpSolver->X           = NULL;
    dsdpSolver->isXComputed = 0;
    
    return retcode;
}

static DSDP_INT DSDPIAlloc( HSDSolver *dsdpSolver ) {
    
    // Allocate memory for the internal solver (level 1)
    // Level 1 allocation only involves pointer/indicator arrays and the rest of memory will be
    // allocated when setting the problem data
    DSDP_INT retcode = DSDP_RETCODE_OK;
    assert( dsdpSolver->insStatus == DSDP_STATUS_INIT_UNSET );
    
    if (dsdpSolver->insStatus != DSDP_STATUS_INIT_UNSET) {
        error(etype, "Level 1 memory cannot be allocated. \n");
        retcode = DSDP_RETCODE_FAILED;
        return retcode;
    }
    
    DSDP_INT nblock = dsdpSolver->nBlock;
    
    dsdpSolver->sdpData   = (sdpMat **)  calloc(nblock, sizeof(sdpMat *));
    dsdpSolver->lpData    = (lpMat   *)  calloc(1,      sizeof(lpMat   ));
    dsdpSolver->clp       = (vec     *)  calloc(1,      sizeof(vec     ));
    dsdpSolver->S         = (spsMat **)  calloc(nblock, sizeof(spsMat *));
    dsdpSolver->dS        = (spsMat **)  calloc(nblock, sizeof(spsMat *));
    dsdpSolver->SinvASinv = (dsMat  **)  calloc(nblock, sizeof(dsMat  *));
    dsdpSolver->asinv     = (vec    **)  calloc(nblock, sizeof(vec    *));
    dsdpSolver->isSDPset  = (DSDP_INT *) calloc(nblock, sizeof(DSDP_INT));
    dsdpSolver->pScaler   = (vec     *)  calloc(1,      sizeof(vec     ));
    
    for (DSDP_INT i = 0; i < nblock; ++i) {
        retcode = sdpMatInit(dsdpSolver->sdpData[i]); checkCode;
    }
    
    if (dsdpSolver->verbosity) {
        printf("Level 1 memory set. \n");
    }
    
    return retcode;
}

static DSDP_INT DSDPICheckData( HSDSolver *dsdpSolver ) {
    
    // Check whether problem data is alreay set up
    DSDP_INT retcode = DSDP_RETCODE_OK;
    DSDP_INT nblock = dsdpSolver->nBlock;
    DSDP_INT nblockSet = 0;
    for (int i = 0; i < nblock; ++i) {
        nblockSet += dsdpSolver->isSDPset[i];
    }
    if (dsdpSolver->isLPset && nblockSet == nblock) {
        dsdpSolver->insStatus = DSDP_STATUS_SET;
    }
    
    if (dsdpSolver->verbosity) {
        printf(ID" out of "ID" blocks are set. \n", nblockSet, nblock);
        if (dsdpSolver->isLPset) {
            printf("LP data is set. \n");
        } else {
            printf("LP data is not set. \n");
        }
    }
    
    return retcode;
}

static DSDP_INT DSDPIFreeLPData ( HSDSolver *dsdpSolver ) {
    
    // Free the internal LP data
    DSDP_INT retcode = DSDP_RETCODE_OK;
    
    if (dsdpSolver->isLPset) {
        retcode = lpMatFree(dsdpSolver->lpData); checkCode;
        dsdpSolver->isLPset = 0;
    }
    
    DSDP_FREE(dsdpSolver->lpData);
    
    return retcode;
}

static DSDP_INT DSDPIFreeSDPData( HSDSolver *dsdpSolver ) {
    
    // Free the internal SDP data. Not responsible for isSDPSet
    DSDP_INT retcode = DSDP_RETCODE_OK;
    
    for (DSDP_INT i = 0; i < dsdpSolver->nBlock + 1; ++i) {
        if (dsdpSolver->isSDPset[i]) {
            retcode = sdpMatFree(dsdpSolver->sdpData[i]); checkCode;
            DSDP_FREE(dsdpSolver->sdpData[i]);
        }
    }
    
    DSDP_FREE(dsdpSolver->sdpData);
    
    return retcode;
}

static DSDP_INT DSDPIFreeAlgIter( HSDSolver *dsdpSolver ) {
    
    // Free the internal iteration data
    DSDP_INT retcode = DSDP_RETCODE_OK;
    
    if (dsdpSolver->insStatus != DSDP_STATUS_SOLVED) {
        return retcode;
    }
    
    DSDP_INT nblock = dsdpSolver->nBlock;
    // Free the internal iteration data
    
    // S
    for (DSDP_INT i = 0; i < nblock; ++i) {
        retcode = spsMatFree(dsdpSolver->S[i]); checkCode;
        DSDP_FREE(dsdpSolver->S[i]);
    }
    
    DSDP_FREE(dsdpSolver->S);
    
    // s
    retcode = vec_free(dsdpSolver->s); checkCode;
    DSDP_FREE(dsdpSolver->s);
    
    // x
    retcode = vec_free(dsdpSolver->x); checkCode;
    DSDP_FREE(dsdpSolver->x);
    
    // SinvASinv
    for (DSDP_INT i = 0; i < nblock; ++i) {
        retcode = denseMatFree(dsdpSolver->SinvASinv[i]); checkCode;
        DSDP_FREE(dsdpSolver->SinvASinv[i]);
    }
    
    DSDP_FREE(dsdpSolver->SinvASinv);
    
    // asinv
    for (DSDP_INT i = 0; i < nblock; ++i) {
        retcode = vec_free(dsdpSolver->asinv[i]); checkCode;
        DSDP_FREE(dsdpSolver->asinv[i]);
    }
    
    DSDP_FREE(dsdpSolver->asinv);
    
    // M
    retcode = denseMatFree(dsdpSolver->M);
    DSDP_FREE(dsdpSolver->M);
    
    // u, d1, d2, yp
    retcode = vec_free(dsdpSolver->u ); checkCode;
    retcode = vec_free(dsdpSolver->d1); checkCode;
    retcode = vec_free(dsdpSolver->d2); checkCode;
    retcode = vec_free(dsdpSolver->yp); checkCode;
    
    DSDP_FREE(dsdpSolver->u );
    DSDP_FREE(dsdpSolver->d1);
    DSDP_FREE(dsdpSolver->d2);
    DSDP_FREE(dsdpSolver->yp);

    // dS
    for (DSDP_INT i = 0; i < nblock; ++i) {
        retcode = spsMatFree(dsdpSolver->dS[i]); checkCode;
        DSDP_FREE(dsdpSolver->dS[i]);
    }
    
    DSDP_FREE(dsdpSolver->dS);
    
    // ds
    retcode = vec_free(dsdpSolver->ds);
    DSDP_FREE(dsdpSolver->ds);
    
    // pScaler
    retcode = vec_free(dsdpSolver->pScaler);
    DSDP_FREE(dsdpSolver->pScaler);
    
    // X
    if (dsdpSolver->isXComputed) {
        for (DSDP_INT i = 0; i < nblock; ++i) {
            retcode = denseMatFree(dsdpSolver->X[i]); checkCode;
            DSDP_FREE(dsdpSolver->X[i]);
        }
        DSDP_FREE(dsdpSolver->X);
    }
    
    return retcode;
}

static DSDP_INT DSDPIFreeCleanUp( HSDSolver *dsdpSolver ) {
    
    // Free the internal indicator arrays and some common data
    DSDP_INT retcode = DSDP_RETCODE_OK;
    
    // isSDPset
    DSDP_FREE(dsdpSolver->isSDPset);
    
    // clp
    if (dsdpSolver->isLPset) {
        retcode = vec_free(dsdpSolver->clp); checkCode;
    }
    
    // dObj
    retcode = vec_free(dsdpSolver->dObj); checkCode;
    DSDP_FREE(dsdpSolver->dObj);
    
    // Other data
    dsdpSolver->m      = 0;
    dsdpSolver->nBlock = 0;
    dsdpSolver->lpDim  = 0;
    dsdpSolver->taup   = 0.0;
    dsdpSolver->alpha  = 0.0;
    dsdpSolver->dtau   = 0.0;
    dsdpSolver->dkappa = 0.0;
    
    return retcode;
}

static DSDP_INT DSDPIPresolve( HSDSolver *dsdpSolver ) {
    
    // Do presolve
    DSDP_INT retcode = DSDP_RETCODE_OK;
    
    assert( dsdpSolver->insStatus == DSDP_STATUS_SET );
    if ( dsdpSolver->insStatus != DSDP_STATUS_SET ) {
        error(etype, "Problem data is not set up. \n");
    }
    
    // Round 1: scale the primal pairs {b_i, A_ip} across different blocks
    retcode = vec_init(dsdpSolver->pScaler); checkCode;
    retcode = vec_alloc(dsdpSolver->pScaler, dsdpSolver->m); checkCode;
    
    double maxNrm     = 0.0;
    double minNrm     = 0.0;
    double tmpnrm     = 0.0;
    double pScalFact  = 0.0;
    DSDP_INT coneSize = 0;
    sdpMat *cone      = NULL;
    
    
    for (DSDP_INT i = 0; i < dsdpSolver->m; ++i) {
        
        maxNrm = 0.0;
        minNrm = 0.0;
        
        for (DSDP_INT j = 0; j < dsdpSolver->nBlock; ++j) {
            
            cone = dsdpSolver->sdpData[j];
            coneSize = cone->dimS;
            
            switch (cone->types[i]) {
                case MAT_TYPE_ZERO:
                    break;
                case MAT_TYPE_DENSE:
                    retcode = denseMatFnorm((dsMat *) cone->sdpData[i], &tmpnrm);
                    checkCode;
                    break;
                case MAT_TYPE_SPARSE:
                    retcode = spsMatFnorm((spsMat *) cone->sdpData[i] , &tmpnrm);
                    checkCode;
                    break;
                case MAT_TYPE_RANK1:
                    retcode = r1MatFnorm((r1Mat *) cone->sdpData[i], &tmpnrm);
                    checkCode;
                    break;
                default:
                    error(etype, "Unknown matrix type. \n");
                    break;
            }
            
            maxNrm = MAX(maxNrm, tmpnrm);
            minNrm = MIN(minNrm, tmpnrm);
        }
        
        pScalFact = maxNrm * minNrm;
        
        if (pScalFact < 3 && pScalFact > 0.8) {
            pScalFact = 1.0;
        } else {
            dsdpSolver->pScaler->x[i] = sqrt(pScalFact);
        }
    }
    
    // Do presolving
    DSDP_INT nblock = dsdpSolver->nBlock;
    for (DSDP_INT j = 0; j < nblock; ++j) {
        cone = dsdpSolver->sdpData[j];
        // Rank-1 detection
        
        
        // Primal coefficient scaling
        retcode = preSDPMatPScale(cone, dsdpSolver->pScaler); checkCode;
        // Dual coefficient scaling
        retcode = preSDPMatDScale(cone); checkCode;
    }
    
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

extern DSDP_INT DSDPSetDim( HSDSolver *dsdpSolver,
                            DSDP_INT  nBlock,
                            DSDP_INT  nConstrs,
                            DSDP_INT  lpDim ) {
    
    /* Set dimension of the DSDP problem instance
       
       nBlock   is the number of SDP varaibles participating in the instance
       nConstrs is the dimension of the dual variable
       lpDim    is the dimension of the LP
     
    */
    
    DSDP_INT retcode = DSDP_RETCODE_OK;
    assert( dsdpSolver->insStatus == DSDP_STATUS_INIT_UNSET );
    if (dsdpSolver->insStatus != DSDP_STATUS_INIT_UNSET) {
        error(etype, "Instance not yet initialized or "
              "dimension is already set. \n");
        retcode = DSDP_RETCODE_FAILED;
        return retcode;
    }
    
    if ((nBlock + lpDim) <= 0 || nConstrs < 0 || lpDim < 0 || nBlock < 0) {
        error(etype, "Invalid dimension. \n");
    }
    
    if (dsdpSolver->verbosity) {
        printf("Dimension is successfully set. \n");
        printf("nSDPCones: "ID" "
               "nConstraints: "ID" "
               "LPDim: "ID". \n", nBlock, nConstrs, lpDim);
    }
    
    dsdpSolver->nBlock = nBlock;
    dsdpSolver->m      = nConstrs;
    dsdpSolver->lpDim  = lpDim;
    
    retcode = DSDPIAlloc(dsdpSolver); checkCode;
    
    return retcode;
}

extern DSDP_INT DSDPSetLPData( HSDSolver *dsdpSolver,
                               DSDP_INT  nCol,
                               DSDP_INT  *Ap,
                               DSDP_INT  *Ai,
                               double    *Ax,
                               double    *lpObj ) {
    
    /*
     LP data interface for the user. DSDP accepts the
     CSC representaion of coefficient A
    */
    
    DSDP_INT retcode = DSDP_RETCODE_OK;
    assert( dsdpSolver->insStatus == DSDP_STATUS_INIT_UNSET );
    assert( nCol > 0 );
    assert( dsdpSolver->m > 0);
    
    if (dsdpSolver->insStatus != DSDP_STATUS_INIT_UNSET) {
        error(etype, "The solver instance is either not initialized or "
              "already set. \n");
    } else if (dsdpSolver->m <= 0) {
        error(etype, "Instance dimension is not set. \n");
    } else if (dsdpSolver->isLPset) {
        error(etype, "LP data is already set. \n");
    } else if (nCol <= 0) {
        error(etype, "Invalid number of columns. \n");
    }
    
    retcode = lpMatInit(dsdpSolver->lpData); checkCode;
    retcode = lpMatSetDim(dsdpSolver->lpData, dsdpSolver->m, nCol); checkCode;
    retcode = lpMatSetData(dsdpSolver->lpData, Ap, Ai, Ax); checkCode;
    
    if (dsdpSolver->verbosity) {
        printf("LP data is set. \n");
    }
    
    dsdpSolver->isLPset = TRUE;
    DSDPICheckData(dsdpSolver);
    
    return retcode;
}

extern DSDP_INT DSDPSetSDPConeData( HSDSolver *dsdpSolver,
                                    DSDP_INT  blockid,
                                    DSDP_INT  coneSize,
                                    DSDP_INT  *typehint,
                                    DSDP_INT  *Asdpp,
                                    DSDP_INT  *Asdpi,
                                    double    *Asdpx ) {
    
    DSDP_INT retcode = DSDP_RETCODE_OK;
    assert( blockid < dsdpSolver->m );
    assert( coneSize > 0 );
    
    if (dsdpSolver->insStatus != DSDP_STATUS_INIT_UNSET) {
        error(etype, "The solver instance is either not initialized or "
              "already set. \n");
    } else if ((blockid >= dsdpSolver->nBlock) || (blockid < 0)) {
        error(etype, "Invalid block id. \n");
    } else if (dsdpSolver->isSDPset[blockid]) {
        error(etype, "SDP block is already set");
    }
    
    retcode = sdpMatSetDim(dsdpSolver->sdpData[blockid],
                           dsdpSolver->m, coneSize, blockid); checkCode;
    if (typehint) {
        retcode = sdpMatSetHint(dsdpSolver->sdpData[blockid], typehint);
        checkCode;
    }
    
    retcode = sdpMatSetData(dsdpSolver->sdpData[blockid],
                            Asdpp, Asdpi, Asdpx); checkCode;
    
    if (dsdpSolver->verbosity) {
        printf("SDP block "ID" is set. \n", blockid);
    }
    
    DSDPICheckData(dsdpSolver);
    
    return retcode;
}

extern DSDP_INT DSDPDestroy( HSDSolver *dsdpSolver ) {
    
    /* Free the internal data structures */
    DSDP_INT retcode = DSDP_RETCODE_OK;
    
    retcode = DSDPIFreeLPData (dsdpSolver); checkCode;
    retcode = DSDPIFreeSDPData(dsdpSolver); checkCode;
    retcode = DSDPIFreeAlgIter(dsdpSolver); checkCode;
    retcode = DSDPIFreeCleanUp(dsdpSolver); checkCode
    
    DSDP_FREE(dsdpSolver);
    
    return retcode;
}
