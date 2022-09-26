#include <stdio.h>
#include <stdlib.h>
#include "dsdppresolve.h"
#include "dsdplog.h"
#include "dsdpoutput.h"
#include "dsdpdinfeas.h"
#include "dsdppfeas.h"
#include "dsdppsol.h"
#include "schurmat.h"
#include "adpcg.h"
#include "lanczos.h"
#include "sparseopts.h"
#include "denseopts.h"
#include "rank1opts.h"
#include "rankkopts.h"
#include "vec.h"

static char etype[] = "DSDP Interface";
#define vec_init_alloc(v, n) vecIter = (vec *) calloc(1, sizeof(vec)); \
                             (v) = vecIter;                            \
                             vec_init(vecIter);                        \
                             retcode = vec_alloc(vecIter, (n));        \
                             checkCode;
#define vec_destroy(v) vec_free(v); DSDP_FREE(v);

/* DSDP internal methods */
static DSDP_INT DSDPIInit( HSDSolver *dsdpSolver ) {
    
    // Allocate memory for the internal solver
    DSDP_INT retcode = DSDP_RETCODE_OK;
    if (dsdpSolver->insStatus != DSDP_STATUS_UNINIT) {
        error(etype, "Instance has been initialized. \n");
    }
    
    // Problem data
    dsdpSolver->sdpData = NULL;
    dsdpSolver->isSDPset = NULL;
    
    dsdpSolver->lpData  = NULL;
    dsdpSolver->lpObj   = NULL;
    dsdpSolver->isLPset = FALSE;
    
    dsdpSolver->dObj = NULL;
    
    // Dimension data
    dsdpSolver->n      = 0;
    dsdpSolver->nall   = 0;
    dsdpSolver->m      = 0;
    dsdpSolver->nBlock = 0;
    dsdpSolver->lpDim  = 0;
    
    // Monitors
    memset(dsdpSolver->eventMonitor, 0, sizeof(DSDP_INT) * nEvent);
    memset(dsdpSolver->iterProgress, 0, sizeof(DSDP_INT) * IterStep);
    
    // Iterator
    dsdpSolver->pObjVal = 0.0;
    dsdpSolver->dObjVal = 0.0;
    dsdpSolver->mu      = 0.0;
    
    dsdpSolver->S    = NULL;
    dsdpSolver->symS = NULL;
    dsdpSolver->s    = NULL;
    dsdpSolver->x    = NULL;
    
    dsdpSolver->asinv       = NULL;
    dsdpSolver->asinvrysinv = NULL;
    dsdpSolver->csinvcsinv  = 0.0;
    dsdpSolver->csinv       = 0.0;
    dsdpSolver->csinvrysinv = 0.0;
    
    dsdpSolver->M        = NULL;
    dsdpSolver->Msdp     = NULL;
    
    dsdpSolver->cg       = NULL;
    dsdpSolver->Mdiag    = NULL;
    dsdpSolver->Mscaler  = 0.0;
    
    dsdpSolver->u      = NULL;
    dsdpSolver->b1     = NULL;
    dsdpSolver->b2     = NULL;
    dsdpSolver->d1     = NULL;
    dsdpSolver->d12    = NULL;
    dsdpSolver->d2     = NULL;
    dsdpSolver->d3     = NULL;
    dsdpSolver->d4     = NULL;
    
    dsdpSolver->y      = NULL;
    dsdpSolver->tau    = 0.0;
    dsdpSolver->kappa  = 0.0;
    dsdpSolver->alpha  = 0.0;
    
    dsdpSolver->ybound = 0.0;
    dsdpSolver->dperturb = 0.0;
    dsdpSolver->Pnrm   = 0.0;
    dsdpSolver->dPotential = 0.0;
    
    // Residuals
    dsdpSolver->rysinv = 0.0;
    dsdpSolver->Ry    = 0.0;
    dsdpSolver->drate = 0.0;

    // Step matrix
    dsdpSolver->dS     = NULL;
    dsdpSolver->lczs   = NULL;
    dsdpSolver->dsaux  = NULL;
    dsdpSolver->rkaux  = NULL;
    dsdpSolver->Scker  = NULL;
    dsdpSolver->scker  = NULL;
    dsdpSolver->ds     = NULL;
    dsdpSolver->dy     = NULL;
    dsdpSolver->dtau   = 0.0;
    dsdpSolver->dkappa = 0.0;
    
    // Primal relaxation
    dsdpSolver->sl = NULL;
    dsdpSolver->slcker = NULL;
    dsdpSolver->su = NULL;
    dsdpSolver->sucker = NULL;
    dsdpSolver->pinfeas = 0.0;
    
    dsdpSolver->param     = &defaultParam;
    dsdpSolver->insStatus = DSDP_STATUS_INIT_UNSET;
    dsdpSolver->solStatus = DSDP_UNKNOWN;
    
    // Verbosity
    dsdpSolver->verbosity = TRUE;
    
    // Primal variable
    dsdpSolver->pScaler = NULL;
    dsdpSolver->cScaler = 0.0;
    dsdpSolver->ymaker  = NULL;
    dsdpSolver->dymaker = NULL;
    dsdpSolver->mumaker  = 0.0;
    dsdpSolver->ymaker2 = NULL;
    dsdpSolver->dymaker2 = NULL;
    dsdpSolver->mumaker2 = 0.0;
    
    // Solution statistic
    DSDPStatInit(&dsdpSolver->dsdpStats);
    dsdpSolver->startTime = my_clock();
    
    return retcode;
}

static DSDP_INT DSDPICheckInputDim( DSDP_INT sdpDim, DSDP_INT nBlock, DSDP_INT nConstr, DSDP_INT lpDim ) {
    
    DSDP_INT retcode = DSDP_RETCODE_OK;
    if (sdpDim <= 0) {
        printf("| Invalid SDP dimension %d <= 0. \n", sdpDim);
        retcode = DSDP_RETCODE_FAILED; return retcode;
    }
    if (lpDim < 0) {
        printf("| Invalid LP dimension %d < 0. \n", lpDim);
        retcode = DSDP_RETCODE_FAILED; return retcode;
    }
    if (nBlock <= 0) {
        printf("| Invalid block number %d <= 0. \n", nBlock);
        retcode = DSDP_RETCODE_FAILED; return retcode;
    }
    if (nConstr <= 0) {
        printf("| Invalid constraint number %d <= 0. \n", nBlock);
        retcode = DSDP_RETCODE_FAILED; return retcode;
    }
    return retcode;
}

static DSDP_INT DSDPIAlloc( HSDSolver *dsdpSolver ) {
    // Allocate memory for the internal solver (level 1)
    // The allocation only involves data/indicator arrays and the rest of memory will be
    // allocated when setting the problem data or starting to optimize
    DSDP_INT retcode = DSDP_RETCODE_OK;
    if (dsdpSolver->insStatus != DSDP_STATUS_INIT_UNSET) {
        error(etype, "| Incorrect solver instance status. \n");
    }
    DSDP_INT nblock = dsdpSolver->nBlock;
    dsdpSolver->sdpData  = (sdpMat  **) calloc(nblock, sizeof(sdpMat *));
    dsdpSolver->lczs     = (lczstep **) calloc(nblock, sizeof(lczstep *));
    dsdpSolver->lpData   = (lpMat    *) calloc(1, sizeof(lpMat));
    dsdpSolver->lpObj    = (vec      *) calloc(1, sizeof(vec));
    dsdpSolver->M        = (symM     *) calloc(1, sizeof(symM));
    
    dsdpSolver->cg       = (adpcg *) calloc(1, sizeof(adpcg));
    
    dsdpSolver->isSDPset = (DSDP_INT *) calloc(nblock, sizeof(DSDP_INT));
    
    cg_init(dsdpSolver->cg);
    
    for (DSDP_INT i = 0; i < nblock; ++i) {
        dsdpSolver->sdpData[i] = (sdpMat *) calloc(1, sizeof(sdpMat));
        sdpMatInit(dsdpSolver->sdpData[i]);
        dsdpSolver->lczs[i] = (lczstep *) calloc(1, sizeof(lczstep));
        lczInit(dsdpSolver->lczs[i]);
    }
    
    return retcode;
}

static DSDP_INT DSDPIAllocIter( HSDSolver *dsdpSolver ) {
    // Invoked after the getIdx routine is called
    // Allocate memory for the iterates
    DSDP_INT retcode = DSDP_RETCODE_OK;
    DSDP_INT nblock = dsdpSolver->nBlock, m = dsdpSolver->m, CGreuse;
    DSDP_INT lpdim = dsdpSolver->lpDim; vec *vecIter = NULL;
    
    DSDPGetIntParam(dsdpSolver, INT_PARAM_CG_REUSE, &CGreuse);
    dsdpSolver->S  = (spsMat **) calloc(nblock, sizeof(spsMat *));
    dsdpSolver->symS = (DSDP_INT **) calloc(nblock, sizeof(DSDP_INT *));
    vec_init_alloc(dsdpSolver->s, lpdim);
    vec_init_alloc(dsdpSolver->x, lpdim);
    vec_init_alloc(dsdpSolver->asinv, m);
    vec_init_alloc(dsdpSolver->asinvrysinv, m);
    
    dsdpSolver->Msdp = (schurMat *) calloc(1, sizeof(schurMat));
    schurMatInit(dsdpSolver->Msdp);
    retcode = schurMatAlloc(dsdpSolver->Msdp, m); checkCode;
    vec_init_alloc(dsdpSolver->Mdiag, m);
    retcode = cg_alloc(dsdpSolver->cg, m, sizeof(vec)); checkCode;
    cg_register(dsdpSolver->cg, dsdpSolver->Msdp,
                dsdpSolver->Mdiag, dsdpSolver->Msdp);
    
    vec_init_alloc(dsdpSolver->u, m);
    vec_init_alloc(dsdpSolver->b1, m);
    vec_init_alloc(dsdpSolver->b2, m);
    vec_init_alloc(dsdpSolver->d1, m);
    vec_init_alloc(dsdpSolver->d12, m);
    vec_init_alloc(dsdpSolver->d2, m);
    vec_init_alloc(dsdpSolver->d3, m);
    vec_init_alloc(dsdpSolver->d4, m);
    vec_init_alloc(dsdpSolver->y, m);
    
    dsdpSolver->dS = (spsMat **) calloc(nblock, sizeof(spsMat *));
    vec_init_alloc(dsdpSolver->sl,       m);
    vec_init_alloc(dsdpSolver->su,       m);
    vec_init_alloc(dsdpSolver->slcker,   m);
    vec_init_alloc(dsdpSolver->sucker,   m);
    vec_init_alloc(dsdpSolver->scker,    lpdim);
    vec_init_alloc(dsdpSolver->ds,       lpdim);
    vec_init_alloc(dsdpSolver->dy,       m);
    vec_init_alloc(dsdpSolver->ymaker,   m);
    vec_init_alloc(dsdpSolver->dymaker,  m);
    vec_init_alloc(dsdpSolver->ymaker2,  m);
    vec_init_alloc(dsdpSolver->dymaker2, m);
    
    dsdpSolver->Scker = (spsMat **) calloc(nblock, sizeof(spsMat *));
    dsdpSolver->dsaux = (dsMat  **) calloc(nblock, sizeof(dsMat *));
    dsdpSolver->rkaux = (rkMat  **) calloc(nblock, sizeof(rkMat *));
    
    for (DSDP_INT i = 0; i < nblock; ++i) {
        DSDP_INT dim = dsdpSolver->sdpData[i]->dimS;
        retcode = lczAlloc(dsdpSolver->lczs[i], dim); checkCode;
        dsdpSolver->Scker[i] = (spsMat *) calloc(1, sizeof(spsMat));
        dsdpSolver->dsaux[i] = (dsMat  *) calloc(1, sizeof(dsMat));
        dsdpSolver->rkaux[i] = (rkMat  *) calloc(1, sizeof(rkMat));
        denseMatInit(dsdpSolver->dsaux[i]);
        retcode = denseMatAlloc(dsdpSolver->dsaux[i], dim, -1); checkCode;
        rkMatInit(dsdpSolver->rkaux[i]);
        retcode = rkMatAllocIter(dsdpSolver->rkaux[i], dim); checkCode;
    }
    
    return retcode;
}

static DSDP_INT DSDPICheckData( HSDSolver *dsdpSolver ) {
    // Check whether problem data is alreay set up
    DSDP_INT retcode = DSDP_RETCODE_OK, nblock = dsdpSolver->nBlock, nblockSet = 0;
    for (DSDP_INT i = 0; i < nblock; ++i) { nblockSet += dsdpSolver->isSDPset[i]; }
    if (nblockSet == nblock) {
        dsdpSolver->insStatus = DSDP_STATUS_SET;
        retcode = DSDPIAllocIter(dsdpSolver);
        if (retcode != DSDP_RETCODE_OK) {
            printf("| Iteration memory allocation failed. \n");
        }
    }
    return retcode;
}

static void DSDPIFreeLPData( HSDSolver *dsdpSolver ) {
    // Free the internal LP data
    if (dsdpSolver->isLPset) {
        lpMatFree(dsdpSolver->lpData); dsdpSolver->isLPset = 0;
    }
    DSDP_FREE(dsdpSolver->lpData);
}

static void DSDPIFreeSDPData( HSDSolver *dsdpSolver ) {
    // Free the internal SDP data. Not responsible for isSDPSet
    for (DSDP_INT i = 0; i < dsdpSolver->nBlock; ++i) {
        if (dsdpSolver->isSDPset[i]) {
            sdpMatFree(dsdpSolver->sdpData[i]);
            DSDP_FREE(dsdpSolver->sdpData[i]);
        }
    }
    DSDP_FREE(dsdpSolver->sdpData);
}

static void DSDPIFreeAlgIter( HSDSolver *dsdpSolver ) {
    // Free the internal iteration data
    
    DSDP_INT nblock = dsdpSolver->nBlock, i;
    // S
    for (i = 0; i < nblock; ++i) {
        if (!dsdpSolver->S) { break; }
        spsMatFree(dsdpSolver->S[i]);
        DSDP_FREE(dsdpSolver->S[i]);
    }
    DSDP_FREE(dsdpSolver->S);
    
    // Ssym
    for (i = 0; i < nblock; ++i) {
        if (!dsdpSolver->symS) { break; }
        DSDP_FREE(dsdpSolver->symS[i]);
    }
    DSDP_FREE(dsdpSolver->symS);
    
    vec_destroy(dsdpSolver->s);
    vec_destroy(dsdpSolver->x);
    vec_destroy(dsdpSolver->asinv);
    vec_destroy(dsdpSolver->asinvrysinv);
    
    schurMatFree(dsdpSolver->Msdp);
    DSDP_FREE(dsdpSolver->Msdp);
    
    cg_free(dsdpSolver->cg);
    DSDP_FREE(dsdpSolver->cg);
    
    vec_destroy(dsdpSolver->Mdiag);
    vec_destroy(dsdpSolver->u);
    vec_destroy(dsdpSolver->b1);
    vec_destroy(dsdpSolver->b2);
    vec_destroy(dsdpSolver->d1);
    vec_destroy(dsdpSolver->d12);
    vec_destroy(dsdpSolver->d2);
    vec_destroy(dsdpSolver->d3);
    vec_destroy(dsdpSolver->d4);
    vec_destroy(dsdpSolver->y);

    // dS
    for (i = 0; i < nblock; ++i) {
        if (!dsdpSolver->dS) { break; }
        spsMatFree(dsdpSolver->dS[i]);
        DSDP_FREE(dsdpSolver->dS[i]);
    }
    DSDP_FREE(dsdpSolver->dS);
    
    // Lanczos
    for (i = 0; i < nblock; ++i) {
        if (!dsdpSolver->lczs) { break; }
        lczFree(dsdpSolver->lczs[i]);
        DSDP_FREE(dsdpSolver->lczs[i]);
    }
    DSDP_FREE(dsdpSolver->lczs);
    
    vec_destroy(dsdpSolver->sl);
    vec_destroy(dsdpSolver->su);
    vec_destroy(dsdpSolver->slcker);
    vec_destroy(dsdpSolver->sucker);
    vec_destroy(dsdpSolver->scker);
    vec_destroy(dsdpSolver->ds);
    vec_destroy(dsdpSolver->dy);
    vec_destroy(dsdpSolver->ymaker);
    vec_destroy(dsdpSolver->dymaker);
    vec_destroy(dsdpSolver->ymaker2);
    vec_destroy(dsdpSolver->dymaker2);
    vec_destroy(dsdpSolver->pScaler);
    
    // Scker
    for (i = 0; i < nblock; ++i) {
        if (!dsdpSolver->Scker) { break; }
        spsMatFree(dsdpSolver->Scker[i]);
        DSDP_FREE(dsdpSolver->Scker[i]);
    }
    DSDP_FREE(dsdpSolver->Scker);
    
    // dsaux
    for (i = 0; i < nblock; ++i) {
        if (!dsdpSolver->dsaux) { break; }
        denseMatFree(dsdpSolver->dsaux[i]);
        DSDP_FREE(dsdpSolver->dsaux[i]);
    }
    DSDP_FREE(dsdpSolver->dsaux);
    
    // rkaux
    for (i = 0; i < nblock; ++i) {
        if (!dsdpSolver->rkaux) { break; }
        rkMatFree(dsdpSolver->rkaux[i]);
        DSDP_FREE(dsdpSolver->rkaux[i]);
    }
    DSDP_FREE(dsdpSolver->rkaux);
    
    symSchurMatFree(dsdpSolver->M);
    DSDP_FREE(dsdpSolver->M);
}

static DSDP_INT DSDPIFreeCleanUp( HSDSolver *dsdpSolver ) {
    // Free the internal indicator arrays and some common data
    DSDP_INT retcode = DSDP_RETCODE_OK;
    
    // isSDPset
    DSDP_FREE(dsdpSolver->isSDPset);
    // lpObj
    if (dsdpSolver->isLPset) {
        vec_destroy(dsdpSolver->lpObj);
    }
    // dObj
    if (dsdpSolver->dObj) {
        vec_destroy(dsdpSolver->dObj);
    }
    // Other data
    dsdpSolver->param = NULL;      dsdpSolver->m = 0;         dsdpSolver->nBlock = 0;
    dsdpSolver->lpDim = 0;         dsdpSolver->mu = 0.0;      dsdpSolver->Mscaler = 0.0;
    dsdpSolver->csinvrysinv = 0.0; dsdpSolver->rysinv = 0.0;  dsdpSolver->Pnrm = 0.0;
    dsdpSolver->dPotential = 0.0;  dsdpSolver->tau = 0.0;     dsdpSolver->kappa = 0.0;
    dsdpSolver->alpha = 0.0;       dsdpSolver->dtau = 0.0;    dsdpSolver->dkappa = 0.0;
    dsdpSolver->mumaker = 0.0;     dsdpSolver->mumaker2 = 0.0;
    dsdpSolver->insStatus = 0;     dsdpSolver->solStatus = DSDP_STATUS_UNINIT;
    dsdpSolver->cScaler = 0.0;     dsdpSolver->ybound = 0.0;
    dsdpSolver->dperturb = 0.0;    dsdpSolver->nall = 0;      dsdpSolver->n = 0;
    dsdpSolver->drate = 0.0;       dsdpSolver->startTime = 0.0;
    
    return retcode;
}

static DSDP_INT DSDPIPresolve( HSDSolver *dsdpSolver ) {
    // Do presolve
    DSDP_INT retcode = DSDP_RETCODE_OK; DSDPStats *stat = &dsdpSolver->dsdpStats;
    double start = my_clock(), tnow, t = 0.0;
    
    if (dsdpSolver->insStatus != DSDP_STATUS_SET) {
        error(etype, "Problem data is not set up. \n");
    }
    
    showBeautifulDashlines();
    printf("| Start presolving \n");
    tnow = my_clock();
    retcode = preRank1Rdc(dsdpSolver); checkCode;
    t = my_clock() - tnow;
    printf("| - Rank one detection completes in %3.3f seconds \n", t);
    DSDPStatUpdate(stat, STAT_RONE_TIME, t);
    
    tnow = my_clock();
    retcode = preRankkRdc(dsdpSolver); checkCode;
    t = my_clock() - tnow;
    printf("| - Eigen decomposition completes in %3.3f seconds \n", t);
    DSDPStatUpdate(stat, STAT_EIG_TIME, t);
    
    tnow = my_clock();
    retcode = getMatIdx(dsdpSolver); checkCode;
    t = my_clock() - tnow;
    printf("| - Matrix statistics ready in %3.3f seconds \n", t);
    DSDPStatUpdate(stat, STAT_MATSTAT_TIME, t);
    
    tnow = my_clock();
    retcode = DSDPPrepareMAssembler(dsdpSolver); checkCode;
    retcode = DSDPCheckSchurType(dsdpSolver->M); checkCode;
    DSDPSchurReorder(dsdpSolver->M);
    t = my_clock() - tnow;
    printf("| - Schur matrix re-ordering completes in %3.3f seconds \n", t);
    DSDPStatUpdate(stat, STAT_SCHURORD_TIME, t);
    
    tnow = my_clock();
    
#if defined(compareMode) || defined(ROTSOLVE)
    dsdpSolver->cScaler = 1.0;
#else
    retcode = preSDPPrimal(dsdpSolver); checkCode;
    retcode = preSDPMatCScale(dsdpSolver); checkCode;
#endif
    t = my_clock() - tnow;
    printf("| - Scaling completes in %3.3f seconds \n", t);
    DSDPStatUpdate(stat, STAT_SCAL_TIME, t);
    
    tnow = my_clock();
    retcode = preSymbolic(dsdpSolver); checkCode;
    t = my_clock() - tnow;
    printf("| - Dual symbolic check completes in %3.3f seconds \n", t);
    DSDPStatUpdate(stat, STAT_SYMBOLIC_TIME, t);
    
    tnow = my_clock();
    printf("| - Detecting special structures \n");
    retcode = preStructureDetect(dsdpSolver); checkCode;
    t = my_clock() - tnow;
    printf("| - Special structure detection completes in %3.3f seconds \n", t);
    DSDPStatUpdate(stat, STAT_SPECIAL_DETECT, t);
    
    dsdpSolver->insStatus = DSDP_STATUS_PRESOLVED;
    t = my_clock() - start;
    DSDPStatUpdate(stat, STAT_PRESOLVE_TIME, t);
    printf("| Presolve Ends. Elapsed Time: %3.3f seconds \n", t);
    
    return retcode;
}

static DSDP_INT DSDPIPostsolve( HSDSolver *dsdpSolver ) {
    
    DSDP_INT retcode = DSDP_RETCODE_OK;
    double start = my_clock();
    
    if (dsdpSolver->pScaler) {
        vec_vdiv(dsdpSolver->y, dsdpSolver->pScaler);
    }
    
    if (dsdpSolver->cScaler) {
        vec_scale(dsdpSolver->y, dsdpSolver->cScaler);
    }
    
    for (DSDP_INT i = 0; i < dsdpSolver->nBlock; ++i) {
        spsMatScale(dsdpSolver->S[i], dsdpSolver->cScaler);
        denseMatScale(dsdpSolver->dsaux[i], dsdpSolver->cScaler);
    }
    
    DSDPStatUpdate(&dsdpSolver->dsdpStats, STAT_POSTSOLVE_TIME, my_clock() - start);
    return retcode;
}

extern void DSDPPrintVersion(void) {
    
    showBeautifulDashlines();
    printf("| Homogeneous Dual Scaling Interior Point Solver. Version %d.%d.%d (Build date %d.%d %d)"
           "                                   \n",
           VERSION_MAJOR, VERSION_MINOR, VERSION_TECHNICAL,
           BUILD_DATE_MONTH, BUILD_DATE_DAY, BUILD_DATE_YEAR);
    showBeautifulDashlines();
}

extern DSDP_INT DSDPCreate( HSDSolver **dsdpSolver, char *modelName ) {
    // Create solver
    DSDP_INT retcode = DSDP_RETCODE_OK;
    HSDSolver *solver = NULL;
    char defaultname[] = "DSDPModel", *name;
    solver = (HSDSolver *) calloc(1, sizeof(HSDSolver));
    retcode = DSDPIInit(solver);
    if (retcode != DSDP_RETCODE_OK) {
        return retcode;
    }
    *dsdpSolver = solver; name = (modelName) ? modelName : defaultname;
    strcpy(solver->SDPModel, name);
    
    return retcode;
}

extern DSDP_INT DSDPSetDim( HSDSolver *dsdpSolver,
                            DSDP_INT  sdpDim,
                            DSDP_INT  nBlock,
                            DSDP_INT  nConstrs,
                            DSDP_INT  lpDim,
                            DSDP_INT  *nNzs ) {

    DSDP_INT retcode = DSDP_RETCODE_OK, nnz;
    
    if (dsdpSolver->insStatus != DSDP_STATUS_INIT_UNSET) {
        error(etype, "| Instance not yet initialized or dimension is already set. \n");
        retcode = DSDP_RETCODE_FAILED; return retcode;
    }
    
    retcode = DSDPICheckInputDim(sdpDim, nBlock, nConstrs, lpDim); checkCode;
    nnz = (nNzs) ? *nNzs : (-1);
    showBeautifulDashlines();
    printf("| nSDPBlock: "ID" | nConstrs: "ID" | LP. Dim: "ID" | SDP. Dim: "ID" | Nonzeros: "ID" \n",
           nBlock, nConstrs, lpDim, sdpDim, nnz);
    dsdpSolver->nBlock = nBlock; dsdpSolver->m = nConstrs;
    dsdpSolver->lpDim = lpDim; dsdpSolver->n = sdpDim;
    
    retcode = DSDPIAlloc(dsdpSolver);
    if (retcode != DSDP_RETCODE_OK) {
        printf("| Failed to allocate memories for dimension-related arrays. \n");
        retcode = DSDP_RETCODE_FAILED;
    }
    return retcode;
}

extern DSDP_INT DSDPSetLPData( HSDSolver *dsdpSolver,
                               DSDP_INT  nCol,
                               DSDP_INT  *Ap,
                               DSDP_INT  *Ai,
                               double    *Ax,
                               double    *lpObj ) {
    // LP data interface for the user. Accept the CSC representaion of coefficient A
    DSDP_INT retcode = DSDP_RETCODE_OK;
    if (dsdpSolver->insStatus != DSDP_STATUS_INIT_UNSET) {
        error(etype, "The solver instance is either uninitialized or already set. \n");
    } else if (dsdpSolver->m <= 0) {
        error(etype, "Instance dimension is not set. \n");
    } else if (dsdpSolver->isLPset) {
        error(etype, "LP data is already set. \n");
    } else if (nCol != dsdpSolver->lpDim) {
        error(etype, "Invalid number of columns. \n");
    }
    lpMatInit(dsdpSolver->lpData);
    lpMatSetDim(dsdpSolver->lpData, dsdpSolver->m, nCol);
    retcode = vec_alloc(dsdpSolver->lpObj, nCol); checkCode;
    retcode = lpMatSetData(dsdpSolver->lpData, Ap, Ai, Ax); checkCode;
    memcpy(dsdpSolver->lpObj->x, lpObj, sizeof(double) * nCol);
    
    dsdpSolver->isLPset = TRUE;
    return retcode;
}

extern DSDP_INT DSDPSetSDPConeData( HSDSolver *dsdpSolver,
                                    DSDP_INT  blockid,
                                    DSDP_INT  coneSize,
                                    DSDP_INT  *Asdpp,
                                    DSDP_INT  *Asdpi,
                                    double    *Asdpx ) {
    
    DSDP_INT retcode = DSDP_RETCODE_OK;
    double cnnz = 0, tmp;
    if (dsdpSolver->insStatus != DSDP_STATUS_INIT_UNSET) {
        error(etype, "The solver instance is either uninitialized or "
              "already set. \n");
    } else if (blockid >= dsdpSolver->nBlock || blockid < 0) {
        error(etype, "Invalid block id. \n");
    } else if (dsdpSolver->isSDPset[blockid]) {
        error(etype, "SDP block is already set. \n");
    }
    
    sdpMat *data = dsdpSolver->sdpData[blockid];
    sdpMatSetDim(data, dsdpSolver->m, coneSize, blockid);
    retcode = sdpMatAlloc(data); checkCode;
    DSDPGetStats(&dsdpSolver->dsdpStats, STAT_NNZ_OBJ, &cnnz);
    retcode = sdpMatSetData(data, Asdpp, Asdpi, Asdpx, &tmp); checkCode;
    cnnz += tmp; DSDPStatUpdate(&dsdpSolver->dsdpStats, STAT_NNZ_OBJ, cnnz);
    dsdpSolver->isSDPset[blockid] = TRUE;
    
    return retcode;
}

extern DSDP_INT DSDPSetObj( HSDSolver *dsdpSolver, double *dObj ) {
    // Set the dual objective
    DSDP_INT retcode = DSDP_RETCODE_OK;
    // Most likely to run out of memory
    retcode = DSDPICheckData(dsdpSolver);
    if (retcode != DSDP_RETCODE_OK) {
        printf("| Failed to allocate internal iteration data. \n");
        return retcode;
    }
    
    dsdpSolver->dObj = (vec *) calloc(1, sizeof(vec));
    vec_init(dsdpSolver->dObj);
    
    retcode = vec_alloc(dsdpSolver->dObj, dsdpSolver->m);
    if (retcode != DSDP_RETCODE_OK) {
        printf("| Failed to allocate memory for the dual objective. \n");
        return retcode;
    }
    
    if (dObj) {
        memcpy(dsdpSolver->dObj->x, dObj, sizeof(double) * dsdpSolver->m);
    } else {
        printf("| Warning: Dual objective not set. Using 0 as dObj \n");
        memset(dsdpSolver->dObj, 0, sizeof(double) * dsdpSolver->m);
    }
    
    return retcode;
}

extern DSDP_INT DSDPOptimize( HSDSolver *dsdpSolver ) {
    // Optimization routine for DSDP
    DSDP_INT retcode = DSDP_RETCODE_OK;
    DSDP_INT gotoB = FALSE;
    
    DSDPStats *stat = &dsdpSolver->dsdpStats;
    if (!dsdpSolver->dObj) {
        retcode = DSDPSetObj(dsdpSolver, NULL);
        checkCode;
    }
    
    assert( dsdpSolver->insStatus == DSDP_STATUS_SET );
    retcode = DSDPIPresolve(dsdpSolver); checkCode;
    
    DSDPDataStatPrint(stat);
    
    assert( dsdpSolver->insStatus == DSDP_STATUS_PRESOLVED );
    retcode = DSDPDInfeasEliminator(dsdpSolver); checkCode;
    printPhaseABConvert(dsdpSolver, &gotoB);
    
    if (gotoB) {
        retcode = DSDPPFeasPhase(dsdpSolver);
    }
    printf("| DSDP Ends. %86s \n", "");
    showBeautifulDashlines();
    
#ifndef ROTSOLVE
    // Compute solution and get DIMACS errors
    computePrimalX(dsdpSolver);
    printf("| Primal solution is extracted.                      "
           "                                              \n");
    computeDIMACS(dsdpSolver);
#endif
    
    // Post-solving
    DSDPIPostsolve(dsdpSolver);
    
    // Summary statistics
    DSDPBProfilerPrint(stat);
    showBeautifulDashlines();
    
#ifndef ROTSOLVE
    DSDPDIMACErrorPrint(stat);
#endif
    
    dsdpSolver->insStatus = DSDP_STATUS_SOLVED;
    
    return retcode;
}

extern DSDP_INT DSDPGetDual( HSDSolver *dsdpSolver, double *y, double **S ) {
    // Extract dual solution y and S
    DSDP_INT retcode = DSDP_RETCODE_OK;
    if (y) {
        memcpy(y, dsdpSolver->y->x, sizeof(double) * dsdpSolver->y->dim);
    }
    if (S) {
        for (DSDP_INT i = 0; i < dsdpSolver->nBlock; ++i) {
            spsMatFill(dsdpSolver->S[i], S[i]);
        }
    }
    return retcode;
}

extern DSDP_INT DSDPGetPrimal( HSDSolver *dsdpSolver, double **X ) {
    // Extract primal solution X
    DSDP_INT retcode = DSDP_RETCODE_OK;
    retcode = computePrimalX(dsdpSolver);
    
    if (retcode != DSDP_RETCODE_OK) {
        printf("| Failed to compute the primal solution.\n");
        return retcode;
    }
    
    for (DSDP_INT i = 0; i < dsdpSolver->nBlock; ++i) {
        denseMatFill(dsdpSolver->dsaux[i], X[i]);
    }
    return retcode;
}

extern DSDP_INT DSDPExport( HSDSolver *dsdpSolver, DSDP_INT output, char *fname ) {
    // Export solver information
    DSDP_INT retcode = DSDP_RETCODE_OK;
    
    if (!dsdpSolver->dObj) {
        retcode = DSDPSetObj(dsdpSolver, NULL); checkCode;
    }
    if (output == DSDP_EXPORT_DSYMBOLIC) {
        retcode = DSDPIPresolve(dsdpSolver);
        if (retcode != DSDP_RETCODE_OK) {
            printf("| Presolve failed. \n"); return retcode;
        }
        dumpDualSymbolic(dsdpSolver, fname);
        
    } else if (output == DSDP_EXPORT_YSOL) {
        retcode = dumpDualSol(dsdpSolver, fname);
    } else {
        error(etype, "Not implemented export type. \n");
    }
    return retcode;
}

extern void DSDPSetDblParam( HSDSolver *dsdpSolver, DSDP_INT pName, double dblVal ) {
    // Set double parameter
    setDblParam(dsdpSolver->param, pName, dblVal);
}

extern void DSDPSetIntParam( HSDSolver *dsdpSolver, DSDP_INT pName, DSDP_INT intVal ) {
    // Set integer parameter
    setIntParam(dsdpSolver->param, pName, intVal);
}

extern void DSDPGetDblParam( HSDSolver *dsdpSolver, DSDP_INT pName, double *dblVal ) {
    // Get double parameter
    getDblParam(dsdpSolver->param, pName, dblVal);
}

extern void DSDPGetIntParam( HSDSolver *dsdpSolver, DSDP_INT pName, DSDP_INT *intVal ) {
    // Get double parameter
    getIntParam(dsdpSolver->param, pName, intVal);
}

extern void DSDPGetStatistic( HSDSolver *dsdpSolver, DSDP_INT sName, double *sVal ) {
    // Get statistic
    DSDPGetStats(&dsdpSolver->dsdpStats, sName, sVal);
}

extern DSDP_INT DSDPDestroy( HSDSolver *dsdpSolver ) {
    /* Free the internal data structures */
    DSDP_INT retcode = DSDP_RETCODE_OK;
    DSDPIFreeLPData (dsdpSolver); DSDPIFreeSDPData(dsdpSolver);
    DSDPIFreeAlgIter(dsdpSolver); DSDPIFreeCleanUp(dsdpSolver);
    return retcode;
}
