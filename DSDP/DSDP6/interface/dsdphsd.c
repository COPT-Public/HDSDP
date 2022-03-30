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
#include "dsdpsolver.h"
#include "hsd.h"
#include "dsdpdinfeas.h"
#include "dsdputils.h"
#include "dsdplog.h"
#include "dsdppfeas.h"
#include "dsdppsol.h"

static char etype[] = "DSDP Interface";


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
    dsdpSolver->lpObj   = NULL;
    dsdpSolver->lpData  = NULL;
    dsdpSolver->dObj    = NULL;
    
    dsdpSolver->isLPset  = FALSE;
    dsdpSolver->isSDPset = NULL;
    
    // Dimension data
    dsdpSolver->n      = 0;
    dsdpSolver->m      = 0;
    dsdpSolver->nBlock = 0;
    dsdpSolver->lpDim  = 0;
    
    // IterProgress monitor
    memset(dsdpSolver->eventMonitor, 0,
           sizeof(DSDP_INT) * nEvent);
    memset(dsdpSolver->iterProgress, 0,
           sizeof(DSDP_INT) * IterStep);
    
    // Residuals
    dsdpSolver->Ry = 0.0;
    dsdpSolver->ry = NULL;
    
    // Iterator
    dsdpSolver->pObjVal = 0.0;
    dsdpSolver->dObjVal = 0.0;
    dsdpSolver->mu      = 0.0;
    
    dsdpSolver->S    = NULL;
    dsdpSolver->symS = NULL;
    dsdpSolver->s    = NULL;
    dsdpSolver->x    = NULL;
    
    dsdpSolver->asinv       = NULL;
    dsdpSolver->csinv       = 0.0;
    dsdpSolver->csinvcsinv  = 0.0;
    dsdpSolver->csinvrysinv = 0.0;
    
    dsdpSolver->M        = NULL;
    dsdpSolver->Msdp     = NULL;
    dsdpSolver->cgSolver = NULL;
    dsdpSolver->Mdiag    = NULL;
    dsdpSolver->Mscaler  = 0.0;
    dsdpSolver->u        = NULL;
    
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
    dsdpSolver->Pnrm   = 0.0;
    dsdpSolver->dPotential = 0.0;
    
    // Step matrix
    dsdpSolver->dS     = NULL;
    dsdpSolver->Scker  = NULL;
    
    dsdpSolver->lczSolver = NULL;
    dsdpSolver->spaux  = NULL;
    dsdpSolver->dsaux  = NULL;
    dsdpSolver->rkaux  = NULL;
    dsdpSolver->ds     = NULL;
    dsdpSolver->dy     = NULL;
    dsdpSolver->dtau   = 0.0;
    dsdpSolver->dkappa = 0.0;
    
    // Primal relaxation
    dsdpSolver->sl = NULL;
    dsdpSolver->su = NULL;
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
    dsdpSolver->mumaker = 0.0;
    
    DSDPStatInit(&dsdpSolver->dsdpStats);
    
    return retcode;
}

static DSDP_INT DSDPIAlloc( HSDSolver *dsdpSolver ) {
    
    // Allocate memory for the internal solver (level 1)
    // The allocation only involves data/indicator arrays and the rest of memory will be
    // allocated when setting the problem data or starting to optimize
    DSDP_INT retcode = DSDP_RETCODE_OK;
    assert( dsdpSolver->insStatus == DSDP_STATUS_INIT_UNSET );
    
    if (dsdpSolver->insStatus != DSDP_STATUS_INIT_UNSET) {
        error(etype, "Level 1 memory cannot be allocated. \n");
        retcode = DSDP_RETCODE_FAILED;
        return retcode;
    }
    
    DSDP_INT nblock = dsdpSolver->nBlock;
    
    dsdpSolver->sdpData   = (sdpMat  **) calloc(nblock, sizeof(sdpMat *));
    dsdpSolver->lczSolver = (DSDPLanczos **) calloc(nblock, sizeof(DSDPLanczos *));
    dsdpSolver->lpData    = (lpMat    *) calloc(1,      sizeof(lpMat   ));
    dsdpSolver->lpObj     = (vec      *) calloc(1,      sizeof(vec     ));
    dsdpSolver->M         = (DSDPSchur*) calloc(1,      sizeof(DSDPSchur));
    dsdpSolver->cgSolver  = (CGSolver *) calloc(1,      sizeof(CGSolver));
    dsdpSolver->isSDPset  = (DSDP_INT *) calloc(nblock, sizeof(DSDP_INT));
    
    retcode = dsdpCGInit(dsdpSolver->cgSolver);
    
    for (DSDP_INT i = 0; i < nblock; ++i) {
        dsdpSolver->sdpData[i] = (sdpMat *) calloc(1, sizeof(sdpMat));
        sdpMatInit(dsdpSolver->sdpData[i]);
        dsdpSolver->lczSolver[i] = (DSDPLanczos *) calloc(1, sizeof(DSDPLanczos));
        dsdpLanczosInit(dsdpSolver->lczSolver[i]);
    }
    
    return retcode;
}

static DSDP_INT DSDPIAllocResi( HSDSolver *dsdpSolver ) {
    
    // Allocate memory for residuals
    DSDP_INT retcode = DSDP_RETCODE_OK;
    DSDP_INT lpdim = dsdpSolver->lpDim;

    // Allocate ry
    dsdpSolver->ry = (vec *) calloc(1, sizeof(vec));
    retcode = vec_init(dsdpSolver->ry);
    retcode = vec_alloc(dsdpSolver->ry, lpdim);
    
    return retcode;
}

static DSDP_INT DSDPIAllocIter( HSDSolver *dsdpSolver ) {
    
    // Invoked after the getIdx routine is called
    // Allocate memory for the iterates
    DSDP_INT retcode = DSDP_RETCODE_OK;
    DSDP_INT nblock = dsdpSolver->nBlock, dim = 0, m = dsdpSolver->m, CGreuse;
    
    retcode = DSDPGetIntParam(dsdpSolver, INT_PARAM_CG_REUSE, &CGreuse);
    
    dsMat *dsIter   = NULL;
    vec *vecIter    = NULL;
    
    // Allocate S and dS
    dsdpSolver->S  = (spsMat **) calloc(nblock, sizeof(spsMat *));
    dsdpSolver->dS = (spsMat **) calloc(nblock, sizeof(spsMat *));
    
    // Allocate symbolic structure
    dsdpSolver->symS = (DSDP_INT **) calloc(nblock, sizeof(DSDP_INT *));

    // Allocate s
    /*
    vecIter = (vec *) calloc(1, sizeof(vec));
    dsdpSolver->s = vecIter;
    retcode = vec_init(vecIter); checkCode;
    retcode = vec_alloc(vecIter, dsdpSolver->lpDim); checkCode;
    
    // Allocate x
    vecIter = (vec *) calloc(1, sizeof(vec));
    dsdpSolver->x = vecIter;
    retcode = vec_init(vecIter); checkCode;
    retcode = vec_alloc(vecIter, dsdpSolver->lpDim); checkCode;
    */
    
    // Allocate asinv
    dsdpSolver->asinv = (vec *) calloc(1, sizeof(vec));
    retcode = vec_init(dsdpSolver->asinv);
    retcode = vec_alloc(dsdpSolver->asinv, m);
        
    // Allocate Msdp
    dsIter = (dsMat *) calloc(1, sizeof(dsMat));
    dsdpSolver->Msdp = dsIter;
    retcode = denseMatInit(dsIter); checkCode;
    retcode = denseMatAlloc(dsIter, m, TRUE); checkCode;
    
    // Allocate CG solver
    retcode = dsdpCGAlloc(dsdpSolver->cgSolver, m);
    retcode = dsdpCGSetTol(dsdpSolver->cgSolver, 1e-05);
    retcode = dsdpCGSetPreReuse(dsdpSolver->cgSolver, CGreuse);
    retcode = dsdpCGSetM(dsdpSolver->cgSolver, dsIter);
    retcode = dsdpCGSetCholPre(dsdpSolver->cgSolver, dsdpSolver->Msdp);
    
    // Allocate Mdiag, u, b1, b2, d1, d12, d2, d3 and d4
    vecIter = (vec *) calloc(1, sizeof(vec)); dsdpSolver->Mdiag = vecIter;
    retcode = vec_init(vecIter); retcode = vec_alloc(vecIter, m); checkCode;
    retcode = dsdpCGSetDPre(dsdpSolver->cgSolver, dsdpSolver->Mdiag);
    vecIter = (vec *) calloc(1, sizeof(vec)); dsdpSolver->u = vecIter;
    retcode = vec_init(vecIter); retcode = vec_alloc(vecIter, m); checkCode;
    vecIter = (vec *) calloc(1, sizeof(vec)); dsdpSolver->b1 = vecIter;
    retcode = vec_init(vecIter); retcode = vec_alloc(vecIter, m); checkCode;
    vecIter = (vec *) calloc(1, sizeof(vec)); dsdpSolver->b2 = vecIter;
    retcode = vec_init(vecIter); retcode = vec_alloc(vecIter, m); checkCode;
    vecIter = (vec *) calloc(1, sizeof(vec)); dsdpSolver->d1 = vecIter;
    retcode = vec_init(vecIter); vec_alloc(vecIter, m); checkCode;
    vecIter = (vec *) calloc(1, sizeof(vec)); dsdpSolver->d12 = vecIter;
    retcode = vec_init(vecIter); retcode = vec_alloc(vecIter, m); checkCode;
    vecIter = (vec *) calloc(1, sizeof(vec)); dsdpSolver->d2 = vecIter;
    retcode = vec_init(vecIter); retcode = vec_alloc(vecIter, m); checkCode;
    vecIter = (vec *) calloc(1, sizeof(vec)); dsdpSolver->d3 = vecIter;
    retcode = vec_init(vecIter); retcode = vec_alloc(vecIter, m); checkCode;
    vecIter = (vec *) calloc(1, sizeof(vec)); dsdpSolver->d4 = vecIter;
    retcode = vec_init(vecIter); retcode = vec_alloc(vecIter, m); checkCode;
    // Allocate y
    vecIter = (vec *) calloc(1, sizeof(vec)); dsdpSolver->y = vecIter;
    retcode = vec_init(vecIter); retcode = vec_alloc(vecIter, m); checkCode;
    // ds, dy, ymaker, dymaker
    vecIter = (vec *) calloc(1, sizeof(vec)); dsdpSolver->ds = vecIter;
    retcode = vec_init(vecIter); retcode = vec_alloc(vecIter, dsdpSolver->lpDim); checkCode;
    vecIter = (vec *) calloc(1, sizeof(vec)); dsdpSolver->dy = vecIter;
    retcode = vec_init(vecIter); retcode = vec_alloc(vecIter, m); checkCode;
    vecIter = (vec *) calloc(1, sizeof(vec)); dsdpSolver->ymaker = vecIter;
    retcode = vec_init(vecIter); retcode = vec_alloc(vecIter, m); checkCode;
    vecIter = (vec *) calloc(1, sizeof(vec)); dsdpSolver->dymaker = vecIter;
    retcode = vec_init(vecIter); retcode = vec_alloc(vecIter, m); checkCode;
    
    // sl, su
    vecIter = (vec *) calloc(1, sizeof(vec)); dsdpSolver->sl = vecIter;
    retcode = vec_init(vecIter); retcode = vec_alloc(vecIter, dsdpSolver->m); checkCode;
    vecIter = (vec *) calloc(1, sizeof(vec)); dsdpSolver->su = vecIter;
    retcode = vec_init(vecIter); retcode = vec_alloc(vecIter, dsdpSolver->m); checkCode;
    
    // Allocate Scker and spaux
    dsdpSolver->Scker = (spsMat **) calloc(nblock, sizeof(spsMat *));
    dsdpSolver->spaux = (spsMat **) calloc(nblock, sizeof(spsMat *));
    dsdpSolver->dsaux = (dsMat  **) calloc(nblock, sizeof(dsMat *));
    dsdpSolver->rkaux = (rkMat  **) calloc(nblock, sizeof(rkMat *));
    
    for (DSDP_INT i = 0; i < nblock; ++i) {
        dim = dsdpSolver->sdpData[i]->dimS;
        retcode = dsdpLanczosAlloc(dsdpSolver->lczSolver[i], dim); checkCode;
        dsdpSolver->spaux[i] = (spsMat *) calloc(1, sizeof(spsMat));
        dsdpSolver->Scker[i] = (spsMat *) calloc(1, sizeof(spsMat));
        dsdpSolver->dsaux[i] = (dsMat  *) calloc(1, sizeof(dsMat));
        dsdpSolver->rkaux[i] = (rkMat  *) calloc(1, sizeof(rkMat));
        retcode = spsMatInit(dsdpSolver->spaux[i]); checkCode;
        retcode = spsMatAlloc(dsdpSolver->spaux[i], dim); checkCode;
        retcode = denseMatInit(dsdpSolver->dsaux[i]); checkCode;
        retcode = denseMatAlloc(dsdpSolver->dsaux[i], dim, -1); checkCode;
        retcode = rkMatInit(dsdpSolver->rkaux[i]); checkCode;
        retcode = rkMatAllocIter(dsdpSolver->rkaux[i], dim); checkCode;
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
    if (nblockSet == nblock) {
        dsdpSolver->insStatus = DSDP_STATUS_SET;
        retcode = DSDPIAllocIter(dsdpSolver);
        retcode = DSDPIAllocResi(dsdpSolver);
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
    
    for (DSDP_INT i = 0; i < dsdpSolver->nBlock; ++i) {
        if (dsdpSolver->isSDPset[i]) {
            retcode = sdpMatFree(dsdpSolver->sdpData[i]); checkCode;
            DSDP_FREE(dsdpSolver->sdpData[i]);
        }
    }
    
    DSDP_FREE(dsdpSolver->sdpData);
    
    return retcode;
}

static DSDP_INT DSDPIFreeResi( HSDSolver *dsdpSolver ) {
    
    // Free the internal residual data
    DSDP_INT retcode = DSDP_RETCODE_OK;
    DSDP_FREE(dsdpSolver->ry);
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
    
    // Ssym
    for (DSDP_INT i = 0; i < nblock; ++i) {
        DSDP_FREE(dsdpSolver->symS[i]);
    }
    DSDP_FREE(dsdpSolver->symS);
    
    retcode = vec_free(dsdpSolver->s);
    DSDP_FREE(dsdpSolver->s);
    retcode = vec_free(dsdpSolver->x);
    DSDP_FREE(dsdpSolver->x);
    retcode = vec_free(dsdpSolver->asinv);
    DSDP_FREE(dsdpSolver->asinv);
    retcode = denseMatFree(dsdpSolver->Msdp);
    DSDP_FREE(dsdpSolver->Msdp);
    retcode = dsdpCGFree(dsdpSolver->cgSolver);
    DSDP_FREE(dsdpSolver->cgSolver);
    retcode = vec_free(dsdpSolver->Mdiag);
    DSDP_FREE(dsdpSolver->Mdiag);
    retcode = vec_free(dsdpSolver->u );
    retcode = vec_free(dsdpSolver->d1);
    retcode = vec_free(dsdpSolver->d12);
    retcode = vec_free(dsdpSolver->d2);
    retcode = vec_free(dsdpSolver->d3);
    retcode = vec_free(dsdpSolver->d4);
    retcode = vec_free(dsdpSolver->y);
    
    DSDP_FREE(dsdpSolver->u ); DSDP_FREE(dsdpSolver->d1);
    DSDP_FREE(dsdpSolver->d12); DSDP_FREE(dsdpSolver->d2);
    DSDP_FREE(dsdpSolver->d3); DSDP_FREE(dsdpSolver->d4);
    DSDP_FREE(dsdpSolver->y);

    // dS
    for (DSDP_INT i = 0; i < nblock; ++i) {
        retcode = spsMatFree(dsdpSolver->dS[i]); checkCode;
        DSDP_FREE(dsdpSolver->dS[i]);
    }
    
    DSDP_FREE(dsdpSolver->dS);
    
    // Scker
    for (DSDP_INT i = 0; i < nblock; ++i) {
        retcode = spsMatFree(dsdpSolver->Scker[i]); checkCode;
        DSDP_FREE(dsdpSolver->Scker[i]);
    }
    DSDP_FREE(dsdpSolver->Scker);
    
    // spaux
    for (DSDP_INT i = 0; i < nblock; ++i) {
        retcode = spsMatFree(dsdpSolver->spaux[i]); checkCode;
        DSDP_FREE(dsdpSolver->spaux[i]);
    }
    DSDP_FREE(dsdpSolver->spaux);
    
    // dsaux
    for (DSDP_INT i = 0; i < nblock; ++i) {
        retcode = denseMatFree(dsdpSolver->dsaux[i]); checkCode;
        DSDP_FREE(dsdpSolver->dsaux[i]);
    }
    DSDP_FREE(dsdpSolver->dsaux);
    
    // rkaux
    for (DSDP_INT i = 0; i < nblock; ++i) {
        retcode = rkMatFree(dsdpSolver->rkaux[i]); checkCode;
        DSDP_FREE(dsdpSolver->rkaux[i]);
    }
    DSDP_FREE(dsdpSolver->rkaux);
    
    retcode = vec_free(dsdpSolver->ds); checkCode;
    retcode = vec_free(dsdpSolver->dy); checkCode;
    retcode = vec_free(dsdpSolver->sl); checkCode;
    retcode = vec_free(dsdpSolver->su); checkCode;
    retcode = vec_free(dsdpSolver->pScaler); checkCode;
    retcode = vec_free(dsdpSolver->ymaker); checkCode;
    retcode = vec_free(dsdpSolver->dymaker); checkCode;
    
    DSDP_FREE(dsdpSolver->ds); DSDP_FREE(dsdpSolver->dy);
    DSDP_FREE(dsdpSolver->sl); DSDP_FREE(dsdpSolver->su);
    DSDP_FREE(dsdpSolver->pScaler); DSDP_FREE(dsdpSolver->ymaker);
    DSDP_FREE(dsdpSolver->dymaker);
    
    SchurMatFree(dsdpSolver->M); DSDP_FREE(dsdpSolver->M);
    
    return retcode;
}

static DSDP_INT DSDPIFreeCleanUp( HSDSolver *dsdpSolver ) {
    
    // Free the internal indicator arrays and some common data
    DSDP_INT retcode = DSDP_RETCODE_OK;
    
    // isSDPset
    DSDP_FREE(dsdpSolver->isSDPset);
    
    // lpObj
    if (dsdpSolver->isLPset) {
        retcode = vec_free(dsdpSolver->lpObj); checkCode;
    }
    DSDP_FREE(dsdpSolver->lpObj);
    
    // dObj
    if (dsdpSolver->dObj) {
        retcode = vec_free(dsdpSolver->dObj); checkCode;
        DSDP_FREE(dsdpSolver->dObj);
    }
    
    // Other data
    dsdpSolver->param = NULL;      dsdpSolver->m = 0;         dsdpSolver->nBlock = 0;
    dsdpSolver->lpDim = 0;         dsdpSolver->mu = 0.0;      dsdpSolver->Mscaler = 0.0;
    dsdpSolver->csinvrysinv = 0.0; dsdpSolver->rysinv = 0.0;     dsdpSolver->Pnrm = 0.0;
    dsdpSolver->dPotential = 0.0;  dsdpSolver->tau = 0.0;     dsdpSolver->kappa = 0.0;
    dsdpSolver->alpha = 0.0;       dsdpSolver->dtau = 0.0;    dsdpSolver->dkappa = 0.0;
    dsdpSolver->mumaker = 0.0;     dsdpSolver->insStatus = 0; dsdpSolver->solStatus = 0;
    dsdpSolver->schurmu = 0.0;     dsdpSolver->cScaler = 0.0;
    
    return retcode;
}

static DSDP_INT DSDPIPresolve( HSDSolver *dsdpSolver ) {
    
    // Do presolve
    DSDP_INT retcode = DSDP_RETCODE_OK;
    DSDPStats *stat = &dsdpSolver->dsdpStats;
    clock_t start = clock(), center;
    double t = 0.0;
    
    assert( dsdpSolver->insStatus == DSDP_STATUS_SET );
    if ( dsdpSolver->insStatus != DSDP_STATUS_SET ) {
        error(etype, "Problem data is not set up. \n");
    }
    
    dsdpshowdash();
    printf("| Start presolving                          "
           "                                                       | \n");
    center = clock();
    retcode = preRank1Rdc(dsdpSolver); checkCode;
    t = (double) (clock() - center) / CLOCKS_PER_SEC;
    printf("| - Rank one detection completes in %g seconds \n", t);
    DSDPStatUpdate(stat, STAT_RONE_TIME, t);
    
    center = clock();
    retcode = preRankkRdc(dsdpSolver); checkCode;
    t = (double) (clock() - center) / CLOCKS_PER_SEC;
    printf("| - Eigen decomposition completes in %g seconds \n", t);
    DSDPStatUpdate(stat, STAT_EIG_TIME, t);
    
    center = clock();
    retcode = preSDPPrimal(dsdpSolver); checkCode;
    retcode = preSDPMatCScale(dsdpSolver);
//    dsdpSolver->cScaler = 1.0;
    // retcode = preSDPDual(dsdpSolver); checkCode;
    t = (double) (clock() - center) / CLOCKS_PER_SEC;
    printf("| - Scaling completes in %g seconds \n", t);
    DSDPStatUpdate(stat, STAT_SCAL_TIME, t);
    
    center = clock();
    retcode = getMatIdx(dsdpSolver); checkCode;
    t = (double) (clock() - center) / CLOCKS_PER_SEC;
    printf("| - Matrix statistics ready in %g seconds \n", t);
    DSDPStatUpdate(stat, STAT_MATSTAT_TIME, t);
    
    center = clock();
    retcode = preSymbolic(dsdpSolver); checkCode;
    t = (double) (clock() - center) / CLOCKS_PER_SEC;
    printf("| - Dual symbolic check completes in %g seconds \n", t);
    DSDPStatUpdate(stat, STAT_SYMBOLIC_TIME, t);
    
    center = clock();
    retcode = DSDPPrepareMAssembler(dsdpSolver); checkCode;
    retcode = DSDPSchurReorder(dsdpSolver->M); checkCode;
    t = (double) (clock() - center) / CLOCKS_PER_SEC;
    printf("| - Schur matrix re-ordering completes in %g seconds \n", t);
    DSDPStatUpdate(stat, STAT_SCHURORD_TIME, t);
    
    dsdpSolver->insStatus = DSDP_STATUS_PRESOLVED;
    t = (double) (clock() - start) / CLOCKS_PER_SEC;
    retcode = DSDPStatUpdate(stat, STAT_PRESOLVE_TIME, t);
    printf("| Presolve Ends. Elapsed Time: %10.3e. \n", t);
    
    return retcode;
}

static DSDP_INT DSDPIPostsolve( HSDSolver *dsdpSolver ) {
    
    DSDP_INT retcode = DSDP_RETCODE_OK;
    clock_t start = clock();
    
    if (dsdpSolver->solStatus != DSDP_OPTIMAL || !dsdpSolver->pScaler) {
        return retcode;
    } else {
        vec_vdiv(dsdpSolver->y, dsdpSolver->pScaler);
    }
    
    if (dsdpSolver->cScaler) {
        vec_scale(dsdpSolver->y, dsdpSolver->cScaler);
    }
    
    for (DSDP_INT i = 0; i < dsdpSolver->nBlock; ++i) {
        spsMatScale(dsdpSolver->S[i], dsdpSolver->cScaler);
        denseMatScale(dsdpSolver->dsaux[i], dsdpSolver->cScaler);
    }
    
    DSDPStatUpdate(&dsdpSolver->dsdpStats, STAT_POSTSOLVE_TIME,
                   (double) (clock() - start) / CLOCKS_PER_SEC);
    
    return retcode;
}

extern void DSDPPrintVersion(void) {
    dsdpshowdash();
    printf("| Homogeneous Dual Scaling Interior Point Solver. Version %d.%d.%d "
           "                                   |\n",
           VERSION_MAJOR, VERSION_MINOR, VERSION_TECHNICAL);
    dsdpshowdash();
}

extern DSDP_INT DSDPCreate( HSDSolver **dsdpSolver ) {
    
    /* Create solver */
    DSDP_INT retcode = DSDP_RETCODE_OK;
    HSDSolver *solver = NULL;
    solver = (HSDSolver *) calloc(1, sizeof(HSDSolver));
    retcode = DSDPIInit(solver);
    *dsdpSolver = solver;
    return retcode;
}

extern DSDP_INT DSDPSetDim( HSDSolver *dsdpSolver,
                            DSDP_INT  sdpDim,
                            DSDP_INT  nBlock,
                            DSDP_INT  nConstrs,
                            DSDP_INT  lpDim,
                            DSDP_INT  *nNzs ) {
    
    /* Set dimension of the DSDP problem instance
       
       nVars    is the total number of variables of the instance
       nBlock   is the number of SDP blocks participating in the instance
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
    
    if ((sdpDim + lpDim) <= 0 || nConstrs < 0 || lpDim < 0 || nBlock < 0 || sdpDim < 0) {
        error(etype, "Invalid dimension. \n");
    }
    
    DSDP_INT nnz = (nNzs) ? *nNzs : (-1);
    if (dsdpSolver->verbosity) {
        // printf("Dimension is successfully set. \n");
        dsdpshowdash();
        printf("| nBlock: "ID" "
               "| nConstrs: "ID" "
               "| nLPVars: "ID" "
               "| nSDPVars: "ID" "
               "| nNonzeros: "ID" |\n", nBlock, nConstrs, lpDim, sdpDim, nnz);
    }
    
    dsdpSolver->nBlock = nBlock; dsdpSolver->m = nConstrs;
    dsdpSolver->lpDim = lpDim; dsdpSolver->n = sdpDim;
    
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
    memcpy(dsdpSolver->lpObj->x, lpObj, sizeof(double) * nCol);
    
    dsdpSolver->isLPset = TRUE; DSDPICheckData(dsdpSolver);
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
    assert( blockid < dsdpSolver->nBlock );
    assert( coneSize > 0 );
    
    if (dsdpSolver->insStatus != DSDP_STATUS_INIT_UNSET) {
        error(etype, "The solver instance is either not initialized or "
              "already set. \n");
    } else if ((blockid >= dsdpSolver->nBlock) || (blockid < 0)) {
        error(etype, "Invalid block id. \n");
    } else if (dsdpSolver->isSDPset[blockid]) {
        error(etype, "SDP block is already set. \n");
    }
    
    sdpMat *data = dsdpSolver->sdpData[blockid];
    retcode = sdpMatSetDim(data, dsdpSolver->m,
                           coneSize, blockid); checkCode;
    retcode = sdpMatAlloc(data); checkCode;
    if (typehint) { retcode = sdpMatSetHint(data, typehint); checkCode; }
    
    retcode = sdpMatSetData(data, Asdpp, Asdpi, Asdpx); checkCode;
    dsdpSolver->isSDPset[blockid] = TRUE;
    DSDPICheckData(dsdpSolver);
    
    return retcode;
}

extern DSDP_INT DSDPSetObj( HSDSolver *dsdpSolver, double *dObj ) {
    // Set the dual objective
    DSDP_INT retcode = DSDP_RETCODE_OK;
    
    assert( dsdpSolver->insStatus == DSDP_STATUS_SET && !dsdpSolver->dObj );
    dsdpSolver->dObj = (vec *) calloc(1, sizeof(vec));
    vec_init(dsdpSolver->dObj); vec_alloc(dsdpSolver->dObj, dsdpSolver->m);
    
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
    DSDPParamPrint(dsdpSolver->param);
    
    assert( dsdpSolver->insStatus == DSDP_STATUS_PRESOLVED );
    retcode = DSDPDInfeasEliminator(dsdpSolver); checkCode;
    printPhaseABConvert(dsdpSolver, &gotoB);
    
    if (gotoB) {
        retcode = DSDPPFeasPhase(dsdpSolver);
    }
    printf("| DSDP Ends. %86s| \n", "");
    dsdpshowdash();
    
    // Compute solution and get DIMACS errors
    computePrimalX(dsdpSolver);
    printf("| Primal solution is extracted.                      "
           "                                              |\n");
    dsdpshowdash();
    computeDIMACS(dsdpSolver);
    
    // Post-solving
    DSDPIPostsolve(dsdpSolver);
    
    // Summary statistics
    DSDPBProfilerPrint(stat);
    dsdpshowdash();
    DSDPDIMACErrorPrint(stat);
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
    
    for (DSDP_INT i = 0; i < dsdpSolver->nBlock; ++i) {
        denseMatFill(dsdpSolver->dsaux[i], X[i]);
    }
    return retcode;
}

extern DSDP_INT DSDPSetDblParam( HSDSolver *dsdpSolver, DSDP_INT pName, double dblVal ) {
    // Set double parameter
    return setDblParam(dsdpSolver->param, pName, dblVal);
}

extern DSDP_INT DSDPSetIntParam( HSDSolver *dsdpSolver, DSDP_INT pName, DSDP_INT intVal ) {
    // Set integer parameter
    return setIntParam(dsdpSolver->param, pName, intVal);
}

extern DSDP_INT DSDPGetDblParam( HSDSolver *dsdpSolver, DSDP_INT pName, double *dblVal ) {
    // Get double parameter
    return getDblParam(dsdpSolver->param, pName, dblVal);
}

extern DSDP_INT DSDPGetIntParam( HSDSolver *dsdpSolver, DSDP_INT pName, DSDP_INT *intVal ) {
    // Get double parameter
    return getIntParam(dsdpSolver->param, pName, intVal);
}

extern DSDP_INT DSDPDestroy( HSDSolver *dsdpSolver ) {
    
    /* Free the internal data structures */
    DSDP_INT retcode = DSDP_RETCODE_OK;
    
    retcode = DSDPIFreeLPData (dsdpSolver); checkCode;
    retcode = DSDPIFreeSDPData(dsdpSolver); checkCode;
    retcode = DSDPIFreeResi(dsdpSolver);    checkCode;
    retcode = DSDPIFreeAlgIter(dsdpSolver); checkCode;
    retcode = DSDPIFreeCleanUp(dsdpSolver); checkCode;
    
    DSDP_FREE(dsdpSolver);
    
    return retcode;
}
