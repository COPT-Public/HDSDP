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
#include "dsdputils.h"

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
    
    dsdpSolver->S = NULL;
    dsdpSolver->s = NULL;
    dsdpSolver->x = NULL;
    
    dsdpSolver->asinv       = NULL;
    dsdpSolver->csinv       = 0.0;
    dsdpSolver->csinvcsinv  = 0.0;
    dsdpSolver->csinvrysinv = 0.0;
    
    dsdpSolver->Msdp   = NULL;
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
    dsdpSolver->Pnrm   = 0.0;
    
    // Step matrix
    dsdpSolver->dS     = NULL;
    dsdpSolver->spaux  = NULL;
    dsdpSolver->ds     = NULL;
    dsdpSolver->dy     = NULL;
    dsdpSolver->dtau   = 0.0;
    dsdpSolver->dkappa = 0.0;
    
    dsdpSolver->param     = &defaultParam;
    dsdpSolver->iterA     = 0;
    dsdpSolver->iterB     = 0;
    dsdpSolver->insStatus = DSDP_STATUS_INIT_UNSET;
    dsdpSolver->solStatus = DSDP_UNKNOWN;
    
    // Verbosity
    dsdpSolver->verbosity = TRUE;
    
    // Primal variable
    dsdpSolver->pScaler     = NULL;
    dsdpSolver->X           = NULL;
    dsdpSolver->isXComputed = FALSE;
    
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
    dsdpSolver->lpData    = (lpMat    *) calloc(1,      sizeof(lpMat   ));
    dsdpSolver->lpObj     = (vec      *) calloc(1,      sizeof(vec     ));
    dsdpSolver->isSDPset  = (DSDP_INT *) calloc(nblock, sizeof(DSDP_INT));
    
    for (DSDP_INT i = 0; i < nblock; ++i) {
        dsdpSolver->sdpData[i] = (sdpMat *) calloc(1, sizeof(sdpMat));
        retcode = sdpMatInit(dsdpSolver->sdpData[i]); checkCode;
    }
    
    if (dsdpSolver->verbosity) {
        printf("Level 1 memory set. \n");
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

static DSDP_INT DSDPIGetBlockSymbolic( HSDSolver *dsdpSolver, DSDP_INT blockid ) {
    
    // Get the symbolic structure for dual matrix S/dS in block k and allocate memory
    // By the time this symbolic phase is called, the problem data should have gone
    // through presolving reduction
    
    DSDP_INT retcode = DSDP_RETCODE_OK;
    
    sdpMat *sdpBlock = dsdpSolver->sdpData[blockid];
    assert( sdpBlock->blockId == blockid );
    
    DSDP_INT dim     = sdpBlock->dimS;
    DSDP_INT nhash   = nsym(dim);
    DSDP_INT nnz     = 0;
    
    dsdpSolver->S[blockid]  = (spsMat *) calloc(1, sizeof(spsMat));
    dsdpSolver->dS[blockid] = (spsMat *) calloc(1, sizeof(spsMat));
    
    retcode = spsMatInit(dsdpSolver->S[blockid]); checkCode;
    retcode = spsMatInit(dsdpSolver->dS[blockid]); checkCode;
    
    spsMat *spsdata  = NULL;
    r1Mat  *r1data   = NULL;
    
    // Get hash table
    DSDP_INT *hash = NULL;
    hash = (DSDP_INT *) calloc(nhash, sizeof(DSDP_INT));
    
    DSDP_INT useDenseS = FALSE;
    DSDP_INT *matIdx = sdpBlock->r1MatIdx;
    
    if (sdpBlock->ndenseMat > 0) {
        useDenseS = TRUE;
    } else {
        for (DSDP_INT i = 0; i < sdpBlock->nr1Mat; ++i) {
            r1data = (r1Mat *) sdpBlock->sdpData[matIdx[i]];
            if (r1data->nnz > 0.5 * dim) {
                useDenseS = TRUE;
            }
        }
    }
    
    // Dense matrix is present
    if (!useDenseS) {
        // Sparse matrix
        matIdx = sdpBlock->spsMatIdx;
        for (DSDP_INT i = 0; i < sdpBlock->nspsMat; ++i) {
            spsdata = (spsMat *) sdpBlock->sdpData[matIdx[i]];
            for (DSDP_INT j = 0; j < dim; ++j) {
                for (DSDP_INT k = spsdata->p[j]; k < spsdata->p[j + 1]; ++k) {
                    if (packIdx(hash, dim, spsdata->i[k], j) == 0) {
                        packIdx(hash, dim, spsdata->i[k], j) = 1;
                        nnz += 1;
                    }
                }
            }
            if (nnz > 0.8 * nhash) {
                useDenseS = TRUE;
            }
        }
        // Rank 1 matrix
        matIdx = sdpBlock->r1MatIdx;
        for (DSDP_INT i = 0; i < sdpBlock->nr1Mat; ++i) {
            r1data = (r1Mat *) sdpBlock->sdpData[matIdx[i]];
            for (DSDP_INT j = 0; j < r1data->nnz; ++j) {
                for (DSDP_INT k = 0; k <= j; ++k) {
                    if (packIdx(hash, dim, r1data->nzIdx[j], r1data->nzIdx[k]) == 0) {
                        packIdx(hash, dim, r1data->nzIdx[j], r1data->nzIdx[k]) = 1;
                        nnz += 1;
                    }
                }
            }
            if (nnz > 0.8 * nhash) {
                useDenseS = TRUE;
            }
        }
    }
    
    if (useDenseS) {
        nnz = nhash;
        DSDP_INT idx = 0;
        for (DSDP_INT i = 0; i < dim; ++i) {
            idx += i;
            for (DSDP_INT j = 0; j <= dim; ++j) {
                hash[j] = idx;
            }
        }
    } else {
        // Get hash table by the symbolic features
        DSDP_INT idx = 0;
        for (DSDP_INT i = 0; i < nhash; ++i) {
            if (hash[i]) {
                hash[i] = idx;
                idx ++;
            }
        }
        assert( idx == nnz );
    }
    
    retcode = spsMatAllocData(dsdpSolver->S[blockid], dim, nnz); checkCode;
    retcode = spsMatAllocData(dsdpSolver->dS[blockid], dim, nnz); checkCode;
    retcode = spsMatAllocSumMat(dsdpSolver->S[blockid]); checkCode;
    retcode = spsMatAllocSumMat(dsdpSolver->dS[blockid]); checkCode;
    
    memcpy(dsdpSolver->S[blockid]->sumHash, hash, nhash);
    memcpy(dsdpSolver->dS[blockid]->sumHash, hash, nhash);
    
#ifdef SHOWALL
        printf("Block "ID" goes through symbolic check. \n", blockid);
#endif
    
    DSDP_FREE(hash);
    return retcode;
}

static DSDP_INT DSDPIAllocIter( HSDSolver *dsdpSolver ) {
    
    // Invoked after the getIdx routine is called
    // Allocate memory for the iterates
    DSDP_INT retcode = DSDP_RETCODE_OK;
    DSDP_INT nblock  = dsdpSolver->nBlock;
    DSDP_INT dim     = 0;
    DSDP_INT m       = dsdpSolver->m;
    
    dsMat *dsIter   = NULL;
    vec *vecIter    = NULL;
    
    // Allocate S and dS with symbolic hash table
    dsdpSolver->S  = (spsMat **) calloc(nblock, sizeof(spsMat *));
    dsdpSolver->dS = (spsMat **) calloc(nblock, sizeof(spsMat *));
    for (DSDP_INT i = 0; i < nblock; ++i) {
        dim = dsdpSolver->sdpData[i]->dimS;
        retcode = DSDPIGetBlockSymbolic(dsdpSolver, i);
    }
    
    // Allocate s
    vecIter = (vec *) calloc(1, sizeof(vec));
    dsdpSolver->s = vecIter;
    retcode = vec_init(vecIter); checkCode;
    retcode = vec_alloc(vecIter, dsdpSolver->lpDim); checkCode;
    
    // Allocate x
    vecIter = (vec *) calloc(1, sizeof(vec));
    dsdpSolver->x = vecIter;
    retcode = vec_init(vecIter); checkCode;
    retcode = vec_alloc(vecIter, dsdpSolver->lpDim); checkCode;
    
    // Allocate asinv
    dsdpSolver->asinv = (vec *) calloc(1, sizeof(vec));
    retcode = vec_init(dsdpSolver->asinv);
    retcode = vec_alloc(dsdpSolver->asinv, m);
        
    // Allocate Msdp
    dsIter = (dsMat *) calloc(1, sizeof(dsMat));
    dsdpSolver->Msdp = dsIter;
    retcode = denseMatInit(dsIter); checkCode;
    retcode = denseMatAlloc(dsIter, m, TRUE); checkCode;
    
    // Allocate u, b1, b2, d1, d12, d2, d3 and d4
    vecIter = (vec *) calloc(1, sizeof(vec));
    dsdpSolver->u = vecIter;
    retcode = vec_init(vecIter); checkCode;
    retcode = vec_alloc(vecIter, m); checkCode;
    
    vecIter = (vec *) calloc(1, sizeof(vec));
    dsdpSolver->b1 = vecIter;
    retcode = vec_init(vecIter); checkCode;
    retcode = vec_alloc(vecIter, m); checkCode;
    
    vecIter = (vec *) calloc(1, sizeof(vec));
    dsdpSolver->b2 = vecIter;
    retcode = vec_init(vecIter); checkCode;
    retcode = vec_alloc(vecIter, m); checkCode;
    
    vecIter = (vec *) calloc(1, sizeof(vec));
    dsdpSolver->d1 = vecIter;
    retcode = vec_init(vecIter); checkCode;
    retcode = vec_alloc(vecIter, m); checkCode;
    
    vecIter = (vec *) calloc(1, sizeof(vec));
    dsdpSolver->d12 = vecIter;
    retcode = vec_init(vecIter); checkCode;
    retcode = vec_alloc(vecIter, m); checkCode;
    
    vecIter = (vec *) calloc(1, sizeof(vec));
    dsdpSolver->d2 = vecIter;
    retcode = vec_init(vecIter); checkCode;
    retcode = vec_alloc(vecIter, m); checkCode;
    
    vecIter = (vec *) calloc(1, sizeof(vec));
    dsdpSolver->d3 = vecIter;
    retcode = vec_init(vecIter); checkCode;
    retcode = vec_alloc(vecIter, m); checkCode;
    
    vecIter = (vec *) calloc(1, sizeof(vec));
    dsdpSolver->d4 = vecIter;
    retcode = vec_init(vecIter); checkCode;
    retcode = vec_alloc(vecIter, m); checkCode;
    
    // Allocate y
    vecIter = (vec *) calloc(1, sizeof(vec));
    dsdpSolver->y = vecIter;
    retcode = vec_init(vecIter); checkCode;
    retcode = vec_alloc(vecIter, m); checkCode;
    
    // Allocate spaux
    dsdpSolver->spaux = (spsMat **) calloc(nblock, sizeof(spsMat *));
    for (DSDP_INT i = 0; i < nblock; ++i) {
        dim = dsdpSolver->sdpData[i]->dimS;
        dsdpSolver->spaux[i] = (spsMat *) calloc(1, sizeof(spsMat));
        retcode = spsMatInit(dsdpSolver->spaux[i]); checkCode;
        retcode = spsMatAlloc(dsdpSolver->spaux[i], dim); checkCode;
    }
    
    // Allocate ds
    vecIter = (vec *) calloc(1, sizeof(vec));
    dsdpSolver->ds = vecIter;
    retcode = vec_init(vecIter); checkCode;
    retcode = vec_alloc(vecIter, dsdpSolver->lpDim); checkCode;
    
    // Allocate dy
    vecIter = (vec *) calloc(1, sizeof(vec));
    dsdpSolver->dy = vecIter;
    retcode = vec_init(vecIter); checkCode;
    retcode = vec_alloc(vecIter, m); checkCode;
    
    return retcode;
}

static DSDP_INT DSDPICheckData( HSDSolver *dsdpSolver ) {
    
    // Check whether problem data is alreay set up
    // Begin presolving after this step
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
    
    // s
    retcode = vec_free(dsdpSolver->s); checkCode;
    DSDP_FREE(dsdpSolver->s);
    
    // x
    retcode = vec_free(dsdpSolver->x); checkCode;
    DSDP_FREE(dsdpSolver->x);
    
    // asinv
    retcode = vec_free(dsdpSolver->asinv); checkCode;
    DSDP_FREE(dsdpSolver->asinv);
    
    // Msdp
    retcode = denseMatFree(dsdpSolver->Msdp);
    DSDP_FREE(dsdpSolver->Msdp);
    
    // u, d1, d2, d3, d4, yp
    retcode = vec_free(dsdpSolver->u ); checkCode;
    retcode = vec_free(dsdpSolver->d1); checkCode;
    retcode = vec_free(dsdpSolver->d12); checkCode;
    retcode = vec_free(dsdpSolver->d2); checkCode;
    retcode = vec_free(dsdpSolver->d3); checkCode;
    retcode = vec_free(dsdpSolver->d4); checkCode;
    retcode = vec_free(dsdpSolver->y); checkCode;
    
    DSDP_FREE(dsdpSolver->u );
    DSDP_FREE(dsdpSolver->d1);
    DSDP_FREE(dsdpSolver->d12);
    DSDP_FREE(dsdpSolver->d2);
    DSDP_FREE(dsdpSolver->d3);
    DSDP_FREE(dsdpSolver->d4);
    DSDP_FREE(dsdpSolver->y);

    // dS
    for (DSDP_INT i = 0; i < nblock; ++i) {
        retcode = spsMatFree(dsdpSolver->dS[i]); checkCode;
        DSDP_FREE(dsdpSolver->dS[i]);
    }
    
    DSDP_FREE(dsdpSolver->dS);
    
    // spaux
    for (DSDP_INT i = 0; i < nblock; ++i) {
        retcode = spsMatFree(dsdpSolver->spaux[i]); checkCode;
        DSDP_FREE(dsdpSolver->spaux[i]);
    }
          DSDP_FREE(dsdpSolver->spaux);
    
    // ds
    retcode = vec_free(dsdpSolver->ds);
    DSDP_FREE(dsdpSolver->ds);
    
    // dy
    retcode = vec_free(dsdpSolver->dy);
    DSDP_FREE(dsdpSolver->dy);
    
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
    dsdpSolver->param  = NULL;
    dsdpSolver->m      = 0;
    dsdpSolver->nBlock = 0;
    dsdpSolver->lpDim  = 0;
    dsdpSolver->mu     = 0.0;
    dsdpSolver->csinvrysinv = 0.0;
    dsdpSolver->rtk    = 0.0;
    dsdpSolver->Pnrm   = 0.0;
    dsdpSolver->tau    = 0.0;
    dsdpSolver->kappa  = 0.0;
    dsdpSolver->alpha  = 0.0;
    dsdpSolver->dtau   = 0.0;
    dsdpSolver->dkappa = 0.0;
    dsdpSolver->isXComputed = FALSE;
    
    dsdpSolver->iterA     = 0;
    dsdpSolver->iterB     = 0;
    dsdpSolver->insStatus = 0;
    dsdpSolver->solStatus = 0;
    
    return retcode;
}

static DSDP_INT DSDPIPresolve( HSDSolver *dsdpSolver ) {
    
    // Do presolve
    DSDP_INT retcode = DSDP_RETCODE_OK;
    
    assert( dsdpSolver->insStatus == DSDP_STATUS_SET );
    if ( dsdpSolver->insStatus != DSDP_STATUS_SET ) {
        error(etype, "Problem data is not set up. \n");
    }
    
    // Round 1: scale the primal pairs {b_i, A_ip} across different constraints
    dsdpSolver->pScaler = (vec *) calloc(1, sizeof(vec));
    retcode = vec_init(dsdpSolver->pScaler); checkCode;
    retcode = vec_alloc(dsdpSolver->pScaler, dsdpSolver->m); checkCode;
    
    double maxNrm     = 0.0;
    double minNrm     = 0.0;
    double tmpnrm     = 0.0;
    double pScalFact  = 0.0;
    sdpMat *cone      = NULL;
    
    for (DSDP_INT i = 0; i < dsdpSolver->m; ++i) {
        
        maxNrm = 0.0;
        minNrm = 0.0;
        
        for (DSDP_INT j = 0; j < dsdpSolver->nBlock; ++j) {
            getMatnrm(dsdpSolver, j, i, &tmpnrm);
            maxNrm = MAX(maxNrm, tmpnrm);
            minNrm = MIN(minNrm, tmpnrm);
        }
        
        pScalFact = maxNrm * minNrm;
        
        if (pScalFact < 1.2 && pScalFact > 0.8) {
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
        retcode = preRank1Rdc(cone); checkCode;
        // Primal coefficient scaling
        retcode = preSDPMatPScale(cone, dsdpSolver->pScaler); checkCode;
        // Dual coefficient scaling
        retcode = preSDPMatDScale(cone); checkCode;
    }
    
    // LP coefficient scaling
    if (dsdpSolver->isLPset) {
        retcode = preLPMatScale(dsdpSolver->lpData,
                                dsdpSolver->lpObj,
                                dsdpSolver->pScaler); checkCode;
    }
    
    // Count index and end presolving
    for (DSDP_INT j = 0; j < nblock; ++j) {
        cone = dsdpSolver->sdpData[j];
        retcode = getMatIdx(cone); checkCode;
    }
    
    return retcode;
}

static DSDP_INT DSDPIPostsolve( HSDSolver *dsdpSolver ) {
    
    DSDP_INT retcode = DSDP_RETCODE_OK;
    // TODO: Post-solver
    
    return retcode;
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
                            DSDP_INT  lpDim ) {
    
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
    
    if (dsdpSolver->verbosity) {
        printf("Dimension is successfully set. \n");
        printf("nSDPCones: "ID" "
               "nConstraints: "ID" "
               "LPDim: "ID" "
               "SDPDim: "ID". \n", nBlock, nConstrs, lpDim, sdpDim);
        
    }
    
    dsdpSolver->nBlock = nBlock;
    dsdpSolver->m      = nConstrs;
    dsdpSolver->lpDim  = lpDim;
    dsdpSolver->n      = sdpDim;
    
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
        error(etype, "SDP block is already set. \n");
    }
    
    sdpMat *data = dsdpSolver->sdpData[blockid];
    retcode = sdpMatSetDim(data, dsdpSolver->m,
                           coneSize, blockid); checkCode;
    retcode = sdpMatAlloc(data); checkCode;
    
    if (typehint) {
        retcode = sdpMatSetHint(data, typehint); checkCode;
    }
    
    retcode = sdpMatSetData(data, Asdpp, Asdpi, Asdpx); checkCode;
    
    if (dsdpSolver->verbosity) {
        printf("SDP block "ID" is set. \n", blockid);
    }
    
    dsdpSolver->isSDPset[blockid] = TRUE;
    DSDPICheckData(dsdpSolver);
    
    return retcode;
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
