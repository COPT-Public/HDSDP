#include "dsdppresolve.h"
#include "dsdplapack.h"
#include "dsdputils.h"
#include "dsdplog.h"
#include "heurpool.h"
#include "sparseopts.h"
#include "denseopts.h"
#include "rank1opts.h"
#include "rankkopts.h"
#include "vec.h"
#include "speigs.h"

static char etype[] = "Presolving operations";

static void preGetRank( DSDP_INT n, double *eigvals, double threshold, DSDP_INT *rank ) {
    // Get low-rank approximation of data
    DSDP_INT r = 0, i;
    for (i = 0; i < n; ++i) { if (fabs(eigvals[i]) > threshold ) { r += 1; } }
    *rank = r;
}

static DSDP_INT preRank1RdcBlock( sdpMat *dataMat ) {
    // Go through the data structure to detect rank-one matrices
    DSDP_INT retcode = DSDP_RETCODE_OK;
    DSDP_INT *types = dataMat->types, isRank1 = FALSE, isDense = FALSE, isSparse = FALSE;
    DSDP_INT m = dataMat->dimy, n = dataMat->dimS;
    void **matdata = dataMat->sdpData;
    if (dataMat->nrkMat == m + 1) { return retcode; }
    dsMat *dsdata = NULL; spsMat *spsdata = NULL; r1Mat *r1data = NULL; rkMat *rkdata = NULL;
    // r1data used as buffer
    r1data = (r1Mat *) calloc(1, sizeof(r1Mat));
    r1MatInit(r1data); retcode = r1MatAlloc(r1data, n); checkCode;
    for (DSDP_INT i = 0; i < m + 1; ++i) {
        isRank1 = FALSE; isDense = FALSE; isSparse = FALSE;
        switch (types[i]) {
            case MAT_TYPE_ZERO: break;
            case MAT_TYPE_DENSE: dsr1check((dsMat *) matdata[i], &isRank1); isDense = TRUE; break;
            case MAT_TYPE_SPARSE: spsMatr1check((spsMat *) matdata[i], &isRank1); isSparse = TRUE; break;
            case MAT_TYPE_RANKK: break;
            default: error(etype, "Unknown matrix type. \n"); break;
        }
        // If rank-one structure is detected
        if (isRank1) {
            rkdata = (rkMat *) calloc(1, sizeof(rkMat));
            rkMatInit(rkdata); r1data->sign = (double) isRank1;
            dataMat->nrkMat += 1; types[i] = MAT_TYPE_RANKK;
            if (isDense) {
                dsdata = matdata[i]; spsMatr1extract(dsdata, r1data->x, isRank1);
                rkMatAllocAndSetData(rkdata, n, 1, &r1data->sign, r1data->x);
                denseMatFree(dsdata); DSDP_FREE(dsdata); dataMat->ndenseMat -= 1;
            }
            if (isSparse) {
                spsdata = matdata[i]; dsr1extract(spsdata, r1data->x, isRank1);
                rkMatAllocAndSetData(rkdata, n, 1, &r1data->sign, r1data->x);
                spsMatFree(spsdata); DSDP_FREE(spsdata); dataMat->nspsMat -= 1;
            }
            matdata[i] = (void *) rkdata;
        }
    }
    r1MatFree(r1data); DSDP_FREE(r1data); return retcode;
}

static DSDP_INT preRankkEvRdcBlock( sdpMat *dataMat, DSDPStats *stat, speigfac *eigfactor ) {
    // Detect rank-k structure in SDP data
    DSDP_INT retcode = DSDP_RETCODE_OK;
    DSDP_INT *types = dataMat->types, isDense = FALSE, isSparse = FALSE, m = dataMat->dimy, n = dataMat->dimS;
    // Matrix number must be correct
    void **matdata = dataMat->sdpData;
    if (dataMat->nrkMat == m + 1) { return retcode; }
    dsMat  *dsdata = NULL; spsMat *spsdata = NULL; rkMat *rkdata = NULL;
    double *eigvals = (double *) calloc(n, sizeof(double));
    double *eigvecs = (double *) calloc(n * n, sizeof(double));
    DSDP_INT rank = 0, status;
    for (DSDP_INT i = 0; i < m + 1; ++i) {
        rank = n + 1; isDense = FALSE; isSparse = FALSE;
        switch (types[i]) {
            case MAT_TYPE_ZERO: break;
            case MAT_TYPE_DENSE: rank = 100000; isDense = TRUE; break; // Do not factorize dense matrices
            case MAT_TYPE_SPARSE: spsdata = (spsMat *) matdata[i]; isSparse = TRUE;
                status = speigSpFactor(eigfactor, spsdata, eigvals, eigvecs);
                if (!status) { break; } preGetRank(n, eigvals, 1e-10, &rank); break;
            case MAT_TYPE_RANKK: break;
            default: error(etype, "Unknown matrix type. \n"); break;
        }
        // Threshold for low-rank matrix
        if (rank <= n) { // rank = n + 1 if the matrix is already rank-one
            rkdata = (rkMat *) calloc(1, sizeof(rkMat)); rkMatInit(rkdata);
            if (isDense) { dsdata = matdata[i]; denseMatStoreFactor(dsdata, rkdata); }
            if (isSparse) { spsdata = matdata[i]; spsMatStoreFactor(spsdata, rkdata); }
            rkMatAllocAndSelectData(rkdata, n, rank, 1e-10, eigvals, eigvecs);
        }
    }
    DSDP_FREE(eigvals); DSDP_FREE(eigvecs); return retcode;
}

static DSDP_INT getBlockStatistic( sdpMat *sdpData ) {
    // Count the matrix index information in sdpData
    DSDP_INT retcode = DSDP_RETCODE_OK, m = sdpData->dimy, *type = sdpData->types;
    DSDP_INT nzeroMat = 0, ndsMat = 0, nspsMat = 0, nr1Mat = 0;
    sdpData->spsMatIdx   = (DSDP_INT *) calloc(sdpData->nspsMat, sizeof(DSDP_INT));
    sdpData->denseMatIdx = (DSDP_INT *) calloc(sdpData->ndenseMat, sizeof(DSDP_INT));
    sdpData->rkMatIdx    = (DSDP_INT *) calloc(sdpData->nrkMat, sizeof(DSDP_INT));
    for (DSDP_INT i = 0; i < m + 1; ++i) {
        switch (type[i]) {
            case MAT_TYPE_ZERO  : nzeroMat += 1; break;
            case MAT_TYPE_DENSE : sdpData->denseMatIdx[ndsMat] = i; ndsMat += 1; break;
            case MAT_TYPE_SPARSE: sdpData->spsMatIdx[nspsMat] = i; nspsMat += 1; break;
            case MAT_TYPE_RANKK : sdpData->rkMatIdx[nr1Mat] = i; nr1Mat += 1; break;
            default             : error(etype, "Unknown matrix type. \n"); break;
        }
    }
    return retcode;
}

static DSDP_INT preSDPMatgetPScaler( HSDSolver *dsdpSolver ) {
    // Compute the primal scaler for DSDP
    DSDP_INT retcode = DSDP_RETCODE_OK;
    DSDPStatUpdate(&dsdpSolver->dsdpStats, STAT_ONE_NORM_B,
                   vec_onenorm(dsdpSolver->dObj));
    return retcode;
}

static DSDP_INT preSDPMatgetDScaler( HSDSolver *dsdpSolver ) {
    // Compute the dual scaler for DSDP
    DSDP_INT retcode = DSDP_RETCODE_OK;
    DSDP_INT m = dsdpSolver->m, i, j;
    double nrm = 0.0, maxnrm = 0.0, minnrm = 0.0;
    for (i = 0; i < dsdpSolver->nBlock; ++i) {
        minnrm = 1.0;
        for (j = 0; j < m; ++j) {
            getMatFnorm(dsdpSolver, i, j, &nrm);
            if (nrm > 0.0) {
                minnrm = MIN(minnrm, nrm); maxnrm = MAX(maxnrm, nrm);
            }
        }
        if (maxnrm == 0.0) { error(etype, "Empty block detected. \n"); }
        if (fabs(sqrt(maxnrm * minnrm) - 1.0) < 0.1) {
            dsdpSolver->sdpData[i]->scaler = 1.0;
        } else {
            dsdpSolver->sdpData[i]->scaler = sqrt(minnrm * maxnrm);
        }
    }
    return retcode;
}

static DSDP_INT preSDPMatDScale( HSDSolver *dsdpSolver ) {
    // Carry out primal scaling of SDP
    DSDP_INT retcode = DSDP_RETCODE_OK;
    double scaler = 0.0;
    for (DSDP_INT i = 0; i < dsdpSolver->nBlock; ++i) {
        scaler = dsdpSolver->sdpData[i]->scaler;
        for (DSDP_INT j = 0; j < dsdpSolver->m; ++j) {
            matRScale(dsdpSolver, i, j, scaler);
        }
    }
    return retcode;
}

static DSDP_INT preSDPgetSymbolic( HSDSolver *dsdpSolver, DSDP_INT blockid ) {
    // Get the symbolic structure for dual matrix S/dS in block k and allocate memory
    // By the time this symbolic phase is called, the problem data should have gone
    // through presolving reduction
    DSDP_INT retcode = DSDP_RETCODE_OK;
    sdpMat *sdpBlock = dsdpSolver->sdpData[blockid];
    DSDP_INT dim = sdpBlock->dimS, nhash = nsym(dim), nnz = 0, i, j, k;
    dsdpSolver->S[blockid]  = (spsMat *) calloc(1, sizeof(spsMat));
    dsdpSolver->dS[blockid] = (spsMat *) calloc(1, sizeof(spsMat));
    spsMatInit(dsdpSolver->S[blockid]); spsMatInit(dsdpSolver->dS[blockid]);
    spsMat *spsdata = NULL; rkMat *rkdata = NULL;
    // Get mapping
    DSDP_INT *hash = NULL;
    hash = (DSDP_INT *) calloc(nhash, sizeof(DSDP_INT));
    DSDP_INT useDenseS = FALSE, isfirstNz = FALSE, tmp = 0, *matIdx = sdpBlock->rkMatIdx;
    if (sdpBlock->ndenseMat > 0) {
        useDenseS = TRUE;
    } else {
        for (i = 0; i < sdpBlock->nrkMat; ++i) {
            rkdata = (rkMat *) sdpBlock->sdpData[matIdx[i]];
            rkMatCheckSparsity(rkdata, &useDenseS, denseThresh);
            if (useDenseS) { break; }
        }
    }
    // Dense matrix is not present
    if (!useDenseS) {
        // Sparse matrix
        matIdx = sdpBlock->spsMatIdx;
        for (i = 0; i < sdpBlock->nspsMat; ++i) {
            spsdata = (spsMat *) sdpBlock->sdpData[matIdx[i]];
            spsMatGetSymbolic(spsdata, hash, &isfirstNz, &nnz);
            if (nnz > denseThresh * nhash) { useDenseS = TRUE; break; }
        }
        
        if (!useDenseS) {
            matIdx = sdpBlock->rkMatIdx;
            for (i = 0; i < sdpBlock->nrkMat; ++i) {
                rkdata = (rkMat *) sdpBlock->sdpData[matIdx[i]];
                rkMatGetSymbolic(rkdata, hash, &isfirstNz, &nnz);
                if (nnz > denseThresh * nhash) {
                    useDenseS = TRUE; break;
                }
            }
        }
    }
    if (useDenseS) {
        nnz = nhash;
    } else {
        // Get hash table by the symbolic features
        for (i = k = 0; i < nhash; ++i) {
            if (hash[i]) { hash[i] = k; k += 1;}
        }
    }
    retcode = spsMatAllocData(dsdpSolver->S[blockid], dim, nnz); checkCodeFree;
    retcode = spsMatAllocData(dsdpSolver->dS[blockid], dim, nnz); checkCodeFree;
    retcode = spsMatAllocData(dsdpSolver->Scker[blockid], dim, nnz); checkCodeFree;
    spsNominalLinkSinv(dsdpSolver->S[blockid], dsdpSolver->M->Sinv[blockid]);
    spsNominalLinkSinv(dsdpSolver->Scker[blockid], dsdpSolver->M->Sinv[blockid]);
    spsNominalLinkSinv(dsdpSolver->dS[blockid], dsdpSolver->M->Sinv[blockid]);
    
    if (useDenseS) {
        DSDP_FREE(hash);
    } else {
        if (isfirstNz) {
            dsdpSolver->S[blockid]->i[tmp] = 0; dsdpSolver->S[blockid]->cidx[tmp] = 0;
            tmp += 1;
        }
        
        for (i = 0; i < dim; ++i) {
            for (j = i; j < dim; ++j) {
                if (packIdx(hash, dim, j, i)) {
                    dsdpSolver->S[blockid]->i[tmp] = j;
                    dsdpSolver->S[blockid]->cidx[tmp] = i;
                    tmp += 1;
                }
            }
            dsdpSolver->S[blockid]->p[i + 1] = tmp;
        }
        dsdpSolver->symS[blockid] = hash;
        memcpy(dsdpSolver->dS[blockid]->p,
               dsdpSolver->S[blockid]->p,
               sizeof(DSDP_INT) * (dim + 1));
        memcpy(dsdpSolver->Scker[blockid]->p,
               dsdpSolver->S[blockid]->p,
               sizeof(DSDP_INT) * (dim + 1));
        memcpy(dsdpSolver->dS[blockid]->i,
               dsdpSolver->S[blockid]->i,
               sizeof(DSDP_INT) * tmp);
        memcpy(dsdpSolver->Scker[blockid]->i,
               dsdpSolver->S[blockid]->i,
               sizeof(DSDP_INT) * tmp);
        memcpy(dsdpSolver->dS[blockid]->cidx,
               dsdpSolver->S[blockid]->cidx,
               sizeof(DSDP_INT) * nnz);
        memcpy(dsdpSolver->Scker[blockid]->cidx,
               dsdpSolver->S[blockid]->cidx,
               sizeof(DSDP_INT) * nnz);
    }
    
    return retcode;
    
clean_up:
    DSDP_FREE(hash);
    return retcode;
}

static void preNoPIntDetect( HSDSolver *dsdpSolver ) {
    // Check if there is no primal interior point
    if (dsdpSolver->nBlock > 1) { return; }
    sdpMat *data = dsdpSolver->sdpData[0]; rkMat *rkdata = NULL;
    DSDP_INT idx, nopint = FALSE;;
    if (data->nrkMat == 0) { return; }
    for (DSDP_INT i = 0; i < data->nrkMat; ++i) {
        idx = data->rkMatIdx[i]; if (idx == data->dimy) { continue; }
        rkdata = data->sdpData[idx];
        if (dsdpSolver->dObj->x[idx] / rkdata->data[0]->sign <= 1e-03) {
            nopint = TRUE; break;
        }
    }
    DSDPStatUpdate(&dsdpSolver->dsdpStats, STAT_NO_PINTERIOR, nopint);
}

static void preNoDIntDetect( HSDSolver *dsdpSolver ) {
    // Check if there is no dual interior implied by linear constraints
    if (!dsdpSolver->isLPset || dsdpSolver->lpDim <= 2) { return; } 
    lpMat *lpdata = dsdpSolver->lpData;
    DSDP_INT *Ap = lpdata->Ap, i, j, nnz, nnz2, m = lpdata->dims, n = lpdata->dimy, m2 = m / 2;
    if (m % 2 != 0) { return; } double *c = dsdpSolver->lpObj->x;
    for (i = 0; i < m2; ++i) { if (c[i] != -c[i + m2]) { return; } }
    double *Ax = lpdata->Ax;
    for (i = 0; i < n; ++i) {
        nnz = Ap[i + 1] - Ap[i]; nnz2 = nnz / 2;
        if (nnz % 2 != 0) { return; }
        for (j = 0; j < nnz2; ++j) {
            if (Ax[Ap[i] + j] != -Ax[Ap[i] + j + nnz2]) { return; }
        }
    }
    DSDPStatUpdate(&dsdpSolver->dsdpStats, STAT_NO_DINTERIOR, TRUE);
}

static void preImpXBoundDetect( HSDSolver *dsdpSolver ) {
    // Check if constraints imply trace(X, I) = T.
    if (dsdpSolver->nBlock > 1 || dsdpSolver->isLPset) { return; }
    sdpMat *data = dsdpSolver->sdpData[0];
    spsMat *spsdata = NULL; r1Mat *r1Mat = NULL;
    DSDP_INT impXbound = FALSE, type = 0, idx; double boundX = 0.0;
    double *aux = dsdpSolver->M->schurAux; memset(aux, 0, sizeof(double) * data->dimS);
    if (data->nspsMat == data->dimy + 1) { type = 1; }
    if (data->nspsMat == data->dimy &&
        data->types[data->dimy] != MAT_TYPE_SPARSE) { type = 1; }
    if (data->nrkMat == data->dimy &&
        data->types[data->dimy] != MAT_TYPE_RANKK) { type = 2; }
    if (!type) { return; }
    if (type == 1) {
        for (DSDP_INT i = 0; i < data->dimy; ++i) {
            spsdata = data->sdpData[i];
            if (spsMatIsDiag(spsdata)) {
                impXbound = TRUE;
                boundX = spsMatGetXbound(spsdata, dsdpSolver->dObj); break;
            }
        }
    } else {
        for (DSDP_INT i = 0; i < data->dimy; ++i) {
            r1Mat = ((rkMat *) data->sdpData[i])->data[0];
            if (r1Mat->nnz == 1) {
                idx = r1Mat->nzIdx[0];
                aux[idx] = dsdpSolver->dObj->x[idx] / r1Mat->x[r1Mat->nzIdx[0]];
            }
        }
        impXbound = TRUE;
        for (DSDP_INT i = 0; i < data->dimS; ++i) {
            if (aux[i]) { boundX += aux[i]; }
            else { impXbound = FALSE; break; }
        }
    }
    if (impXbound) {
        DSDPStatUpdate(&dsdpSolver->dsdpStats, STAT_IMP_BOUNDX, boundX);
    }
}

static void preImpYBoundDetect( HSDSolver *dsdpSolver ) {
    //  Check if LP constraints imply u <= y <= l
    if (!dsdpSolver->isLPset) { return; }
    
    lpMat *lpdata = dsdpSolver->lpData;
    DSDP_INT *Ap = lpdata->Ap, *Ai = lpdata->Ai, i, j, impy = TRUE, impylb = FALSE, impyub = FALSE;
    // If y is only restricted by a' * y = c
    if (dsdpSolver->lpDim <= 2) { return; }
    
    double *Ax = lpdata->Ax, yubound = 0.0, ylbound = 0.0, *c = dsdpSolver->lpObj->x;
    double *ylb = dsdpSolver->d1->x, *yub = dsdpSolver->d2->x, tmp;
    for (i = 0; i < lpdata->dimy; ++i) {
        for (j = Ap[i]; j < Ap[i + 1]; ++j) {
            if (Ap[i + 1] - Ap[i] > 2) { impy = FALSE; break; }
            if (Ax[j] > 0.0) {
                if (yub[i]) { impy = FALSE; break; } impyub = TRUE;
                tmp = c[Ai[j]] / Ax[j]; yub[i] = MAX(yub[i], tmp);
            } else {
                if (ylb[i]) { impy = FALSE; break; } impylb = TRUE;
                tmp = c[Ai[j]] / Ax[j]; ylb[i] = MIN(ylb[i], tmp);
            }
        }
        if (!impy) { break; }
    }
    if (impy) {
        if (impyub) {
            for (i = 0; i < lpdata->dimy; ++i) { yubound = MAX(yubound, yub[i]); }
            yubound = (yubound == 0.0) ? 1.0 : yubound;
            DSDPStatUpdate(&dsdpSolver->dsdpStats, STAT_IMP_UBOUNDY, yubound);
        }
        if (impylb) {
            for (i = 0; i < lpdata->dimy; ++i) { ylbound = MIN(ylbound, ylb[i]); }
            ylbound = (ylbound == 0.0) ? -1.0 : ylbound;
            DSDPStatUpdate(&dsdpSolver->dsdpStats, STAT_IMP_LBOUNDY, ylbound);
        }
    }
}

extern DSDP_INT DSDPPrepareMAssembler( HSDSolver *dsdpSolver ) {
    // Initialize the internal Schur matrix structure
    DSDP_INT retcode = DSDP_RETCODE_OK;
    
    symM *M = dsdpSolver->M;
    symSchurMatInit(M); symSchurMatSetDim(M, dsdpSolver->m, dsdpSolver->nBlock);
    retcode = symSchurMatAlloc(M);
    
    if (retcode != DSDP_RETCODE_OK) {
        printf("| Failed to allocate memory for symbolic schur matrix. \n");
        return retcode;
    }
    
    symSchurMatRegister(M, dsdpSolver->S,
                        dsdpSolver->dsaux,
                        dsdpSolver->sdpData,
                        dsdpSolver->Msdp,
                        dsdpSolver->asinv,
                        dsdpSolver->asinvrysinv,
                        dsdpSolver->u,
                        &dsdpSolver->csinvrysinv,
                        &dsdpSolver->csinv,
                        &dsdpSolver->csinvcsinv,
                        &dsdpSolver->rysinv,
                        &dsdpSolver->Ry,
                        dsdpSolver->rkaux,
                        &dsdpSolver->eventMonitor[EVENT_IN_PHASE_A],
                        &dsdpSolver->eventMonitor[EVENT_HSD_UPDATE],
                        &dsdpSolver->param->intParams[INT_PARAM_GOLDSEARCH],
                        &dsdpSolver->isLPset);
    return retcode;
}

extern DSDP_INT preRank1Rdc( HSDSolver *dsdpSolver ) {
    // Do rank 1 reduction
    DSDP_INT retcode = DSDP_RETCODE_OK;
    
    DSDP_INT sps = 0, ds = 0, r1 = 0, zero = 0;
    DSDPStats *stat = &dsdpSolver->dsdpStats;
    
    for (DSDP_INT i = 0; i < dsdpSolver->nBlock; ++i) {
        retcode = preRank1RdcBlock(dsdpSolver->sdpData[i]);
        sps  += dsdpSolver->sdpData[i]->nspsMat;
        ds   += dsdpSolver->sdpData[i]->ndenseMat;
        r1   += dsdpSolver->sdpData[i]->nrkMat;
        zero += dsdpSolver->sdpData[i]->nzeroMat;
        checkCode;
    }
    
    DSDPStatUpdate(stat, STAT_NUM_SPARSE_MAT, sps);
    DSDPStatUpdate(stat, STAT_NUM_DENSE_MAT, ds);
    DSDPStatUpdate(stat, STAT_NUM_RONE_MAT, r1);
    DSDPStatUpdate(stat, STAT_NUM_ZERO_MAT, zero);

    return retcode;
}

extern DSDP_INT preRankkRdc( HSDSolver *dsdpSolver ) {
    // Do rank k reduction
    DSDP_INT retcode = DSDP_RETCODE_OK;
    
    DSDP_INT largestblock = 0;
    for (DSDP_INT i = 0; i < dsdpSolver->nBlock; ++i) {
        largestblock = MAX(largestblock, dsdpSolver->sdpData[i]->dimS);
    }
    
    DSDPStatUpdate(&dsdpSolver->dsdpStats, STAT_LARGEST_BLOCK, largestblock);
    speigfac *eigfactor = (speigfac *) calloc(1, sizeof(speigfac));
    speigInit(eigfactor);
    retcode = speigAlloc(eigfactor, largestblock);
    
    for (DSDP_INT i = 0; i < dsdpSolver->nBlock; ++i) {
        retcode = preRankkEvRdcBlock(dsdpSolver->sdpData[i], &dsdpSolver->dsdpStats, eigfactor);
        checkCode;
    }
    
    speigFree(eigfactor); DSDP_FREE(eigfactor);
    return retcode;
}

extern DSDP_INT preSDPPrimal( HSDSolver *dsdpSolver ) {
    // Do matrix coefficient scaling given preScaler for the primal
    DSDP_INT retcode = DSDP_RETCODE_OK;
    retcode = preSDPMatgetPScaler(dsdpSolver); checkCode;
    return retcode;
}

extern DSDP_INT preSDPMatCScale( HSDSolver *dsdpSolver ) {
    // Scale C by its one norm or decide that the problem is a feasibility problem
    DSDP_INT retcode = DSDP_RETCODE_OK;
    DSDPStatUpdate(&dsdpSolver->dsdpStats, STAT_ONE_NORM_A,
                   DSDPConic( COPS_GET_A_ONE_NORM ) (dsdpSolver));
    DSDPStatUpdate(&dsdpSolver->dsdpStats, STAT_ONE_NORM_C,
                   DSDPConic( COPS_GET_C_ONE_NORM ) (dsdpSolver));
    DSDPGetStats(&dsdpSolver->dsdpStats, STAT_ONE_NORM_C,
                 &dsdpSolver->cScaler);
    DSDP_HEURS( adjustCScaler ) (dsdpSolver);
    DSDPConic( COPS_DO_C_SCALE ) (dsdpSolver);
    return retcode;
}

extern DSDP_INT preSDPDual( HSDSolver *dsdpSolver ) {
    // Dual coefficient scaling. Not in use
    DSDP_INT retcode = DSDP_RETCODE_OK;
    retcode = preSDPMatgetDScaler(dsdpSolver); checkCode;
    retcode = preSDPMatDScale(dsdpSolver); checkCode;
    return retcode;
}

extern DSDP_INT getMatIdx( HSDSolver *dsdpSolver ) {
    
    DSDP_INT retcode = DSDP_RETCODE_OK;
    for (DSDP_INT i = 0; i < dsdpSolver->nBlock; ++i) {
        retcode = getBlockStatistic(dsdpSolver->sdpData[i]);
        checkCode;
    }
    return retcode;
}

extern DSDP_INT preSymbolic( HSDSolver *dsdpSolver ) {
    // Get the symbolic structure of S
    DSDP_INT retcode = DSDP_RETCODE_OK;
    for (DSDP_INT i = 0; i < dsdpSolver->nBlock; ++i) {
        retcode = preSDPgetSymbolic(dsdpSolver, i);
        if (retcode != DSDP_RETCODE_OK) {
            printf("| Failed to get symbolic structure. \n");
            break;
        }
    }
    return retcode;
}

extern DSDP_INT preStructureDetect( HSDSolver *dsdpSolver ) {
    // Detect four speccial structures in the dual formulation
    DSDP_INT retcode = DSDP_RETCODE_OK;
    
    if (dsdpSolver->nBlock < 10) {
        preNoPIntDetect(dsdpSolver);
        preNoDIntDetect(dsdpSolver);
        preImpXBoundDetect(dsdpSolver);
        preImpYBoundDetect(dsdpSolver);
    }
    
    double nopint, nodint, impX, impyub, impylb;
    DSDPGetStats(&dsdpSolver->dsdpStats, STAT_NO_PINTERIOR, &nopint);
    DSDPGetStats(&dsdpSolver->dsdpStats, STAT_NO_DINTERIOR, &nodint);
    DSDPGetStats(&dsdpSolver->dsdpStats, STAT_IMP_BOUNDX, &impX);
    DSDPGetStats(&dsdpSolver->dsdpStats, STAT_IMP_UBOUNDY, &impyub);
    DSDPGetStats(&dsdpSolver->dsdpStats, STAT_IMP_LBOUNDY, &impylb);
    
    if (nopint || nodint || impX || impyub || impylb) {
        printf("| - Special structures found \n");
        if (nopint) {
            printf("|    %s : %s \n", "No primal interior: tr(X * aa') -> 0", "PRelax penalty tightened");
        }
        if (nodint) {
            printf("|    %s : %s \n", "No dual interior: A' * y = c", "PRelax penalty tightened");
        }
        if (impX) {
            printf("|    %s = %6.2e : %s \n", "tr(X)", impX, "Bound of X fixed");
        }
        if (impyub) {
            printf("|    %s <= u : %s \n", "y", "PRelax penalty adjusted");
        }
        if (impylb) {
            printf("|    %s >= l : %s \n", "y", "PRelax penalty adjusted");
        }
    } else {
        printf("| - No special structure is available \n");
    }
    
    return retcode;
}
