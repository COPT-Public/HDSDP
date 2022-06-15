#include "dsdppresolve.h"
#include "dsdputils.h"
#include "symschur.h"
#include "heurpool.h"
#include "dsdpeigfact.h"
#include "speigs.h"

#ifdef SHOWALL
#undef SHOWALL
#endif

static char etype[] = "Presolving operations";

static DSDP_INT isDenseRank1Acc( dsMat *dataMat, DSDP_INT *isRank1 ) {
    // Detect if a dense matrix is rank one by directly computing the outer product
    // Slower but accurate
    DSDP_INT retcode = DSDP_RETCODE_OK;
    
    double *A = dataMat->array, *a = NULL;
    DSDP_INT n = dataMat->dim, i, j, r1 = TRUE, col = 0, isNeg = FALSE;
    
    // Get the first column that contains non-zero elements
    for (i = 0; i < n; ++i) {
        if (packIdx(A, n, i, i) != 0) {
            break;
        }
    }
    
    if (i == n) {*isRank1 = FALSE; return retcode;}
    if (i >= n - 1 && !packIdx(A, n, i, i)) {
        *isRank1 = FALSE;
        return retcode;
    }
    
    col = i;
    
    a = (double *) calloc(n, sizeof(double));
    
    double adiag = packIdx(A, n, col, col);
    
    if (adiag < 0) {
        isNeg = TRUE;
        adiag = sqrt(- adiag);
    } else {
        adiag = sqrt(adiag);
    }
    
    for (i = col; i < n; ++i) {
        a[i] = packIdx(A, n, i, col) / adiag;
    }
    
    // Check if A = a * a' by computing ||A - a * a'||_F
    double *start = NULL;
    double err    = 0.0;
    double diff   = 0.0;
    DSDP_INT idx  = 0;
    
    if (isNeg) {
        for (i = 0; i < n; ++i) {
            start = &A[idx];
            for (j = 0; j < n - i; ++j) {
                diff = start[j] + a[i] * a[i + j];
                err += diff * diff;
            }
            idx += n - i;
            if (err > 1e-08) {
                r1 = FALSE;
                break;
            }
        }
    } else {
        for (i = 0; i < n; ++i) {
            start = &A[idx];
            for (j = 0; j < n - i; ++j) {
                diff = start[j] - a[i] * a[i + j];
                err += diff * diff;
            }
            idx += n - i;
            if (err > 1e-08) {
                r1 = FALSE;
                break;
            }
        }
    }
    
    if (r1) {
        *isRank1 = (DSDP_INT) (1 - 2 * isNeg);
    } else {
        *isRank1 = FALSE;
    }
    
    DSDP_FREE(a);
    
    return retcode;
}

static DSDP_INT isSparseRank1( spsMat *dataMat, DSDP_INT *isRank1 ) {
    // Check if a sparse matrix is rank-one
    DSDP_INT retcode = DSDP_RETCODE_OK;
    DSDP_INT isR1 = TRUE;
    
    DSDP_INT *Ap = dataMat->p, *Ai = dataMat->i;
    double *Ax = dataMat->x, err = 0.0, diff = 0.0;
    DSDP_INT n = dataMat->dim, i, j, col = 0, isNeg = FALSE, nnz = 0;
    // First detect the first column containing nonzeros
    for (DSDP_INT i = 0; i < n; ++i) {
        col = i;
        if (Ap[i + 1] - Ap[i] > 0) {
            break;
        }
    }
    
    assert( col <= n - 1 ); // Otherwise the matrix is empty
    
    if (Ai[0] != col) {
        isR1 = FALSE; *isRank1 = isR1;
        return retcode;
    }
    
    double *a = NULL, adiag = 0.0;
    a = (double *) calloc(n, sizeof(double));
    
    adiag = Ax[0];
    if (adiag < 0) {
        isNeg = TRUE; adiag = - sqrt(-adiag);
    } else {
        adiag = sqrt(adiag);
    }
    
    // Get the sparse rank 1 matrix
    for (i = Ap[col]; i < Ap[col + 1]; ++i) {
        // If the diagonal is zero but other rows contain non-zeros
        a[Ai[i]] = Ax[i] / adiag;
        nnz += (Ax[i] != 0);
    }
    
    if (dataMat->nnz != nsym(nnz)) {
        isR1 = FALSE;
    }
    
    if (isR1) {
        for (i = col + 1; i < n; ++i) {
            if (Ap[i] > Ap[col + 1] && Ap[i] < dataMat->nnz) {
                if (Ai[Ap[i]] < col) {
                    isR1 = FALSE; break;
                }
            }
        }
    }

    if (isR1) {
        // Ready to check rank-one property
        if (isNeg) {
            for (i = 0; i < n; ++i) {
                for (j = Ap[i]; j < Ap[i + 1]; ++j) {
                    diff = Ax[j] + a[i] * a[Ai[j]];
                    err += fabs(diff);
                }
                if (err > 1e-10) {
                    isR1 = FALSE; break;
                }
            }
        } else {
            for (i = 0; i < n; ++i) {
                for (j = Ap[i]; j < Ap[i + 1]; ++j) {
                    diff = Ax[j] - a[i] * a[Ai[j]];
                    err += fabs(diff);
                }
                
                if (err > 1e-10) {
                    isR1 = FALSE; break;
                }
            }
        }
    }
    
    DSDP_FREE(a);
    
    if (isR1) {
        *isRank1 = (DSDP_INT) (1 - 2 * isNeg);
    } else {
        *isRank1 = FALSE;
    }
    
    return retcode;
}

static DSDP_INT extractR1fromDs( dsMat *dataMat, double *a, DSDP_INT isNeg ) {
    // Extract the rank 1 data from dense data structure
    DSDP_INT retcode = DSDP_RETCODE_OK;
    
    double *A    = dataMat->array;
    DSDP_INT n   = dataMat->dim;
    DSDP_INT col = 0;
    
    // Get the first column that contains non-zero elements
    for (DSDP_INT i = 0; i < n; ++i) {
        col = i;
        if (packIdx(A, n, i, i) != 0) {
            break;
        }
    }
    
    double adiag = packIdx(A, n, col, col);
    
    if (isNeg == -1) {
        adiag = sqrt(- adiag);
    } else {
        adiag = sqrt(adiag);
    }
    
    for (DSDP_INT i = col; i < n; ++i) {
        a[i] = packIdx(A, n, i, col) / adiag;
    }
    
    return retcode;
}

static DSDP_INT extractR1fromSps( spsMat *dataMat, double *a, DSDP_INT isNeg ) {
    // Extract the rank 1 data from sparse data structure
    DSDP_INT retcode = DSDP_RETCODE_OK;
    
    DSDP_INT n = dataMat->dim, col = 0, *Ap = dataMat->p, *Ai = dataMat->i;
    double *Ax  = dataMat->x, adiag = 0.0;
    
    memset(a, 0, sizeof(double) * n);
    
    for (DSDP_INT i = 0; i < n; ++i) {
        col = i;
        if (Ap[i + 1] - Ap[i] > 0) {
            break;
        }
    }
    
    assert( col <= n - 1 ); // Otherwise the matrix is empty
    
    if (isNeg == -1) {
        adiag = sqrt(- Ax[0]);
    } else {
        adiag = sqrt(Ax[0]);
    }
    
    if (adiag != adiag) {
        error(etype, "NAN encountered when extracting rank-1 vector. \n");
    }
    
    // Get the sparse rank 1 matrix
    for (DSDP_INT j = Ap[col]; j < Ap[col + 1]; ++j) {
        // If the diagonal is zero but other rows contain non-zeros
        a[Ai[j]] = Ax[j] / adiag;
    }
    
    return retcode;
}

static DSDP_INT preGetRank( DSDP_INT n, double *eigvals, double threshold, DSDP_INT *rank ) {
    // Get low-rank approximation of data
    DSDP_INT r = 0;
    for (DSDP_INT i = 0; i < n; ++i) {
        if (fabs(eigvals[i]) > threshold ) { r += 1; }
    }
    *rank = r;
    return DSDP_RETCODE_OK;
}

static DSDP_INT preRank1RdcBlock( sdpMat *dataMat ) {
    
    // Detect rank-one structure in SDP data
    DSDP_INT retcode  = DSDP_RETCODE_OK;
    DSDP_INT *types   = dataMat->types;
    DSDP_INT isRank1  = FALSE,
             isDense  = FALSE,
             isSparse = FALSE;
    
    DSDP_INT m = dataMat->dimy, n = dataMat->dimS;
    void **matdata = dataMat->sdpData;
    // If all the matrices are already rank-one
    if (dataMat->nrkMat == m + 1) { return retcode; }
    dsMat  *dsdata  = NULL; spsMat *spsdata = NULL;
    r1Mat  *r1data  = NULL; rkMat  *rkdata  = NULL;
    
    // r1data used as buffer
    r1data  = (r1Mat *) calloc(1, sizeof(r1Mat)); checkCode;
    r1MatInit(r1data); r1MatAlloc(r1data, n); checkCode;
    
    // Recall that C is located at the end of the array
    for (DSDP_INT i = 0; i < m + 1; ++i) {
        isRank1 = FALSE; isDense = FALSE; isSparse = FALSE;
        switch (types[i]) {
            case MAT_TYPE_ZERO: break;
            case MAT_TYPE_DENSE:
                isDenseRank1Acc((dsMat *) matdata[i], &isRank1); checkCode;
                isDense = TRUE; break;
            case MAT_TYPE_SPARSE:
                isSparseRank1((spsMat *) matdata[i], &isRank1); checkCode;
                isSparse = TRUE; break;
            case MAT_TYPE_RANKK: break;
            default: error(etype, "Unknown matrix type. \n"); break;
        }
        
        // If rank-one structure is detected,
        // free current structure and use rank-one structure
        if (isRank1) {
            rkdata = (rkMat *) calloc(1, sizeof(rkMat));
            rkMatInit(rkdata); r1data->sign = (double) isRank1;
            dataMat->nrkMat += 1; types[i] = MAT_TYPE_RANKK;
            if (isDense) {
                dsdata = matdata[i]; extractR1fromDs(dsdata, r1data->x, isRank1);
                rkMatAllocAndSetData(rkdata, n, 1, &r1data->sign, r1data->x);
                retcode = denseMatFree(dsdata);
                DSDP_FREE(dsdata); dataMat->ndenseMat -= 1;
            }
            if (isSparse) {
                spsdata = matdata[i]; extractR1fromSps(spsdata, r1data->x, isRank1);
                rkMatAllocAndSetData(rkdata, n, 1, &r1data->sign, r1data->x);
                spsMatFree(spsdata); DSDP_FREE(spsdata);
                dataMat->nspsMat -= 1;
            }
            matdata[i] = (void *) rkdata;
        }
    }
    
    r1MatFree(r1data); DSDP_FREE(r1data); return retcode;
}

static DSDP_INT preRankkEvRdcBlock( sdpMat *dataMat, DSDPStats *stat, speigfac *eigfactor ) {
    
    // Detect rank-k structure in SDP data
    DSDP_INT retcode = DSDP_RETCODE_OK;
    DSDP_INT *types = dataMat->types;
    DSDP_INT isDense  = FALSE, isSparse = FALSE, m = dataMat->dimy, n = dataMat->dimS;
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
            case MAT_TYPE_DENSE:
                rank = 100000; dsdata = (dsMat *) matdata[i];
                isDense = TRUE; break;
            case MAT_TYPE_SPARSE:
                spsdata = (spsMat *) matdata[i]; isSparse = TRUE;
                status = speigSfac(eigfactor, spsdata, eigvals, eigvecs);
                if (!status) {
                    break;
                }
                preGetRank(n, eigvals, 1e-10, &rank); break;
            case MAT_TYPE_RANKK: break;
            default: error(etype, "Unknown matrix type. \n"); break;
        }
        // Threshold for low-rank matrix
        if (rank <= n) { // rank = n + 1 if the matrix is already rank-one
            rkdata = (rkMat *) calloc(1, sizeof(rkMat)); checkCode;
            retcode = rkMatInit(rkdata);
            if (isDense) { dsdata = matdata[i]; denseMatStoreFactor(dsdata, rkdata); }
            if (isSparse) { spsdata = matdata[i]; spsMatStoreFactor(spsdata, rkdata); }
            rkMatAllocAndSelectData(rkdata, n, rank, 1e-10, eigvals, eigvecs);
        }
    }
    DSDP_FREE(eigvals); DSDP_FREE(eigvecs); return retcode;
}

static DSDP_INT getBlockStatistic( sdpMat *sdpData ) {
    // Count the matrix index information in sdpData
    // This routine also checks validity of data
    DSDP_INT retcode  = DSDP_RETCODE_OK;
    DSDP_INT m = sdpData->dimy, *type = sdpData->types;
    DSDP_INT nzeroMat = 0, ndsMat = 0, nspsMat = 0, nr1Mat = 0;
    
    // Ready to give index
    assert( (!sdpData->spsMatIdx) && (!sdpData->denseMatIdx) && (!sdpData->rkMatIdx) );
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
    DSDP_INT m = dsdpSolver->m;
    DSDPStatUpdate(&dsdpSolver->dsdpStats, STAT_ONE_NORM_B,
                   vec_onenorm(dsdpSolver->dObj));
    return retcode;
    
    dsdpSolver->pScaler = (vec *) calloc(1, sizeof(vec));
    vec_init(dsdpSolver->pScaler);
    vec_alloc(dsdpSolver->pScaler, m);
    
    double nrm = 0.0, maxnrm = 0.0, minnrm = 0.0, bnrm = 0.0, bval;
    
    for (DSDP_INT i = 0; i < m; ++i) {
        
        maxnrm = 0.0;
        bval = fabs(dsdpSolver->dObj->x[i]);
        
        if (bval > 0.0) {
            minnrm = bval;
            maxnrm = minnrm;
        } else {
            dsdpSolver->pScaler->x[i] = 1.0;
            continue;
        }
                
        for (DSDP_INT j = 0; j < dsdpSolver->nBlock; ++j) {
            getMatFnorm(dsdpSolver, j, i, &nrm);
            if (nrm > 0) {
                minnrm = MIN(minnrm, nrm);
                maxnrm = MAX(maxnrm, nrm);
            }
        }
        
        if (maxnrm == 0.0) {
            error(etype, "Empty row detected. \n");
        } else {
            if (fabs(maxnrm * minnrm - 1.0) < 0.1) {
                dsdpSolver->pScaler->x[i] = 1.0;
            } else {
                dsdpSolver->pScaler->x[i] = sqrt(maxnrm * minnrm);
            }
        }
        
        bnrm += bval; // dsdpSolver->pScaler->x[i];
    }
    
    DSDPStatUpdate(&dsdpSolver->dsdpStats, STAT_ONE_NORM_B, bnrm);
    
    return retcode;
}

static DSDP_INT preSDPMatgetDScaler( HSDSolver *dsdpSolver ) {
    // Compute the dual scaler for DSDP
    DSDP_INT retcode = DSDP_RETCODE_OK;
    
    DSDP_INT m = dsdpSolver->m;
    double nrm    = 0.0;
    double maxnrm = 0.0;
    double minnrm = 0.0;
    
    for (DSDP_INT i = 0; i < dsdpSolver->nBlock; ++i) {
        minnrm = 1.0;
        for (DSDP_INT j = 0; j < m; ++j) {
            getMatFnorm(dsdpSolver, i, j, &nrm);
            if (nrm > 0.0) {
                minnrm = MIN(minnrm, nrm);
                maxnrm = MAX(maxnrm, nrm);
            }
        }
        
        if (maxnrm == 0.0) {
            error(etype, "Empty block detected. \n");
        }
        
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
        assert( scaler > 0 );
        for (DSDP_INT j = 0; j < dsdpSolver->m; ++j) {
            matRScale(dsdpSolver, i, j, scaler);
        }
        checkCode;
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
    
    // Get hash table
    DSDP_INT *hash = NULL;
    hash = (DSDP_INT *) calloc(nhash, sizeof(DSDP_INT));
    DSDP_INT useDenseS = FALSE, isfirstNz = FALSE, tmp = 0;
    DSDP_INT *matIdx = sdpBlock->rkMatIdx;
    
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
            dsdpSolver->S[blockid]->i[tmp] = 0;
            dsdpSolver->S[blockid]->cidx[tmp] = 0;
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

extern DSDP_INT DSDPPrepareMAssembler( HSDSolver *dsdpSolver ) {
    // Initialize the internal Schur matrix structure
    DSDP_INT retcode = DSDP_RETCODE_OK;
    
    DSDPSymSchur *M = dsdpSolver->M;
    retcode = symSchurMatInit(M);
    retcode = symSchurMatSetDim(M, dsdpSolver->m, dsdpSolver->nBlock);
    retcode = symSchurMatAlloc(M);
    retcode = symSchurMatRegister(M, dsdpSolver->S,
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
    DSDPStatUpdate(stat, STAT_NUM_RANKONE_MAT, r1);
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
    DSDP_INT nblock = dsdpSolver->nBlock;
    
    for (DSDP_INT i = 0; i < nblock; ++i) {
        preSDPgetSymbolic(dsdpSolver, i);
    }
 
    return retcode;
}
