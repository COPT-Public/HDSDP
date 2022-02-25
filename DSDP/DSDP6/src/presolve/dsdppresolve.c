#include "dsdppresolve.h"
#include "dsdputils.h"
#include "dsdpeigfact.h"

#ifdef SHOWALL
#undef SHOWALL
#endif

static char etype[] = "Presolving operations";

static DSDP_INT isDenseRank1Acc( dsMat *dataMat, DSDP_INT *isRank1 ) {
    // Detect if a dense matrix is rank one by directly computing the outer product
    // Slower but accurate
    DSDP_INT retcode = DSDP_RETCODE_OK;
    
    double *A    = dataMat->array;
    double *a    = NULL;
    DSDP_INT n   = dataMat->dim;
    DSDP_INT r1  = TRUE;
    DSDP_INT col = 0;
    DSDP_INT isNeg = FALSE;
    
    // Get the first column that contains non-zero elements
    for (DSDP_INT i = 0; i < n; ++i) {
        col = i;
        if (packIdx(A, n, i, i) != 0) {
            break;
        }
    }
    
    if (col == n - 1) {
        *isRank1 = FALSE;
        return retcode;
    }
    
    a = (double *) calloc(n, sizeof(double));
    
    double adiag = packIdx(A, n, col, col);
    
    if (adiag < 0) {
        isNeg = TRUE;
        adiag = sqrt(- adiag);
    } else {
        adiag = sqrt(adiag);
    }
    
    for (DSDP_INT i = col; i < n; ++i) {
        a[i] = packIdx(A, n, i, col) / adiag;
    }
    
    // Check if A = a * a' by computing ||A - a * a'||_F
    double *start = NULL;
    double err    = 0.0;
    double diff   = 0.0;
    DSDP_INT idx  = 0;
    
    if (isNeg) {
        for (DSDP_INT i = 0; i < n; ++i) {
            start = &A[idx];
            for (DSDP_INT j = 0; j < n - i; ++j) {
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
        for (DSDP_INT i = 0; i < n; ++i) {
            start = &A[idx];
            for (DSDP_INT j = 0; j < n - i; ++j) {
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
    
    DSDP_INT *Ap   = dataMat->p;
    DSDP_INT *Ai   = dataMat->i;
    double   *Ax   = dataMat->x;
    DSDP_INT n     = dataMat->dim;
    DSDP_INT col   = 0;
    DSDP_INT isNeg = FALSE;
    DSDP_INT nnz   = 0;
    double   err   = 0.0;
    double   diff  = 0.0;
    
    // First detect the first column containing nonzeros
    for (DSDP_INT i = 0; i < n; ++i) {
        col = i;
        if (Ap[i + 1] - Ap[i] > 0) {
            break;
        }
    }
    
    assert( col <= n - 1 ); // Otherwise the matrix is empty
    
    if (Ai[0] != col) {
        isR1 = FALSE;
        *isRank1 = isR1;
        return retcode;
    }
    
    double *a = NULL;
    double adiag = 0.0;
    a = (double *) calloc(n, sizeof(double));
    
    adiag = Ax[0];
        
    if (adiag < 0) {
        isNeg = TRUE;
        adiag = - sqrt(-adiag);
    } else {
        adiag = sqrt(adiag);
    }
    
    nnz = 0;
    // Get the sparse rank 1 matrix
    for (DSDP_INT j = Ap[col]; j < Ap[col + 1]; ++j) {
        // If the diagonal is zero but other rows contain non-zeros
        a[Ai[j]] = Ax[j] / adiag;
        nnz += (Ax[j] != 0);
    }
    
    if (dataMat->nnz != (DSDP_INT) (nnz + 1) * nnz / 2) {
        // Currently only nnz is checked and there is still chance of mismatch
        isR1 = FALSE;
    } else {
        // Ready to check rank-one property
        if (isNeg) {
            for (DSDP_INT i = 0; i < n; ++i) {
                for (DSDP_INT j = Ap[i]; j < Ap[i + 1]; ++j) {
                    diff = Ax[j] + a[i] * a[Ai[j]];
                    err += diff * diff;
                }
                if (err > 1e-06) {
                    isR1 = FALSE;
                    break;
                }
            }
        } else {
            for (DSDP_INT i = 0; i < n; ++i) {
                for (DSDP_INT j = Ap[i]; j < Ap[i + 1]; ++j) {
                    diff = Ax[j] - a[i] * a[Ai[j]];
                    err += diff * diff;
                }
                if (err > 1e-06) {
                    isR1 = FALSE;
                    break;
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
    
    assert( col != n - 1 ); // or it is a zero matrix
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
    
    DSDP_INT n     = dataMat->dim;
    DSDP_INT *Ap   = dataMat->p;
    DSDP_INT *Ai   = dataMat->i;
    double   *Ax   = dataMat->x;
    DSDP_INT col   = 0;
    double adiag   = 0.0;
    
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
        if (fabs(eigvals[i]) > threshold ) {
            r += 1;
        }
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
    
    // Matrix number must be correct
    assert( dataMat->nzeroMat + dataMat->ndenseMat \
           + dataMat->nrkMat + dataMat->nspsMat == m + 1 );
    
    void **matdata = dataMat->sdpData;
    
#ifdef SHOWALL
    printf("| Block %d before reduction: \n", blockid);
    printf("| Sparse: %-3d | Dense: %-3d | Rankk: %-3d | Zero: %-3d \n",
           nspsMat, ndsMat, nr1Mat, nzeroMat);
#endif
    // If all the matrices are already rank-one
    if (dataMat->nrkMat == m + 1) {
        return retcode;
    }
    
    dsMat  *dsdata  = NULL;
    spsMat *spsdata = NULL;
    r1Mat  *r1data  = NULL;
    rkMat  *rkdata  = NULL;
    
    // r1data used as buffer
    r1data  = (r1Mat *) calloc(1, sizeof(r1Mat)); checkCode;
    retcode = r1MatInit(r1data); checkCode;
    retcode = r1MatAlloc(r1data, n); checkCode;
    
    // Recall that C is located at the end of the array
    for (DSDP_INT i = 0; i < m + 1; ++i) {
        
        isRank1  = FALSE;
        isDense  = FALSE;
        isSparse = FALSE;
        
        switch (types[i]) {
            case MAT_TYPE_ZERO:
                break;
            case MAT_TYPE_DENSE:
                retcode = isDenseRank1Acc((dsMat *) matdata[i], &isRank1); checkCode;
                isDense = TRUE;
                break;
            case MAT_TYPE_SPARSE:
                retcode = isSparseRank1((spsMat *) matdata[i], &isRank1); checkCode;
                isSparse = TRUE;
                break;
            case MAT_TYPE_RANKK:
                break;
            default:
                error(etype, "Unknown matrix type. \n");
                break;
        }
        
        // If rank-one structure is detected,
        // free current structure and use rank-one structure
        if (isRank1) {
            
            rkdata = (rkMat *) calloc(1, sizeof(rkMat)); checkCode;
            retcode = rkMatInit(rkdata);
        
            r1data->sign = (double) isRank1;
            dataMat->nrkMat += 1;
            types[i] = MAT_TYPE_RANKK;
            
            if (isDense) {
                dsdata = matdata[i];
                retcode = extractR1fromDs(dsdata, r1data->x, isRank1);
                checkCode;
                rkMatAllocAndSetData(rkdata, n, 1, &r1data->sign, r1data->x);
                // Free and re-allocate
                retcode = denseMatFree(dsdata); checkCode;
                DSDP_FREE(dsdata);
                dataMat->ndenseMat -= 1;
                assert(dataMat->ndenseMat >= 0);
            }
            
            if (isSparse) {
                spsdata = matdata[i];
                retcode = extractR1fromSps(spsdata, r1data->x, isRank1);
                checkCode;
                rkMatAllocAndSetData(rkdata, n, 1, &r1data->sign, r1data->x);
                // Free and re-allocate
                retcode = spsMatFree(spsdata); checkCode;
                DSDP_FREE(spsdata);
                dataMat->nspsMat -= 1;
                assert(dataMat->nspsMat >= 0);
            }
            
            rkMatStoreOriginalData(rkdata, MAT_TYPE_RANKK, NULL);
            matdata[i] = (void *) rkdata;
        }
        
        assert( dataMat->nzeroMat + dataMat->ndenseMat \
               + dataMat->nrkMat + dataMat->nspsMat == m + 1 );
    }
    
#ifdef SHOWALL
    printf("| Block %d after reduction: \n", blockid);
    printf("| Sparse: %-3d | Dense: %-3d | Rankk(1): %-3d | Zero: %-3d \n",
           dataMat->nspsMat, dataMat->ndenseMat, dataMat->nrkMat,
           dataMat->nzeroMat);
#endif
    
    r1MatFree(r1data);
    DSDP_FREE(r1data);
    
    return retcode;
}

// TODO: Add special reduction for diagonal / tridiagonal / matrix of two elements
static DSDP_INT preRankkEvRdcBlock( sdpMat *dataMat, DSDPStats *stat ) {
    
    // Detect rank-k structure in SDP data
    DSDP_INT retcode = DSDP_RETCODE_OK;
    DSDP_INT *types = dataMat->types;
    DSDP_INT isDense  = FALSE, isSparse = FALSE, m = dataMat->dimy, n = dataMat->dimS;
    
    // Matrix number must be correct
    assert( dataMat->nzeroMat + dataMat->ndenseMat \
           + dataMat->nrkMat + dataMat->nspsMat == m + 1 );
    
    void **matdata = dataMat->sdpData;
    
#ifdef SHOWALL
    printf("| Block %d before reduction: \n", blockid);
    printf("| Sparse: %-3d | Dense: %-3d | Rankk: %-3d | Zero: %-3d \n",
           nspsMat, ndsMat, nrkMat, nzeroMat);
#endif
    
    if (dataMat->nrkMat == m + 1) {
        return retcode;
    }
    
    dsMat  *dsdata  = NULL;
    spsMat *spsdata = NULL;
    rkMat  *rkdata  = NULL;
    
    double *eigvals = (double *) calloc(n, sizeof(double));
    double *eigvecs = (double *) calloc(n * n, sizeof(double));
    double onenrm, Cnrm, Anrm;
    DSDP_INT rank = 0;
    
    DSDPGetStats(stat, STAT_ONE_NORM_C, &Cnrm);
    DSDPGetStats(stat, STAT_ONE_NORM_A, &Anrm);
    
    for (DSDP_INT i = 0; i < m + 1; ++i) {
        
        rank = n + 1; onenrm = 0.0;
        isDense  = FALSE; isSparse = FALSE;
        
        switch (types[i]) {
            case MAT_TYPE_ZERO:
                break;
            case MAT_TYPE_DENSE:
                dsdata = (dsMat *) matdata[i];
                retcode = denseMatOneNorm(dsdata, &onenrm);
                retcode = factorizeDenseData(dsdata, -(onenrm + 1e-03), eigvals, eigvecs);
                isDense = TRUE;
                // TODO: Add 1e-10 as a controllable parameter in the solver
                retcode = preGetRank(n, eigvals, 1e-10, &rank);
                break;
            case MAT_TYPE_SPARSE:
                spsdata = (spsMat *) matdata[i];
                if (spsdata->nnz > n * n) continue; // Why ?
                retcode = spsMatOneNorm(spsdata, &onenrm);
                retcode = factorizeSparseData(spsdata, -onenrm, eigvals, eigvecs);
                isSparse = TRUE;
                retcode = preGetRank(n, eigvals, 1e-10, &rank);
                break;
            case MAT_TYPE_RANKK:
                rkdata = (rkMat *) matdata[i];
                retcode = r1MatOneNorm(rkdata->data[0], &onenrm);
                break;
            default:
                error(etype, "Unknown matrix type. \n");
                break;
        }
        
        if (i < m) { Anrm += onenrm; } else { Cnrm += onenrm; }
        
        // Threshold for low-rank matrix
        if (rank <= n) { // rank = n + 1 if the matrix is already rank-one
            rkdata = (rkMat *) calloc(1, sizeof(rkMat)); checkCode;
            retcode = rkMatInit(rkdata);
            if (isDense) { dsdata = matdata[i]; denseMatStoreFactor(dsdata, rkdata); }
            if (isSparse) { spsdata = matdata[i]; spsMatStoreFactor(spsdata, rkdata); }
            rkMatAllocAndSelectData(rkdata, n, rank, 1e-10, eigvals, eigvecs);
        }
        
        assert( dataMat->nzeroMat + dataMat->ndenseMat \
               + dataMat->nrkMat + dataMat->nspsMat == m + 1 );

    }
    
    DSDPStatUpdate(stat, STAT_ONE_NORM_A, Anrm);
    DSDPStatUpdate(stat, STAT_ONE_NORM_C, Cnrm);
    
#ifdef SHOWALL
    printf("| Block %d after reduction: \n", blockid);
    printf("| Sparse: %-3d | Dense: %-3d | Rankk: %-3d | Zero: %-3d \n",
           dataMat->nspsMat, dataMat->ndenseMat, dataMat->nrkMat,
           dataMat->nzeroMat);
#endif
    
    DSDP_FREE(eigvals);
    DSDP_FREE(eigvecs);
    return retcode;
}

static DSDP_INT getBlockStatistic( sdpMat *sdpData ) {
    // Count the matrix index information in sdpData
    // This routine also checks validity of data
    
    DSDP_INT retcode  = DSDP_RETCODE_OK;
    DSDP_INT m        = sdpData->dimy;
    DSDP_INT *type    = sdpData->types;
    DSDP_INT nzeroMat = 0;
    DSDP_INT ndsMat   = 0;
    DSDP_INT nspsMat  = 0;
    DSDP_INT nr1Mat   = 0;
    
    // Ready to give index
    assert( (!sdpData->spsMatIdx) && (!sdpData->denseMatIdx) && (!sdpData->rkMatIdx) );
    sdpData->spsMatIdx   = (DSDP_INT *) calloc(sdpData->nspsMat, sizeof(DSDP_INT));
    sdpData->denseMatIdx = (DSDP_INT *) calloc(sdpData->ndenseMat, sizeof(DSDP_INT));
    sdpData->rkMatIdx    = (DSDP_INT *) calloc(sdpData->nrkMat, sizeof(DSDP_INT));
    
    for (DSDP_INT i = 0; i < m + 1; ++i) {
        switch (type[i]) {
            case MAT_TYPE_ZERO:
                nzeroMat += 1; break;
            case MAT_TYPE_DENSE:
                sdpData->denseMatIdx[ndsMat] = i; ndsMat += 1; break;
            case MAT_TYPE_SPARSE:
                sdpData->spsMatIdx[nspsMat] = i; nspsMat += 1; break;
            case MAT_TYPE_RANKK:
                sdpData->rkMatIdx[nr1Mat] = i; nr1Mat += 1; break;
            default:
                error(etype, "Unknown matrix type. \n"); break;
        }
    }
    
    return retcode;
}

static DSDP_INT preSDPMatgetPScaler( HSDSolver *dsdpSolver ) {
    
    // Compute the primal scaler for DSDP
    DSDP_INT retcode = DSDP_RETCODE_OK;
    DSDP_INT m = dsdpSolver->m;
    
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
            retcode = getMatFnorm(dsdpSolver, j, i, &nrm);
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
        
        bnrm += bval / dsdpSolver->pScaler->x[i];
    }
    
    DSDPStatUpdate(&dsdpSolver->dsdpStats, STAT_ONE_NORM_B, bnrm);
    
    return retcode;
}

static DSDP_INT preSDPMatPScale( HSDSolver *dsdpSolver ) {
    
    // Carry out primal scaling of SDP
    DSDP_INT retcode = DSDP_RETCODE_OK;
    assert( dsdpSolver->pScaler->x );
    
    double scaler = 0.0;
    for (DSDP_INT i = 0; i < dsdpSolver->m; ++i) {
        scaler = dsdpSolver->pScaler->x[i];
        assert( scaler > 0 );
        dsdpSolver->dObj->x[i] = dsdpSolver->dObj->x[i] / scaler;
        for (DSDP_INT j = 0; j < dsdpSolver->nBlock; ++j) {
            retcode = matRScale(dsdpSolver, j, i, scaler);
        }
        checkCode;
    }
    
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
            retcode = getMatFnorm(dsdpSolver, i, j, &nrm);
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
            retcode = matRScale(dsdpSolver, i, j, scaler);
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
    assert( sdpBlock->blockId == blockid );
    
    DSDP_INT dim     = sdpBlock->dimS;
    DSDP_INT nhash   = nsym(dim);
    DSDP_INT nnz     = 0;
    
    dsdpSolver->S[blockid]  = (spsMat *) calloc(1, sizeof(spsMat));
    dsdpSolver->dS[blockid] = (spsMat *) calloc(1, sizeof(spsMat));
    
    retcode = spsMatInit(dsdpSolver->S[blockid]); checkCode;
    retcode = spsMatInit(dsdpSolver->dS[blockid]); checkCode;
    
    spsMat *spsdata  = NULL;
    rkMat  *rkdata   = NULL;
    r1Mat  *r1data   = NULL;
    
    // Get hash table
    DSDP_INT *hash = NULL;
    hash = (DSDP_INT *) calloc(nhash, sizeof(DSDP_INT));
    
    DSDP_INT useDenseS = FALSE, isfirstNz = FALSE, tmp = 0;
    DSDP_INT *matIdx = sdpBlock->rkMatIdx;
    
    if (sdpBlock->ndenseMat > 0) {
        useDenseS = TRUE;
    } else {
        for (DSDP_INT i = 0; i < sdpBlock->nrkMat; ++i) {
            rkdata = (rkMat *) sdpBlock->sdpData[matIdx[i]];
            for (DSDP_INT j = 0; j < rkdata->rank; ++j) {
                r1data = rkdata->data[j];
                if (r1data->nnz > denseThresh * dim) {
                    useDenseS = TRUE;
                    break;
                }
            }
            if (useDenseS) {
                break;
            }
        }
    }
    
    // Dense matrix is not present
    if (!useDenseS) {
        // Sparse matrix
        matIdx = sdpBlock->spsMatIdx;
        for (DSDP_INT i = 0; i < sdpBlock->nspsMat; ++i) {
            spsdata = (spsMat *) sdpBlock->sdpData[matIdx[i]];
            
            if ((spsdata->i[0] == 0) && (spsdata->p[1] > 0)) {
                isfirstNz = TRUE;
            }
            
            for (DSDP_INT j = 0; j < dim; ++j) {
                for (DSDP_INT k = spsdata->p[j]; k < spsdata->p[j + 1]; ++k) {
                    
                    if (packIdx(hash, dim, spsdata->i[k], j) == 0) {
                        packIdx(hash, dim, spsdata->i[k], j) = 1;
                        nnz += 1;
                    }
                }
            }
            
            if (nnz > denseThresh * nhash) {
                useDenseS = TRUE;
                break;
            }
        }
        
        if (!useDenseS) {
            // Rank 1 matrix
            matIdx = sdpBlock->rkMatIdx;
            for (DSDP_INT i = 0; i < sdpBlock->nrkMat; ++i) {
                rkdata = (rkMat *) sdpBlock->sdpData[matIdx[i]];
                for (DSDP_INT r = 0; r < rkdata->rank; ++r) {
                    r1data = rkdata->data[r];
                    if (r1data->nzIdx[0] == 0) {
                        isfirstNz = TRUE;
                    }
                    
                    for (DSDP_INT j = 0; j < r1data->nnz; ++j) {
                        for (DSDP_INT k = 0; k <= j; ++k) {
                            if (packIdx(hash, dim, r1data->nzIdx[j], r1data->nzIdx[k]) == 0) {
                                packIdx(hash, dim, r1data->nzIdx[j], r1data->nzIdx[k]) = 1;
                                nnz += 1;
                            }
                        }
                        if (nnz > denseThresh * nhash) {
                            useDenseS = TRUE;
                            break;
                        }
                    }
                }
            }
        }
    }
    
    if (useDenseS) {
        nnz = nhash;
    } else {
        // Get hash table by the symbolic features
        DSDP_INT idx = 0;
        
        for (DSDP_INT i = 0; i < nhash; ++i) {
            if (hash[i]) {
                hash[i] = idx;
                idx += 1;
            }
        }
        
        assert( idx == nnz );
    }
    
    retcode = spsMatAllocData(dsdpSolver->S[blockid], dim, nnz); checkCodeFree;
    retcode = spsMatAllocData(dsdpSolver->dS[blockid], dim, nnz); checkCodeFree;
    retcode = spsMatAllocData(dsdpSolver->Scker[blockid], dim, nnz); checkCodeFree;
    
    if (isfirstNz && !useDenseS) {
        dsdpSolver->S[blockid]->i[tmp] = 0;
        dsdpSolver->S[blockid]->nzHash[tmp] = 0;
        tmp += 1;
    }
    
    for (DSDP_INT i = 0; i < dim; ++i) {
        for (DSDP_INT j = i; j < dim; ++j) {
            if (useDenseS) {
                dsdpSolver->S[blockid]->i[tmp] = j;
                tmp += 1;
            } else if (packIdx(hash, dim, j, i)) {
                dsdpSolver->S[blockid]->i[tmp] = j;
                dsdpSolver->S[blockid]->nzHash[tmp] = i;
                tmp += 1;
            }
        }
        dsdpSolver->S[blockid]->p[i + 1] = tmp;
    }
    
    if (useDenseS) {
        DSDP_FREE(hash);
    } else {
        dsdpSolver->symS[blockid] = hash;
    }
    
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
    memcpy(dsdpSolver->dS[blockid]->x,
           dsdpSolver->S[blockid]->x,
           sizeof(double) * tmp);
    memcpy(dsdpSolver->Scker[blockid]->x,
           dsdpSolver->S[blockid]->x,
           sizeof(double) * tmp);
    
    if (useDenseS) {
        DSDP_FREE(dsdpSolver->dS[blockid]->nzHash);
        DSDP_FREE(dsdpSolver->Scker[blockid]->nzHash);
        DSDP_FREE(dsdpSolver->S[blockid]->nzHash);
    } else {
        memcpy(dsdpSolver->dS[blockid]->nzHash,
               dsdpSolver->S[blockid]->nzHash,
               sizeof(DSDP_INT) * nnz);
        memcpy(dsdpSolver->Scker[blockid]->nzHash,
               dsdpSolver->S[blockid]->nzHash,
               sizeof(DSDP_INT) * nnz);
    }
    
#ifdef SHOWALL
        printf("Block "ID" goes through symbolic check. \n", blockid);
#endif
    
    return retcode;
    
clean_up:
    DSDP_FREE(hash);
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
    
    for (DSDP_INT i = 0; i < dsdpSolver->nBlock; ++i) {
        retcode = preRankkEvRdcBlock(dsdpSolver->sdpData[i], &dsdpSolver->dsdpStats);
        checkCode;
    }
    
    return retcode;
}

extern DSDP_INT preSDPPrimal( HSDSolver *dsdpSolver ) {
    // Do matrix coefficient scaling given preScaler for the primal
    DSDP_INT retcode = DSDP_RETCODE_OK;
    
    retcode = preSDPMatgetPScaler(dsdpSolver); checkCode;
    retcode = preSDPMatPScale(dsdpSolver); checkCode;
    
    return retcode;
}

extern DSDP_INT preSDPDual( HSDSolver *dsdpSolver ) {
    // Dual coefficient scaling
    DSDP_INT retcode = DSDP_RETCODE_OK;
    
    retcode = preSDPMatgetDScaler(dsdpSolver); checkCode;
    retcode = preSDPMatDScale(dsdpSolver); checkCode;
    
    return retcode;
}

extern DSDP_INT preLPMatScale( lpMat *lpData, vec *lpObj, vec *pScaler ) {
    
    // Do matrix coefficient scaling for LP
    DSDP_INT retcode = DSDP_RETCODE_OK;
    
    DSDP_INT m = lpData->dimy;
    DSDP_INT n = lpData->dims;
    
    assert( m == pScaler->dim );
    assert( n == lpObj->dim );
    
    if ( (m != pScaler->dim) || (n != pScaler->dim) ) {
        error(etype, "Presolver and problem dimension mismatch. \n");
    }
    
    DSDP_INT *Ap = lpData->lpdata->p;
    DSDP_INT *Ai = lpData->lpdata->i;
    double   *Ax = lpData->lpdata->x;
    double *scal = pScaler->x;
    
    // Dual coefficient scaling
    assert ( !lpData->xscale );
    lpData->xscale = (double *) calloc(n, sizeof(double));
    double *dscal  = lpData->xscale;
    
    double maxNrm = 0.0;
    double minNrm = 0.0;
    double tmp    = 0.0;
    
    // Primal scaling w.r.t. pScaler
    for (DSDP_INT i = 0; i < n; ++i) {
        
        for (DSDP_INT j = Ap[i]; j < Ap[i + 1]; ++j) {
            tmp = Ax[j] / scal[Ai[j]];
            Ax[j] = tmp;
            maxNrm = MAX(maxNrm, tmp);
            minNrm = MIN(minNrm, tmp);
        }
        
        tmp = maxNrm * minNrm;
        if ((tmp < 1.1) && (tmp > 0.9)) {
            dscal[i] = 1.0;
            break;
        }
        
        dscal[i] = sqrt(tmp);
        
        for (DSDP_INT j = Ap[i]; j < Ap[i + 1]; ++j) {
            Ax[j] = Ax[j] / dscal[i];
        }
    }

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
        retcode = preSDPgetSymbolic(dsdpSolver, i);
        checkCode;
    }
 
    return retcode;
}
