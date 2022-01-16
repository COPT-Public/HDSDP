#include "dsdppresolve.h"

/*
 TODO: Matrix rank-one structure detector
 TODO: Matrix coefficient scaling
*/

static char *etype = "Presolving operations";

static DSDP_INT isDenseRank1InAcc( dsMat *dataMat, DSDP_INT *isRank1 ) {
    // Check if a dense matrix is rank-one
    // This version is fast but not accurate due to potential numerical error
    DSDP_INT retcode = DSDP_RETCODE_OK;
    double *A    = dataMat->array;
    DSDP_INT n   = dataMat->dim;
    DSDP_INT r1  = TRUE;
    DSDP_INT col = 0;
    
    double benchCol  = 0.0;
    double scaleCol  = 0.0;
    double benchCol2 = 0.0;
    double scaleCol2 = 0.0;
    double diff      = 0.0;
    
    // Get the first column that contains non-zero elements
    for (DSDP_INT i = 0; i < n; ++i) {
        col = i;
        if (packIdx(A, n, i, i) != 0) {
            break;
        }
    }
    
    assert( col != n - 1 ); // or it is a zero matrix
    
    // Check the scaling coefficient
    for (DSDP_INT i = col + 1; i < n; ++i) {
        scaleCol = packIdx(A, n, i, col);
        for (DSDP_INT j = col; j < n; ++j) {
            benchCol2 = packIdx(A, n, j, col);
            if (i <= j) {
                scaleCol2 = packIdx(A, n, j, i);
            } else {
                scaleCol2 = packIdx(A, n, i, j);
            }
            diff = benchCol * scaleCol2 - benchCol2 * scaleCol;
            if (fabs(diff) > 1e-04 * MAX(1, benchCol)) {
                r1 = FALSE;
                break;
            }
        }
    }
    
    *isRank1 = r1;
    
    return retcode;
}

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
    
    assert( col != n - 1 ); // or it is a zero matrix
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
    double   err   = 0.0;
    double   diff  = 0.0;
    
    // First detect the first column containing nonzeros
    for (DSDP_INT i = 0; i < n; ++i) {
        col = i;
        if (Ap[i + 1] - Ap[i] > 0) {
            break;
        }
    }
    
    assert( col < n - 1 ); // Otherwise the matrix is empty
    
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
    
    // Get the sparse rank 1 matrix
    for (DSDP_INT j = Ap[col]; j < Ap[col + 1]; ++j) {
        // If the diagonal is zero but other rows contain non-zeros
        a[Ai[j]] = Ax[j] / adiag;
    }
    
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
    DSDP_FREE(a);
    
    if (isR1) {
        *isRank1 = (DSDP_INT) (1 - 2 * isNeg);
    } else {
        *isRank1 = FALSE;
    }
    
    return retcode;
}

static DSDP_INT extractR1fromDs ( dsMat *dataMat, double *a, DSDP_INT isNeg ) {
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
    
    if (isNeg) {
        adiag = sqrt(- adiag);
    } else {
        adiag = sqrt(adiag);
    }
    
    for (DSDP_INT i = col; i < n; ++i) {
        a[i] = packIdx(A, n, i, col) / adiag;
    }
    
    return retcode;
}

static DSDP_INT extractR1fromSps ( spsMat *dataMat, double *a, DSDP_INT isNeg ) {
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
    
    assert( col < n - 1 ); // Otherwise the matrix is empty
    
    if (isNeg) {
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

extern DSDP_INT preRank1Rdc( sdpMat *dataMat ) {
    // Detect rank-one structure in SDP data, some block
    
    DSDP_INT retcode  = DSDP_RETCODE_OK;
    DSDP_INT *types   = dataMat->types;
    DSDP_INT ndsMat   = dataMat->ndenseMat;
    DSDP_INT nspsMat  = dataMat->nspsMat;
    DSDP_INT nr1Mat   = dataMat->nr1Mat;
    DSDP_INT nzeroMat = dataMat->nzeroMat;
    DSDP_INT blockid  = dataMat->blockId;
    DSDP_INT isRank1  = FALSE;
    DSDP_INT isDense  = FALSE;
    DSDP_INT isSparse = FALSE;
    
    DSDP_INT m = dataMat->dimy;
    DSDP_INT n = dataMat->dimS;
    
    // Matrix number must be correct
    assert( dataMat->nzeroMat + dataMat->ndenseMat \
           + dataMat->nr1Mat + dataMat->nspsMat == m + 1 );
    
    void **matdata = dataMat->sdpData;
    
    printf("Block "ID" Before reduction: \n", blockid);
    printf("nSparse: "ID"  nDense: "ID"  nRank1: "ID" nZero: "ID" \n",
           nspsMat, ndsMat, nr1Mat, nzeroMat);
    
    // If all the matrices are already rank-one
    if (nr1Mat == m) {
        return retcode;
    }
    
    dsMat  *dsdata  = NULL;
    spsMat *spsdata = NULL;
    r1Mat  *r1data  = NULL;
    
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
            case MAT_TYPE_RANK1:
                break;
            default:
                error(etype, "Unknown matrix type. \n");
                break;
        }
        
        // If rank-one structure is detected,
        // free current structure and use rank-one structure
        if (isRank1) {
            
            r1data  = (r1Mat *) calloc(1, sizeof(r1Mat)); checkCode;
            retcode = r1MatInit(r1data); checkCode;
            retcode = r1MatAlloc(r1data, n); checkCode;
            
            r1data->sign = (double) isRank1;
            dataMat->nr1Mat += 1;
            types[i] = MAT_TYPE_RANK1;
            
            if (isDense) {
                dsdata = matdata[i];
                retcode = extractR1fromDs(dsdata, r1data->x, isRank1);
                checkCode;
                retcode = r1MatCountNnz(r1data); checkCode;
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
                retcode = r1MatCountNnz(r1data); checkCode;
                retcode = spsMatFree(spsdata); checkCode;
                DSDP_FREE(spsdata);
                dataMat->nspsMat -= 1;
                assert(dataMat->nspsMat >= 0);
            }
            
            matdata[i] = (void *) r1data;
        }
        
        assert( dataMat->nzeroMat + dataMat->ndenseMat \
               + dataMat->nr1Mat + dataMat->nspsMat == m + 1 );
    }
    
    printf("Block "ID" after reduction: \n", blockid);
    printf("nSparse: "ID"  nDense: "ID"  nRank1: "ID" nZero: "ID" \n",
           dataMat->nspsMat, dataMat->ndenseMat, dataMat->nr1Mat,
           dataMat->nzeroMat);
    
    return retcode;
}

extern DSDP_INT preSDPMatPScale( sdpMat *dataMat, vec *pScaler ) {
    // Do matrix coefficient scaling given preScaler for the primal
    DSDP_INT retcode = DSDP_RETCODE_OK;
    
    assert( dataMat->dimy == pScaler->dim );
    if (dataMat->dimy != pScaler->dim) {
        error(etype, "Presolving vector does not match "
              "the number of matrices in the block. \n");
    }
    
    // Here the scaling does not consider C
    for (DSDP_INT i = 0; i < dataMat->dimy; ++i) {
        
        if (pScaler->x[i] == 1.0) {
            continue;
        }
        
        switch (dataMat->types[i]) {
            case MAT_TYPE_ZERO:
                break;
            case MAT_TYPE_DENSE:
                retcode = denseMatRscale((dsMat *) dataMat->sdpData[i],
                                         pScaler->x[i]); checkCode;
                break;
            case MAT_TYPE_SPARSE:
                retcode = spsMatRscale((spsMat *) dataMat->sdpData[i],
                                       pScaler->x[i]); checkCode;
                break;
            case MAT_TYPE_RANK1:
                retcode = r1MatRscale((r1Mat *) dataMat->sdpData[i],
                                      pScaler->x[i]); checkCode;
                checkCode;
                break;
            default:
                error(etype, "Unknown matrix type. \n");
                break;
        }
    }
    
    return retcode;
}

extern DSDP_INT preSDPMatDScale( sdpMat *dataMat ) {
    // Do matrix coefficient scaling for one SDP block
    DSDP_INT retcode = DSDP_RETCODE_OK;
    DSDP_INT m = dataMat->dimy;
    
    double maxNrm     = 0.0;
    double minNrm     = 0.0;
    double tmpnrm     = 0.0;
    
    assert ( dataMat->scaler == 0.0 );
    
    // Dual normalization considers C, while primal does not
    for (DSDP_INT i = 0; i < m + 1; ++i) {
        
        switch (dataMat->types[i]) {
            case MAT_TYPE_ZERO:
                break;
            case MAT_TYPE_DENSE:
                retcode = denseMatFnorm((dsMat *) dataMat->sdpData[i], &tmpnrm);
                checkCode;
                break;
            case MAT_TYPE_SPARSE:
                retcode = spsMatFnorm((spsMat *) dataMat->sdpData[i] , &tmpnrm);
                checkCode;
                break;
            case MAT_TYPE_RANK1:
                retcode = r1MatFnorm((r1Mat *) dataMat->sdpData[i], &tmpnrm);
                checkCode;
                break;
            default:
                error(etype, "Unknown matrix type. \n");
                break;
        }
        
        maxNrm = MAX(maxNrm, tmpnrm);
        minNrm = MIN(minNrm, tmpnrm);
    }
    
    dataMat->scaler = sqrt(maxNrm * minNrm);
    
    if (dataMat->scaler < 1.05 && dataMat->scaler > 0.95) {
        dataMat->scaler = 1.0;
    } else {
        
        // Do scaling
        for (DSDP_INT i = 0; i < m; ++i) {
            
            tmpnrm = dataMat->scaler;
            switch (dataMat->types[i]) {
                case MAT_TYPE_ZERO:
                    break;
                case MAT_TYPE_DENSE:
                    retcode = denseMatRscale((dsMat *) dataMat->sdpData[i],
                                             tmpnrm); checkCode;
                    break;
                case MAT_TYPE_SPARSE:
                    retcode = spsMatRscale((spsMat *) dataMat->sdpData[i],
                                            tmpnrm); checkCode;
                    break;
                case MAT_TYPE_RANK1:
                    retcode = r1MatRscale((r1Mat *) dataMat->sdpData[i],
                                          tmpnrm); checkCode;
                    checkCode;
                    break;
                default:
                    error(etype, "Unknown matrix type. \n");
                    break;
            }
        }
    }
    
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

extern DSDP_INT getMatIdx( sdpMat *sdpData ) {
    // Count the matrix index information in sdpData
    // This routine also checks validity of data
    
    DSDP_INT retcode  = DSDP_RETCODE_OK;
    DSDP_INT blockid  = sdpData->blockId;
    DSDP_INT m        = sdpData->dimy;
    DSDP_INT *type    = sdpData->types;
    DSDP_INT nzeroMat = 0;
    DSDP_INT ndsMat   = 0;
    DSDP_INT nspsMat  = 0;
    DSDP_INT nr1Mat   = 0;
    
    // Ready to give index
    assert( (!sdpData->spsMatIdx) && (!sdpData->denseMatIdx) && (!sdpData->r1MatIdx) );
    sdpData->spsMatIdx   = (DSDP_INT *) calloc(nspsMat, sizeof(DSDP_INT));
    sdpData->denseMatIdx = (DSDP_INT *) calloc(ndsMat, sizeof(DSDP_INT));
    sdpData->r1MatIdx    = (DSDP_INT *) calloc(nr1Mat, sizeof(DSDP_INT));
    
    printf("- Cone "ID"\n", blockid);
    
    for (DSDP_INT i = 0; i < m + 1; ++i) {
        switch (type[i]) {
            case MAT_TYPE_ZERO:
                nzeroMat += 1;
                break;
            case MAT_TYPE_DENSE:
                sdpData->denseMatIdx[ndsMat] = i;
                ndsMat += 1;
                break;
            case MAT_TYPE_SPARSE:
                sdpData->spsMatIdx[nspsMat] = i;
                nspsMat += 1;
                break;
            case MAT_TYPE_RANK1:
                sdpData->r1MatIdx[nr1Mat] = i;
                nr1Mat += 1;
                break;
            default:
                error(etype, "Unknown matrix type. \n");
                break;
        }
    }
    
    assert( nzeroMat == sdpData->nzeroMat );
    assert( ndsMat   == sdpData->ndenseMat );
    assert( nspsMat  == sdpData->nspsMat );
    assert( nr1Mat   == sdpData->nr1Mat );
    
    return retcode;
}
