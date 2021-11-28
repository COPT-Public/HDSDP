#include "dsdppresolve.h"

/*
 TODO: Matrix rank-one structure detector
 TODO: Matrix coefficient scaling
*/

static char *etype = "Presolving operations";

static DSDP_INT isDenseRank1InAcc( dsMat *dataMat, DSDP_INT *isRank1 ) {
    // Check if a dense matrix is rank-one
    DSDP_INT retcode = DSDP_RETCODE_OK;
    double *A    = dataMat->array;
    DSDP_INT n   = dataMat->dim;
    DSDP_INT r1  = TRUE;
    DSDP_INT col = 0;
    double benchCol = 0.0;
    double scaleCol = 0.0;
    double benchCol2 = 0.0;
    double scaleCol2 = 0.0;
    
    // Get the first column that contains non-zero elements
    for (DSDP_INT i = 0; i < n; ++i) {
        if (packIdx(A, n, i, i) > 0) {
            col = i;
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
            if (fabs(benchCol * scaleCol2 - benchCol2 * scaleCol) > 1e-04 * MAX(1, benchCol)) {
                r1 = FALSE;
                break;
            }
        }
    }
    
    *isRank1 = r1;
    
    return retcode;
}

static DSDP_INT isDenseRank1Acc( dsMat *dataMat, DSDP_INT *isRank1 ) {
    // Detect if a dense matrix is rank one by directly computing the  outer product
    DSDP_INT retcode = DSDP_RETCODE_OK;
    
    double *A    = dataMat->array;
    double *a    = NULL;
    DSDP_INT n   = dataMat->dim;
    DSDP_INT r1  = TRUE;
    DSDP_INT col = 0;
    
    // Get the first column that contains non-zero elements
    for (DSDP_INT i = 0; i < n; ++i) {
        if (packIdx(A, n, i, i) > 0) {
            col = i;
            break;
        }
    }
    
    assert( col != n - 1 ); // or it is a zero matrix
    a = (double *) calloc(n, sizeof(double));
    double adiag = sqrt(packIdx(A, n, col, col));
    
    for (DSDP_INT i = col; i < n; ++i) {
        a[i] = packIdx(A, n, i, col) / adiag;
    }
    
    // Check if A = a * a' by computing ||A - a * a'||_F
    double *start = NULL;
    double err    = 0.0;
    double diff   = 0.0;
    DSDP_INT idx  = 0;
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
    
    *isRank1 = r1;
    DSDP_FREE(a);
    
    return retcode;
}

static DSDP_INT isSparseRank1( spsMat *dataMat, DSDP_INT *isRank1 ) {
    // TODO: Check if a sparse matrix is rank-one
    DSDP_INT retcode = DSDP_RETCODE_OK;
    
    
    return retcode;
}

extern DSDP_INT preRank1Rdc( sdpMat *dataMat ) {
    // Detect rank-one structure in SDP data
    DSDP_INT retcode = DSDP_RETCODE_OK;
    
    
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
    
    if (dataMat->scaler < 1.2 && dataMat->scaler > 0.8) {
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

extern DSDP_INT preLPMatScale ( lpMat *lpData, vec *lpObj ) {
    // Do matrix coefficient scaling for LP
    DSDP_INT retcode = DSDP_RETCODE_OK;
    
    
    return retcode;
}
