#include <string.h>
#include "dsdplapack.h"
#include "dsdpdata.h"
#include "dsdpsort.h"
#include "dsdpcs.h"
#include "sparseopts.h"
#include "denseopts.h"
#include "rank1opts.h"
#include "rankkopts.h"
#include "vec.h"

// Implement the data interface for DSDP-HSD Solver

static char *etype = "Data interface";

/* Index sorting */
static DSDP_INT idxsort( DSDP_INT *Ai, double *Ax, DSDP_INT nnz ) {
    // Sort the Ai and Ax by Ai in ascending order
    dsdpSort(Ax, Ai, 0, nnz - 1);
    return DSDP_RETCODE_OK;
}

/* LP internal methods */
static DSDP_INT lpMatIAlloc( lpMat *lpData, DSDP_INT nnz ) {
    // Allocate memory for internal LP matrix
    DSDP_INT retcode = DSDP_RETCODE_OK;
    assert( (lpData->dimy > 0) && (lpData->dims > 0) );
    
    if ((lpData->dimy <= 0) || (lpData->dims <= 0)) {
        error(etype, "Incorrect size for LP data matrix for allocation. \n");
    }
    lpData->Ap = (DSDP_INT *) calloc(lpData->dimy + 1, sizeof(DSDP_INT));
    lpData->Ai = (DSDP_INT *) calloc(nnz, sizeof(DSDP_INT));
    lpData->Ax = (double   *) calloc(nnz, sizeof(double));
    
    if (!lpData->Ap || !lpData->Ai || !lpData->Ax) {
        printf("| Failed to allocate memory for internal LP data. \n");
        retcode = DSDP_RETCODE_FAILED; return retcode;
    }
    
    return retcode;
}

/* SDP internal methods */
static DSDP_INT sdpMatIAllocByType( sdpMat *sdpData, DSDP_INT k, DSDP_INT *Ai,
                                    double *Ax, DSDP_INT nnz ) {
    // Automatically detect the type of a matrix and allocate memory for it
    // Ai identifies row indices of the lower-triangular packed format of the k th matrix
    // Ax stores the data, nnz specifies the number of nonzeros
    DSDP_INT retcode = DSDP_RETCODE_OK;
    DSDP_INT n = sdpData->dimS;
    void *userdata = NULL;
    
    // Check sparsity
    if (nnz == 0) {
        sdpData->nzeroMat += 1;
        sdpData->types[k] = MAT_TYPE_ZERO;
        
    } else if (((nnz <= denseThresh * nsym(n)) && (sdpData->types[k] == MAT_TYPE_UNKNOWN)) ||
               (sdpData->types[k] == MAT_TYPE_SPARSE)) {
        
        DSDP_INT ordered = (nnz < 1000) ? FALSE : checkIsOrdered(Ai, nnz);
        if (!ordered) {
            idxsort(Ai, Ax, nnz);
        }
        
        // Sparse
        sdpData->types[k] = MAT_TYPE_SPARSE; sdpData->nspsMat += 1;
        spsMat *data = (spsMat *) calloc(1, sizeof(spsMat)); spsMatInit(data);
        retcode = spsMatAllocData(data, n, nnz); checkCode;
        
        data->p[n] = nnz; data->nnz = nnz;
        
        if (n > 10000) {
            DSDP_INT i;
            for (i = 0; i < n - n % 4; i+=4) {
                data->p[i + 1] = nnz; data->p[i + 2] = nnz;
                data->p[i + 3] = nnz; data->p[i + 4] = nnz;
            }
            for (i = n - n % 4; i < n; ++i) {
                data->p[i + 1] = nnz;
            }
        } else {
            for (DSDP_INT i = 0; i < n; ++i) {
                data->p[i + 1] = nnz;
            }
        }
        
        memcpy(data->x, Ax, sizeof(double) * nnz);
        DSDP_INT rowidx = 0, where = 0, colnnz = 0, idxthresh = n;
        data->p[0] = 0;
        
        for (DSDP_INT i = 0; i < nnz; ++i) {
            rowidx = Ai[i];
            while (rowidx >= idxthresh) {
                where += 1; data->p[where] = colnnz;
                idxthresh += n - where;
            }
            colnnz += 1; data->i[i] = rowidx - idxthresh + n;
            data->cidx[i] = where; // Record column index
        }
        userdata = (void *) data;
        
    } else {
        // Dense
        sdpData->types[k] = MAT_TYPE_DENSE;
        sdpData->ndenseMat += 1;
        dsMat *data = NULL;
        data = (dsMat *) calloc(1, sizeof(dsMat));
        
        denseMatInit(data); checkCode;
        retcode = denseMatAlloc(data, n, FALSE); checkCode;
        
        for (DSDP_INT i = 0; i < nnz - nnz % 4; i+=4) {
            data->array[Ai[i    ]] = Ax[i    ];
            data->array[Ai[i + 1]] = Ax[i + 1];
            data->array[Ai[i + 2]] = Ax[i + 2];
            data->array[Ai[i + 3]] = Ax[i + 3];
        }
        
        for (DSDP_INT i = nnz - nnz % 4; i < nnz; ++i) {
            data->array[Ai[i]] = Ax[i];
        }
        userdata = (void *) data;
    }
    
    sdpData->sdpData[k] = userdata;
    
    return retcode;
}

/* LP public methods */
extern void lpMatInit( lpMat *lpData ) {
    // Initialize LP data
    lpData->dims = 0; lpData->dimy = 0; lpData->nnz = 0;
    lpData->Ap = NULL; lpData->Ai = NULL; lpData->Ax = NULL;
    lpData->xscale = NULL;
}

extern void lpMatSetDim( lpMat *lpData, DSDP_INT dimy, DSDP_INT dims ) {
    // Set problem dimension
    lpData->dimy = dimy; lpData->dims = dims;
}

extern DSDP_INT lpMatSetData( lpMat *lpData, DSDP_INT *Ap, DSDP_INT *Ai, double *Ax ) {
    // Set LP data
    DSDP_INT retcode = DSDP_RETCODE_OK, m = lpData->dimy, nnz = Ap[m];
    retcode = lpMatIAlloc(lpData, nnz); checkCode;
    memcpy(lpData->Ap, Ap, sizeof(DSDP_INT) * (m + 1));
    memcpy(lpData->Ai, Ai, sizeof(DSDP_INT) * nnz);
    memcpy(lpData->Ax, Ax, sizeof(double) * nnz);
    lpData->nnz = nnz; return retcode;
}

extern void lpMataATy( double alpha, lpMat *lpData, vec *y, double *ATy ) {
    // Compute A' * y
    DSDP_INT n = lpData->dims, i, j;
    DSDP_INT *Ap = lpData->Ap, *Ai = lpData->Ai;
    double *Ax = lpData->Ax, *ydata = y->x; memset(ATy, 0, sizeof(double) * n);
    for (i = 0; i < n; ++i) {
        for (j = Ap[i]; j < Ap[i + 1]; ++i) {
            ATy[i] += alpha * Ax[Ai[j]] * ydata[Ai[j]];
        }
    }
}

extern void lpMatAx( lpMat *lpData, vec *x, vec *Ax ) {
    // Compute Ax = Ax + A * x
    DSDP_INT *lpp = lpData->Ap, *lpi = lpData->Ai, i, j;
    double *lpx = lpData->Ax, *xdata = x->x, *Axres = Ax->x;
    for (i = 0; i < lpData->dimy; ++i) {
        for (j = lpp[i]; j < lpp[i + 1]; ++j) {
            Axres[i] += lpx[j] * xdata[lpi[j]];
        }
    }
}

extern double lpMatcx( vec *c, vec *x ) {
    
    double cx; vec_dot(c, x, &cx); return cx;
}

extern void lpMatFree( lpMat *lpData ) {
    // Free the lpData structure
    if (!lpData) { return; }
    lpData->dimy = 0; lpData->dims = 0;
    DSDP_FREE(lpData->Ap); DSDP_FREE(lpData->Ai); DSDP_FREE(lpData->Ax);
    DSDP_FREE(lpData->xscale);
}

extern void lpMatView( lpMat *lpData ) {
    
    dcs mat; mat.p = lpData->Ap; mat.i = lpData->Ai; mat.x = lpData->Ax;
    mat.nz = -1; mat.nzmax = lpData->nnz; mat.n = lpData->dimy;
    mat.m = lpData->dims; dcs_print(&mat, FALSE);
}

/* SDP public methods */
extern DSDP_INT sdpMatAlloc( sdpMat *sdpData ) {
    // Allocate memory for internal SDP matrix
    DSDP_INT retcode = DSDP_RETCODE_OK;
    if ((sdpData->dimy <= 0) || (sdpData->dimS <= 0)) {
        error(etype, "Incorrect size of SDP data matrix for allocation. \n");
    }
    sdpData->types   = (DSDP_INT *) calloc(sdpData->dimy + 1, sizeof(DSDP_INT));
    sdpData->sdpData = (void **) calloc(sdpData->dimy + 1, sizeof(void *));
    if (!sdpData->types || !sdpData->sdpData) {
        printf("| Failed to allocate internal SDP data. \n");
        retcode = DSDP_RETCODE_FAILED; return retcode;
    }
    memset(sdpData->types, 0, sizeof(DSDP_INT) * sdpData->dimy);
    return retcode;
}

extern void sdpMatInit( sdpMat *sdpData ) {
    // Initialize SDP data
    sdpData->blockId     = 0;
    sdpData->dimy        = 0;
    sdpData->dimS        = 0;
    
    sdpData->nzeroMat    = 0;
    sdpData->nnzAmat     = 0;
    sdpData->nspsMat     = 0;
    sdpData->spsMatIdx   = NULL;
    sdpData->ndenseMat   = 0;
    sdpData->denseMatIdx = NULL;
    sdpData->nrkMat      = 0;
    sdpData->rkMatIdx    = NULL;
    
    sdpData->types       = NULL;
    sdpData->sdpData     = NULL;
    sdpData->scaler      = 0.0;
    
    sdpData->nzIdx       = NULL;
    sdpData->schurspIdx  = NULL;
}

extern void sdpMatSetDim( sdpMat *sdpData, DSDP_INT dimy, DSDP_INT dimS, DSDP_INT blockId ) {
    // Set dimension of the corresponding block
    sdpData->blockId = blockId; sdpData->dimy = dimy; sdpData->dimS = dimS;
}

extern DSDP_INT sdpMatSetData( sdpMat *sdpData, DSDP_INT *Ap, DSDP_INT *Ai, double *Ax, double *cnnz ) {
    // Set SDP data
    DSDP_INT retcode = DSDP_RETCODE_OK;
    
    for (DSDP_INT i = 0; i < sdpData->dimy + 1; ++i) {
        retcode = sdpMatIAllocByType(sdpData, i, &Ai[Ap[i]],
                                     &Ax[Ap[i]], Ap[i + 1] - Ap[i]);
        checkCode;
    }
    
    *cnnz = Ap[sdpData->dimy + 1] - Ap[sdpData->dimy];
        
    return retcode;
}

extern DSDP_INT sdpMatScatterNnz( sdpMat *sdpData, DSDP_INT start, DSDP_INT col, DSDP_INT *colNnz ) {
    // Scatter nonzero pattern of a block
    DSDP_INT i, j, *nzIdx = sdpData->nzIdx;
    for (i = start; i < sdpData->nzeroMat; ++i) {
        if (nzIdx[i] >= col) break;
    }
    
    if (i >= sdpData->nzeroMat) { return i; }
    if (sdpData->nzIdx[i] == col) {
        for (j = i; j < sdpData->nnzAmat; ++j) {
            colNnz[nzIdx[j]] = 1;
        }
    }
    
    return i;
}

extern void sdpMatSetSchurIndex( sdpMat *sdpData, DSDP_INT start, DSDP_INT col, DSDP_INT *csum, DSDP_INT ishift ) {
    // Associate the index for a block
    if (start >= sdpData->nzeroMat ) return;
    if (col != sdpData->nzIdx[start]) return;
    DSDP_INT i, nnzAmats = sdpData->nnzAmat;
    DSDP_INT *schurIdx = sdpData->schurspIdx;
    DSDP_INT idxshift = (2 * nnzAmats - start + 1) * start / 2;
    for (i = start; i < nnzAmats; ++i) {
        schurIdx[idxshift + i - start] = ishift + csum[sdpData->nzIdx[i]] - 1;
    }
}

extern double sdpMatGetCFnorm( sdpMat *sdpData ) {
    
    double fnrm = 0.0;
    if (sdpData->schurspIdx) {
        if (sdpData->nnzAmat == sdpData->nzeroMat) {
            fnrm = 0.0;
        } else {
            switch (sdpData->types[sdpData->nzeroMat - 1]) {
                case MAT_TYPE_SPARSE: spsMatFnorm(sdpData->sdpData[sdpData->nzeroMat - 1], &fnrm); break;
                case MAT_TYPE_DENSE : denseMatFnorm(sdpData->sdpData[sdpData->nzeroMat - 1], &fnrm); break;
                case MAT_TYPE_RANKK : rkMatFnorm(sdpData->sdpData[sdpData->nzeroMat - 1], &fnrm); break;
                default: assert(FALSE); break;
            }
        }
    } else {
        switch (sdpData->types[sdpData->dimy]) {
            case MAT_TYPE_ZERO  : fnrm = 0.0; break;
            case MAT_TYPE_SPARSE: spsMatFnorm(sdpData->sdpData[sdpData->dimy], &fnrm); break;
            case MAT_TYPE_DENSE : denseMatFnorm(sdpData->sdpData[sdpData->dimy], &fnrm); break;
            case MAT_TYPE_RANKK : rkMatFnorm(sdpData->sdpData[sdpData->dimy], &fnrm); break;
            default: assert(FALSE); break;
        }
    }
    return fnrm;
}

extern double sdpMatGetAOneNorm( sdpMat *sdpData ) {
    
    double onenorm = 0.0; DSDP_INT i, j, m = sdpData->dimy;
    if (sdpData->schurspIdx) {
        if (sdpData->nspsMat) {
            for (i = 0; i < sdpData->nspsMat - 1; ++i) {
                j = sdpData->spsMatIdx[i];
                onenorm += spsMatOneNorm(sdpData->sdpData[j]);
            }
            j = sdpData->spsMatIdx[sdpData->nspsMat - 1];
            if (sdpData->nzIdx[j] < m) onenorm += spsMatOneNorm(sdpData->sdpData[j]);
        }
        if (sdpData->ndenseMat) {
            for (i = 0; i < sdpData->ndenseMat - 1; ++i) {
                j = sdpData->denseMatIdx[i];
                onenorm += denseMatOneNorm(sdpData->sdpData[j]);
            }
            j = sdpData->denseMatIdx[sdpData->ndenseMat - 1];
            if (sdpData->nzIdx[j] < m) onenorm += denseMatOneNorm(sdpData->sdpData[j]);
        }
        if (sdpData->nrkMat) {
            for (i = 0; i < sdpData->nrkMat - 1; ++i) {
                j = sdpData->rkMatIdx[i];
                onenorm += r1MatOneNorm(((rkMat *) sdpData->sdpData[j])->data[0]);
            }
            j = sdpData->rkMatIdx[sdpData->nrkMat - 1];
            if (sdpData->nzIdx[j] < m) onenorm += \
                r1MatOneNorm(((rkMat *) sdpData->sdpData[j])->data[0]);
        }
        
    } else {
        if (sdpData->nspsMat) {
            for (i = 0; i < sdpData->nspsMat - 1; ++i) {
                j = sdpData->spsMatIdx[i];
                onenorm += spsMatOneNorm(sdpData->sdpData[j]);
            }
            j = sdpData->spsMatIdx[sdpData->nspsMat - 1];
            if (j < m) onenorm += spsMatOneNorm(sdpData->sdpData[j]);
        }
        if (sdpData->ndenseMat) {
            for (i = 0; i < sdpData->ndenseMat - 1; ++i) {
                j = sdpData->denseMatIdx[i];
                onenorm += denseMatOneNorm(sdpData->sdpData[j]);
            }
            j = sdpData->denseMatIdx[sdpData->ndenseMat - 1];
            if (j < m) onenorm += denseMatOneNorm(sdpData->sdpData[j]);
        }
        if (sdpData->nrkMat) {
            for (i = 0; i < sdpData->nrkMat - 1; ++i) {
                j = sdpData->rkMatIdx[i];
                onenorm += r1MatOneNorm(((rkMat *) sdpData->sdpData[j])->data[0]);
            }
            j = sdpData->rkMatIdx[sdpData->nrkMat - 1];
            if (j < m) onenorm += \
                r1MatOneNorm(((rkMat *) sdpData->sdpData[j])->data[0]);
        }
    }
    return onenorm;
}

extern double sdpMatGetCOneNorm( sdpMat *sdpData ) {
    if (sdpData->schurspIdx) {
        if (sdpData->nnzAmat == sdpData->nzeroMat) {
            return 0.0;
        } else {
            switch (sdpData->types[sdpData->nzeroMat - 1]) {
                case MAT_TYPE_SPARSE: return spsMatOneNorm(sdpData->sdpData[sdpData->nzeroMat - 1]); break;
                case MAT_TYPE_DENSE : return denseMatOneNorm(sdpData->sdpData[sdpData->nzeroMat - 1]); break;
                case MAT_TYPE_RANKK : return r1MatOneNorm(((rkMat *) sdpData->sdpData[sdpData->nzeroMat - 1])->data[0]); break;
                default: assert(FALSE); break;
            }
        }
    } else {
        switch (sdpData->types[sdpData->dimy]) {
            case MAT_TYPE_ZERO  : return 0.0;
            case MAT_TYPE_SPARSE: return spsMatOneNorm(sdpData->sdpData[sdpData->dimy]); break;
            case MAT_TYPE_DENSE : return denseMatOneNorm(sdpData->sdpData[sdpData->dimy]); break;
            case MAT_TYPE_RANKK : return r1MatOneNorm(((rkMat *) sdpData->sdpData[sdpData->dimy])->data[0]); break;
            default: assert(FALSE); break;
        }
    }
    return 0.0;
}

extern void sdpMatRScaleC( sdpMat *sdpData, double r ) {
    if (sdpData->schurspIdx) {
        if (sdpData->nnzAmat == sdpData->nzeroMat) {
            return;
        } else {
            switch (sdpData->types[sdpData->nzeroMat - 1]) {
                case MAT_TYPE_SPARSE: spsMatRscale(sdpData->sdpData[sdpData->nzeroMat - 1], r); break;
                case MAT_TYPE_DENSE : denseMatRscale(sdpData->sdpData[sdpData->nzeroMat - 1], r); break;
                case MAT_TYPE_RANKK : rkMatRscale(sdpData->sdpData[sdpData->nzeroMat - 1], r); break;
                default: assert(FALSE); break;
            }
        }
    } else {
        switch (sdpData->types[sdpData->dimy]) {
            case MAT_TYPE_ZERO  : return; break;
            case MAT_TYPE_SPARSE: spsMatRscale(sdpData->sdpData[sdpData->dimy], r); break;
            case MAT_TYPE_DENSE : denseMatRscale(sdpData->sdpData[sdpData->dimy], r); break;
            case MAT_TYPE_RANKK : rkMatRscale(sdpData->sdpData[sdpData->dimy], r); break;
            default: assert(FALSE); break;
        }
    }
}

extern void sdpMatAX( sdpMat *sdpData, dsMat *X, vec *AX ) {
    // AX = AX + <A, X>
    DSDP_INT i, j, issparse = (sdpData->schurspIdx) ? TRUE : FALSE, \
             nmats = (issparse) ? sdpData->nnzAmat : sdpData->dimy;
    double tmp;
    for (i = 0; i < nmats; ++i) {
        switch (sdpData->types[i]) {
            case MAT_TYPE_ZERO  : tmp = 0.0; break;
            case MAT_TYPE_DENSE : denseDsTrace(X, sdpData->sdpData[i], &tmp); break;
            case MAT_TYPE_SPARSE: denseSpsTrace(X, sdpData->sdpData[i], &tmp); break;
            case MAT_TYPE_RANKK : rkMatdenseTrace(sdpData->sdpData[i], X, &tmp); break;
            default: assert( FALSE ); break;
        }
        j = (issparse) ? sdpData->nzIdx[i] : i; AX->x[j] += tmp;
    }
}

extern double sdpMatCX( sdpMat *sdpData, dsMat *X ) {
    // Compute <C, X>
    double tmp; void *data = NULL;
    if (sdpData->schurspIdx) {
        if (sdpData->nnzAmat == sdpData->nzeroMat) {
            return 0.0;
        } else {
            data = sdpData->sdpData[sdpData->nzeroMat - 1];
            switch (sdpData->types[sdpData->nzeroMat - 1]) {
                case MAT_TYPE_SPARSE: denseSpsTrace(X, data, &tmp); break;
                case MAT_TYPE_DENSE : denseDsTrace(X, data, &tmp); break;
                case MAT_TYPE_RANKK : rkMatdenseTrace(data, X, &tmp); break;
                default: assert( FALSE ); break;
            }
        }
    } else {
        data = sdpData->sdpData[sdpData->dimy];
        switch (sdpData->types[sdpData->dimy]) {
            case MAT_TYPE_ZERO  : tmp = 0.0; break;
            case MAT_TYPE_SPARSE: denseSpsTrace(X, data, &tmp); break;
            case MAT_TYPE_DENSE : denseDsTrace(X, data, &tmp); break;
            case MAT_TYPE_RANKK : rkMatdenseTrace(data, X, &tmp); break;
            default: assert( FALSE ); break;
        }
    }
    return tmp;
}

extern void sdpMatATy( sdpMat *sdpData, double ycoef, vec *y, double tau, spsMat *S, DSDP_INT *sumHash ) {
    // Compute conic combination S = tau * C + S + ycoef * ATy
    DSDP_INT i, j, m = sdpData->dimy, *nzidx = sdpData->nzIdx; double coef = 0.0;
    DSDP_INT nsps = sdpData->nspsMat, nds = sdpData->ndenseMat, nr1 = sdpData->nrkMat;
    
    if (sdpData->schurspIdx) {
        // Add sparse
        if (nsps) {
            for (i = 0; i < nsps - 1; ++i) {
                j = sdpData->spsMatIdx[i]; coef = ycoef * y->x[nzidx[j]];
                spsMataXpbY(coef, sdpData->sdpData[j], 1.0, S, sumHash);
            }
            j = sdpData->spsMatIdx[nsps - 1];
            coef = (nzidx[j] == m) ? tau : ycoef * y->x[nzidx[j]];
            spsMataXpbY(coef, sdpData->sdpData[j], 1.0, S, sumHash);
        }
        // Add dense
        if (nds) {
            for (i = 0; i < nds - 1; ++i) {
                j = sdpData->denseMatIdx[i]; coef = ycoef * y->x[nzidx[j]];
                spsMatAddds(S, coef, sdpData->sdpData[j]);
            }
            j = sdpData->denseMatIdx[nds - 1];
            coef = (nzidx[j] == m) ? tau : ycoef * y->x[nzidx[j]];
            spsMatAddds(S, coef, sdpData->sdpData[j]);
        }
        // Add rank-1
        if (nr1) {
            for (i = 0; i < nr1 - 1; ++i) {
                j = sdpData->rkMatIdx[i]; coef = ycoef * y->x[nzidx[j]];
                spsMatAddrk(S, coef, sdpData->sdpData[j], sumHash);
            }
            j = sdpData->rkMatIdx[nr1 - 1];
            coef = (nzidx[j] == m) ? tau : ycoef * y->x[nzidx[j]];
            spsMatAddrk(S, coef, sdpData->sdpData[j], sumHash);
        }
    } else {
        // Add sparse
        if (nsps) {
            for (i = 0; i < nsps - 1; ++i) {
                j = sdpData->spsMatIdx[i]; coef = ycoef * y->x[j];
                spsMataXpbY(coef, sdpData->sdpData[j], 1.0, S, sumHash);
            }
            j = sdpData->spsMatIdx[nsps - 1];
            coef = (j == m) ? tau : ycoef * y->x[j];
            spsMataXpbY(coef, sdpData->sdpData[j], 1.0, S, sumHash);
        }
        // Add dense
        if (nds) {
            for (i = 0; i < nds - 1; ++i) {
                j = sdpData->denseMatIdx[i]; coef = ycoef * y->x[j];
                spsMatAddds(S, coef, sdpData->sdpData[j]);
            }
            j = sdpData->denseMatIdx[nds - 1];
            coef = (j == m) ? tau : ycoef * y->x[j];
            spsMatAddds(S, coef, sdpData->sdpData[j]);
        }
        // Add rank-1
        if (nr1) {
            for (i = 0; i < nr1 - 1; ++i) {
                j = sdpData->rkMatIdx[i]; coef = ycoef * y->x[j];
                spsMatAddrk(S, coef, sdpData->sdpData[j], sumHash);
            }
            j = sdpData->rkMatIdx[nr1 - 1];
            coef = (j == m) ? tau : ycoef * y->x[j];
            spsMatAddrk(S, coef, sdpData->sdpData[j], sumHash);
        }
    }
}

extern void sdpMatFree( sdpMat *sdpData ) {
    // Free SDP data
    if (!sdpData) { return; }
    void *data = NULL;
    DSDP_INT ndatatofree = (sdpData->schurspIdx) ? sdpData->nzeroMat : sdpData->dimy + 1;
    
    for (DSDP_INT i = 0; i < ndatatofree; ++i) {
        data = sdpData->sdpData[i];
        switch (sdpData->types[i]) {
            case MAT_TYPE_ZERO: break;
            case MAT_TYPE_DENSE: denseMatFree((dsMat *) data); break;
            case MAT_TYPE_SPARSE: spsMatFree((spsMat *) data); break;
            case MAT_TYPE_RANKK: rkMatFree((rkMat *) data); break;
            default: assert( FALSE );
        }
    }
    DSDP_FREE(sdpData->sdpData);
    DSDP_FREE(sdpData->types); DSDP_FREE(sdpData->schurspIdx);
    DSDP_FREE(sdpData->spsMatIdx); DSDP_FREE(sdpData->denseMatIdx);
    DSDP_FREE(sdpData->rkMatIdx); DSDP_FREE(sdpData->nzIdx);
    sdpData->nnzAmat = 0; sdpData->nspsMat = 0; sdpData->ndenseMat = 0;
    sdpData->nrkMat  = 0; sdpData->nzeroMat = 0;
    sdpData->blockId = 0; sdpData->dimy = 0; sdpData->dimS = 0;
}
