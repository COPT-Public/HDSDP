#include <string.h>
#include "dsdpdata.h"
#include "dsdpsort.h"
// Implement the data interface for DSDP-HSD Solver
/*
 Problem data is presented using two structures
 
 sdpMat
 lpMat
 
 for SDP data and LP data respectively.
*/

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
    return retcode;
}

/* SDP internal methods */
static DSDP_INT sdpMatIAllocByType( sdpMat *sdpData, DSDP_INT k, DSDP_INT *Ai,
                                    double *Ax, DSDP_INT nnz ) {
    
    // Automatically detect the type of a matrix and allocate memory for it
    // Ai identifies row indices of the lower-triangular packed format of the k th matrix
    // Ax stores the data, nnz specifies the number of nonzeros
    DSDP_INT retcode = DSDP_RETCODE_OK;
    
    DSDP_INT m = sdpData->dimy, n = sdpData->dimS;
    // Sort the arrays
    assert( k < m + 1 );
    void *userdata = NULL;
    
    // Check sparsity
    if (nnz == 0) {
        sdpData->nzeroMat += 1;
        sdpData->types[k] = MAT_TYPE_ZERO;
        
    } else if (sdpData->types[k] == MAT_TYPE_RANKK) {
        // Rank 1
        sdpData->nrkMat += 1;
        rkMat *data = NULL;
        data = (rkMat *) calloc(1, sizeof(rkMat)); rkMatInit(data);
        r1Mat *r1data = (r1Mat *) calloc(1, sizeof(r1Mat)); r1MatInit(r1data);
        retcode = r1MatAlloc(r1data, n); checkCode;
        
        // The first non-zero element must be a diagonal element
        double adiag = Ax[0];
        DSDP_INT isNeg = FALSE, rowidx = Ai[0];
        
        if (adiag > 0) {
            adiag = sqrt(adiag); r1data->sign = 1.0;
        } else {
            r1data->sign = - 1.0; adiag = - sqrt(-adiag);
            isNeg = TRUE;
        }
        
        // Seek which column is correct
        DSDP_INT where = 0, diagidx = 0;
        for (DSDP_INT i = 0; i < n; ++i) {
            if (rowidx == diagidx) { break; }
            diagidx += n - where; where += 1;
        }
        
        if (where == n) {
            error(etype, "The matrix is not rank-one. \n");
        }
        
        for (DSDP_INT i = 0; i < n - where; ++i) {
            rowidx = Ai[i];
            // The compressed matrix has n * (n + 1) / 2 rows and one column
            if (rowidx >= diagidx + n - where) {
                break;
            }
            r1data->x[rowidx - diagidx + where] = Ax[i] / adiag;
        }
        rkMatAllocAndSetData(data, n, 1, &r1data->sign, r1data->x);
        r1MatFree(r1data);
        DSDP_FREE(r1data);
        userdata = (void *) data;
        
    } else if (((nnz <= denseThresh * nsym(n)) && (sdpData->types[k] == MAT_TYPE_UNKNOWN)) ||
               (sdpData->types[k] == MAT_TYPE_SPARSE)) {
        
        // May be put in earlier parts
        idxsort(Ai, Ax, nnz);
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
        
        retcode = denseMatInit(data); checkCode;
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
    DSDP_INT retcode = DSDP_RETCODE_OK;
    DSDP_INT m = lpData->dimy, nnz = Ap[m];
    retcode = lpMatIAlloc(lpData, nnz);
    
    // Set problem data
    memcpy(lpData->Ap, Ap, sizeof(DSDP_INT) * (m + 1));
    memcpy(lpData->Ai, Ai, sizeof(DSDP_INT) * nnz);
    memcpy(lpData->Ax, Ax, sizeof(double)   * nnz);
    
    // Some extra parameters
    lpData->nnz = nnz;
    
    return retcode;
}

extern void lpMataATy( double alpha, lpMat *lpData, vec *y, double *ATy ) {
    // Compute A' * y
    DSDP_INT n = lpData->dims;
    DSDP_INT *Ap = lpData->Ap, *Ai = lpData->Ai;
    double *Ax = lpData->Ax, *ydata = y->x;
        memset(ATy, 0, sizeof(double) * n);
    for (DSDP_INT i = 0; i < n; ++i) {
        for (DSDP_INT j = Ap[i]; j < Ap[i + 1]; ++i) {
            ATy[i] += alpha * Ax[Ai[j]] * ydata[Ai[j]];
        }
    }
}

extern void lpMatFree( lpMat *lpData ) {

    // Free the lpData structure
    lpData->dimy = 0;
    lpData->dims = 0;
    DSDP_FREE(lpData->Ap); DSDP_FREE(lpData->Ai); DSDP_FREE(lpData->Ax);
    DSDP_FREE(lpData->xscale);
}

extern void lpMatView( lpMat *lpData ) {
    
    cs mat; mat.p = lpData->Ap; mat.i = lpData->Ai; mat.x = lpData->Ax;
    mat.nz = -1; mat.nzmax = lpData->nnz; mat.n = lpData->dimy;
    mat.m = lpData->dims; cs_print(&mat, FALSE);
}

/* SDP public methods */
extern DSDP_INT sdpMatAlloc( sdpMat *sdpData ) {
    
    // Allocate memory for internal SDP matrix
    DSDP_INT retcode = DSDP_RETCODE_OK;
    assert( (sdpData->dimy > 0) && (sdpData->dimS > 0) );
    
    if ((sdpData->dimy <= 0) || (sdpData->dimS <= 0)) {
        error(etype, "Incorrect size for SDP data matrix for allocation. \n");
    }
    
    sdpData->types   = (DSDP_INT *) calloc(sdpData->dimy + 1, sizeof(DSDP_INT));
    sdpData->sdpData = (void **) calloc(sdpData->dimy + 1, sizeof(void *));
    
    for (DSDP_INT i = 0; i < sdpData->dimy; ++i) {
        sdpData->types[i] = MAT_TYPE_UNKNOWN;
    }
    
    return retcode;
}

extern void sdpMatInit( sdpMat *sdpData ) {
    
    // Initialize SDP data
    sdpData->blockId     = 0;
    sdpData->dimy        = 0;
    sdpData->dimS        = 0;
    
    sdpData->nzeroMat    = 0;
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

extern void sdpMatSetHint( sdpMat *sdpData, DSDP_INT *hint ) {
    
    // Set type hint for matrix type
    for (DSDP_INT i = 0; i < sdpData->dimy + 1; ++i) {
        switch (hint[i]) {
            case MAT_TYPE_RANKK  : sdpData->types[i] = MAT_TYPE_RANKK; break;
            case MAT_TYPE_DENSE  : sdpData->types[i] = MAT_TYPE_DENSE; break;
            case MAT_TYPE_SPARSE : sdpData->types[i] = MAT_TYPE_SPARSE; break;
            case MAT_TYPE_UNKNOWN: sdpData->types[i] = MAT_TYPE_ZERO; break;
            default: sdpData->types[i] = MAT_TYPE_UNKNOWN; break;
        }
    }
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
        for (j = i; j < sdpData->nzeroMat; ++j) {
            colNnz[nzIdx[j]] = 1;
        }
    }
    return i;
}

extern void sdpMatSetSchurIndex( sdpMat *sdpData, DSDP_INT start, DSDP_INT col, DSDP_INT *csum, DSDP_INT ishift ) {
    // Associate the index for a block
    if (start >= sdpData->nzeroMat ) return;
    if (col != sdpData->nzIdx[start]) return;
    DSDP_INT i, nnzmat = sdpData->nzeroMat;
    DSDP_INT *schurIdx = sdpData->schurspIdx;
    DSDP_INT idxshift = (2 * nnzmat - start + 1) * start / 2;
    for (i = start; i < nnzmat; ++i) {
        schurIdx[idxshift + i - start] = ishift + csum[sdpData->nzIdx[i]];
    }
}

extern void sdpMatFree( sdpMat *sdpData ) {
    
    // Free SDP data
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
        
    DSDP_FREE(sdpData->types); DSDP_FREE(sdpData->schurspIdx);
    DSDP_FREE(sdpData->spsMatIdx); DSDP_FREE(sdpData->denseMatIdx);
    DSDP_FREE(sdpData->rkMatIdx); DSDP_FREE(sdpData->nzIdx);
    sdpData->nspsMat   = 0; sdpData->ndenseMat = 0;
    sdpData->nrkMat    = 0; sdpData->nzeroMat  = 0;
    sdpData->blockId = 0; sdpData->dimy = 0; sdpData->dimS = 0;
}
