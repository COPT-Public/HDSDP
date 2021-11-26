#include "dsdpdata.h"
// Implement the data interface for DSDP-HSD Solver
/*
 Problem data is presented using two structures
 
 sdpMat
 lpMat
 
 for SDP data and LP data respectively.
*/

static char *etype = "Data interface";

/* LP internal methods */
static DSDP_INT lpMatIAlloc( lpMat *lpData, DSDP_INT nnz ) {
    
    // Allocate memory for internal LP matrix
    DSDP_INT retcode = DSDP_RETCODE_OK;
    assert( (lpData->dimy > 0) && (lpData->dims > 0) );
    
    if ((lpData->dimy <= 0) || (lpData->dims <= 0)) {
        error(etype, "Incorrect size for LP data matrix for allocation. \n");
    }
    
    lpData->lpdata = cs_spalloc(lpData->dimy, lpData->dims, nnz, TRUE, FALSE);
    return retcode;
}

/* SDP internal methods */
static DSDP_INT sdpMatIAlloc( sdpMat *sdpData ) {
    
    // Allocate memory for internal SDP matrix
    DSDP_INT retcode = DSDP_RETCODE_OK;
    assert( (sdpData->dimy > 0) && (sdpData->dimS > 0) );
    
    if ((sdpData->dimy <= 0) || (sdpData->dimS <= 0)) {
        error(etype, "Incorrect size for SDP data matrix for allocation. \n");
    }
    
    sdpData->types   = (DSDP_INT *) calloc(sdpData->dimy, sizeof(DSDP_INT));
    sdpData->sdpData = (void **) calloc(sdpData->dimy, sizeof(void *));
    
    for (DSDP_INT i = 0; i < sdpData->dimy; ++i) {
        sdpData->types = MAT_TYPE_UNKNOWN;
    }
    
    return retcode;
}

static DSDP_INT sdpMatIAllocByType( sdpMat *sdpData, DSDP_INT k, DSDP_INT *Ai,
                                    double *Ax, DSDP_INT nnz ) {
    
    // Automatically detect the type of a matrix and allocate memory for it
    // Ai identifies row indices of the lower-triangular packed format of the k th matrix
    // Ax stores the data, nnz specifies the number of nonzeros
    DSDP_INT retcode = DSDP_RETCODE_OK;
    
    DSDP_INT m = sdpData->dimy;
    DSDP_INT n = sdpData->dimS;
    assert( k < m );
    
    void *userdata = NULL;
    
    // Check sparsity
    if (nnz == 0) {
        sdpData->types[k] = MAT_TYPE_ZERO;
    } else if (sdpData->types[k] == MAT_TYPE_RANK1) {
        // Rank 1
        sdpData->nr1Mat += 1;
        r1Mat *data = NULL;
        data = (r1Mat *) calloc(1, sizeof(r1Mat));
        retcode = vec_init(data); checkCode;
        retcode = vec_alloc(data, m); checkCode;
        
        for (DSDP_INT i = 0; i < nnz; ++i) {
            if (Ai[i] >= m) {
                break;
            }
            data->x[Ai[i]] = Ax[i];
        }
        
        userdata = (void *) data;
        
    } else if (((nnz <= nsym(n) / 10) && (sdpData->types[k] == MAT_TYPE_UNKNOWN)) ||
               (sdpData->types[k] == MAT_TYPE_SPARSE)) {
        // Sparse
        sdpData->types[k] = MAT_TYPE_SPARSE;
        sdpData->nspsMat += 1;
        spsMat *data      = NULL;
        data    = (spsMat *) calloc(1, sizeof(spsMat));
        retcode = spsMatInit(data); checkCode;
        retcode = spsMatAllocData(data, n, nnz); checkCode;
        
        for (DSDP_INT i = 0; i < n; ++i) {
            data->cscMat->p[i + 1] = nnz;
        }
        
        memcpy(data->cscMat->x, Ax, sizeof(double) * nnz);
        
        DSDP_INT rowidx    = 0;
        DSDP_INT where     = 0;
        DSDP_INT colnnz    = 0;
        DSDP_INT idxthresh = n;
        DSDP_INT diff      = n;
        data->cscMat->p[0] = 0;
        
        for (DSDP_INT i = 0; i < nnz; ++i) {
            rowidx = Ai[i];
            if (rowidx >= idxthresh) {
                where += 1;
                data->cscMat->p[where] = colnnz;
                diff = n - where;
                idxthresh += diff;
                colnnz = 0;
            }
            colnnz += 1;
            data->cscMat->i[i] = rowidx - idxthresh + diff;
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
        
        DSDP_INT rowidx    = 0;
        DSDP_INT where     = 0;
        DSDP_INT idxthresh = n;
        DSDP_INT diff      = n;
        
        for (DSDP_INT i = 0; i < nnz; ++i) {
            rowidx = Ai[i];
            if (rowidx >= idxthresh) {
                where += 1;
                diff = n - where;
                idxthresh += diff;
            }
            data->array[rowidx - idxthresh + diff] = Ax[i];
        }
        
        userdata = (void *) data;
    }
    
    sdpData->sdpData[k] = userdata;
    
    return retcode;
}

/* LP public methods */
extern DSDP_INT lpMatInit( lpMat *lpData ) {
    
    // Initialize LP data
    DSDP_INT retcode = DSDP_RETCODE_OK;
    
    lpData->dims = 0;
    lpData->dimy = 0;
    lpData->lpdata = NULL;
    lpData->xscale = NULL;
    
    return retcode;
}

extern DSDP_INT lpMatSetDim( lpMat *lpData, DSDP_INT dimy, DSDP_INT dims ) {
    
    // Set problem dimension
    DSDP_INT retcode = DSDP_RETCODE_OK;
    assert( (lpData->dimy > 0) && (lpData->dims > 0) );
    
    lpData->dimy = dimy;
    lpData->dims = dims;
    
    return retcode;
}

extern DSDP_INT lpMatSetData( lpMat *lpData, DSDP_INT *Ap, DSDP_INT *Ai, double *Ax ) {
    
    // Set LP data
    DSDP_INT retcode = DSDP_RETCODE_OK;
    DSDP_INT n = lpData->dims;
    
    DSDP_INT nnz = Ap[n];
    retcode = lpMatIAlloc(lpData, nnz);
    
    // Set problem data
    memcpy(lpData->lpdata->p, Ap, sizeof(DSDP_INT) * (n + 1));
    memcpy(lpData->lpdata->i, Ai, sizeof(DSDP_INT) * nnz);
    memcpy(lpData->lpdata->x, Ax, sizeof(double)   * nnz);
    
    // Some extra parameters
    lpData->lpdata->nz = nnz;
    
    return retcode;
}

/* SDP public methods */
extern DSDP_INT sdpMatInit( sdpMat *sdpData ) {
    
    // Initialize SDP data
    DSDP_INT retcode = DSDP_RETCODE_OK;
    
    sdpData->blockId     = 0;
    sdpData->dimy        = 0;
    sdpData->dimS        = 0;
    
    sdpData->nspsMat     = 0;
    sdpData->spsMatIdx   = NULL;
    sdpData->ndenseMat   = 0;
    sdpData->denseMatIdx = NULL;
    sdpData->nr1Mat      = 0;
    sdpData->r1MatIdx    = NULL;
    
    sdpData->types       = NULL;
    sdpData->sdpData     = NULL;
    sdpData->scaler      = 0.0;
    
    return retcode;
}

extern DSDP_INT sdpMatSetDim( sdpMat *sdpData, DSDP_INT dimy, DSDP_INT dimS, DSDP_INT blockId ) {
    
    // Set dimension of the corresponding block
    DSDP_INT retcode = DSDP_RETCODE_OK;
    assert( (dimy > 0) && (dimS > 0) );
    
    sdpData->blockId = blockId;
    sdpData->dimy    = dimy;
    sdpData->dimS    = dimS;
    
    return retcode;
}

extern DSDP_INT sdpMatSetHint( sdpMat *sdpData, DSDP_INT *hint ) {
    
    // Set type hint for matrix type
    DSDP_INT retcode = DSDP_RETCODE_OK;
    
    for (DSDP_INT i = 0; i < sdpData->dimy; ++i) {
        
        switch (hint[i]) {
            case MAT_TYPE_RANK1:
                sdpData->types[i] = MAT_TYPE_RANK1;
                break;
            case MAT_TYPE_DENSE:
                sdpData->types[i] = MAT_TYPE_DENSE;
                break;
            case MAT_TYPE_SPARSE:
                sdpData->types[i] = MAT_TYPE_SPARSE;
                break;
            default:
                sdpData->types[i] = MAT_TYPE_UNKNOWN;
                break;
        }
    }
    
    return retcode;
}

extern DSDP_INT sdpMatSetData( sdpMat *sdpData, DSDP_INT *Ap, DSDP_INT *Ai, double *Ax ) {
    
    // Set SDP data
    DSDP_INT retcode = DSDP_RETCODE_OK;
    retcode = sdpMatIAlloc(sdpData); checkCode;
    
    for (DSDP_INT i = 0; i < sdpData->dimy; ++i) {
        retcode = sdpMatIAllocByType(sdpData, i, &Ai[Ap[i]],
                                     &Ax[Ap[i]], Ap[i + 1] - Ap[i]);
        checkCode;
    }
        
    return retcode;
}
