#include "schurmat.h"
#include "residualsetup.h"
#include "dsdpdata.h"
#include "structs.h"
#include "dsdpsolver.h"
#include "hsd.h"

static char etype[] = "Schur matrix setup";
/* Setup the auxiliary arrays as well as the schur matrix
 
 In DSDP, we set up the schur matrix by
 
 
 for k = 1 : nblock
 
    Sk = S[k]
    Ak = A[k]
 
    for i = 1 : m + 1
 
        X = Skinv * Ak[i] * Skinv
        asinv[k][i] = trace(A[i] * Skinv)
 
        for j = 1 : i
            M[i, j] = trace(X * A[j])
 
        u[i] = trace(X * A[m + 1])
     
*/

// A helping routine to build SinvASinv
static DSDP_INT SinvASinv( spsMat *S, DSDP_INT typeA, void *A, double *asinv, void *SinvASinv ) {
    
    // Given S and A, the routine computes A, asinv and trace(S, Sinv A Sinv)
    DSDP_INT retcode = DSDP_RETCODE_OK;
    
    if (typeA == MAT_TYPE_RANK1) {
        r1Mat *dataA = (r1Mat *) A;
        r1Mat *dataSinvASinv = (r1Mat *) SinvASinv;
        retcode = spsSinvR1SinvSolve(S, dataA, dataSinvASinv, asinv); checkCode;
    } else if (typeA == MAT_TYPE_SPARSE) {
        spsMat *dataA = (spsMat *) A;
        dsMat *dataSinvASinv = (dsMat *) SinvASinv;
        retcode = spsSinvSpSinvSolve(S, dataA, dataSinvASinv, asinv); checkCode;
    } else if (typeA == MAT_TYPE_DENSE) {
        dsMat *dataA = (dsMat *) A;
        dsMat *dataSinvASinv = (dsMat *) SinvASinv;
        retcode = spsSinvDsSinvSolve(S, dataA, dataSinvASinv, asinv); checkCode;
    } else {
        error(etype, "Invalid matrix type. \n");
    }
    
    return retcode;
}

static DSDP_INT setupSDPSchur( HSDSolver *dsdpSolver ) {
    // Set up the schur matrix for SDP
    // After calling this routine, Msdp, asinv, u for SDP and csinv will be filled
    // TODO: Also, d3 = ASinvRySinv will be set in this routine
    DSDP_INT retcode = DSDP_RETCODE_OK;
    DSDP_INT mattype = MAT_TYPE_UNKNOWN;
    
    assert( !dsdpSolver->iterProgress[ITER_SDP_SCHUR] );
    if (!dsdpSolver->iterProgress[ITER_SDP_SCHUR]) {
        error(etype, "Schur matrix is already setup. \n");
    }
    
    DSDP_INT nblock = dsdpSolver->nBlock;
    DSDP_INT dim = 0;
    DSDP_INT m   = dsdpSolver->m;
    
    sdpMat **sdpAllBlock = dsdpSolver->sdpData;
    sdpMat *sdpData = NULL;
    
    vec   **asinvOfnBlock     = dsdpSolver->asinv;
    double *csinv             = dsdpSolver->csinv;
    vec *u = dsdpSolver->u;
    dsMat *Msdp = dsdpSolver->Msdp;
    double *Mdata = Msdp->array;
    
    // Clear the Schur matrix and other arrays
    memset(Msdp->array, 0, sizeof(double) * nsym(m));
    retcode = vec_reset(u);
    
    r1Mat *r1data = NULL;
    dsMat *dsdata = NULL;
    spsMat *S     = NULL;
    vec    *asinv = NULL;
    
    void **blockdata = NULL;
    double trace = 0.0;
    
    // Start setting up the Schur matrix
    for (DSDP_INT k = 0; k < nblock; ++k) {
        
        // Setting up data for block k
        S = dsdpSolver->S[k];
        sdpData = sdpAllBlock[k];
        blockdata = sdpData->sdpData;
        asinv   = asinvOfnBlock[k];
        retcode = vec_reset(asinv);
        dim = S->dim;
        r1data = (r1Mat *) calloc(1, sizeof(r1Mat));
        dsdata = (dsMat *) calloc(1, sizeof(dsMat));
        
        // Matrices to store SinvASinv
        retcode = r1MatInit(r1data);
        retcode = r1MatAlloc(r1data, dim);
        retcode = denseMatInit(dsdata);
        retcode = denseMatAlloc(dsdata, dim, FALSE);
        
        for (DSDP_INT i = 0; i < m + 1; ++i) {
            
            // Compute SinvASinv
            mattype = sdpData->types[i];
            if (mattype == MAT_TYPE_ZERO) {
                continue;
            } else if (mattype == MAT_TYPE_RANK1) {
                retcode = SinvASinv(S, mattype, sdpData->sdpData[i],
                                    &asinv->x[i], (void *) r1data);
                checkCode;
                // retcode = r1MatCountNnz(r1data);
            } else {
                retcode = SinvASinv(S, mattype, sdpData->sdpData[i],
                                    &asinv->x[i], (void *) dsdata);
                checkCode;
            }
            
            if (mattype == MAT_TYPE_RANK1) {
                for (DSDP_INT j = 0; j <= i; ++i) {
                    switch (sdpData->types[i]) {
                        case MAT_TYPE_ZERO:
                            break;
                        case MAT_TYPE_RANK1:
                            retcode = r1Matr1Trace(r1data,
                                                   (r1Mat *) blockdata[j], &trace);
                            break;
                        case MAT_TYPE_SPARSE:
                            retcode = r1MatspsTrace(r1data,
                                                    (spsMat *) blockdata[j], &trace);
                            break;
                        case MAT_TYPE_DENSE:
                            retcode = r1MatdenseTrace(r1data,
                                                      (dsMat *) blockdata[j], &trace);
                            break;
                        default:
                            error(etype, "Unknown matrix type. \n");
                            break;
                    }
                    if (i < m) {
                        packIdx(Mdata, m, i, j) += trace;
                    } else {
                        if (j == m) {
                            csinv[k] += trace;
                        } else {
                            u->x[j] += trace;
                        }
                    }
                }
            } else {
                for (DSDP_INT j = 0; j <= i; ++i) {
                    switch (sdpData->types[i]) {
                        case MAT_TYPE_ZERO:
                            break;
                        case MAT_TYPE_RANK1:
                            retcode = r1MatdenseTrace((r1Mat *) blockdata[j], dsdata, &trace);
                            break;
                        case MAT_TYPE_SPARSE:
                            retcode = denseSpsTrace(dsdata, (spsMat *) blockdata[j], &trace);
                            break;
                        case MAT_TYPE_DENSE:
                            retcode = denseDsTrace(dsdata, (dsMat *) blockdata[j], &trace);
                            break;
                        default:
                            error(etype, "Unknown matrix type. \n");
                            break;
                    }
                    if (i < m) {
                        packIdx(Mdata, m, i, j) += trace;
                    } else {
                        if (j == m) {
                            csinv[k] += trace;
                        } else {
                            u->x[j] += trace;
                        }
                    }
                }
            }
        }

        retcode = r1MatFree(r1data); checkCode;
        retcode = denseMatFree(dsdata); checkCode;
    }
    
    dsdpSolver->iterProgress[ITER_SDP_SCHUR] = TRUE;
    
    return retcode;
}

static DSDP_INT setupLPSchur( HSDSolver *dsdpSolver ) {
    // Set up the schur matrix for LP
    // After calling this routine, Mlp will be filled. asinv, u, csinv and d3 will be updated
    
    DSDP_INT retcode = DSDP_RETCODE_OK;
    
    
    
    return retcode;
}

static DSDP_INT setupRM( HSDSolver *dsdpSolver, vec *RM ) {
    // Set up the auxiliary vector RM
    DSDP_INT retcode = DSDP_RETCODE_OK;
    
    
    return retcode;
}


static DSDP_INT setuprM( HSDSolver *dsdpSolver, vec *rM ) {
    // Set up the auxiliary vector rM
    DSDP_INT retcode = DSDP_RETCODE_OK;
    
    
    return retcode;
}

