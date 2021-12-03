#include "residualsetup.h"
#include "dsdpdata.h"
#include "structs.h"
#include "dsdpsolver.h"
#include "hsd.h"

// Setup all the residuals for the dual scaling algorithm to work
// We also note that the dual residual is not contained in the solver since
// we will eliminate it as soon as possible
static char etype[] = "Residual setup";

static DSDP_INT getSDPBlockResidualRy( HSDSolver *dsdpSolver, spsMat *Ry, DSDP_INT blockid ) {
    /* There are three types of matrices: sparse/dense packed/rank 1
        
     We consider the following three cases
     
     A. There exists dense matrix
     B. There are no dense matrices but the NNZ of some rank1 Matrix > 0.8 * n
     
     In the above two cases we treat Ry as dense and directly operate over
     the x array of the CSC matrix
     
     C. There is no dense or r1 matrix with more than 0.5 * n nonzeros in the array
     
     In this case we do addition column-wise by
     
     1. Add all the sparse matrices to the residual
     2. Add rank-1 matrices to the residual
     
     Note that this routine is also used for recovering dual step dS
    */
    
    // Setup the dual infeasibility of some block
    DSDP_INT retcode = DSDP_RETCODE_OK;
    
    DSDP_INT nblock = dsdpSolver->nBlock;
    DSDP_INT m = dsdpSolver->m;
    
    assert( blockid < nblock );
    sdpMat *cone = dsdpSolver->sdpData[blockid];
    assert( cone->dimy == m );
    
    vec *y = dsdpSolver->y;
    DSDP_INT n = cone->dimS;
    
    assert( Ry->cscMat->nzmax == nsym(n) );
    assert( Ry->cscMat->n == n );
    
    DSDP_INT nspsMat  = cone->nspsMat;
    DSDP_INT ndsMat   = cone->ndenseMat;
    DSDP_INT nr1Mat   = cone->nr1Mat;
    
    DSDP_INT *spsMatIdx = cone->spsMatIdx;
    DSDP_INT *dsMatIdx  = cone->denseMatIdx;
    DSDP_INT *r1MatIdx  = cone->r1MatIdx;
    
    void **conedata = cone->sdpData;
    
    DSDP_INT *Ryp = Ry->cscMat->p;
    DSDP_INT *Ryi = Ry->cscMat->i;
    double   *Ryx = Ry->cscMat->x;
    
    dsMat  *dsdata = NULL;
    spsMat *spdata = NULL;
    r1Mat  *r1data = NULL;
    
    // Clean the residual
    memset(Ryx, 0, sizeof(double) * nsym(n));
    
    // Setup dual residual
    DSDP_INT planA    = FALSE;
    DSDP_INT maxR1nnz = 0;
    double coeff      = 0.0;
    double sign       = 0.0;
    
    // TODO: save the plan to avoid computing again
    for (DSDP_INT i = 0; i < ndsMat; ++i) {
        if (fabs(y->x[dsMatIdx[i]]) > 1e-12) {
            planA = TRUE;
        }
    }
    
    if (!planA) {
        for (DSDP_INT i = 0; i < nr1Mat; ++i) {
            if (fabs(y->x[r1MatIdx[i]]) > 1e-12) {
                continue;
            }
            maxR1nnz = MAX(maxR1nnz,
                           ((r1Mat *) conedata[r1MatIdx[i]])->nnz);
        }
        if (maxR1nnz > 0.8 * n) {
            planA = TRUE;
        }
    }
    
    if (planA) {
        
        DSDP_INT row  = 0;
        DSDP_INT *spp = NULL;
        DSDP_INT *spi = NULL;
        double   *spx = NULL;
        
        // Directly sum the arrays as if they are all dense
        // First sum all the dense matrices
        for (DSDP_INT i = 0; i < ndsMat; ++i) {
            
            // In case of C
            if (dsMatIdx[i] == m) {
                coeff = dsdpSolver->tau;
            } else {
                coeff = - y->x[dsMatIdx[i]];
            }
            
            // Do not sum if the dimension is empty
            if (fabs(coeff) <= 1e-10) {
                continue;
            }
            
            dsdata = (dsMat *) conedata[dsMatIdx[i]];
            for (DSDP_INT j = 0; j < (DSDP_INT) (nsym(n) / 4); j+=4) {
                Ryx[j    ] += coeff * dsdata->array[j    ];
                Ryx[j + 1] += coeff * dsdata->array[j + 1];
                Ryx[j + 2] += coeff * dsdata->array[j + 2];
                Ryx[j + 3] += coeff * dsdata->array[j + 3];
            }
            
            for (DSDP_INT j = 4 * (DSDP_INT) (nsym(n) / 4); j < nsym(n); ++j) {
                Ryx[j] += coeff * dsdata->array[j];
            }
        }
        
        // Then sum the rank 1 matrices
        for (DSDP_INT i = 0; i < nr1Mat; ++i) {
            
            if (r1MatIdx[i] == m) {
                coeff = dsdpSolver->tau;
            } else {
                coeff = - y->x[r1MatIdx[i]];
            }
            
            // Do not sum if the dimension is empty
            if (fabs(coeff) <= 1e-10) {
                continue;
            }
            
            r1data = (r1Mat *) conedata[r1MatIdx[i]];
            r1denseSpsUpdate(Ry, coeff, r1data);
        }
        
        // Then collect the sparse data
        for (DSDP_INT k = 0; k < nr1Mat; ++k) {
            
            if (spsMatIdx[k] == m) {
                coeff = dsdpSolver->tau;
            } else {
                coeff = - y->x[spsMatIdx[k]];
            }
            
            if (fabs(coeff) <= 1e-10) {
                continue;
            }
            
            spdata = (spsMat *) conedata[spsMatIdx[k]];
            row    = 0;
            spp    = spdata->cscMat->p;
            spi    = spdata->cscMat->i;
            spx    = spdata->cscMat->x;
            
            // We take summation column-wise
            for (DSDP_INT i = 0; i < n; ++i) {
                for (DSDP_INT j = spp[i]; j < spp[i + 1]; ++j) {
                    row = spi[j];
                    packIdx(Ryx, n, row, i) += coeff * spx[j];
                }
            }
        }
        
        // Last we subtract S off to finish the construction
        coeff  = - 1.0;
        spdata = dsdpSolver->S[blockid];
        row    = 0;
        spp    = spdata->cscMat->p;
        spi    = spdata->cscMat->i;
        spx    = spdata->cscMat->x;
        
        for (DSDP_INT i = 0; i < n; ++i) {
            for (DSDP_INT j = spp[i]; j < spp[i + 1]; ++j) {
                row = spi[j];
                packIdx(Ryx, n, row, i) += coeff * spx[j];
            }
        }
        
        // Change the Ry data (p and i)
        DSDP_INT counter = 0;
        Ryp[0] = 0;
        for (DSDP_INT i = 0; i < n; ++i) {
            for (DSDP_INT j = 0; j < n - i; ++j) {
                Ryi[i] = j;
            }
            counter += n - i;
            Ryp[i + 1] = counter;
        }
        
        // Clear the zero elements
        cs_dropzeros(Ry->cscMat);
        
    } else {
        // Plan B
        // There is no dense matrix
        // Sparse matrix
        double *colNnzs   = NULL;
        DSDP_INT *spp     = NULL;
        DSDP_INT *spi     = NULL;
        double   *spx     = NULL;
        DSDP_INT *r1nzIdx = NULL;
        double   *r1a     = NULL;
        DSDP_INT nnz      = 0;
        double   anz      = 0.0;
        
        colNnzs = (double *) calloc(n, sizeof(double));
        
        for (DSDP_INT i = 0; i < n; ++i) {
            memset(colNnzs, 0, sizeof(double) * n);
            
            // Sparse columns from A
            for (DSDP_INT k = 0; k < nspsMat; ++k) {
                
                if (spsMatIdx[k] == m) {
                    coeff = dsdpSolver->tau;
                } else {
                    coeff = - y->x[spsMatIdx[k]];
                }
                
                if (fabs(coeff) < 1e-10) {
                    continue;
                }
                
                spdata = (spsMat *) conedata[spsMatIdx[k]];
                spp = spdata->cscMat->p;
                spi = spdata->cscMat->i;
                spx = spdata->cscMat->x;
                
                for (DSDP_INT j = spdata->cscMat->p[i]; j < spdata->cscMat->p[i + 1]; ++j) {
                    colNnzs[spi[j]] += coeff * spx[j];
                }
            }
            
            // Sparse column from S
            coeff  = - 1.0;
            spdata = dsdpSolver->S[blockid];
            spp    = spdata->cscMat->p;
            spi    = spdata->cscMat->i;
            spx    = spdata->cscMat->x;
            
            for (DSDP_INT j = spdata->cscMat->p[i]; j < spdata->cscMat->p[i + 1]; ++j) {
                colNnzs[spi[j]] += coeff * spx[j];
            }
            
            // Rank 1 column
            for (DSDP_INT k = 0; k < nr1Mat; ++k) {
                
                if (r1MatIdx[k] == m) {
                    coeff = dsdpSolver->tau;
                } else {
                    coeff = - y->x[r1MatIdx[k]];
                }
                
                if (fabs(coeff) < 1e-10) {
                    continue;
                }
                r1data  = (r1Mat *) conedata[r1MatIdx[k]];
                r1nzIdx = r1data->nzIdx;
                r1a     = r1data->x;
                nnz     = r1data->nnz;
                sign    = (double) r1data->sign;
                anz     = coeff * r1a[i] * sign;
                
                if (anz == 0) {
                    continue;
                } else {
                    for (DSDP_INT j = 0; j < nnz; ++j) {
                        colNnzs[r1nzIdx[j]] += r1a[r1nzIdx[j]] * anz;
                    }
                }
            }
            
            nnz = 0;
            for (DSDP_INT j = 0; j < n; ++j) {
                if (colNnzs[j] > 0) {
                    Ryi[nnz] = j;
                    Ryx[nnz] = colNnzs[j];
                    nnz += 1;
                }
            }
            Ryp[i + 1] = nnz;
        }
        
        DSDP_FREE(colNnzs);
    }
    
    return retcode;
}

static DSDP_INT getLPResidualry( HSDSolver *dsdpSolver, vec *ry ) {
    // Compute ry = - A' * y + c * tau - s
    DSDP_INT retcode = DSDP_RETCODE_OK;
    
    lpMat *Adata = dsdpSolver->lpData;
    DSDP_INT m   = Adata->dimy;
    
    vec *y = dsdpSolver->y;
    vec *s = dsdpSolver->s;
    vec *c = dsdpSolver->lpObj;
    double tau = dsdpSolver->tau;
    double *ydata = ry->x;
    
    assert( ry->dim = m );
    
    // Setup - A' * y
    double alpha = - 1.0;
    retcode = lpMataATy(alpha, Adata, y, ydata);
    
    // Get - A' * y + c * tau
    retcode = vec_axpy(tau, c, ry);
    
    // Get - A' * y + c * tau - s
    retcode = vec_axpy(alpha, s, ry);
    
    return retcode;
}

static DSDP_INT getRkappaTau( HSDSolver *dsdpSolver, double *rtk ) {
    // Setup the complementarity residual rtk
    DSDP_INT retcode = DSDP_RETCODE_OK;
    
    *rtk = dsdpSolver->tau * dsdpSolver->kappa - dsdpSolver->mu;
    
    return retcode;
}

extern DSDP_INT setupRes( HSDSolver *dsdpSolver ) {
    
    // Setup residuals used for LP and SDP
    DSDP_INT retcode = DSDP_RETCODE_OK;
    
    spsMat **Rys = dsdpSolver->Rys;
    vec *ry = dsdpSolver->ry;
    
    assert( !dsdpSolver->iterProgress[ITER_LP_RESIDUAL] );
    if (dsdpSolver->iterProgress[ITER_LP_RESIDUAL]) {
        error(etype, "LP residual has been setup. \n");
    }
    
    retcode = checkIterProgress(dsdpSolver, ITER_LP_RESIDUAL);
    checkCode;
    
    if (!dsdpSolver->eventMonitor[EVENT_LP_NO_RY]) {
        retcode = getLPResidualry(dsdpSolver, ry); checkCode;
    }
    
    dsdpSolver->iterProgress[ITER_LP_RESIDUAL] = TRUE;
    
    // Set up SDP residual
    assert( !dsdpSolver->iterProgress[ITER_SDP_RESIDUAL] );
    if (dsdpSolver->iterProgress[ITER_SDP_RESIDUAL]) {
        error(etype, "SDP residual has been setup for the current block. \n");
    }
    
    retcode = checkIterProgress(dsdpSolver, ITER_SDP_RESIDUAL);
    checkCode;
    
    DSDP_INT nblock = dsdpSolver->nBlock;
    
    if (!dsdpSolver->eventMonitor[EVENT_SDP_NO_RY]) {
        for (DSDP_INT i = 0; i < nblock; ++i) {
            retcode = getSDPBlockResidualRy(dsdpSolver, Rys[i], i);
        }
    }
    
    // Setup rtk
    retcode = getRkappaTau(dsdpSolver, &dsdpSolver->rtk);
    dsdpSolver->iterProgress[ITER_SDP_RESIDUAL] = TRUE;
    
    return retcode;
}
