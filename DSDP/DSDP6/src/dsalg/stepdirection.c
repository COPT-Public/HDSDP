#include "stepdirection.h"
/* Recover the steps to take after solving the Schur system */

static char etype[] = "Stepsize recovery";

static DSDP_INT getdTau( HSDSolver *dsdpSolver, double rM, vec *b2 ) {
    
    DSDP_INT retcode = DSDP_RETCODE_OK;
    
    vec *d1 = dsdpSolver->d1;
    vec *d2 = dsdpSolver->d2;
    
    DSDP_INT m = d1->dim;
    DSDP_INT one = 1;
    
    double mu = dsdpSolver->mu;
    double T = 0.0;
    double b2dotd1 = 0.0;
    double b2dotd2 = 0.0;
    
    b2dotd1 = ddot(&m, b2->x, &one, d1->x, &one);
    T = b2dotd1 + mu * dsdpSolver->csinvcsinv + dsdpSolver->kappa / dsdpSolver->tau;
    if (fabs(T) < 1e-12) {
        dsdpSolver->dtau = 0.0;
    } else {
        b2dotd2 = ddot(&m, b2->x, &one, d2->x, &one);
        dsdpSolver->dtau = (rM - b2dotd2) / T;
    }
    
    return retcode;
}

static DSDP_INT getdy( HSDSolver *dsdpSolver ) {
    // dy = d1 * dtau + d2
    DSDP_INT retcode = DSDP_RETCODE_OK;
    
    vec *dy = dsdpSolver->dy;
    vec *d1 = dsdpSolver->d1;
    vec *d2 = dsdpSolver->d2;
    
    retcode = vec_zaxpby(dy, dsdpSolver->dtau, d1, 1.0, d2);
    
    return retcode;
}

static DSDP_INT getBlockdS( HSDSolver *dsdpSolver, DSDP_INT blockid ) {
    // Compute dS for some block (dS = Ry - ATdy + C * dtau)
    // dS is computed exactly in the same way as Ry
    DSDP_INT retcode = DSDP_RETCODE_OK;
    
    DSDP_INT nblock = dsdpSolver->nBlock;
    DSDP_INT m = dsdpSolver->m;
    
    assert( blockid < nblock );
    sdpMat *cone = dsdpSolver->sdpData[blockid];
    spsMat *Ry = dsdpSolver->Rys[blockid];
    spsMat *dS = dsdpSolver->dS[blockid];
    assert( cone->dimy == m );
    
    vec *dy = dsdpSolver->dy;
    DSDP_INT n = cone->dimS;
    
    assert( dS->cscMat->nzmax == nsym(n) );
    assert( dS->cscMat->n == n );
    
    DSDP_INT nspsMat  = cone->nspsMat;
    DSDP_INT ndsMat   = cone->ndenseMat;
    DSDP_INT nr1Mat   = cone->nr1Mat;
    
    DSDP_INT *spsMatIdx = cone->spsMatIdx;
    DSDP_INT *dsMatIdx  = cone->denseMatIdx;
    DSDP_INT *r1MatIdx  = cone->r1MatIdx;
    
    void **conedata = cone->sdpData;
    
    DSDP_INT *dSp = dS->cscMat->p;
    DSDP_INT *dSi = dS->cscMat->i;
    double   *dSx = dS->cscMat->x;
    
    dsMat  *dsdata = NULL;
    spsMat *spdata = NULL;
    r1Mat  *r1data = NULL;
    
    // Clean the dS array
    memset(dSx, 0, sizeof(double) * nsym(n));
    
    // Set up the direction
    DSDP_INT planA    = FALSE;
    DSDP_INT maxR1nnz = 0;
    double coeff      = 0.0;
    double sign       = 0.0;
    
    for (DSDP_INT i = 0; i < ndsMat; ++i) {
        if (fabs(dy->x[dsMatIdx[i]]) > 1e-12) {
            planA = TRUE;
        }
    }
    
    if (!planA) {
        for (DSDP_INT i = 0; i < nr1Mat; ++i) {
            if (fabs(dy->x[r1MatIdx[i]]) > 1e-12) {
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
                coeff = dsdpSolver->dtau;
            } else {
                coeff = - dy->x[dsMatIdx[i]];
            }
            
            // Do not sum if the dimension is empty
            if (fabs(coeff) <= 1e-10) {
                continue;
            }
            
            dsdata = (dsMat *) conedata[dsMatIdx[i]];
            for (DSDP_INT j = 0; j < (DSDP_INT) (nsym(n) / 4); j+=4) {
                dSx[j    ] += coeff * dsdata->array[j    ];
                dSx[j + 1] += coeff * dsdata->array[j + 1];
                dSx[j + 2] += coeff * dsdata->array[j + 2];
                dSx[j + 3] += coeff * dsdata->array[j + 3];
            }
            
            for (DSDP_INT j = 4 * (DSDP_INT) (nsym(n) / 4); j < nsym(n); ++j) {
                dSx[j] += coeff * dsdata->array[j];
            }
        }
        
        // Then sum the rank 1 matrices
        for (DSDP_INT i = 0; i < nr1Mat; ++i) {
            
            if (r1MatIdx[i] == m) {
                coeff = dsdpSolver->dtau;
            } else {
                coeff = - dy->x[r1MatIdx[i]];
            }
            
            // Do not sum if the dimension is empty
            if (fabs(coeff) <= 1e-10) {
                continue;
            }
            
            r1data = (r1Mat *) conedata[r1MatIdx[i]];
            r1denseSpsUpdate(dS, coeff, r1data);
        }
        
        // Then collect the sparse data
        for (DSDP_INT k = 0; k < nr1Mat; ++k) {
            
            if (spsMatIdx[k] == m) {
                coeff = dsdpSolver->dtau;
            } else {
                coeff = - dy->x[spsMatIdx[k]];
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
                    packIdx(dSx, n, row, i) += coeff * spx[j];
                }
            }
        }
        
        // Last we add Rd to finish the construction
        coeff  = 1.0;
        spdata = Ry;
        row    = 0;
        spp    = spdata->cscMat->p;
        spi    = spdata->cscMat->i;
        spx    = spdata->cscMat->x;
        
        for (DSDP_INT i = 0; i < n; ++i) {
            for (DSDP_INT j = spp[i]; j < spp[i + 1]; ++j) {
                row = spi[j];
                packIdx(dSx, n, row, i) += spx[j]; // 1.0 * spx[j]
            }
        }
        
        // Change the Ry data (p and i)
        DSDP_INT counter = 0;
        dSp[0] = 0;
        for (DSDP_INT i = 0; i < n; ++i) {
            for (DSDP_INT j = 0; j < n - i; ++j) {
                dSi[i] = j;
            }
            counter += n - i;
            dSp[i + 1] = counter;
        }
        
        // Clear the zero elements
        cs_dropzeros(dS->cscMat);
        
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
                    coeff = - dy->x[spsMatIdx[k]];
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
            coeff  = 1.0;
            spdata = Ry;
            spp    = spdata->cscMat->p;
            spi    = spdata->cscMat->i;
            spx    = spdata->cscMat->x;
            
            for (DSDP_INT j = spdata->cscMat->p[i]; j < spdata->cscMat->p[i + 1]; ++j) {
                colNnzs[spi[j]] += spx[j]; // 1.0 * spx[j]
            }
            
            // Rank 1 column
            for (DSDP_INT k = 0; k < nr1Mat; ++k) {
                
                if (r1MatIdx[k] == m) {
                    coeff = dsdpSolver->tau;
                } else {
                    coeff = - dy->x[r1MatIdx[k]];
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
                if (fabs(colNnzs[j]) > 1e-12) {
                    dSi[nnz] = j;
                    dSx[nnz] = colNnzs[j];
                    nnz += 1;
                }
            }
            dSp[i + 1] = nnz;
        }
        
        DSDP_FREE(colNnzs);
    }
    
    
    return retcode;
}

static DSDP_INT getdS( HSDSolver *dsdpSolver ) {
    
    DSDP_INT retcode = DSDP_RETCODE_OK;
    DSDP_INT nblock = dsdpSolver->nBlock;
    
    for (DSDP_INT i = 0; i < nblock; ++i) {
        retcode = getBlockdS(dsdpSolver, i); checkCode;
    }
    
    return retcode;
}

static DSDP_INT getdKappa( HSDSolver *dsdpSolver ) {
    
    DSDP_INT retcode = DSDP_RETCODE_OK;
    
    dsdpSolver->dkappa = - dsdpSolver->kappa
                         + dsdpSolver->mu / dsdpSolver->tau
                         - (dsdpSolver->kappa * dsdpSolver->dtau) / dsdpSolver->tau;
    return retcode;
}

static DSDP_INT getdsLP( HSDSolver *dsdpSolver ) {
    // Compute the dual direction for LP ry - A' * dy + c * dtau
    DSDP_INT retcode = DSDP_RETCODE_OK;
    
    retcode = checkIterProgress(dsdpSolver, ITER_RECOVER_LP_DIR);
    assert( !dsdpSolver->iterProgress[ITER_RECOVER_LP_DIR] );
    
    if (dsdpSolver->iterProgress[ITER_RECOVER_LP_DIR]) {
        error(etype, "LP directions have been set up. \n");
    }
    
    vec *ds = dsdpSolver->ds;
    double *dydata = dsdpSolver->dy->x;
    double *dsdata = ds->x;
    DSDP_INT *Ap = dsdpSolver->lpData->lpdata->p;
    DSDP_INT *Ai = dsdpSolver->lpData->lpdata->i;
    double   *Ax = dsdpSolver->lpData->lpdata->x;
    DSDP_INT n = ds->dim;
    
    retcode = vec_copy(dsdpSolver->ry, ds);
    retcode = vec_axpy(dsdpSolver->dtau, dsdpSolver->lpObj, ds);
    
    double tmp = 0.0;
    
    for (DSDP_INT i = 0; i < n; ++i) {
        tmp = 0.0;
        for (DSDP_INT j = Ap[i]; j < Ap[i + 1]; ++i) {
            tmp += dydata[Ai[j]] * Ax[j];
        }
        dsdata[i] -= tmp;
    }
    
    dsdpSolver->iterProgress[ITER_RECOVER_LP_DIR] = TRUE;
    
    return retcode;
}

static DSDP_INT getSDPDirs( HSDSolver *dsdpSolver, double rM, vec *b2 ) {
    
    DSDP_INT retcode = DSDP_RETCODE_OK;
    retcode = checkIterProgress(dsdpSolver, ITER_RECOVER_SDP_DIR);
    assert( !dsdpSolver->iterProgress[ITER_RECOVER_SDP_DIR] );
    
    if (dsdpSolver->iterProgress[ITER_RECOVER_SDP_DIR]) {
        error(etype, "SDP directions have been set up. \n");
    }
    
    retcode = getdTau(dsdpSolver, rM, b2); checkCode;
    retcode = getdy(dsdpSolver); checkCode;
    retcode = getdS(dsdpSolver); checkCode;
    retcode = getdKappa(dsdpSolver); checkCode;
    
    dsdpSolver->iterProgress[ITER_RECOVER_SDP_DIR] = TRUE;
    
    return retcode;
}

extern DSDP_INT getStepDirs( HSDSolver *dsdpSolver, double rM, vec *b2 ) {
    
    DSDP_INT retcode = DSDP_RETCODE_OK;
    
    retcode = getSDPDirs(dsdpSolver, rM, b2); checkCode;
    retcode = getdsLP(dsdpSolver); checkCode;
    
    return retcode;
}
