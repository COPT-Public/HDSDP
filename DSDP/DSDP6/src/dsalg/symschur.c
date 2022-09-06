#include "symschur.h"
#include "dsdplapack.h"
#include "dsdpsolver.h"
#include "schurmat.h"
#include "sparseopts.h"
#include "denseopts.h"
#include "rank1opts.h"
#include "rankkopts.h"
#include "vec.h"
#include "dsdpsort.h"

// Implement the advanced Schur matrix setup in DSDP using M1, M2, M3, M4, M5 techniques


#ifndef continue_if_not_in_A_or_not_building_hsd
#define continue_if_not_in_A_or_not_building_hsd if (!(*M->phaseA) || !(*M->buildhsd) ) { continue; }
#define return_if_not_in_A_or_not_building_hsd   if ((!(*M->phaseA) || !(*M->buildhsd)) && computeC) { return retcode; }
#endif

static char etype[] = "Advanced Schur matrix setup";

static DSDP_INT schurBlockContract( sdpMat *sdpData ) {
    // Contract a SDP block if it only contains few entries for the Schur matrix
    DSDP_INT retcode = DSDP_RETCODE_OK;
    DSDP_INT nzero = sdpData->nzeroMat, m = sdpData->dimy, nnzmats = m + 1 - nzero, nnzAmats, i, j;
    nnzAmats = (sdpData->types[m] == MAT_TYPE_ZERO) ? nnzmats : nnzmats - 1;
    DSDP_INT *types = sdpData->types; void **data = sdpData->sdpData;
    void **newdata = (void **) calloc(nnzmats, sizeof(void *));
    DSDP_INT *newtypes = (DSDP_INT *) calloc(nnzmats, sizeof(DSDP_INT));
    DSDP_INT *nzidx = (DSDP_INT *) calloc(nnzmats, sizeof(DSDP_INT));
    DSDP_INT spidx = 0, dsidx = 0, r1idx = 0;
    
    for (i = j = 0; i < m + 1; ++i) {
        switch (types[i]) {
            case MAT_TYPE_ZERO  : break;
            case MAT_TYPE_SPARSE:
                newdata[j] = data[i]; nzidx[j] = i; newtypes[j] = types[i];
                sdpData->spsMatIdx[spidx] = j; ++j; ++spidx; break;
            case MAT_TYPE_DENSE :
                newdata[j] = data[i]; nzidx[j] = i; newtypes[j] = types[i];
                sdpData->denseMatIdx[dsidx] = j; ++j; ++dsidx; break;
            case MAT_TYPE_RANKK :
                newdata[j] = data[i]; nzidx[j] = i; newtypes[j] = types[i];
                sdpData->rkMatIdx[r1idx] = j; ++j; ++r1idx; break;
            default:
                assert( FALSE );
        }
    }
    
    assert( j = nnzmats && r1idx == sdpData->nrkMat &&
            spidx == sdpData->nspsMat && dsidx == sdpData->ndenseMat);
    DSDP_FREE(sdpData->types); DSDP_FREE(sdpData->sdpData);
    sdpData->types = newtypes; sdpData->sdpData = newdata;
    sdpData->nzIdx = nzidx; sdpData->nzeroMat = nnzmats; sdpData->nnzAmat = nnzAmats; // nzero now stands for nnz
    sdpData->schurspIdx = (DSDP_INT *) calloc(nsym(nnzAmats), sizeof(DSDP_INT));
    return retcode;
}

static DSDP_INT getScore( DSDP_INT *ranks, DSDP_INT *nnzs, DSDP_INT *perm,
                        DSDP_INT m, DSDP_INT n, DSDP_INT idx, DSDP_INT *MX,
                        double *M1M2, double *M1M5 ) {
    // Compute the estimated number of multiplications for technique M1 to M5
    
    DSDP_INT i, j, k = perm[idx]; //, npack = nsym(n);
    double d1, d2, d3, d4, d5, nsqr = n * n, ncbe = n * nsqr,\
           sumnnz = 0.0, rsigma = ranks[k], best = SCHUR_M1;
    double m1m5 = 0.0, m1m2 = 0.0;
    // On exit, MX is set 10 * MX12 + M12345
    for (i = idx; i < m + 1; ++i) { sumnnz += nnzs[i]; }
    // TODO: Change the way to compute d1
    d1 = rsigma * (2.0 * nsqr) + KAPPA * sumnnz;
    // d2 = rsigma * (nsqr + 2 * KAPPA * sumnnz);
    d2 = rsigma * ((double) nnzs[idx] * n + 3 * KAPPA * sumnnz);
    // d1 = rsigma * (nsqr + 3 * npack) + 2 * KAPPA * sumnnz;
    // d2 = rsigma * (nsqr + 3 * KAPPA * sumnnz);
    d3 = (double) n * KAPPA * nnzs[idx] + ncbe + KAPPA * sumnnz + ncbe / m;
    d4 = (double) n * KAPPA * nnzs[idx] + KAPPA * (n + 1) * sumnnz + ncbe / m;
    d5 = KAPPA * (2.0 * KAPPA * nnzs[idx] + 1) * sumnnz + ncbe / m;
    // Get better of d1 and d2
    j = 0;
    if (d1 < d2) {
        best = SCHUR_M1; m1m2 = d1; m1m5 = d1;
    } else {
        best = SCHUR_M2; m1m2 = d2; m1m5 = d2;
    }
    j += best * 10;
    if (d3 < m1m5) { best = SCHUR_M3; m1m5 = d3; }
    if (d4 < m1m5) { best = SCHUR_M4; m1m5 = d4; }
    if (d5 < m1m5) { best = SCHUR_M5; m1m5 = d5; }
    *M1M2 = m1m2; *M1M5 = m1m5; j += best; *MX = j;

    return DSDP_RETCODE_OK;
}

static DSDP_INT schurBlockAnalysis( sdpMat *Adata, DSDP_INT *permk, DSDP_INT *MXk, DSDP_INT *useM1M2 ) {
    // Reorder block
    DSDP_INT retcode = DSDP_RETCODE_OK;
    
    DSDP_INT *dsidx = Adata->denseMatIdx, ndsMat = Adata->ndenseMat,
             *spsidx = Adata->spsMatIdx, nspsMat = Adata->nspsMat,
             *rkidx = Adata->rkMatIdx, nrkMat = Adata->nrkMat,
             m = Adata->dimy, n = Adata->dimS;
    
    memset(permk, 0, sizeof(DSDP_INT) * (m + 1));
    memset(MXk, 0, sizeof(DSDP_INT) * (m + 1));
    double *scoreM1M2 = (double *) calloc(m + 1, sizeof(double)), ksiM1M2 = 0.0;
    double *scoreM1M5 = (double *) calloc(m + 1, sizeof(double)), ksiM1M5 = 0.0;
    DSDP_INT *ranks   = (DSDP_INT *) calloc(m + 1, sizeof(DSDP_INT));
    DSDP_INT *nnzs = (DSDP_INT *) calloc(m + 1, sizeof(DSDP_INT));
    DSDP_INT i, j, k, M1M2 = FALSE;
    
    // Get f_1, ..., f_m in MX
    for (i = 0; i < m + 1; ++i) { permk[i] = i; }
    for (i = 0, j = n * n; i < ndsMat; ++i) {
        k = dsidx[i];
        nnzs[k] = j;
        denseMatGetRank(Adata->sdpData[k], &ranks[k]);
    }
    for (i = 0; i < nspsMat; ++i) {
        k = spsidx[i];
        nnzs[k] = ((spsMat *) Adata->sdpData[k])->nnz;
        ranks[k] = spsMatGetRank(Adata->sdpData[k]);
    }
    for (i = 0; i < nrkMat; ++i) {
        k = rkidx[i];
        nnzs[k] = ((rkMat *) (Adata->sdpData[k]))->data[0]->nnz;
        nnzs[k] = nsym(nnzs[k]); ranks[k] = 1;
    }
    
    // Sort f_1, ..., f_m in descending order
    dsdpSort2(permk, nnzs, 0, m);
    // Determine strategies for setting up the Schur matrix
    for (i = 0; i < m + 1; ++i) {
        getScore(ranks, nnzs, permk, m, // m here is correct
                 n, i, &MXk[i], &scoreM1M2[i], &scoreM1M5[i]);
    }
    
    // Determine whether to use M3 to M5
    for (i = 0; i < m + 1; ++i) {
        ksiM1M2 += scoreM1M2[i]; ksiM1M5 += scoreM1M5[i];
    }
    
    double ncbe = n; ncbe *= n;
    M1M2 = FALSE;// (ksiM1M2 <= ksiM1M5 + ncbe) ? TRUE: FALSE;
    
    // Determine which strategy to use
    if (M1M2) {
        for (i = 0; i < m + 1; ++i) {
            MXk[i] = (DSDP_INT) (MXk[i] / 10);
        }
    } else {
        for (i = 0; i < m + 1; ++i) {
            MXk[i] = MXk[i] % 10;
        }
    }
    
    *useM1M2 = M1M2;
    DSDP_FREE(scoreM1M2); DSDP_FREE(scoreM1M5); DSDP_FREE(ranks); DSDP_FREE(nnzs);
    return retcode;
}

static void schurMatSinvCleanup( symM *M ) {
    // Set up the inverse of matrices when necessary
    for (DSDP_INT i = 0, j, n; i < M->nblock; ++i) {
        n = M->Adata[i]->dimS; memset(M->Sinv[i], 0, sizeof(double) * n * n);
        for (j = 0; j < n; ++j) { M->Sinv[i][j * n + j] = 1.0; }
    }
}

static void schurMatCleanup( symM *M ) {
    // Clean up the arrays related to the Schur matrix
    schurMatReset(M->M, TRUE); vec_reset(M->asinv);
    if (*M->phaseA) {
        if (*M->buildhsd) {
            vec_reset(M->asinvcsinv);
            *M->csinvrysinv = 0.0; *M->csinvcsinv = 0.0;
        }
        vec_reset(M->asinvrysinv); *M->rysinv = 0.0;
    }
    *M->csinv = 0.0;
    schurMatSinvCleanup(M);
}

static void schurMatGetSinv( symM *M ) {
    // Compute inverse of the dual matrix when M3, M4 or M5 techniques are used
    M->scaler = 1.0;
    
    for (DSDP_INT i = 0; i < M->nblock; ++i) {
        spsMatInverse(M->S[i], M->Sinv[i], M->schurAux);
    }
    
    if (*M->phaseA) {
        DSDP_INT n; double tracesinv = 0.0;
        for (DSDP_INT i = 0, j; i < M->nblock; ++i) {
            for (j = 0, n = M->S[i]->dim; j < n; ++j) {
                tracesinv += M->Sinv[i][j * n + j];
            }
        }
        *M->rysinv = tracesinv * (*M->Ry);
    }
}

static DSDP_INT schurM1rowSetup( symM *M, DSDP_INT blockid, DSDP_INT row ) {
    // Apply M1 technique to setup a row of the Schur complement matrix
    /*
     M1 technique first applies rank-1 update to get B_i = <S^-1 * A_i * S^-1>
     and then computes M_{ij} = <B_i, A_j>
     
     When each row is set up, not only M is computed, but also
        
     csinvrysinv, asinvrysinv, csinvcsinv, asinvcsinv, csinv, asinv
     
     when available
    */
    
    DSDP_INT retcode = DSDP_RETCODE_OK;
    DSDP_INT i, j, k, m = M->m, *perm = M->perms[blockid], computeC = FALSE, gold = *M->gold;
    dsMat *B = M->B[blockid]; rkMat *factor, *rkaux = M->rkaux[blockid];
    double *val, res, Ry = *M->Ry; k = perm[row];
    
    switch (M->Adata[blockid]->types[k]) {
        case MAT_TYPE_ZERO  : return retcode;
        case MAT_TYPE_SPARSE: factor = spsMatGetFactor(M->Adata[blockid]->sdpData[k]); break;
        case MAT_TYPE_DENSE : factor = denseMatGetFactor(M->Adata[blockid]->sdpData[k]); break;
        case MAT_TYPE_RANKK : factor = M->Adata[blockid]->sdpData[k]; break;
        default             : error(etype, "Invalid matrix type. \n"); break;
    }
    
    if (LIKELY(k < m)) {
        val = &M->asinv->x[k];
    } else {
        val = M->csinv; if (!gold) { return retcode; }
        computeC = TRUE;
    }
    
    res = SinvRkSinv(M->S[blockid], factor, rkaux);
    *val += res; return_if_not_in_A_or_not_building_hsd
    
    // Get B through rank-one update
    denseMatReset(B); rkMatdenseUpdate(B, rkaux);
    
    /* Start M1 */
    if (computeC) {
        assert( *M->phaseA );
        val = M->csinvrysinv; *val += denseDiagTrace(B, Ry);
        for (i = row; i < m + 1; ++i) {
            j = perm[i];
            if (j == m) {
                if (!*M->buildhsd) { continue; } val = M->csinvcsinv;
            } else {
                val = &M->asinvcsinv->x[j];
            }
            switch (M->Adata[blockid]->types[j]) {
                case MAT_TYPE_ZERO  : continue; break;
                case MAT_TYPE_DENSE : denseDsTrace(B, M->Adata[blockid]->sdpData[j], &res); break;
                case MAT_TYPE_SPARSE: denseSpsTrace(B, M->Adata[blockid]->sdpData[j], &res); break;
                case MAT_TYPE_RANKK : rkMatdenseTrace(M->Adata[blockid]->sdpData[j], B, &res); break;
                default             : error(etype, "Invalid matrix type. \n"); break;
            }
            *val += res;
        }
    } else {
        if (*M->phaseA) { val = &M->asinvrysinv->x[k]; *val += denseDiagTrace(B, Ry); }
        double *array = M->M->denseM->array;
        for (i = row; i < m + 1; ++i) {
            j = perm[i]; if (M->Adata[blockid]->types[j] == MAT_TYPE_ZERO) { continue; }
            if (j == m) {
                continue_if_not_in_A_or_not_building_hsd
                val = &M->asinvcsinv->x[k];
            } else {
                val = (j > k) ? (&fullIdx(array, m, j, k)) : (&fullIdx(array, m, k, j));
            }
            switch (M->Adata[blockid]->types[j]) {
                case MAT_TYPE_ZERO  : continue; break;
                case MAT_TYPE_DENSE : denseDsTrace(B, M->Adata[blockid]->sdpData[j], &res); break;
                case MAT_TYPE_SPARSE: denseSpsTrace(B, M->Adata[blockid]->sdpData[j], &res); break;
                case MAT_TYPE_RANKK : rkMatdenseTrace(M->Adata[blockid]->sdpData[j], B, &res); break;
                default             : error(etype, "Invalid matrix type. \n"); break;
            }
            *val += res;
        }
    }
    /* End M1 */
    
    return retcode;
}

static DSDP_INT schurM2rowSetup( symM *M, DSDP_INT blockid, DSDP_INT row ) {
    // Apply M2 technique to setup a row of the Schur complement matrix
    /*
     M2 technique does not form <S^-1 * A_i * S^-1> but compute the value of quadratic forms
     directly.
     
     When each row is set up, now only M is computed, but also
     
     csinvrysinv, asinvrysinv, csinvcsinv, asinvcsinv, csinv, asinv
     
     when available
    */
    
    DSDP_INT retcode = DSDP_RETCODE_OK;
    
    DSDP_INT i, j, k, r, m = M->m, *perm = M->perms[blockid], computeC = FALSE, rank, gold = *M->gold;
    rkMat *factor, *rkaux = M->rkaux[blockid]; r1Mat *r1aux = NULL;
    double *val, res, coeff, Ry = *M->Ry; k = perm[row];
    
    switch (M->Adata[blockid]->types[k]) {
        case MAT_TYPE_ZERO  : return retcode;
        case MAT_TYPE_SPARSE: factor = spsMatGetFactor(M->Adata[blockid]->sdpData[k]); break;
        case MAT_TYPE_DENSE : factor = denseMatGetFactor(M->Adata[blockid]->sdpData[k]); break;
        case MAT_TYPE_RANKK : factor = M->Adata[blockid]->sdpData[k]; break;
        default             : error(etype, "Invalid matrix type. \n"); break;
    }
    
    if (k < m) {
        val = &M->asinv->x[k];
    } else {
        val = M->csinv; if (!gold) { return retcode; }
        computeC = TRUE;
    }
    
    res = SinvRkSinv(M->S[blockid], factor, rkaux);
    *val += res; return_if_not_in_A_or_not_building_hsd
    rank = rkMatGetRank(factor);
    
    /* Start M2 */
    if (computeC) {
        assert(*M->phaseA);
        val = M->csinvrysinv; rkMatdiagTrace(rkaux, Ry, &res); *val += res;
        for (r = 0; r < rank; ++r) {
            r1aux = rkMatGetBase(rkaux, r); coeff = r1aux->sign;
            for (i = row; i < m + 1; ++i) {
                j = perm[i]; if (M->Adata[blockid]->types[j] == MAT_TYPE_ZERO) { continue; }
                if (j == m) {
                    if (!*M->buildhsd) { continue; } val = M->csinvcsinv;
                } else {
                    val = &M->asinvcsinv->x[j];
                }
                switch (M->Adata[blockid]->types[j]) {
                    case MAT_TYPE_DENSE : res = denseMatxTAx(M->Adata[blockid]->sdpData[j], M->schurAux, r1aux->x);
                        res *= coeff; break;
                    case MAT_TYPE_SPARSE: res = spsMatxTAx(M->Adata[blockid]->sdpData[j], r1aux->x);
                        res *= coeff; break;
                    case MAT_TYPE_RANKK :
                        res = r1Matr1Trace(((rkMat *) M->Adata[blockid]->sdpData[j])->data[0], r1aux);
                        break;
                    default             : error(etype, "Invalid matrix type. \n"); break;
                }
                *val += res;
            }
        }
    } else {
        if (*M->phaseA) { val = &M->asinvrysinv->x[k]; rkMatdiagTrace(rkaux, Ry, &res); *val += res; }
        double *array = M->M->denseM->array;
        for (r = 0; r < rank; ++r) {
            r1aux = rkMatGetBase(rkaux, r); coeff = r1aux->sign;
            for (i = row; i < m + 1; ++i) {
                j = perm[i]; if (M->Adata[blockid]->types[j] == MAT_TYPE_ZERO) { continue; }
                if (j == m) {
                    continue_if_not_in_A_or_not_building_hsd
                    val = &M->asinvcsinv->x[k];
                } else {
                    val = (j > k) ? (&fullIdx(array, m, j, k)) : (&fullIdx(array, m, k, j));
                }
                
                switch (M->Adata[blockid]->types[j]) {
                    case MAT_TYPE_DENSE : res = denseMatxTAx(M->Adata[blockid]->sdpData[j], M->schurAux, r1aux->x);
                        res *= coeff; break;
                    case MAT_TYPE_SPARSE: res = spsMatxTAx(M->Adata[blockid]->sdpData[j], r1aux->x);
                        res *= coeff; break;
                    case MAT_TYPE_RANKK :
                        res = r1Matr1Trace(((rkMat *) M->Adata[blockid]->sdpData[j])->data[0], r1aux);
                        break;
                    default             : error(etype, "Invalid matrix type. \n"); break;
                }

                *val += res;
            }
        }
    }
    /* End M2 */
    
    return retcode;
}

static DSDP_INT schurM3rowSetup( symM *M, DSDP_INT blockid, DSDP_INT row ) {
    // Apply M3 technique to setup a row of the Schur complement matrix
    
    /*
     M3 Technique is like M1 but it computes B_i = <S^-1 * A_i * S^-1> explicitly
     by inverting S at the beginning of each iteration
     
     When each row is set up, not only M is computed, but also
        
     csinvrysinv, asinvrysinv, csinvcsinv, asinvcsinv, csinv, asinv
     
     when available
    */
    DSDP_INT retcode = DSDP_RETCODE_OK;
    DSDP_INT i, j, k, m = M->m, *perm = M->perms[blockid], computeC = FALSE, gold = *M->gold;
    dsMat *B = M->B[blockid]; spsMat *S; rkMat *rkaux;
    double *val, res, Ry = *M->Ry, *aux = M->schurAux, *Sinv = M->Sinv[blockid]; k = perm[row];
    
    if (k < m) {
        val = &M->asinv->x[k];
    } else {
        val = M->csinv; if (!gold) { return retcode; }
        computeC = TRUE;
    }
    
    // Get B through direct matrix multiplication
    denseMatReset(B);
    // Compute B_i = <S^-1 * A_i * S^-1>
    switch (M->Adata[blockid]->types[k]) {
        case MAT_TYPE_ZERO: return retcode;
        case MAT_TYPE_SPARSE:
            res = SinvSpSinv(Sinv, aux, M->Adata[blockid]->sdpData[k], B); break;
        case MAT_TYPE_DENSE:
            res = SinvDsSinv(Sinv, aux, M->Adata[blockid]->sdpData[k], B); break;
        case MAT_TYPE_RANKK:
            S = M->S[blockid]; rkaux = M->rkaux[blockid];
            res = SinvRkSinv(S, M->Adata[blockid]->sdpData[k], rkaux);
            rkMatdenseUpdate(B, rkaux); break;
        default: error(etype, "Invalid matrix type. \n"); break;
    }
    
    *val += res; return_if_not_in_A_or_not_building_hsd
    
    /* Start M3 */
    if (computeC) {
        assert( *M->phaseA ); val = M->csinvrysinv; *val += denseDiagTrace(B, Ry);
        for (i = row; i < m + 1; ++i) {
            j = perm[i]; if (M->Adata[blockid]->types[j] == MAT_TYPE_ZERO) { continue; }
            if (j == m) {
                if (!*M->buildhsd) { continue; } val = M->csinvcsinv;
            } else {
                val = &M->asinvcsinv->x[j];
            }
            switch (M->Adata[blockid]->types[j]) {
                case MAT_TYPE_ZERO  : continue; break;
                case MAT_TYPE_DENSE : denseDsTrace(B, M->Adata[blockid]->sdpData[j], &res); break;
                case MAT_TYPE_SPARSE: denseSpsTrace(B, M->Adata[blockid]->sdpData[j], &res); break;
                case MAT_TYPE_RANKK : rkMatdenseTrace(M->Adata[blockid]->sdpData[j], B, &res); break;
                default             : error(etype, "Invalid matrix type. \n"); break;
            }
            *val += res;
        }
    } else {
        if (*M->phaseA) { val = &M->asinvrysinv->x[k]; *val += denseDiagTrace(B, Ry); }
        double *array = M->M->denseM->array;
        for (i = row; i < m + 1; ++i) {
            j = perm[i];
            if (j == m) {
                continue_if_not_in_A_or_not_building_hsd
                val = &M->asinvcsinv->x[k];
            } else {
                val = (j > k) ? (&fullIdx(array, m, j, k)) : (&fullIdx(array, m, k, j));
            }
            switch (M->Adata[blockid]->types[j]) {
                case MAT_TYPE_ZERO  : continue; break;
                case MAT_TYPE_DENSE : denseDsTrace(B, M->Adata[blockid]->sdpData[j], &res); break;
                case MAT_TYPE_SPARSE: denseSpsTrace(B, M->Adata[blockid]->sdpData[j], &res); break;
                case MAT_TYPE_RANKK : rkMatdenseTrace(M->Adata[blockid]->sdpData[j], B, &res); break;
                default             : error(etype, "Invalid matrix type. \n"); break;
            }
            *val += res;
        }
    }
    
    /* End M3 */
    return retcode;
}

static DSDP_INT schurM4rowSetup( symM *M, DSDP_INT blockid, DSDP_INT row ) {
    // Apply M4 technique to setup a row of the Schur complement matrix
    DSDP_INT retcode = DSDP_RETCODE_OK;
    
    DSDP_INT i, j, k, m = M->m, *perm = M->perms[blockid], computeC = FALSE, gold = *M->gold;
    double *val, res, tmp, Ry = *M->Ry, *ASinv = M->schurAux, *Sinv = M->Sinv[blockid]; k = perm[row];
    double *aux = M->rkaux[blockid]->data[0]->x;
    
    if (k < m) {
        val = &M->asinv->x[k];
    } else {
        val = M->csinv; if (!gold) { return retcode; }
        computeC = TRUE;
    }
    
    // Prepare ASinv for later computation. asinv and asinvrysinv will also be ready
    switch (M->Adata[blockid]->types[k]) {
        case MAT_TYPE_ZERO  : return retcode;
        case MAT_TYPE_SPARSE:
            res = spsSinvSolve(Sinv, M->Adata[blockid]->sdpData[k], ASinv, &tmp, Ry); break;
        case MAT_TYPE_DENSE :
            res = denseSinvSolve(Sinv, M->Adata[blockid]->sdpData[k], ASinv, &tmp, Ry); break;
        case MAT_TYPE_RANKK :
            res = r1MatSinvSolve(Sinv, ((rkMat *) M->Adata[blockid]->sdpData[k])->data[0], ASinv, aux, &tmp, Ry); break;
        default: error(etype, "Invalid matrix type. \n"); break;
    }
    
    *val += tmp; return_if_not_in_A_or_not_building_hsd
    
    /* Start M4 */
    if (computeC) {
        assert( *M->phaseA ); val = M->csinvrysinv; *val += res;
        for (i = row; i < m + 1; ++i) {
            j = perm[i]; if (M->Adata[blockid]->types[j] == MAT_TYPE_ZERO) { continue; }
            if (j == m) {
                if (!*M->buildhsd) { continue; } val = M->csinvcsinv;
            } else {
                val = &M->asinvcsinv->x[j];
            }
            switch (M->Adata[blockid]->types[j]) {
                case MAT_TYPE_ZERO  : continue; break;
                case MAT_TYPE_SPARSE:
                    res = spsSinvASinv(Sinv, M->Adata[blockid]->sdpData[j], ASinv); break;
                case MAT_TYPE_DENSE :
                    res = denseSinvASinv(Sinv, M->Adata[blockid]->sdpData[j], ASinv); break;
                case MAT_TYPE_RANKK :
                    res = r1MatSinvASinv(Sinv, ((rkMat *) M->Adata[blockid]->sdpData[j])->data[0], ASinv); break;
                default: error(etype, "Invalid matrix type. \n"); break;
            }
            *val += res;
        }
    } else {
        val = &M->asinvrysinv->x[k]; *val += res;
        double *array = M->M->denseM->array;
        for (i = row; i < m + 1; ++i) {
            j = perm[i];
            if (j == m) {
                continue_if_not_in_A_or_not_building_hsd
                val = &M->asinvcsinv->x[k];
            } else {
                val = (j > k) ? (&fullIdx(array, m, j, k)) : (&fullIdx(array, m, k, j));
            }
            
            switch (M->Adata[blockid]->types[j]) {
                case MAT_TYPE_ZERO: continue; break;
                case MAT_TYPE_SPARSE:
                    res = spsSinvASinv(Sinv, M->Adata[blockid]->sdpData[j], ASinv); break;
                case MAT_TYPE_DENSE:
                    res = denseSinvASinv(Sinv, M->Adata[blockid]->sdpData[j], ASinv); break;
                case MAT_TYPE_RANKK:
                    res = r1MatSinvASinv(Sinv, ((rkMat *) M->Adata[blockid]->sdpData[j])->data[0], ASinv); break;
                default: error(etype, "Invalid matrix type. \n"); break;
            }
            *val += res;
        }
    }
    
    /* End M4 */
    return retcode;
}

static double schurM5MatAux( DSDP_INT type1, void *A1, DSDP_INT type2, void *A2, double *Sinv ) {
    if (type1 == MAT_TYPE_SPARSE) {
        switch (type2) {
            // TODO: Consider changing the order of the two matrices by sparsity
            case MAT_TYPE_SPARSE:
                return spsSinvspsSinv(A2, A1, Sinv);
            case MAT_TYPE_RANKK : return spsSinvr1Sinv(A1, ((rkMat *) A2)->data[0], Sinv);
            case MAT_TYPE_ZERO  : return 0.0;
            default: assert( FALSE ); break;
        }
    } else {
        switch (type2) {
            case MAT_TYPE_SPARSE: return r1Sinvsps(A2, ((rkMat *) A1)->data[0], Sinv);
            case MAT_TYPE_RANKK : return r1Sinvr1(((rkMat *) A1)->data[0], ((rkMat *) A2)->data[0], Sinv);
            case MAT_TYPE_ZERO  : return 0.0;
            default: assert( FALSE ); break;
        }
    }
    assert( FALSE ); return 0.0;
}

static DSDP_INT schurM5rowSetup( symM *M, DSDP_INT blockid, DSDP_INT row ) {
    // Apply M5 technique to setup a row of the Schur complement matrix
    /*
     M5 Technique directly computes elements of the Schur matrix without any intermediate variables
    */
    DSDP_INT retcode = DSDP_RETCODE_OK;
    DSDP_INT i, j, k, m = M->m, *perm = M->perms[blockid], computeC = FALSE, gold = *M->gold;
    DSDP_INT *types = M->Adata[blockid]->types;
    double *Sinv = M->Sinv[blockid];
    double *val, res, tmp, Ry = *M->Ry; k = perm[row];
    void **data = M->Adata[blockid]->sdpData;
    
    /* Start M5 */
    if (k < m) {
        val = &M->asinv->x[k];
    } else {
        val = M->csinv; if (!gold) { return retcode; }
        computeC = TRUE;
    }
    
    switch (types[k]) {
        case MAT_TYPE_ZERO: return retcode;
        case MAT_TYPE_SPARSE:
            res = spsRySinv(data[k], Sinv, &tmp, Ry); break;
        case MAT_TYPE_RANKK:
            res = r1RySinv(((rkMat *) data[k])->data[0], Sinv, &tmp, Ry, M->schurAux); break;
        default:
            error(etype, "Invalid matrix type. Dense matrix should not appear in M5. \n");
            break;
    }
    
    *val += tmp; return_if_not_in_A_or_not_building_hsd
    
    if (computeC) {
        *M->csinvrysinv += res;
        for (i = row; i < m + 1; ++i) {
            j = perm[i]; if (M->Adata[blockid]->types[j] == MAT_TYPE_ZERO) { continue; }
            if (j == m) {
                if (!*M->buildhsd) { continue; } val = M->csinvcsinv;
            } else {
                val = &M->asinvcsinv->x[j];
            }
            *val += schurM5MatAux(types[k], data[k], types[j], data[j], Sinv);
        }
    } else {
        M->asinvrysinv->x[k] += res;
        double *array = M->M->denseM->array;
        for (i = row; i < m + 1; ++i) {
            j = perm[i];
            if (j == m) {
                continue_if_not_in_A_or_not_building_hsd
                val = &M->asinvcsinv->x[k];
            } else {
                val = (j > k) ? (&fullIdx(array, m, j, k)) : (&fullIdx(array, m, k, j));
            }
            *val += schurM5MatAux(types[k], data[k], types[j], data[j], Sinv);
        }
    }
    
    /* End M5 */
    return retcode;
}

// #define SCHUR_PROFILER
#ifdef SCHUR_PROFILER
static double schurTime[5] = {0.0};
#endif

static DSDP_INT schurMatSetupBlock( symM *M, DSDP_INT blockid ) {
    // Set up a block of the Schur matrix
    DSDP_INT retcode = DSDP_RETCODE_OK;
#ifdef SCHUR_PROFILER
    clock_t start;
    memset(schurTime, 0, sizeof(double) * 5);
#endif
    for (DSDP_INT i = 0; i < M->m + 1; ++i) {
#ifdef SCHUR_PROFILER
        start = clock();
        switch (M->MX[blockid][i]) {
            case SCHUR_M1: schurM1rowSetup(M, blockid, i);
                schurTime[0] += (double) (clock() - start) / CLOCKS_PER_SEC; break;
            case SCHUR_M2: schurM2rowSetup(M, blockid, i);
                schurTime[1] += (double) (clock() - start) / CLOCKS_PER_SEC; break;
            case SCHUR_M3: schurM3rowSetup(M, blockid, i);
                schurTime[2] += (double) (clock() - start) / CLOCKS_PER_SEC; break;
            case SCHUR_M4: schurM4rowSetup(M, blockid, i);
                schurTime[3] += (double) (clock() - start) / CLOCKS_PER_SEC; break;
            case SCHUR_M5: schurM5rowSetup(M, blockid, i);
                schurTime[4] += (double) (clock() - start) / CLOCKS_PER_SEC; break;
            default: error(etype, "Invalid technique type. \n"); break;
        }
#else
        switch (M->MX[blockid][i]) {
            case SCHUR_M1: schurM1rowSetup(M, blockid, i); break;
            case SCHUR_M2: schurM2rowSetup(M, blockid, i); break;
            case SCHUR_M3: schurM3rowSetup(M, blockid, i); break;
            case SCHUR_M4: schurM4rowSetup(M, blockid, i); break;
            case SCHUR_M5: schurM5rowSetup(M, blockid, i); break;
            default: error(etype, "Invalid technique type. \n"); break;
        }
#endif
    }
#ifdef SCHUR_PROFILER
        printf("|1: %f 2: %f 3: %f 4: %f 5: %f \n", schurTime[0], schurTime[1], schurTime[2], schurTime[3], schurTime[4]);
#endif
    return retcode;
}

static DSDP_INT schurspM2Setup( symM *M, DSDP_INT blockid, DSDP_INT idx ) {
    // Apply M2 technique to set up a sub-matrix of the sparse schur matrix. Sparse version of schurM2Setup.
    DSDP_INT retcode = DSDP_RETCODE_OK;
    DSDP_INT i, j, k, r, m = M->m, nzmats = M->Adata[blockid]->nzeroMat, computeC = FALSE, rank, \
    gold = *M->gold, *nzidx = M->Adata[blockid]->nzIdx, nzAmats = M->Adata[blockid]->nnzAmat;
    rkMat *factor, *rkaux = M->rkaux[blockid]; r1Mat *r1aux = NULL;
    double *val, res, coeff, Ry = *M->Ry;
    
    k = nzidx[idx];
    switch (M->Adata[blockid]->types[idx]) {
        case MAT_TYPE_SPARSE: factor = spsMatGetFactor(M->Adata[blockid]->sdpData[idx]); break;
        case MAT_TYPE_RANKK : factor = M->Adata[blockid]->sdpData[idx]; break;
        default: error(etype, "Invalid matrix type. \n"); break;
    }
    
    if (k < m) {
        val = &M->asinv->x[k];
    } else {
        val = M->csinv; if (!gold) { return retcode; }
        computeC = TRUE;
    }
    
    res = SinvRkSinv(M->S[blockid], factor, rkaux);
    *val += res; return_if_not_in_A_or_not_building_hsd
    rank = rkMatGetRank(factor);
    
    /* Start sparse M2 */
    if (computeC) {
        assert( *M->phaseA );
        val = M->csinvrysinv; rkMatdiagTrace(rkaux, Ry, &res); *val += res;
        for (r = 0; r < rank; ++i) {
            r1aux = rkMatGetBase(rkaux, r); coeff = r1aux->sign;
            for (i = idx; i < nzmats; ++i) {
                j = nzidx[i];
                if (j == m) {
                    if (!*M->buildhsd) { continue; }
                    val = M->csinvcsinv;
                } else {
                    val = &M->asinvcsinv->x[j];
                }
                switch (M->Adata[blockid]->types[i]) {
                    case MAT_TYPE_DENSE: res = denseMatxTAx(M->Adata[blockid]->sdpData[i], M->schurAux, r1aux->x);
                        res *= coeff; break;
                    case MAT_TYPE_SPARSE: res = spsMatxTAx(M->Adata[blockid]->sdpData[i], r1aux->x);
                        res *= coeff; break;
                    case MAT_TYPE_RANKK:
                        res = r1Matr1Trace(((rkMat *) M->Adata[blockid]->sdpData[i])->data[0], r1aux); break;
                    default: error(etype, "Invalid matrix type. \n"); break;
                }
                *val += res;
            }
        }
    } else {
        if (*M->phaseA) { val = &M->asinvrysinv->x[k]; rkMatdiagTrace(rkaux, Ry, &res); *val += res; }
        double *array = M->M->spsM->x;
        for (r = 0; r < rank; ++r) {
            r1aux = rkMatGetBase(rkaux, r); coeff = r1aux->sign;
            for (i = idx; i < nzmats; ++i) {
                j = nzidx[i];
                if (j == m) {
                    continue_if_not_in_A_or_not_building_hsd
                    val = &M->asinvcsinv->x[k];
                } else {
                    val = &array[packIdx(M->Adata[blockid]->schurspIdx, nzAmats, i, idx)];
                }
                
                switch (M->Adata[blockid]->types[i]) {
                    case MAT_TYPE_DENSE : res = denseMatxTAx(M->Adata[blockid]->sdpData[i], M->schurAux, r1aux->x);
                        res *= coeff; break;
                    case MAT_TYPE_SPARSE: res = spsMatxTAx(M->Adata[blockid]->sdpData[i], r1aux->x);
                        res *= coeff; break;
                    case MAT_TYPE_RANKK :
                        res = r1Matr1Trace(((rkMat *) M->Adata[blockid]->sdpData[i])->data[0], r1aux);
                        break;
                    default             : error(etype, "Invalid matrix type. \n"); break;
                }
                *val += res;
            }
        }
    }
    
    /* End sparse M2 */
    return retcode;
}

static DSDP_INT schurspM3Setup( symM *M, DSDP_INT blockid, DSDP_INT idx ) {
    // Apply M3 technique to set up a sub-matrix of the sparse schur matrix. Sparse version of schurM3Setup.
    DSDP_INT retcode = DSDP_RETCODE_OK;
    DSDP_INT i, j, k, m = M->m, computeC = FALSE, nzmats = M->Adata[blockid]->nzeroMat, \
             gold = *M->gold, *nzidx = M->Adata[blockid]->nzIdx, nzAmats = M->Adata[blockid]->nnzAmat;
    dsMat *B = M->B[blockid];
    double *val, res, Ry = *M->Ry, *aux = M->schurAux, *Sinv = M->Sinv[blockid]; k = nzidx[idx];
    
    if (k < m) {
        val = &M->asinv->x[k];
    } else {
        val = M->csinv; if (!gold) { return retcode; }
        computeC = TRUE;
    }
    
    denseMatReset(B);
    switch (M->Adata[blockid]->types[idx]) {
        case MAT_TYPE_DENSE: res = SinvDsSinv(Sinv, aux, M->Adata[blockid]->sdpData[idx], B); break;
        default: error(etype, "Invalid matrix type. \n"); break;
    }
    
    *val += res; return_if_not_in_A_or_not_building_hsd
    
    /* Start sparse M3 */
    if (computeC) {
        val = M->csinvrysinv; *val += denseDiagTrace(B, Ry);
        for (i = idx; i < nzmats; ++i) {
            j = nzidx[i]; if (M->Adata[blockid]->types[j] == MAT_TYPE_ZERO) { continue; }
            if (j == m) {
                if (!*M->buildhsd) { continue; } val = M->csinvcsinv;
            } else {
                val = &M->asinvcsinv->x[j];
            }
            switch (M->Adata[blockid]->types[i]) {
                case MAT_TYPE_DENSE : denseDsTrace(B, M->Adata[blockid]->sdpData[i], &res); break;
                case MAT_TYPE_SPARSE: denseSpsTrace(B, M->Adata[blockid]->sdpData[i], &res); break;
                case MAT_TYPE_RANKK : rkMatdenseTrace(M->Adata[blockid]->sdpData[i], B, &res); break;
                default             : error(etype, "Invalid matrix type. \n"); break;
            }
            *val += res;
        }
    } else {
        if (*M->phaseA) { val = &M->asinvrysinv->x[k]; *val += denseDiagTrace(B, Ry); }
        double *array = M->M->spsM->x;
        for (i = idx; i < nzmats; ++i) {
            j = nzidx[i];
            if (j == m) {
                continue_if_not_in_A_or_not_building_hsd
                val = &M->asinvcsinv->x[k];
            } else {
                val = &array[packIdx(M->Adata[blockid]->schurspIdx, nzAmats, i, idx)];
            }
            switch (M->Adata[blockid]->types[i]) {
                case MAT_TYPE_DENSE : denseDsTrace(B, M->Adata[blockid]->sdpData[i], &res); break;
                case MAT_TYPE_SPARSE: denseSpsTrace(B, M->Adata[blockid]->sdpData[i], &res); break;
                case MAT_TYPE_RANKK : rkMatdenseTrace(M->Adata[blockid]->sdpData[i], B, &res); break;
                default             : error(etype, "Invalid matrix type. \n"); break;
            }
            *val += res;
        }
    }
    
    /* End sparse M3 */
    return retcode;
}

static DSDP_INT schurspM5Setup( symM *M, DSDP_INT blockid, DSDP_INT idx ) {
    // Apply M5 technique to set up a sub-matrix of the sparse schur matrix. Sparse version of schurM5Setup.
    DSDP_INT retcode = DSDP_RETCODE_OK;
    DSDP_INT i, j, k, m = M->m, computeC = FALSE, gold = *M->gold, *nzidx = M->Adata[blockid]->nzIdx, \
             nzmats = M->Adata[blockid]->nzeroMat, nzAmats = M->Adata[blockid]->nnzAmat;
    
    DSDP_INT *types = M->Adata[blockid]->types;
    double *Sinv = M->Sinv[blockid], *val, res, tmp, Ry = *M->Ry; k = nzidx[idx];
    void **data = M->Adata[blockid]->sdpData;
    
    /* Start sparse M5 */
    if (k < m) {
        val = &M->asinv->x[k];
    } else {
        val = M->csinv; if (!gold) { return retcode; }
        computeC = TRUE;
    }
    
    switch (types[idx]) {
        case MAT_TYPE_SPARSE: res = spsRySinv(data[idx], Sinv, &tmp, Ry); break;
        default: assert(FALSE); break;
    }
    
    *val += tmp; return_if_not_in_A_or_not_building_hsd
    if (computeC) {
        *M->csinvrysinv += res;
        for (i = idx; i < nzmats; ++i) {
            j = nzidx[i]; if (M->Adata[blockid]->types[j] == MAT_TYPE_ZERO) { continue; }
            if (j == m) {
                if (!*M->buildhsd) { continue; } val = M->csinvcsinv;
            } else {
                val = &M->asinvcsinv->x[j];
            }
            *val += schurM5MatAux(MAT_TYPE_SPARSE, data[idx], types[i], data[i], Sinv);
        }
    } else {
        M->asinvrysinv->x[k] += res;
        double *array = M->M->spsM->x;
        for (i = idx; i < nzmats; ++i) {
            j = nzidx[i];
            if (j == m) {
                continue_if_not_in_A_or_not_building_hsd
                val = &M->asinvcsinv->x[k];
            } else {
                val = &array[packIdx(M->Adata[blockid]->schurspIdx, nzAmats, i, idx)];
            }
            *val += schurM5MatAux(MAT_TYPE_SPARSE, data[idx], types[i], data[i], Sinv);
        }
    }
    
    /* End sparse M5 */
    return retcode;
}

static DSDP_INT schurMatspsSetupBlock( symM *M, DSDP_INT blockid ) {
    // Set up a block of the Schur matrix. Way of computing is heuristically chosen from M2 M3 and M5
    DSDP_INT retcode = DSDP_RETCODE_OK;
    assert( M->M->stype = SCHUR_TYPE_SPARSE );
    DSDP_INT Mk = SCHUR_M5, nnzmats, i;
    sdpMat *Adata = M->Adata[blockid]; nnzmats = Adata->nzeroMat;
    spsMat *spdata = NULL;
    
    // Begin setup
    for (i = 0; i < nnzmats; ++i) {
        // Choose schur technique
        if (Adata->types[i] == MAT_TYPE_SPARSE) {
            spdata = Adata->sdpData[i];
            Mk = (spsMatGetRank(spdata) <= 10) ? SCHUR_M2 : SCHUR_M5;
        } else if (Adata->types[i] == MAT_TYPE_DENSE) {
            Mk = SCHUR_M3;
        } else {
            Mk = SCHUR_M2; assert( Adata->types[i] == MAT_TYPE_RANKK );
        }
        
        switch (Mk) {
            case SCHUR_M2: schurspM2Setup(M, blockid, i); break;
            case SCHUR_M3: schurspM3Setup(M, blockid, i); break;
            case SCHUR_M5: schurspM5Setup(M, blockid, i); break;
            default: assert(FALSE); break;
        }
    }

    // Finish setup
    return retcode;
}

static DSDP_INT schurMatReorder( symM *M ) {
    
    DSDP_INT retcode = DSDP_RETCODE_OK;
    DSDP_INT maxblock = 0;
    // Information of all the matrices must have been supplied
    for (DSDP_INT i = 0, dim; i < M->nblock; ++i) {
        retcode = schurBlockAnalysis(M->Adata[i], M->perms[i],
                                     M->MX[i], &M->useTwo[i]);
        dim = M->Adata[i]->dimS; maxblock = MAX(dim, maxblock);
        M->Sinv[i] = (double *) calloc(dim * dim, sizeof(double));
    }
    
    M->schurAux = (double *) calloc(MAX(M->m, maxblock * maxblock), sizeof(double));
    M->Mready = TRUE;
    
    DSDP_INT m1 = 0, m2 = 0, m3 = 0, m4 = 0, m5 = 0;
    for (DSDP_INT i = 0, j; i < M->nblock; ++i) {
        for (j = 0; j < M->m + 1; ++j) {
            switch (M->MX[i][j]) {
                case SCHUR_M1: m1 += 1; break;
                case SCHUR_M2: m2 += 1; break;
                case SCHUR_M3: m3 += 1; break;
                case SCHUR_M4: m4 += 1; break;
                case SCHUR_M5: m5 += 1; break;
                default: assert( FALSE ); break;
            }
        }
    }
    printf("|    Schur Re-ordering: M1: %d  M2: %d  M3: %d  M4: %d  M5: %d \n", m1, m2, m3, m4, m5);
    return retcode;
}

#define SCHUR_BUFFER 100000
static DSDP_INT schurMatSpsReorder( symM *M ) {
    // Implement the sparse Schur matrix re-ordering technique
    DSDP_INT retcode = DSDP_RETCODE_OK;
    
    assert( M->M->stype == SCHUR_TYPE_SPARSE );
    DSDP_INT m = M->m, nblock = M->nblock;
    DSDP_INT *start = M->useTwo, i, j, k, maxblock = 0;
    sdpMat **Adata = M->Adata;
    memset(start, 0, sizeof(DSDP_INT) * nblock);
    DSDP_INT *Mp, *Mi, *colNnz, schurNnz = 0, memNnz = MAX(2 * m, SCHUR_BUFFER);
    
    colNnz = (DSDP_INT *) calloc(m + 1, sizeof(DSDP_INT));
    Mp = (DSDP_INT *) calloc(m + 1, sizeof(DSDP_INT));
    Mi = (DSDP_INT *) calloc(memNnz, sizeof(DSDP_INT));
    
    // Allocate memory for Sinv and auxiliary
    for (i = 0; i < M->nblock; ++i) {
        j = M->Adata[i]->dimS; maxblock = MAX(maxblock, j);
        M->Sinv[i] = (double *) calloc(j * j, sizeof(double));
    }
        
    M->schurAux = (double *) calloc(MAX(M->m, maxblock * maxblock), sizeof(double));
    
    // Go through columns
    for (i = 0; i < m; ++i) {
        memset(colNnz, 0, sizeof(DSDP_INT) * (m + 1));
        for (k = 0; k < nblock; ++k) {
            start[k] = sdpMatScatterNnz(Adata[k], start[k], i, colNnz);
        }
        // Accumulate the index
        j = k = 0; if (colNnz[j]) { Mi[schurNnz + k] = j; k += 1;}
        for (j = 1; j < m; ++j) {
            if (colNnz[j]) { Mi[schurNnz + k] = j; k += 1;}
            colNnz[j] = colNnz[j] + colNnz[j - 1];
        }
        // Associate the index
        for (k = 0; k < nblock; ++k) {
            sdpMatSetSchurIndex(Adata[k], start[k], i, colNnz, schurNnz);
        }
        // Count nnzs
        schurNnz += colNnz[m - 1];
        if (memNnz <= schurNnz + m) {
            memNnz += m; Mi = (DSDP_INT *) realloc(Mi, sizeof(DSDP_INT) * memNnz);
        }
        Mp[i + 1] = schurNnz;
    }
    
    // Finish setup
    printf("| - Schur sparsity: %d out of %d. \n", schurNnz, nsym(m));
    Mi = (DSDP_INT *) realloc(Mi, sizeof(DSDP_INT) * schurNnz);
    retcode = spsMatAllocData(M->M->spsM, m, schurNnz);
    memcpy(M->M->spsM->p, Mp, sizeof(DSDP_INT) * (m + 1));
    memcpy(M->M->spsM->i, Mi, sizeof(DSDP_INT) * schurNnz);
    spsMatSymbolic(M->M->spsM);
    
    // Associate the diagonal elements
    for (i = 0; i < m; ++i) {
        M->M->diag[i] = &M->M->spsM->x[M->M->spsM->p[i]];
    }
    
    M->Mready = TRUE;
    DSDP_FREE(colNnz); DSDP_FREE(Mp); DSDP_FREE(Mi);
    return retcode;
}

static void asinvdsSetup( symM *M, vec *asinv ) {
    // Set up dense asinv for the corrector step
    DSDP_INT m = M->m; double res = 0.0, *Sinv; void *data;
    for (DSDP_INT blockid = 0, i = 0; blockid < M->nblock; ++blockid) {
        Sinv = M->Sinv[blockid];
        for (i = 0; i < m; ++i) {
            data = M->Adata[blockid]->sdpData[i];
            switch (M->Adata[blockid]->types[i]) {
                case MAT_TYPE_ZERO  : continue; break;
                case MAT_TYPE_SPARSE: res = spsFullTrace((spsMat *) data, Sinv); break;
                case MAT_TYPE_DENSE : res = denseFullTrace((dsMat *) data, Sinv); break;
                case MAT_TYPE_RANKK : res = r1MatFullTrace(((rkMat *) data)->data[0], Sinv, M->schurAux); break;
                default             : assert( FALSE ); break;
            }
            asinv->x[i] += res;
        }
    }
}

static void asinvspSetup( symM *M, vec *asinv ) {
    // Set up sparse asinv for the corrector step
    double res = 0.0, *Sinv; void *data;
    for (DSDP_INT blockid = 0, i = 0, j; blockid < M->nblock; ++blockid) {
        Sinv = M->Sinv[blockid];
        for (i = 0; i < M->Adata[blockid]->nnzAmat; ++i) {
            data = M->Adata[blockid]->sdpData[i]; j = M->Adata[blockid]->nzIdx[i];
            switch (M->Adata[blockid]->types[i]) {
                case MAT_TYPE_SPARSE: res = spsFullTrace((spsMat *) data, Sinv); break;
                case MAT_TYPE_DENSE : res = denseFullTrace((dsMat *) data, Sinv); break;
                case MAT_TYPE_RANKK : res = r1MatFullTrace(((rkMat *) data)->data[0], Sinv, M->schurAux); break;
                default             : assert( FALSE ); break;
            }
            asinv->x[j] += res;
        }
    }
}

static void arysinvdsSetup( symM *M, vec *asinv, vec *asinvrysinv ) {
    // Set up both asinv and asinvrysinv for the corrector step, M5 routines are used
    DSDP_INT m = M->m; double res1 = 0.0, res2 = 0.0, *Sinv, Ry = *M->Ry; void *data;
    for (DSDP_INT blockid = 0, i = 0; blockid < M->nblock; ++blockid) {
        Sinv = M->Sinv[blockid];
        for (i = 0; i < m; ++i) {
            data = M->Adata[blockid]->sdpData[i];
            switch (M->Adata[blockid]->types[i]) {
                case MAT_TYPE_ZERO  : continue; break;
                case MAT_TYPE_SPARSE: res2 = spsRySinv((spsMat *) data, Sinv, &res1, Ry); break;
                case MAT_TYPE_DENSE : res2 = denseSinvSolve2(Sinv, (dsMat *) data, &res1, Ry); break;
                case MAT_TYPE_RANKK : res2 = r1RySinv(((rkMat *) data)->data[0], Sinv, &res1, Ry, M->schurAux); break;
                default             : assert( FALSE ); break;
            }
            asinv->x[i] += res1; asinvrysinv->x[i] += res2;
        }
    }
}

static void arysinvspSetup( symM *M, vec *asinv, vec *asinvrysinv ) {
    // Set up both asinv and asinvrysinv for the corrector step, M5 routines are used
    double res1 = 0.0, res2 = 0.0, *Sinv, Ry = *M->Ry; void *data;
    for (DSDP_INT blockid = 0, i = 0, j; blockid < M->nblock; ++blockid) {
        Sinv = M->Sinv[blockid];
        for (i = 0; i < M->Adata[blockid]->nnzAmat; ++i) {
            data = M->Adata[blockid]->sdpData[i]; j = M->Adata[blockid]->nzIdx[i];
            switch (M->Adata[blockid]->types[i]) {
                case MAT_TYPE_SPARSE: res2 = spsRySinv((spsMat *) data, Sinv, &res1, Ry); break;
                case MAT_TYPE_DENSE :  res2 = denseSinvSolve2(Sinv, (dsMat *) data, &res1, Ry); break;
                case MAT_TYPE_RANKK :  res2 = r1RySinv(((rkMat *) data)->data[0], Sinv, &res1, Ry, M->schurAux); break;
                default             : assert( FALSE ); break;
            }
            asinv->x[j] += res1; asinvrysinv->x[j] += res2;
        }
    }
}

extern DSDP_INT symSchurMatInit( symM *M ) {
    
    DSDP_INT retcode = DSDP_RETCODE_OK;
    
    M->m = 0; M->nblock = 0; M->Mready = FALSE; M->useTwo = NULL;
    M->perms = NULL; M->MX = NULL; M->S = NULL; M->B = NULL;
    M->Adata = NULL; M->M = NULL; M->asinv = NULL; M->Ry = NULL;
    M->asinvrysinv = NULL; M->asinvcsinv = NULL; M->csinvrysinv = NULL;
    M->csinv = NULL; M->csinvcsinv = NULL; M->Sinv = NULL; M->rkaux = NULL;
    M->scaler = 0.0; M->buildhsd = NULL; M->rysinv = NULL; M->gold = NULL;
    M->lpSet = NULL;
    
    return retcode;
}

extern void symSchurMatSetDim( symM *M, DSDP_INT m, DSDP_INT nblock ) {
    
    M->m = m; M->nblock = nblock;
}

extern DSDP_INT symSchurMatAlloc( symM *M ) {
    // Allocate memory for internal perms
    DSDP_INT retcode = DSDP_RETCODE_OK;
    
    M->perms  = (DSDP_INT **) calloc(M->nblock, sizeof(DSDP_INT *));
    M->MX     = (DSDP_INT **) calloc(M->nblock, sizeof(DSDP_INT *));
    M->useTwo = (DSDP_INT  *) calloc(M->nblock, sizeof(DSDP_INT));
    M->Sinv   = (double   **) calloc(M->nblock, sizeof(double *));
    M->scaler = 0.0;
    
    if (!M->perms || !M->MX || !M->useTwo || !M->Sinv) {
        retcode = DSDP_RETCODE_FAILED;
        return retcode;
    }
    
    for (DSDP_INT i = 0; i < M->nblock; ++i) {
        M->perms[i] = NULL; M->MX[i] = NULL; M->Sinv[i] = NULL;
    }

    return retcode;
}

extern void symSchurMatRegister( symM *M, spsMat **S, dsMat **B, sdpMat **Adata, schurMat *Msdp,
                                     vec *asinv, vec *asinvrysinv, vec *asinvcsinv, double *csinvrysinv,
                                     double *csinv, double *csinvcsinv, double *rysinv, double *Ry,
                                     rkMat **rkaux, DSDP_INT *phaseA, DSDP_INT *buildHSD, DSDP_INT *gold,
                                     DSDP_INT *lpSet ) {
    // Register information into M
    M->S = S; M->B = B; M->Adata = Adata; M->M = Msdp;
    M->asinv = asinv; M->asinvrysinv = asinvrysinv; M->Ry = Ry;
    M->asinvcsinv = asinvcsinv; M->csinvrysinv = csinvrysinv;
    M->csinv = csinv; M->csinvcsinv = csinvcsinv; M->rkaux = rkaux;
    M->phaseA = phaseA; M->buildhsd = buildHSD; M->rysinv = rysinv;
    M->gold = gold; M->lpSet = lpSet;
}

extern DSDP_INT symSchurMatFree( symM *M ) {
    
    DSDP_INT retcode = DSDP_RETCODE_OK;
    for (DSDP_INT i = 0; i < M->nblock; ++i) {
        DSDP_FREE(M->perms[i]); DSDP_FREE(M->MX[i]); DSDP_FREE(M->Sinv[i]);
    }
    
    DSDP_FREE(M->perms); DSDP_FREE(M->MX); DSDP_FREE(M->Sinv); DSDP_FREE(M->schurAux);
    M->Mready     = FALSE; M->m           = 0;    M->nblock = 0; M->perms = NULL;
    M->MX         = NULL;  M->S           = NULL; M->B = NULL;   M->Adata = NULL;
    M->M          = NULL;  M->asinv       = NULL; M->asinvrysinv = NULL;
    M->asinvcsinv = NULL;  M->csinvrysinv = NULL; M->useTwo = FALSE;
    M->csinv      = NULL;  M->csinvcsinv  = NULL; M->scaler = 0.0;
    M->buildhsd   = NULL;  M->rysinv      = NULL; M->gold = NULL;
    M->lpSet      = NULL;
    
    return retcode;
}

// Schur matrix re-ordering heuristic
extern DSDP_INT DSDPSchurReorder( symM *M ) {
    // Invoke the re-ordering heuristic for dense Schur matrix setup
    return (M->M->stype == SCHUR_TYPE_SPARSE) ? schurMatSpsReorder(M) : schurMatReorder(M);
}

extern DSDP_INT DSDPCheckSchurType( symM *M ) {
    // Determine the sparsity pattern of the Schur matrix
    DSDP_INT retcode = DSDP_RETCODE_OK;
    
    if (M->nblock <= 100 || *M->lpSet) {
        for (DSDP_INT i = 0; i < M->nblock; ++i) {
            M->perms[i] = (DSDP_INT *) calloc(M->m + 1, sizeof(DSDP_INT));
            M->MX[i] = (DSDP_INT *) calloc(M->m + 1, sizeof(DSDP_INT));
            M->Sinv[i] = NULL;
        }
        retcode = schurMatselectType(M->M, SCHUR_TYPE_DENSE);
    } else {
        printf("| - Found many blocks. Trying sparse Schur Matrix \n");
        DSDP_INT minblknnz = M->m, i; sdpMat *data;
        for (i = 0; i < M->nblock; ++i) {
            data = (sdpMat *) M->Adata[i];
            minblknnz = MIN(minblknnz, data->nzeroMat);
        }
        if (minblknnz >= M->m * 0.8) {
            printf("| - Maximum block density %3.2f%%."
                   " Using sparse Schur matrix \n",
                   100.0 - (double) 100 * minblknnz / M->m);
        }
        retcode = schurMatselectType(M->M, SCHUR_TYPE_SPARSE);
        // Contract data
        for (i = 0; i < M->nblock; ++i) {
            schurBlockContract(M->Adata[i]);
        }
    }
    return retcode;
}

// static double schurtime = 0.0;

extern DSDP_INT DSDPSchurSetup( symM *M ) {
    // Routine for setting up the Schur complement matrix M
    /*
     In phase A where HSD is to get solved, we need
     
     1. Conventional components <S^-1 * A_i * S^-1, A_j>
     2. HSD components <S^-1 * A_i * S^-1, C>
     3. Residual components <S^-1 * A_i * S^-1, I> * Ry
     
     We set up the Schur matrix as if there are (m + 1) {A_i} matrices and the
     residual related terms are set up during the set up (in different parts of the techniques
     */
    DSDP_INT retcode = DSDP_RETCODE_OK;
    
    // Clean up current values and invert blocks
    schurMatCleanup(M); schurMatGetSinv(M);
    
    if (M->M->stype == SCHUR_TYPE_DENSE) {
        for (DSDP_INT i = 0; i < M->nblock; ++i) {
            schurMatSetupBlock(M, i);
        }
    } else {
        for (DSDP_INT i = 0; i < M->nblock; ++i) {
            schurMatspsSetupBlock(M, i);
        }
    }
    
    return retcode;
}

extern void asinvSetup( symM *M ) {
    // Set up asinv for the corrector step
    vec *asinv = M->asinv; vec_reset(asinv);
    schurMatSinvCleanup(M); schurMatGetSinv(M);
    if (M->M->stype == SCHUR_TYPE_SPARSE) {
        asinvspSetup(M, asinv);
    } else {
        asinvdsSetup(M, asinv);
    }
}

extern void arysinvSetup( symM *M ) {
    // Set up both asinv and asinvrysinv for the corrector step, M5 routines are used
    vec *asinv = M->asinv, *asinvrysinv = M->asinvrysinv;
    vec_reset(asinv); vec_reset(asinvrysinv);
    schurMatSinvCleanup(M); schurMatGetSinv(M);
    if (M->M->stype == SCHUR_TYPE_SPARSE) {
        arysinvspSetup(M, asinv, asinvrysinv);
    } else {
        arysinvdsSetup(M, asinv, asinvrysinv);
    }
}
