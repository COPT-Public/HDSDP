#include "symschur.h"
#include "dsdputils.h"
#include "dsdpsolver.h"
#include "dsdpSort.h"

// Implement the advanced Schur matrix setup in DSDP using M1, M2, M3, M4 techniques
static char etype[] = "Advanced Schur matrix setup";

static DSDP_INT getScore( DSDP_INT *ranks, DSDP_INT *nnzs, DSDP_INT *perm,
                        DSDP_INT m, DSDP_INT n, DSDP_INT idx, DSDP_INT *MX,
                        double *M1M2, double *M1M5 ) {
    // Compute the estimated number of multiplications for technique M1 to M5
    
    DSDP_INT i, j, k = perm[idx], npack = nsym(n);
    double d1, d2, d3, d4, d5, nsqr = n * n, ncbe = n * nsqr, sumnnz = 0.0, rsigma = ranks[k], best = SCHUR_M1;
    double m1m5 = 0.0, m1m2 = 0.0;
    // On exit, MX is set 10 * MX12 + M12345

    for (i = idx; i < m + 1; ++i) {
        sumnnz += nnzs[i];
    }
  
    d1 = rsigma * (2 * nsqr) + KAPPA * sumnnz;
    d2 = rsigma * (nsqr + 2 * KAPPA * sumnnz);
    // d1 = rsigma * (nsqr + 3 * npack) + 2 * KAPPA * sumnnz;
    // d2 = rsigma * (nsqr + 3 * KAPPA * sumnnz);
    d3 = n * KAPPA * nnzs[idx] + ncbe + KAPPA * sumnnz + ncbe / m;
    d5 = KAPPA * (2 * KAPPA * nnzs[idx] + 1) * sumnnz + ncbe / m;
    
    // Turn off d4 since not implemented
    d4 = n * KAPPA * nnzs[idx] + KAPPA * (n + 1) * sumnnz + ncbe / m;
    
        
    // printf("d1: %f   d2: %f \n", d1, d2);
    
    // Get better of d1 and d2
    j = 0;
    if (d1 < d2) {
        best = SCHUR_M1; m1m2 = d1; m1m5 = d1;
    } else {
        best = SCHUR_M2; m1m2 = d2; m1m5 = d2;
    }
    
    j += best * 10;
    
    if (d3 < m1m5) { best = SCHUR_M3; m1m5 = d3; }
    // if (d4 < m1m5) { best = SCHUR_M4; m1m5 = d4; }
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
             nzeroMat = Adata->nzeroMat, m = Adata->dimy, n = Adata->dimS;
    
    assert( ndsMat + nspsMat + nrkMat + nzeroMat == m + 1 );
    
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
        nnzs[k] = nsym(nnzs[k]);
        ranks[k] = 1;
    }
    
    // Sort f_1, ..., f_m in descending order
    dsdpSort2(permk, nnzs, 0, m);
    
    // Determine strategies for setting up the Schur matrix
    for (i = 0; i < m + 1; ++i) {
        getScore(ranks, nnzs, permk, m, // m here is correct
                 n, i,
                 &MXk[i], &scoreM1M2[i], &scoreM1M5[i]);
    }
    
    // Determine whether to use M3 to M5
    for (i = 0; i < m + 1; ++i) {
        ksiM1M2 += scoreM1M2[i];
        ksiM1M5 += scoreM1M5[i];
    }
    
    double ncbe = n;
    ncbe *= n; ncbe *= n;
    if (ksiM1M2 <= ksiM1M5 + ncbe) {
        M1M2 = TRUE;
    } else {
        M1M2 = FALSE;
    }
    
    /* For Debugging purpose only */
    M1M2 = FALSE;
    
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
    
    DSDP_FREE(scoreM1M2); DSDP_FREE(scoreM1M5);
    DSDP_FREE(ranks); DSDP_FREE(nnzs);
    
    return retcode;
}

static void schurMatCleanup( DSDPSchur *M ) {
    // Clean up the arrays related to the Schur matrix
    denseMatReset(M->M); vec_reset(M->asinv);
    if (*M->phaseA) {
        vec_reset(M->asinvcsinv); vec_reset(M->asinvrysinv);
        *M->csinvrysinv = 0.0; *M->csinv = 0.0; *M->csinvcsinv = 0.0;
    }
    // Set up the inverse of matrices when necessary
    for (DSDP_INT i = 0, j, n; i < M->nblock; ++i) {
        if (M->useTwo[i]) {
            continue;
        } else {
            n = M->Adata[i]->dimS; memset(M->Sinv[i], 0, sizeof(double) * n * n);
            for (j = 0; j < n; ++j) { M->Sinv[i][j * n + j] = 1.0; }
        }
    }
}

static void schurMatGetSinv( DSDPSchur *M ) {
    // Compute inverse of the dual matrix when M3, M4 or M5 techniques are used
    double maxdiag = 0.0; M->scaler = 1.0;
    DSDP_INT n, nsqr, one = 1;
    
    // Compute S^-1 and scale it to improve numerical stability
    for (DSDP_INT i = 0; i < M->nblock; ++i) {
        if (!M->useTwo[i]) {
            spsMatInverse(M->S[i], M->Sinv[i], M->schurAux);
            n = M->S[i]->dim;
            /* for (DSDP_INT j = 0; j < n; ++j) {
                    maxdiag = MAX(M->Sinv[i][j * n + j], maxdiag);
               } */
        }
    }
    
    if (maxdiag > 1e+03 && (FALSE)) {
        M->scaler = sqrt(maxdiag);
        for (DSDP_INT i = 0; i < M->nblock; ++i) {
            if (!M->useTwo[i]) {
                n = M->S[i]->dim; nsqr = n * n;
                drscl(&nsqr, &M->scaler, M->Sinv[i], &one);
            }
        }
    }
}

static DSDP_INT schurM1rowSetup( DSDPSchur *M, DSDP_INT blockid, DSDP_INT row ) {
    // Apply M1 technique to setup a row of the Schur complement matrix
    /*
     M1 technique first applies rank-1 update to get B_i = <S^-1 * A_i * S^-1>
     and then computes M_{ij} = <B_i, A_j>
     
     When each row is set up, not only M is computed, but also
        
     csinvrysinv, asinvrysinv, csinvcsinv, asinvcsinv, csinv, asinv
     
     when available
    */
    
    DSDP_INT retcode = DSDP_RETCODE_OK;
    DSDP_INT i, j, k, m = M->m, *perm = M->perms[blockid], computeC = FALSE;
    dsMat *B = M->B[blockid]; rkMat *factor, *rkaux = M->rkaux[blockid];
    double *val, res, Ry = *M->Ry; k = perm[row];
    
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
        if (!*M->phaseA) { return retcode; }
        val = M->csinv; computeC = TRUE;
    }
    res = spsSinvRkSinvSolve(M->S[blockid], factor, rkaux); *val += res;
    
    // Get B through rank-one update
    denseMatReset(B); rkMatdenseUpdate(B, rkaux);
    
    /* Start M1 */
    if (computeC) {
        assert( *M->phaseA );
        val = M->csinvrysinv; denseDiagTrace(B, Ry, &res); *val += res;
        for (i = row; i < m + 1; ++i) {
            j = perm[i]; val = (j == m) ? M->csinvcsinv : &M->asinvcsinv->x[j];
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
        if (*M->phaseA) { val = &M->asinvrysinv->x[k]; denseDiagTrace(B, Ry, &res); *val += res; }
        double *array = M->M->array;
        for (i = row; i < m + 1; ++i) {
            j = perm[i];
            if (j == m) {
                if (!*M->phaseA) { continue; }
                val = &M->asinvcsinv->x[k];
            } else {
                val = (j > k) ? (&packIdx(array, m, j, k)) : (&packIdx(array, m, k, j));
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

static DSDP_INT schurM2rowSetup( DSDPSchur *M, DSDP_INT blockid, DSDP_INT row ) {
    // Apply M2 technique to setup a row of the Schur complement matrix
    /*
     M2 technique does not form <S^-1 * A_i * S^-1> but compute the value of quadratic forms
     directly.
     
     When each row is set up, now only M is computed, but also
     
     csinvrysinv, asinvrysinv, csinvcsinv, asinvcsinv, csinv, asinv
     
     when available
    */
    
    DSDP_INT retcode = DSDP_RETCODE_OK;
    
    DSDP_INT i, j, k, r, m = M->m, *perm = M->perms[blockid], computeC = FALSE, rank;
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
        if (!*M->phaseA) { return retcode; }
        val = M->csinv; computeC = TRUE;
    }
    res = spsSinvRkSinvSolve(M->S[blockid], factor, rkaux); *val += res;
    rank = rkMatGetRank(factor);
    
    /* Start M2 */
    if (computeC) {
        assert(*M->phaseA);
        val = M->csinvrysinv; rkMatdiagTrace(rkaux, Ry, &res); *val += res;
        for (r = 0; r < rank; ++r) {
            r1aux = rkMatGetBase(rkaux, r); coeff = r1aux->sign;
            for (i = row; i < m + 1; ++i) {
                j = perm[i]; if (M->Adata[blockid]->types[j] == MAT_TYPE_ZERO) { continue; }
                val = (j == m) ? M->csinvcsinv : &M->asinvcsinv->x[j];
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
        double *array = M->M->array;
        for (r = 0; r < rank; ++r) {
            r1aux = rkMatGetBase(rkaux, r); coeff = r1aux->sign;
            for (i = row; i < m + 1; ++i) {
                j = perm[i]; if (M->Adata[blockid]->types[j] == MAT_TYPE_ZERO) { continue; }
                if (j == m) {
                    if (!*M->phaseA) { continue; }
                    val = &M->asinvcsinv->x[k];
                } else {
                    val = (j > k) ? (&packIdx(array, m, j, k)) : (&packIdx(array, m, k, j));
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

static DSDP_INT schurM3rowSetup( DSDPSchur *M, DSDP_INT blockid, DSDP_INT row ) {
    // Apply M3 technique to setup a row of the Schur complement matrix
    
    /*
     M3 Technique is like M1 but it computes B_i = <S^-1 * A_i * S^-1> explicitly
     by inverting S at the beginning of each iteration
     
     When each row is set up, not only M is computed, but also
        
     csinvrysinv, asinvrysinv, csinvcsinv, asinvcsinv, csinv, asinv
     
     when available
    */
    DSDP_INT retcode = DSDP_RETCODE_OK;
    DSDP_INT i, j, k, m = M->m, *perm = M->perms[blockid], computeC = FALSE;
    dsMat *B = M->B[blockid]; spsMat *S; rkMat *rkaux;
    double *val, res, Ry = *M->Ry, *aux = M->schurAux, *Sinv = M->Sinv[blockid]; k = perm[row];
    
    if (k < m) {
        val = &M->asinv->x[k];
    } else {
        if (!*M->phaseA) { return retcode; }
        val = M->csinv; computeC = TRUE;
    }
    
    // Get B through direct matrix multiplication
    denseMatReset(B);
    // Compute B_i = <S^-1 * A_i * S^-1>
    switch (M->Adata[blockid]->types[k]) {
        case MAT_TYPE_ZERO: return retcode;
        case MAT_TYPE_SPARSE:
            res = spsSinvSpSinvSolve(Sinv, aux, M->Adata[blockid]->sdpData[k], B); break;
        case MAT_TYPE_DENSE:
            res = spsSinvDsSinvSolve(Sinv, aux, M->Adata[blockid]->sdpData[k], B); break;
        case MAT_TYPE_RANKK:
            S = M->S[blockid]; rkaux = M->rkaux[blockid];
            res = spsSinvRkSinvSolve(S, M->Adata[blockid]->sdpData[k], rkaux);
            rkMatdenseUpdate(B, rkaux); break;
        default: error(etype, "Invalid matrix type. \n"); break;
    }
    
    *val += res;
    
    /* Start M3 */
    if (computeC) {
        assert( *M->phaseA );
        val = M->csinvrysinv; denseDiagTrace(B, Ry, &res); *val += res;
        for (i = row; i < m + 1; ++i) {
            j = perm[i]; val = (j == m) ? M->csinvcsinv : &M->asinvcsinv->x[j];
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
        if (*M->phaseA) { val = &M->asinvrysinv->x[k]; denseDiagTrace(B, Ry, &res); *val += res; }
        double *array = M->M->array;
        for (i = row; i < m + 1; ++i) {
            j = perm[i];
            if (j == m) {
                if (!*M->phaseA) { continue; }
                val = &M->asinvcsinv->x[k];
            } else {
                val = (j > k) ? (&packIdx(array, m, j, k)) : (&packIdx(array, m, k, j));
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

static DSDP_INT schurM4rowSetup( DSDPSchur *M, DSDP_INT blockid, DSDP_INT row ) {
    // Apply M4 technique to setup a row of the Schur complement matrix
    DSDP_INT retcode = DSDP_RETCODE_OK;
    
    
    return retcode;
}

static double schurM5MatAux( DSDP_INT type1, void *A1, DSDP_INT type2, void *A2, double *Sinv ) {
    if (type1 == MAT_TYPE_SPARSE) {
        switch (type2) {
            // A1 is more sparse
            case MAT_TYPE_SPARSE:
                return spsSinvspsSinv(A2, A1, Sinv);
            case MAT_TYPE_RANKK : return spsSinvr1Sinv(A1, ((rkMat *) A2)->data[0], Sinv);
            case MAT_TYPE_ZERO  : return 0.0;
            default: assert( FALSE ); break;
        }
    } else if (type1 == MAT_TYPE_RANKK) {
        switch (type2) {
            case MAT_TYPE_SPARSE: return r1Sinvsps(A2, ((rkMat *) A1)->data[0], Sinv);
            case MAT_TYPE_RANKK : return r1Sinvr1(((rkMat *) A1)->data[0], ((rkMat *) A2)->data[0], Sinv);
            case MAT_TYPE_ZERO  : return 0.0;
            default: assert( FALSE ); break;
        }
    }
    assert( FALSE ); return 0.0;
}

static DSDP_INT schurM5rowSetup( DSDPSchur *M, DSDP_INT blockid, DSDP_INT row ) {
    // Apply M5 technique to setup a row of the Schur complement matrix
    /*
     M5 Technique directly computes elements of the Schur matrix without any intermediate variables
    */
    DSDP_INT retcode = DSDP_RETCODE_OK;
    DSDP_INT i, j, k, m = M->m, *perm = M->perms[blockid], computeC = FALSE;
    DSDP_INT *types = M->Adata[blockid]->types;
    double *Sinv = M->Sinv[blockid];
    double *val, res, tmp, Ry = *M->Ry; k = perm[row];
    void **data = M->Adata[blockid]->sdpData;
    
    /* Start M5 */
    if (k < m) {
        val = &M->asinv->x[k];
    } else {
        if (!*M->phaseA) { return retcode; }
        val = M->csinv; computeC = TRUE;
    }
    
    switch (types[k]) {
        case MAT_TYPE_ZERO: return retcode;
        case MAT_TYPE_SPARSE:
            res = spsRySinv(data[k], Sinv, &tmp, Ry);
            break;
        case MAT_TYPE_RANKK:
            res = r1RySinv(((rkMat *) data[k])->data[0], Sinv, &tmp, Ry);
            break;
        default:
            error(etype, "Invalid matrix type. Dense matrix should not appear in M5. \n");
            break;
    }
    
    *val += tmp;
    
    if (computeC) {
        *M->csinvrysinv += res;
        assert( *M->phaseA );
        for (i = row; i < m + 1; ++i) {
            j = perm[i]; val = (j == m) ? M->csinvcsinv : &M->asinvcsinv->x[j];
            *val += schurM5MatAux(types[k], data[k], types[j], data[j], Sinv);
        }
    } else {
        M->asinvrysinv->x[k] += res;
        double *array = M->M->array;
        for (i = row; i < m + 1; ++i) {
            j = perm[i];
            if (j == m) {
                if (!*M->phaseA) { continue; }
                val = &M->asinvcsinv->x[k];
            } else {
                val = (j > k) ? (&packIdx(array, m, j, k)) : (&packIdx(array, m, k, j));
            }
            *val += schurM5MatAux(types[k], data[k], types[j], data[j], Sinv);
        }
    }
    
    /* End M5 */
    return retcode;
}

static DSDP_INT schurMatSetupBlock( DSDPSchur *M, DSDP_INT blockid ) {
    // Set up a block of the Schur matrix
    DSDP_INT retcode = DSDP_RETCODE_OK;
    
    for (DSDP_INT i = 0; i < M->m + 1; ++i) {
        switch (M->MX[blockid][i]) {
            case SCHUR_M1: schurM1rowSetup(M, blockid, i); break;
            case SCHUR_M2: schurM2rowSetup(M, blockid, i); break;
            case SCHUR_M3: schurM3rowSetup(M, blockid, i); break;
            case SCHUR_M4: schurM4rowSetup(M, blockid, i); break;
            case SCHUR_M5: schurM5rowSetup(M, blockid, i); break;
            default: error(etype, "Invalid technique type. \n"); break;
        }
    }
    return retcode;
}

extern DSDP_INT SchurMatInit( DSDPSchur *M ) {
    
    DSDP_INT retcode = DSDP_RETCODE_OK;
    
    M->m = 0; M->nblock = 0; M->Mready = FALSE; M->useTwo = NULL;
    M->perms = NULL; M->MX = NULL; M->S = NULL; M->B = NULL;
    M->Adata = NULL; M->M = NULL; M->asinv = NULL; M->Ry = NULL;
    M->asinvrysinv = NULL; M->asinvcsinv = NULL; M->csinvrysinv = NULL;
    M->csinv = NULL; M->csinvcsinv = NULL; M->Sinv = NULL; M->rkaux = NULL;
    M->scaler = 0.0;
    
    return retcode;
}

extern DSDP_INT SchurMatSetDim( DSDPSchur *M, DSDP_INT m, DSDP_INT nblock ) {
    assert( m > 0 && nblock > 0);
    M->m = m; M->nblock = nblock;
    return DSDP_RETCODE_OK;
}

extern DSDP_INT SchurMatAlloc( DSDPSchur *M ) {
    // Allocate memory for internal perms
    DSDP_INT retcode = DSDP_RETCODE_OK;
    
    M->perms  = (DSDP_INT **) calloc(M->nblock, sizeof(DSDP_INT *));
    M->MX     = (DSDP_INT **) calloc(M->nblock, sizeof(DSDP_INT *));
    M->useTwo = (DSDP_INT  *) calloc(M->nblock, sizeof(DSDP_INT));
    M->Sinv   = (double   **) calloc(M->nblock, sizeof(double *));
    M->scaler = 0.0;
    
    for (DSDP_INT i = 0; i < M->nblock; ++i) {
        M->perms[i] = (DSDP_INT *) calloc(M->m + 1, sizeof(DSDP_INT));
        M->MX[i] = (DSDP_INT *) calloc(M->m + 1, sizeof(DSDP_INT));
        M->Sinv[i] = NULL;
    }

    return retcode;
}

extern DSDP_INT SchurMatRegister( DSDPSchur *M, spsMat **S, dsMat **B, sdpMat **Adata, dsMat *Msdp,
                                  vec *asinv, vec *asinvrysinv, vec *asinvcsinv, double *csinvrysinv,
                                  double *csinv, double *csinvcsinv, double *Ry, rkMat **rkaux, DSDP_INT *phaseA ) {
    // Register information into M
    DSDP_INT retcode = DSDP_RETCODE_OK;
    assert( S && B && Adata && Msdp && !M->Mready);
    M->S = S; M->B = B; M->Adata = Adata; M->M = Msdp;
    
    assert( asinv && asinvrysinv && asinvcsinv && csinvrysinv );
    M->asinv = asinv; M->asinvrysinv = asinvrysinv; M->Ry = Ry;
    M->asinvcsinv = asinvcsinv; M->csinvrysinv = csinvrysinv;
    M->csinv = csinv; M->csinvcsinv = csinvcsinv; M->rkaux = rkaux;
    M->phaseA = phaseA;
    
    return retcode;
}

extern DSDP_INT SchurMatFree( DSDPSchur *M ) {
    
    DSDP_INT retcode = DSDP_RETCODE_OK;
    
    for (DSDP_INT i = 0; i < M->nblock; ++i) {
        DSDP_FREE(M->perms[i]); DSDP_FREE(M->MX[i]); DSDP_FREE(M->Sinv[i]);
    }
    
    DSDP_FREE(M->perms); DSDP_FREE(M->MX); DSDP_FREE(M->Sinv); DSDP_FREE(M->schurAux);
    M->Mready = FALSE; M->m = 0; M->nblock = 0; M->perms = NULL;
    M->MX= NULL; M->S = NULL; M->B = NULL; M->Adata = NULL;
    M->M = NULL; M->asinv = NULL; M->asinvrysinv = NULL;
    M->asinvcsinv = NULL; M->csinvrysinv = NULL; M->useTwo = FALSE;
    M->csinv = NULL; M->csinvcsinv = NULL; M->scaler = 0.0;
    
    return retcode;
}

// Schur matrix re-ordering heuristic
extern DSDP_INT DSDPSchurReorder( DSDPSchur *M ) {
    // Invoke the re-ordering heuristic for Schur matrix setup
    DSDP_INT retcode = DSDP_RETCODE_OK;
    
    DSDP_INT maxblock = 0;
    // Information of all the matrices must have been supplied
    for (DSDP_INT i = 0, dim; i < M->nblock; ++i) {
        retcode = schurBlockAnalysis(M->Adata[i], M->perms[i],
                                     M->MX[i], &M->useTwo[i]);
        if (!M->useTwo[i]) {
            dim = M->Adata[i]->dimS; maxblock = MAX(dim, maxblock);
            M->Sinv[i] = (double *) calloc(dim * dim, sizeof(double));
        }
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

extern DSDP_INT DSDPSchurSetup( DSDPSchur *M ) {
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
    
    assert( M->Mready );
    
//#define superDebug
#ifdef superDebug
    DSDP_INT m = M->m, mpack = nsym(m);
    
    double *MM1 = (double *) calloc(mpack, sizeof(double));
    double *MM2 = (double *) calloc(mpack, sizeof(double));
    double *MM3 = (double *) calloc(mpack, sizeof(double));
    double *MM5 = (double *) calloc(mpack, sizeof(double));
    double *Mref = M->M->array;
    
    schurMatCleanup(M); schurMatGetSinv(M);
    for (DSDP_INT i = 0; i < M->nblock; ++i) {
        for (DSDP_INT k = 0; k < M->m + 1; ++k) {
            M->MX[i][k] = SCHUR_M1;
        }
    }
    for (DSDP_INT i = 0; i < M->nblock; ++i) {
        schurMatSetupBlock(M, i);
    }
    memcpy(MM1, Mref, sizeof(double) * mpack);
    
    schurMatCleanup(M); schurMatGetSinv(M);
    for (DSDP_INT i = 0; i < M->nblock; ++i) {
        for (DSDP_INT k = 0; k < M->m + 1; ++k) {
            M->MX[i][k] = SCHUR_M2;
        }
    }
    for (DSDP_INT i = 0; i < M->nblock; ++i) {
        schurMatSetupBlock(M, i);
    }
    
    memcpy(MM2, Mref, sizeof(double) * mpack);
    schurMatCleanup(M); schurMatGetSinv(M);
    for (DSDP_INT i = 0; i < M->nblock; ++i) {
        for (DSDP_INT k = 0; k < M->m + 1; ++k) {
            M->MX[i][k] = SCHUR_M3;
        }
    }
    for (DSDP_INT i = 0; i < M->nblock; ++i) {
        schurMatSetupBlock(M, i);
    }
    
    memcpy(MM3, Mref, sizeof(double) * mpack);
    
    schurMatCleanup(M); schurMatGetSinv(M);
    for (DSDP_INT i = 0; i < M->nblock; ++i) {
        for (DSDP_INT k = 0; k < M->m + 1; ++k) {
            M->MX[i][k] = SCHUR_M5;
        }
    }
    for (DSDP_INT i = 0; i < M->nblock; ++i) {
        schurMatSetupBlock(M, i);
    }
    
    memcpy(MM5, Mref, sizeof(double) * mpack);
    double diff12 = 0.0, diff13 = 0.0, diff23 = 0.0, diff25 = 0.0, diff35 = 0.0, tmp;
    
    for (DSDP_INT i = 0; i < mpack; ++i) {
        tmp = MM1[i] - MM2[i];
        diff12 += tmp * tmp;
        tmp = MM1[i] - MM3[i];
        diff13 += tmp * tmp;
        tmp = MM2[i] - MM3[i];
        diff23 += tmp * tmp;
        tmp = MM2[i] - MM5[i];
        diff25 += tmp * tmp;
        tmp = MM3[i] - MM5[i];
        diff35 += tmp * tmp;
        
        if (fabs(tmp) * M->scaler * M->scaler > 1) {
             assert( FALSE );
        }
        
    }
    
    if (MAX(MAX(diff12, diff23), diff13) > 0.1) {
        // assert(FALSE);
    }
    printf("| Difference: 12 %e  13 %e  23 %e  25 %e  35 %e \n", diff12, diff13, diff23, diff25, diff35);
    DSDP_FREE(MM1); DSDP_FREE(MM2); DSDP_FREE(MM3); DSDP_FREE(MM5);
    
#else
    // Clean up current values and invert blocks
    schurMatCleanup(M);
    schurMatGetSinv(M);
    
    for (DSDP_INT i = 0; i < M->nblock; ++i) {
        schurMatSetupBlock(M, i);
    }
    
#endif
    
    // Scale asinv and
    vec_rscale(M->asinv, M->scaler);
    *M->csinv /= M->scaler;
    
    return retcode;
}
