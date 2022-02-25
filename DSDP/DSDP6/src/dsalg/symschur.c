#include "symschur.h"
#include "dsdputils.h"
#include "dsdpSort.h"

// Implement the advanced Schur matrix setup in DSDP using M1, M2, M3, M4 techniques
static char etype[] = "Advanced Schur matrix setup";

static DSDP_INT SchurMatInit( DSDPSchur *M ) {
    
    DSDP_INT retcode = DSDP_RETCODE_OK;
    
    M->m = 0; M->nblock = 0; M->Mready = FALSE; M->useTwo = NULL;
    M->perms = NULL; M->MX = NULL; M->S = NULL; M->B = NULL;
    M->Adata = NULL; M->M = NULL; M->asinv = NULL; M->Ry = NULL;
    M->asinvrysinv = NULL; M->asinvcsinv = NULL; M->csinvrysinv = NULL;
    M->csinv = NULL; M->csinvcsinv = NULL; M->Sinv = NULL; M->rkaux = NULL;
    
    return retcode;
}

static DSDP_INT SchurMatSetDim( DSDPSchur *M, DSDP_INT m, DSDP_INT nblock ) {
    assert( m > 0 && nblock > 0);
    M->m = m; M->nblock = nblock;
    return DSDP_RETCODE_OK;
}

static DSDP_INT SchurMatAlloc( DSDPSchur *M ) {
    // Allocate memory for internal perms
    DSDP_INT retcode = DSDP_RETCODE_OK;
    
    M->perms = (DSDP_INT **) calloc(M->nblock, sizeof(DSDP_INT));
    M->MX = (DSDP_INT **) calloc(M->nblock, sizeof(DSDP_INT));
    M->useTwo = (DSDP_INT *) calloc(M->nblock, sizeof(DSDP_INT));
    M->Sinv = (double **) calloc(M->nblock, sizeof(double));
    
    for (DSDP_INT i = 0; i < M->nblock; ++i) {
        M->perms[i] = (DSDP_INT *) calloc(M->m, sizeof(DSDP_INT));
        M->MX[i] = (DSDP_INT *) calloc(M->m, sizeof(DSDP_INT));
        M->Sinv = NULL;
    }

    return retcode;
}

static DSDP_INT SchurMatRegister( DSDPSchur *M, spsMat **S, dsMat **B, sdpMat **Adata, dsMat *Msdp,
                                  vec *asinv, vec *asinvrysinv, vec *asinvcsinv, double *csinvrysinv,
                                  double *csinv, double *csinvcsinv, double *Ry, rkMat **rkaux ) {
    // Register information into M
    DSDP_INT retcode = DSDP_RETCODE_OK;
    assert( S && B && Adata && Msdp && !M->Mready);
    M->S = S; M->B = B; M->Adata = Adata; M->M = Msdp;
    
    assert( asinv && asinvrysinv && asinvcsinv && csinvrysinv );
    M->asinv = asinv; M->asinvrysinv = asinvrysinv; M->Ry = Ry;
    M->asinvcsinv = asinvcsinv; M->csinvrysinv = csinvrysinv;
    M->csinv = csinv; M->csinvcsinv = csinvcsinv; M->rkaux = rkaux;
    
    return retcode;
}

static DSDP_INT SchurMatFree( DSDPSchur *M ) {
    
    DSDP_INT retcode = DSDP_RETCODE_OK;
    
    for (DSDP_INT i = 0; i < M->nblock; ++i) {
        DSDP_FREE(M->perms[i]); DSDP_FREE(M->MX[i]); DSDP_FREE(M->Sinv[i]);
    }
    
    DSDP_FREE(M->perms); DSDP_FREE(M->MX); DSDP_FREE(M->Sinv);
    M->Mready = FALSE; M->m = 0; M->nblock = 0; M->perms = NULL;
    M->MX= NULL; M->S = NULL; M->B = NULL; M->Adata = NULL;
    M->M = NULL; M->asinv = NULL; M->asinvrysinv = NULL;
    M->asinvcsinv = NULL; M->csinvrysinv = NULL; M->useTwo = FALSE;
    M->csinv = NULL; M->csinvcsinv = NULL;
    
    return retcode;
}

static DSDP_INT getScore( DSDP_INT *ranks, DSDP_INT *nnzs, DSDP_INT *perm,
                        DSDP_INT m, DSDP_INT n, DSDP_INT idx, DSDP_INT *MX,
                        double *M1M2, double *M1M5 ) {
    // Compute the estimated number of multiplications for technique M1 to M5
    
    DSDP_INT i, j, k = perm[idx], nsqr = n * n, ncbe = n * n * n;
    double d1, d2, d3, d4, d5, sumnnz = 0.0, rsigma = ranks[k], best = SCHUR_M1;
    double m1m5 = 0.0, m1m2 = 0.0;
    // On exit, MX is set 10 * MX12 + M12345

    for (i = idx; i < m + 1; ++i) {
        sumnnz += nnzs[perm[i]];
    }
  
    d1 = rsigma * (1 + 2 * KAPPA) * nsqr + KAPPA * sumnnz;
    d2 = rsigma * (nsqr + KAPPA * sumnnz);
    d3 = n * KAPPA * nnzs[k] + ncbe + KAPPA * sumnnz;
    d4 = n * KAPPA * nnzs[k] + KAPPA * (n + 1) * sumnnz;
    d5 = KAPPA * (2 * KAPPA * nnzs[k] + 1) * sumnnz;
    
    // Get better of d1 and d2
    j = 0;
    if (d1 < d2) {
        best = SCHUR_M1;
        m1m2 = d1; m1m5 = d1;
    } else {
        best = SCHUR_M2;
        m1m2 = d2; m1m5 = d2;
    }
    
    j += best * 10;
    
    if (d3 < m1m5) { best = SCHUR_M3; m1m5 = d3; }
    if (d4 < m1m5) { best = SCHUR_M4; m1m5 = d4; }
    if (d5 < m1m5) { best = SCHUR_M5; m1m5 = d5; }
    
    *M1M2 = m1m2; *M1M5 = m1m5;
    
    j += best;

    return DSDP_RETCODE_OK;
}

static DSDP_INT schurBlockAnalysis( sdpMat *Adata, DSDP_INT *permk, DSDP_INT *MXk, DSDP_INT *useM1M2 ) {
    // Reorder block
    DSDP_INT retcode = DSDP_RETCODE_OK;
    
    DSDP_INT *dsidx = Adata->denseMatIdx, ndsMat = Adata->ndenseMat,
             *spsidx = Adata->spsMatIdx, nspsMat = Adata->nspsMat,
             *rkidx = Adata->rkMatIdx, nrkMat = Adata->nrkMat,
             nzeroMat = Adata->nzeroMat, m = Adata->dimy, n = Adata->dimy;
    
    assert( ndsMat + nspsMat + nrkMat + nzeroMat == m + 1 );
    
    memset(permk, 0, sizeof(DSDP_INT) * m);
    memset(MXk, 0, sizeof(DSDP_INT) * m);
    
    double *scoreM1M2 = (double *) calloc(m, sizeof(double)), ksiM1M2 = 0.0;
    double *scoreM1M5 = (double *) calloc(m, sizeof(double)), ksiM1M5 = 0.0;
    DSDP_INT *ranks   = (DSDP_INT *) calloc(m, sizeof(DSDP_INT));
    DSDP_INT *nnzs = (DSDP_INT *) calloc(m, sizeof(DSDP_INT));
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
        spsMatGetRank(Adata->sdpData[k], &ranks[k]);
    }
    for (i = 0; i < nrkMat; ++i) {
        k = rkidx[i];
        nnzs[k] = ((r1Mat *) (Adata->sdpData[i]))->nnz;
        nnzs[k] *= MXk[k];
        ranks[k] = 1;
    }
    
    // Sort f_1, ..., f_m in descending order
    dsdpSort2(permk, nnzs, 0, m);
    
    // Determine strategies for setting up the Schur matrix
    for (i = 0; i < m + 1; ++i) {
        getScore(ranks, nnzs, permk, m + 1, n, i,
                 &MXk[i], &scoreM1M2[i], &scoreM1M5[i]);
    }
    
    // Determine whether to use M3 to M5
    for (i = 0; i < m + 1; ++i) {
        ksiM1M2 += scoreM1M2[i];
        ksiM1M5 += scoreM1M5[i];
    }
    
    if (ksiM1M2 <= ksiM1M5 + n * n * n) {
        M1M2 = TRUE;
    } else {
        M1M2 = FALSE;
    }
    
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
    
    DSDP_FREE(scoreM1M2); DSDP_FREE(scoreM1M5);
    DSDP_FREE(ranks); DSDP_FREE(nnzs);
    
    return retcode;
}

static DSDP_INT schurM1rowSetup( DSDPSchur *M, DSDP_INT blockid, DSDP_INT row ) {
    // Apply M1 technique to setup a row of the Schur complement matrix
    /*
     M1 technique first applies rank-1 update to get B_i = <S^-1 * A_i * S^-1>
     and then computes M_{ij} = <B_i, A_j>
    */
    
    DSDP_INT retcode = DSDP_RETCODE_OK;
    
    DSDP_INT i, j, k, m = M->m, *perm = M->perms[blockid];
    dsMat *B = M->B[blockid]; rkMat *factor, *rkaux = M->rkaux[blockid];
    k = perm[row];
    
    switch (M->Adata[blockid]->types[k]) {
        case MAT_TYPE_ZERO: return retcode;
        case MAT_TYPE_SPARSE: factor = spsMatGetFactor(M->Adata[blockid]->sdpData[k]);
            break;
        case MAT_TYPE_DENSE: factor = denseMatGetFactor(M->Adata[blockid]->sdpData[k]);
            break;
        case MAT_TYPE_RANKK: factor = M->Adata[blockid]->sdpData[k];
            break;
        default: error(etype, "Invalid matrix type. \n");
            break;
    }
    
    if (k < m) {
        spsSinvRkSinvSolve(M->S[blockid], factor, rkaux, &M->asinv->x[k]);
    } else {
        spsSinvRkSinvSolve(M->S[blockid], factor, rkaux, M->csinv);
    }
    
    // Get B through rank-one update
    denseMatReset(B); rkMatdenseUpdate(B, rkaux);
    
    // Compute B_i = <S^-1 * A_i * S^-1>
    for (i = row; i < m + 1; ++i) {
        j = perm[i];
        
        
        
    }

    return retcode;
}

static DSDP_INT schurM2rowSetup( DSDPSchur *M, DSDP_INT row ) {
    // Apply M2 technique to setup a row of the Schur complement matrix
    DSDP_INT retcode = DSDP_RETCODE_OK;
    
    
    
    return retcode;
}

static DSDP_INT schurM3rowSetup( DSDPSchur *M, DSDP_INT row ) {
    // Apply M3 technique to setup a row of the Schur complement matrix
    DSDP_INT retcode = DSDP_RETCODE_OK;
    
    
    
    return retcode;
}

static DSDP_INT schurM4rowSetup( DSDPSchur *M, DSDP_INT row ) {
    // Apply M4 technique to setup a row of the Schur complement matrix
    DSDP_INT retcode = DSDP_RETCODE_OK;
    
    
    
    return retcode;
}

static DSDP_INT schurM5rowSetup( DSDPSchur *M, DSDP_INT row ) {
    // Apply M5 technique to setup a row of the Schur complement matrix
    DSDP_INT retcode = DSDP_RETCODE_OK;
    
    
    
    return retcode;
}

static DSDP_INT schurMatSetupBlock( DSDPSchur *M, DSDP_INT blockid ) {
    // Set up a block of the Schur matrix
    DSDP_INT retcode = DSDP_RETCODE_OK;
    
    
    return retcode;
}

extern DSDP_INT DSDPPrepareMAssembler( HSDSolver *dsdpSolver ) {
    // Initialize the internal Schur matrix structure
    DSDP_INT retcode = DSDP_RETCODE_OK;
    
    DSDPSchur *M = dsdpSolver->M;
    retcode = SchurMatInit(M);
    retcode = SchurMatSetDim(M, dsdpSolver->m, dsdpSolver->nBlock);
    retcode = SchurMatAlloc(M);
    retcode = SchurMatRegister(M, dsdpSolver->S,
                               dsdpSolver->dsaux,
                               dsdpSolver->sdpData,
                               dsdpSolver->Msdp,
                               dsdpSolver->asinv,
                               dsdpSolver->d4,
                               dsdpSolver->u,
                               &dsdpSolver->csinvrysinv,
                               &dsdpSolver->csinv,
                               &dsdpSolver->csinvcsinv,
                               &dsdpSolver->Ry,
                               dsdpSolver->rkaux);
    return retcode;
}

// Schur matrix re-ordering heuristic
extern DSDP_INT DSDPSchurReorder( DSDPSchur *M ) {
    // Invoke the re-ordering heuristic for Schur matrix setup
    DSDP_INT retcode = DSDP_RETCODE_OK;
    
    // Information of all the matrices must have been supplied
    for (DSDP_INT i = 0, dim; i < M->nblock; ++i) {
        retcode = schurBlockAnalysis(M->Adata[i], M->perms[i],
                                     M->MX[i], &M->useTwo[i]);
        
        if (!M->useTwo[i]) {
            dim = M->Adata[i]->dimS;
            M->Sinv[i] = (double *) calloc(dim * dim, sizeof(double));
        }
    }
    
    M->Mready = TRUE;
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
    
    // Clean up current values
    denseMatReset(M->M);
    vec_reset(M->asinv); vec_reset(M->asinvcsinv); vec_reset(M->asinvrysinv);
    *M->csinvrysinv = 0.0; *M->csinv = 0.0; *M->csinvcsinv = 0.0;
    
    for (DSDP_INT i = 0; i < M->nblock; ++i) {
        schurMatSetupBlock(M, i);
    }
    
    return retcode;
}
