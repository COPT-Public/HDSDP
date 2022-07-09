#include "sparsemat.h"
#include "dsdplanczos.h"

static char etype[] = "Lanczos iteration";

#define LWORK   30  // >= MAXITER + 6
#define IWORK   12  // >= 10

static DSDP_INT lanczosAlloc( vec **v,  vec **w,
                              vec **z1, vec **z2,
                              vec **vecaux,
                              double **V, double **H,
                              double **Y, double **d,
                              double **mataux, double **eigaux,
                              DSDP_INT **eigintaux, DSDP_INT n,
                              DSDP_INT maxiter ) {
    
    DSDP_INT retcode = DSDP_RETCODE_OK;
    vec *vecdata = NULL;
    vecdata = (vec *) calloc(1, sizeof(vec));
    vec_init(vecdata); retcode = vec_alloc(vecdata, n);
    checkCode; *v = vecdata;
    vecdata = (vec *) calloc(1, sizeof(vec));
    vec_init(vecdata); retcode = vec_alloc(vecdata, n);
    checkCode; *w = vecdata;
    vecdata = (vec *) calloc(1, sizeof(vec));
    vec_init(vecdata); retcode = vec_alloc(vecdata, n);
    checkCode; *z1 = vecdata;
    vecdata = (vec *) calloc(1, sizeof(vec));
    vec_init(vecdata); retcode = vec_alloc(vecdata, n);
    checkCode; *z2 = vecdata;
    vecdata = (vec *) calloc(1, sizeof(vec));
    vec_init(vecdata); retcode = vec_alloc(vecdata, n);
    checkCode; *vecaux = vecdata;
    
    *V = (double *) calloc(n * (maxiter + 1), sizeof(double));
    *H = (double *) calloc((maxiter + 1) * (maxiter + 1), sizeof(double));
    *Y = (double *) calloc(2 * maxiter, sizeof(double));
    *d = (double *) calloc(maxiter, sizeof(double));
    *mataux = (double *) calloc(maxiter * maxiter, sizeof(double));
    *eigaux = (double *) calloc(maxiter * LWORK, sizeof(double));
    *eigintaux = (DSDP_INT *) calloc(maxiter * IWORK, sizeof(DSDP_INT));
    
    if (!V || !H || !Y || !d || !mataux || !eigaux || !eigintaux) {
        printf("| Failed to allocate memory for Lanczos iterator. \n");
        retcode = DSDP_RETCODE_FAILED; return retcode;
    }
    
    return retcode;
}

static DSDP_INT lanczosFree( vec **v,  vec **w,
                             vec **z1, vec **z2,
                             vec **vecaux,
                             double **V, double **H,
                             double **Y, double **d,
                             double **mataux,
                             double **eigaux,
                             DSDP_INT **eigintaux ) {
    
    vec_free(*v);      DSDP_FREE(*v);
    vec_free(*w);      DSDP_FREE(*w);
    vec_free(*z1);     DSDP_FREE(*z1);
    vec_free(*z2);     DSDP_FREE(*z2);
    vec_free(*vecaux); DSDP_FREE(*vecaux);
    
    DSDP_FREE(*V); DSDP_FREE(*H); DSDP_FREE(*Y);
    DSDP_FREE(*d); DSDP_FREE(*mataux); DSDP_FREE(*eigaux);
    DSDP_FREE(*eigintaux);
    
    return DSDP_RETCODE_OK;
}

static DSDP_INT lanczositerInitalize( vec *v, double *V ) {
    // v = v / norm(v); V(:, 1) = v;
    DSDP_INT retcode = DSDP_RETCODE_OK;
    srand(v->dim); double nrm = 0.0;
    for (DSDP_INT i = 0; i < v->dim; ++i) {
        srand(rand());
        v->x[i] = sqrt(sqrt((rand() % 1627)));
    }
    vec_norm(v, &nrm); vec_rscale(v, nrm);
    // V(:, 1) = v;
    memcpy(V, v->x, sizeof(double) * v->dim);
    return retcode;
}

extern DSDP_INT dsdpLanczosInit( DSDPLanczos *lczSolver ) {
    
    lczSolver->n = 0; lczSolver->v = NULL; lczSolver->w = NULL;
    lczSolver->z1 = NULL; lczSolver->z2 = NULL; lczSolver->vecaux = NULL;
    lczSolver->V = NULL; lczSolver->H = NULL; lczSolver->Y = NULL;
    lczSolver->d = NULL; lczSolver->mataux = NULL; lczSolver->eigaux = NULL;
    lczSolver->eigintaux = NULL; lczSolver->lwork = LWORK * MAXITER;
    lczSolver->iwork = IWORK * MAXITER;
    return DSDP_RETCODE_OK;
}

extern DSDP_INT dsdpLanczosAlloc( DSDPLanczos *lczSolver, DSDP_INT n ) {
    
    DSDPLanczos *l = lczSolver;
    lczSolver->n = n;
    return lanczosAlloc(&l->v, &l->w, &l->z1, &l->z2, &l->vecaux,
                        &l->V, &l->H, &l->Y, &l->d, &l->mataux,
                        &l->eigaux, &l->eigintaux, l->n, MAXITER);
}

extern void dsdpLanczosFree( DSDPLanczos *lczSolver ) {
    
    DSDPLanczos *l = lczSolver;
    lanczosFree(&l->v, &l->w, &l->z1, &l->z2, &l->vecaux,
                &l->V, &l->H, &l->Y, &l->d, &l->mataux,
                &l->eigaux, &l->eigintaux);
    l->n = 0; l->lwork = 0; l->iwork = 0;
}

extern DSDP_INT dsdpLanczosStep( DSDPLanczos *lczSolver, spsMat *S, spsMat *dS, double *lbd, double *delta ) {
    // Lanczos algorithm that computes lambda_max(L^-1 dS L^-T)
    DSDP_INT retcode = DSDP_RETCODE_OK;    
    DSDPLanczos *l = lczSolver;
    double fnrm = 0.0; spsMatFnorm(dS, &fnrm);
    if (fnrm < 1e-13) { *lbd = 0.0; *delta = 0.0; return retcode;}
    
    /* Prepare working array */
    DSDP_INT n = S->dim, maxiter = MAXITER, one = 1, mH = maxiter + 1,
             il, iu, lwork = lczSolver->lwork, iwork = lczSolver->iwork, failed = FALSE, \
            idx, neigs, info, *eigintaux = lczSolver->eigintaux, *isuppz = lczSolver->isuppz;
    
    vec *v = l->v, *w = l->w, *z1 = l->z1, *z2 = l->z2, *vecaux = l->vecaux;
    
    double *V = l->V, *H = l->H, *Y = l->Y, *d = l->d,
           *mataux = l->mataux, *eigaux = l->eigaux,
            nrm, nrm2, alp, alpha, beta = 0.0, res;
    char notrans = DSDP_MAT_NOTRANSPOSE, jobz = 'V',
         range = 'I', uplo = DSDP_MAT_UP;
    
    /* Initialize */
    retcode = lanczositerInitalize(v, V);
    /* Start Iterating */
    for (DSDP_INT k = 0; k < maxiter; ++k) {
        spsMatVecBSolve(S, v, w); spsMatAx(dS, w, vecaux); // w = dS * (LT \ v);
        spsMatVecFSolve(S, vecaux, w); vec_norm(w, &nrm); // w = L \ w;
        
        if (k > 0) {
            idx = (k - 1) * n;
            memcpy(vecaux->x, &V[idx], sizeof(double) * n); // w = w - H(k, k - 1) * V(:, k - 1);
            vec_axpy(- H[(k - 1) * mH + k], vecaux, w);
        }
        
        idx = k * n;
        alp = ddot(&n, w->x, &one, &V[idx], &one); // alp = w' * V(:, k);
        memcpy(vecaux->x, &V[idx], sizeof(double) * n);
        vec_axpy(- alp, vecaux, w); // w = w - alp * V(:, k);
        vec_norm(w, &nrm2);
        H[k * mH + k] = alp; // H(k, k) = alp;
        if (nrm2 < 0.8 * nrm && (FALSE)) {
            alpha = 1.0; beta  = 0.0; // s = (w' * V(:, 1:k))';
            idx = k + 1; memset(mataux, 0, sizeof(double) * idx);
            for (DSDP_INT p = 0, q; p < idx; ++p) {
                for (q = 0; q < n; ++q) {
                    mataux[p] += V[p * n + q] * w->x[q];
                }
            }
            
            dgemv(&notrans, &n, &idx, &alpha, V, &n, mataux,
                  &one, &alpha, w->x, &one); // w = w - V(:, 1:k) * s;
            daxpy(&idx, &alpha, mataux, &one, &H[mH * k], &one); // H(1:k, k) = H(1:k, k) + s;
            vec_norm(w, &nrm2);
        }
        
        if (nrm2 > 0.0) {
            vec_rscale(w, nrm2);
            vec_copy(w, v);
        
            // V(:, k + 1) = v;
            memcpy(&V[n * k + n], v->x, sizeof(double) * n);
            H[mH *  k      + (k + 1)] = nrm2;
            H[mH * (k + 1) +  k     ] = nrm2;
        }
        
        // Check convergence
        if ( ((k + 1) % 5 == 0 || k >= maxiter - 1 ) || nrm2 == 0 ) {
            DSDP_INT kp1 = k + 1;
            // Hk = H(1:k, 1:k);
            for (DSDP_INT i = 0; i < kp1; ++i) {
                memcpy(&mataux[kp1 * i], &H[mH * i], sizeof(double) * kp1);
            }
            alpha = 0; il = k; iu = kp1;
            dsyevr(&jobz, &range, &uplo, &kp1, mataux, &kp1,
                   NULL, NULL, &il, &iu, &alpha, &neigs, d,
                   Y, &kp1, isuppz, eigaux, &lwork, eigintaux,
                   &iwork, &info);
            
            if (info) { failed = TRUE; }
            
            res = fabs(H[k * mH + k + 1] * Y[kp1 + k]);
            if (res <= 1e-04 || k >= maxiter - 1 || nrm == 0.0) {
                //lambda = eigH(idx(k)); lambda2 = eigH(idx(k - 1));
                double lambda1 = d[1], lambda2 = d[0], res1, res2, tmp, gamma, alpha = 1.0;
                // z = V(:, 1:k) * Y(:, idx(k));
                dgemv(&notrans, &n, &kp1, &alpha, V, &n,
                      &Y[kp1], &one, &beta, z1->x, &one);
                // tmp  = dS * (LT \ z);
                spsMatVecBSolve(S, z1, z2); spsMatAx(dS, z2, vecaux);
                spsMatVecFSolve(S, vecaux, z2);
                // res1 = norm(L \ tmp - lambda * z);
                vec_axpby(1.0, z2, -lambda1, z1);
                vec_norm(z1, &res1);
                // z2 = V(:, 1:k) * Y(:, idx(k - 1));
                dgemv(&notrans, &n, &kp1, &alpha, V, &n,
                      Y, &one, &beta, z2->x, &one);
                // tmp  = dS * (LT \ z2);
                spsMatVecBSolve(S, z2, z1);
                spsMatAx(dS, z1, vecaux);
                spsMatVecFSolve(S, vecaux, z1);
                // res2 = norm(L \ tmp - lambda * z);
                vec_axpby(1.0, z1, -lambda1, z2);
                vec_norm(z2, &res2);
                tmp = lambda1 - lambda2 - res2;
                gamma = (tmp > 0) ? tmp : 1e-15;
                tmp = res1 * res1 / gamma;
                gamma = MIN(res1, tmp);
                
                if (gamma < 1e-03 || gamma + lambda1 <= 0.8 || failed) {
                    if (failed) { lambda1 = MAX(1000.0, lambda1); }
                    *delta = gamma; *lbd = lambda1; break;
                } else {
                    if (nrm2 == 0.0) {
                        error(etype, "Bad lanczos initialization \n");
                    }
                    *delta = gamma; *lbd = lambda1;
                }
            }
        }
    }
    
    return retcode;
}
