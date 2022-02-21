#include <assert.h>
#include "dsdplanczos.h"
#include "sparsemat.h"

static char etype[] = "Lanczos iteration";

#define MAXITER 10
#define LWORK   30  // >= MAXITER + 6
#define IWORK   12  // >= 10

static double firstNnzSign( double *Yk, DSDP_INT k ) {
    // Return the sign of the first non-zero element
    for (DSDP_INT i = 0; i < k; ++i) {
        if (Yk[i] > 0) {
            return (-1.0);
        } else if (Yk[i] < 0) {
            return ( 1.0);
        }
    }
    // We should never be here
    assert( 0 );
    return 0.0;
}

// TODO: Fix potential bugs in Lanczos
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
    vec_init(vecdata);
    vec_alloc(vecdata, n);
    *v = vecdata;
    
    vecdata = (vec *) calloc(1, sizeof(vec));
    vec_init(vecdata);
    vec_alloc(vecdata, n);
    *w = vecdata;
    
    vecdata = (vec *) calloc(1, sizeof(vec));
    vec_init(vecdata);
    vec_alloc(vecdata, n);
    *z1 = vecdata;
    
    vecdata = (vec *) calloc(1, sizeof(vec));
    vec_init(vecdata);
    vec_alloc(vecdata, n);
    *z2 = vecdata;
    
    vecdata = (vec *) calloc(1, sizeof(vec));
    vec_init(vecdata);
    vec_alloc(vecdata, n);
    *vecaux = vecdata;
    
    *V = (double *) calloc(n * (maxiter + 1), sizeof(double));
    *H = (double *) calloc((maxiter + 1) * (maxiter + 1), sizeof(double));
    *Y = (double *) calloc(2 * maxiter, sizeof(double));
    *d = (double *) calloc(maxiter, sizeof(double));
    *mataux = (double *) calloc(maxiter * maxiter, sizeof(double));
    *eigaux = (double *) calloc(maxiter * LWORK, sizeof(double));
    *eigintaux = (DSDP_INT *) calloc(maxiter * IWORK, sizeof(DSDP_INT));
    
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

static DSDP_INT lanczosInitialize( vec *v, double *V ) {
    /*
     v = v / norm(v);
     V(:, 1) = v;
    */
    DSDP_INT retcode = DSDP_RETCODE_OK;
    srand(123);
    double nrm = 0.0;
    
    /*
     for (DSDP_INT i = 0; i < v->dim; ++i) {
         v->x[i] = rand();
         srand(v->x[i]);
         v->x[i] /= RAND_MAX;
     }
    */
    
    for (DSDP_INT i = 0; i < v->dim; ++i) {
        v->x[i] = 1;
    }
    
    vec_norm(v, &nrm);
    vec_rscale(v, nrm);
    
    // V(:, 1) = v;
    memcpy(V, v->x, sizeof(double) * v->dim);
    return retcode;
}

extern DSDP_INT dsdpLanczos( spsMat *S, spsMat *dS, double *lbd, double *delta ) {
    
    /* Lanczos algorithm that computes lambda_max(L^-1 dS L^-T) */
    DSDP_INT retcode = DSDP_RETCODE_OK;
    assert( S->dim == dS->dim && S->isFactorized );
    
    // spsMatLinvView(S);
    
    double fnrm = 0.0;
    spsMatFnorm(dS, &fnrm);
    
    if (fnrm < 1e-15) {
        *lbd = 0.0;
        *delta = 0.0;
    }
    
    /* Prepare working array */
    DSDP_INT n = S->dim, maxiter = MAXITER, one = 1, mH = maxiter + 1,
            il, iu, lwork = LWORK * maxiter, iwork = IWORK * maxiter,
            idx, neigs, info, *eigintaux, *isuppz;
    
    // 2 * 2
    isuppz = (DSDP_INT *) calloc(2 * 2, sizeof(DSDP_INT));
    
    vec *v = NULL, *w = NULL, *z1 = NULL, *z2 = NULL, *vecaux = NULL;
    
    double *V = NULL, *H      = NULL, *Y = NULL,
           *d = NULL, *mataux = NULL, *eigaux = NULL,
            nrm, nrm2, alp, alpha, beta, res, seig1, seig2;
    
    char trans = DSDP_MAT_TRANSPOSE, notrans = DSDP_MAT_NOTRANSPOSE,
         jobz = 'V', range = 'I', uplo = DSDP_MAT_UP;
    
    retcode = lanczosAlloc(&v, &w, &z1, &z2, &vecaux,
                           &V, &H, &Y, &d, &mataux, &eigaux, &eigintaux,
                           n, maxiter);
    
    /* Initialize */
    retcode = lanczosInitialize(v, V);
    
    /* Start Iterating */
    for (DSDP_INT k = 0; k < maxiter; ++k) {
        
        // w = dS * (LT \ v);
        // w = L \ w;
        spsMatVecBSolve(S, v, w);
        spsMatAx(dS, w, vecaux);
        spsMatVecFSolve(S, vecaux, w);
        
        // norm(wold)
        vec_norm(w, &nrm);
        
        if (k > 0) {
            idx = (k - 1) * n;
            // w = w - H(k, k - 1) * V(:, k - 1);
            memcpy(vecaux->x, &V[idx], sizeof(double) * n);
            vec_axpy(- H[(k - 1) * mH + k], vecaux, w);
        }
        
        // alp = w' * V(:, k);
        idx = k * n;
        alp = ddot(&n, w->x, &one, &V[idx], &one);
        
        // w = w - alp * V(:, k);
        memcpy(vecaux->x, &V[idx], sizeof(double) * n);
        vec_axpy(- alp, vecaux, w);
        // H(k, k) = alp;
        H[k * mH + k] = alp;
        vec_norm(w, &nrm2);
            
        if (nrm2 < 0.8 * nrm) {
            
            alpha = 1.0;
            beta  = 0.0;
            // s = (w' * V(:, 1:k))';
            idx = k + 1;
            memset(mataux, 0, sizeof(double) * idx);
            for (DSDP_INT p = 0, q; p < idx; ++p) {
                for (q = 0; q < n; ++q) {
                    mataux[p] += V[p * n + q] * w->x[q];
                }
            }
            /*
            dgemv(&trans, &n, &idx, &alpha, V, &n, w->x,
                  &one, &beta, mataux, &one);
            */
            // w = w - V(:, 1:k) * s;
            dgemv(&notrans, &n, &idx, &alpha, V, &n, mataux,
                  &one, &alpha, w->x, &one);
            // H(1:k, k) = H(1:k, k) + s;
            daxpy(&idx, &alpha, mataux, &one, &H[mH * k], &one);
            vec_norm(w, &nrm2);
        }
        
        // v = w / nrm;
        
        assert( nrm2 > 0 );
        
        vec_rscale(w, nrm2);
        vec_copy(w, v);
        
        // V(:, k + 1) = v;
        memcpy(&V[n * k + n], v->x, sizeof(double) * n);
        H[mH *  k      + (k + 1)] = nrm2;
        H[mH * (k + 1) +  k     ] = nrm2;
        
        // Check convergence
        if ( ((k + 1) % 5 == 0 || k >= maxiter - 1) ) {
            DSDP_INT kp1 = k + 1;
            // Hk = H(1:k, 1:k);
            for (DSDP_INT i = 0; i < kp1; ++i) {
                memcpy(&mataux[kp1 * i], &H[mH * i], sizeof(double) * kp1);
            }
            
            alpha = -1;
            il = k;
            iu = kp1;
            
            dsyevr(&jobz, &range, &uplo, &kp1, mataux, &kp1,
                   NULL, NULL, &il, &iu, &alpha, &neigs, d,
                   Y, &kp1, isuppz, eigaux, &lwork, eigintaux,
                   &iwork, &info);
            
            
            seig1 = firstNnzSign(&Y[0], kp1);
            seig2 = firstNnzSign(&Y[kp1], kp1);
            
//            if (Y[0] > 0) {
//                seig2 = -1.0;
//            } else {
//                seig2 = 1.0;
//            }
//
//            if (Y[kp1] > 0) {
//                seig1 = -1.0;
//            } else {
//                seig1 = 1.0;
//            }

            if (info) {
                error(etype, "Lanczos Eigenvalue decomposition failed. \n");
            }
            
            res = fabs(H[k * mH + k + 1] * Y[kp1 + k]);
            
            if (res <= 1e-04 || k >= maxiter - 1) {
                //lambda = eigH(idx(k)); lambda2 = eigH(idx(k - 1));
                double lambda1 = d[1], lambda2 = d[0], res1, res2, tmp, gamma;
                alpha = 1.0;
                
                // z = V(:, 1:k) * Y(:, idx(k));
                dgemv(&notrans, &n, &kp1, &alpha, V, &n,
                      &Y[kp1], &one, &beta, z1->x, &one);

                // tmp  = dS * (LT \ z);
                spsMatVecBSolve(S, z1, z2);
                spsMatAx(dS, z2, vecaux);
                spsMatVecFSolve(S, vecaux, z2);
                
                // res1 = norm(L \ tmp - lambda * z);
                vec_axpby(1.0, z2, seig1 * lambda1, z1);
                vec_norm(z1, &res1);
                
                // z2 = V(:, 1:k) * Y(:, idx(k - 1));
                dgemv(&notrans, &n, &kp1, &alpha, V, &n,
                      Y, &one, &beta, z2->x, &one);
                
                // tmp  = dS * (LT \ z2);
                spsMatVecBSolve(S, z2, z1);
                spsMatAx(dS, z1, vecaux);
                spsMatVecFSolve(S, vecaux, z1);
                
                // res2 = norm(L \ tmp - lambda * z);
                vec_axpby(1.0, z1, seig2 * lambda1, z2);
                vec_norm(z2, &res2);
                
                tmp = lambda1 - lambda2 - res2;
                
                if (tmp > 0) {
                    gamma = tmp;
                } else {
                    gamma = 1e-15;
                }
                
                gamma = MIN(res1, res1 * res1 / gamma);
                
                if (gamma < 1e-03 || gamma + lambda1 <= 0.8) {
                    *delta = gamma;
                    *lbd = lambda1;
                    break;
                } else {
                    *delta = gamma;
                    *lbd = lambda1;
                }
            }
        }
    }
    
clean_up:

    DSDP_FREE(isuppz);
    retcode = lanczosFree(&v, &w, &z1, &z2, &vecaux,
                          &V, &H, &Y, &d, &mataux,
                          &eigaux, &eigintaux);
    
    return retcode;
}
