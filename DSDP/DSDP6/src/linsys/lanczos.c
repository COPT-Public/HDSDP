/** @file lanczos.c
 *  @brief Implement the Lanczos method to compute step-length in SDP
 *
 * Implement the Lanczos routine to compute the maximum \f$ \alpha \f$ such that
 * \f$ S + \alpha \Delta S \f$ is still in the cone.
 *
 *  @author Wenzhi Gao, Shanghai University of Finance and Economics
 *  @date Aug, 28th, 2022
 *  
 */

/* Include headers */
#include "lanczos.h"
#include "sparseopts.h"
#include "vec.h"
#include "dsdplapack.h"

static char etype[] = "Lanczos iteration";

#define MAXITER (30)

/** @brief Initialize the Lanczos iteration using random starting point
 *  @param[out] v Starting vector of Lanczos
 *  @param[out] V Matrix of past "v" s
 *
 *  Initialize v to start the Lanczos iteration. v = randn(n, 1) / norm(v);
 *  Sometimes the random initialization causes numerical difficulties if matrix is small.
 *  So be careful
 */
static void lczIterInit( vec *v, double *V ) {
    
    srand(v->dim); double nrm = 0.0;
    for (DSDP_INT i = 0; i < v->dim; ++i) {
        srand(rand());
        v->x[i] = sqrt(sqrt((rand() % 1627))) * (rand() % 2 - 0.5);
    }
    
    vec_norm(v, &nrm); vec_rscale(v, nrm);
    memcpy(V, v->x, sizeof(double) * v->dim); // V(:, 1) = v;
    
    return;
}

/** @brief Initialize the Lanczos solver
 *  @param[out] lcz The Lanczos solver struct
 *
 *  All the sizes are set 0 and pointers are set NULL
 */
extern void lczInit( lczstep *lcz ) {
    
    lcz->n = 0; lcz->v = NULL; lcz->w = NULL; lcz->z1 = NULL;
    lcz->z2 = NULL; lcz->va = NULL; lcz->V = NULL;
    lcz->H = NULL; lcz->Y = NULL; lcz->d = NULL;
    lcz->mataux = NULL; lcz->eigaux = NULL;
    lcz->eigintaux = NULL; lcz->lwork = LWORK * MAXITER;
    lcz->iwork = IWORK * MAXITER;
    
    return;
}

/** @brief Allocate the internal space for Lancozs struct
 *  @param[out] lcz Lanczos solver struct
 *  @param[in] n Size of the matrix
 *  @return DSDP_RETCODE_OK if memory is successfully allocated
 *
 * The memory is related to the MAXITER macro that defines how many iterations Lanczos would run.
 */
extern DSDP_INT lczAlloc( lczstep *lcz, DSDP_INT n ) {
    
    DSDP_INT retcode = DSDP_RETCODE_OK, maxiter = MAXITER; lcz->n = n;
    vec *v = (vec *) calloc(1, sizeof(vec));
    vec_init(v); retcode = vec_alloc(v, n); checkCode; lcz->v = v;
    v = (vec *) calloc(1, sizeof(vec));
    vec_init(v); retcode = vec_alloc(v, n); checkCode; lcz->w = v;
    v = (vec *) calloc(1, sizeof(vec));
    vec_init(v); retcode = vec_alloc(v, n); checkCode; lcz->z1 = v;
    v = (vec *) calloc(1, sizeof(vec));
    vec_init(v); retcode = vec_alloc(v, n); checkCode; lcz->z2 = v;
    v = (vec *) calloc(1, sizeof(vec));
    vec_init(v); retcode = vec_alloc(v, n); checkCode; lcz->va = v;
    
    lcz->V = (double *) calloc(n * (maxiter + 1), sizeof(double));
    lcz->H = (double *) calloc((maxiter + 1) * (maxiter + 1), sizeof(double));
    lcz->Y = (double *) calloc(2 * maxiter, sizeof(double));
    lcz->d = (double *) calloc(maxiter, sizeof(double));
    lcz->mataux = (double *) calloc(maxiter * maxiter, sizeof(double));
    lcz->eigaux = (double *) calloc(maxiter * LWORK, sizeof(double));
    lcz->eigintaux = (DSDP_INT *) calloc(maxiter * IWORK, sizeof(DSDP_INT));
    
    if (!lcz->V || !lcz->H || !lcz->Y || !lcz->d ||
        !lcz->mataux || !lcz->eigaux || !lcz->eigintaux) {
        printf("| Failed to allocate memory for Lanczos iterator. \n");
        retcode = DSDP_RETCODE_FAILED;
    }
    
    return retcode;
}

/** @brief Free the memory allocated by Lanczos struct
 *  @param[out] lcz The Lanczos solver struct
 *
 *  The routine frees the internal memory allocated by the Lanczos solver
 *  The solver pointer itself should be freed by the user
 */
extern void lczFree( lczstep *lcz ) {
    
    lcz->n = 0; lcz->lwork = 0; lcz->iwork = 0;
    vec_free(lcz->v);  DSDP_FREE(lcz->v);
    vec_free(lcz->w);  DSDP_FREE(lcz->w);
    vec_free(lcz->z1); DSDP_FREE(lcz->z1);
    vec_free(lcz->z2); DSDP_FREE(lcz->z2);
    vec_free(lcz->va); DSDP_FREE(lcz->va);
    
    DSDP_FREE(lcz->V); DSDP_FREE(lcz->H); DSDP_FREE(lcz->Y);
    DSDP_FREE(lcz->d); DSDP_FREE(lcz->mataux); DSDP_FREE(lcz->eigaux);
    DSDP_FREE(lcz->eigintaux);
    
    return;
}

/** @brief Compute the maximum step length using Lanczos algorithm
 *  @param[in] lcz Lanczos struct
 *  @param[in] S Dual iteration matrix \f$ S \f$
 *  @param[in] dS Dual step matrix \f$ \Delta S \f$
 *  @param[out] lbd  Estimated \f$ \lambda_1(L^{-1} \Delta S L^{-T}) \f$
 *  @param[out] delta Error in the Lanczos estimation
 *
 * The algorithm implements a classic Lanczos method that estimates the maximum stepsize the
 * current dual iterate can take without getting out of the cone
 */
extern void lczStepLength( lczstep *lcz, spsMat *S, spsMat *dS, double *lbd, double *delta ) {
    // Lanczos algorithm that computes lambda_max(L^-1 dS L^-T)
    lczstep *l = lcz; double fnrm = 0.0; spsMatFnorm(dS, &fnrm);
    
    if (fnrm < 1e-13) {
        *lbd = 0.0; *delta = 0.0; return;
    }
    
    /* Prepare working array */
    DSDP_INT n = S->dim, maxiter = MAXITER, one = 1, mH = maxiter + 1, i, k, kp1,
             il, iu, lwork = lcz->lwork, iwork = lcz->iwork, failed = FALSE, \
             neigs, info, *eigintaux = lcz->eigintaux, *isuppz = lcz->isuppz;
    vec *v = l->v, *w = l->w, *z1 = l->z1, *z2 = l->z2, *va = l->va;
    double *V = l->V, *H = l->H, *Y = l->Y, *d = l->d, *mataux = l->mataux, \
           *eigaux = l->eigaux, *Vp, lbd1, lbd2, nrm, nrm2, alp, alpha, \
           beta = 0.0, res, res1, res2, tmp, gm;
    char notrans = DSDP_MAT_NOTRANSPOSE, jobz = 'V', range = 'I', uplo = DSDP_MAT_UP;
    
    /* Initialize */
    lczIterInit(v, V);
    
    /* Start iterating */
    for (k = 0; k < maxiter; ++k) {
        
        spsMatVecBSolve(S, v, w); // w = (LT \ v)
        spsMatAx(dS, w, va); // va = dS * (LT \ v);
        spsMatVecFSolve(S, va, w); vec_norm(w, &nrm); // w = L \ va;
        
        if (k > 0) {
            memcpy(va->x, V + (k - 1) * n, sizeof(double) * n);
            vec_axpy(- H[(k - 1) * mH + k], va, w); // w = w - H(k, k - 1) * V(:, k - 1);
        }
        
        Vp = V + k * n;
        alp = -ddot(&n, w->x, &one, Vp, &one); // alp = - w' * V(:, k);
        daxpy(&n, &alp, Vp, &one, w->x, &one);  // w = w + alp * V(:, k);
        vec_norm(w, &nrm2); // nrm2 = norm(w)
        H[k * mH + k] = - alp; // H(k, k) = - alp;
        
        if (nrm2 < 0.8 * nrm && FALSE) {
            for (i = 0; i < k + 1; ++i) {
                Vp = V + i * n;
                alpha = - ddot(&n, Vp, &one, w->x, &one);
                daxpy(&n, &alpha, Vp, &one, w->x, &one); // w = w - V(:, 1:k) * s;
                H[mH * k + i] -= alpha; // H(1:k, k) = H(1:k, k) + s;
            }
            vec_norm(w, &nrm2);
        }
        
        if (nrm2 > 0.0) {
            vec_rscale(w, nrm2); vec_copy(w, v);
            memcpy(V + n * k + n, w->x, sizeof(double) * n); // V(:, k + 1) = w / norm(w);
            H[mH * k + (k + 1)] = H[mH * (k + 1) + k] = nrm2;
        }
        
        /* Check convergence */
        if ( ((k + 1) % 5 == 0 || k >= maxiter - 1 ) || nrm2 == 0 ) {
            kp1 = k + 1;
            // Hk = H(1:k, 1:k);
            for (i = 0; i < kp1; ++i) {
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
                
                lbd1 = d[1]; lbd2 = d[0]; alpha = 1.0; // lambda = eigH(idx(k)); lambda2 = eigH(idx(k - 1));
                
                dgemv(&notrans, &n, &kp1, &alpha, V, &n,
                      &Y[kp1], &one, &beta, z1->x, &one); // z = V(:, 1:k) * Y(:, idx(k));
                
                spsMatVecBSolve(S, z1, z2);
                spsMatAx(dS, z2, va);
                spsMatVecFSolve(S, va, z2); // tmp = dS * (LT \ z);
                
                vec_axpby(1.0, z2, -lbd1, z1);
                vec_norm(z1, &res1); // res1 = norm(L \ tmp - lambda * z);
                
                dgemv(&notrans, &n, &kp1, &alpha, V, &n,
                      Y, &one, &beta, z2->x, &one); // z2 = V(:, 1:k) * Y(:, idx(k - 1));
                
                spsMatVecBSolve(S, z2, z1);
                spsMatAx(dS, z1, va);
                spsMatVecFSolve(S, va, z1); // tmp  = dS * (LT \ z2);
                
                
                vec_axpby(1.0, z1, -lbd1, z2);
                vec_norm(z2, &res2); // res2 = norm(L \ tmp - lambda * z);
                
                tmp = lbd1 - lbd2 - res2;
                gm = (tmp > 0) ? tmp : 1e-15;
                tmp = res1 * res1 / gm;
                gm = MIN(res1, tmp);
                
                if (gm < 1e-03 || gm + lbd1 <= 0.8 || failed) {
                    if (failed) {
                        lbd1 = MAX(1000.0, lbd1);
                    }
                    *delta = gm; *lbd = lbd1; break;
                } else {
                    if (nrm2 == 0.0) {
                        fatal_error_msg("Bad Lanczos initialization. \n");
                    }
                    *delta = gm; *lbd = lbd1;
                }
            }
        }
    }
    
    return;
}
