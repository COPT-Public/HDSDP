/** @file adpcg.c
 *  @brief Header for basic types and routine list
 *
 * Solve a sequence of linear systems by either pre-conditioned conjugate gradient or direct solver, which
 * is chosen heuristically based on problem conditioning. The routine also implements an adaptive pre-conditioning
 * mechanism that updates the pre-conditioner automatically. Diagonal and Cholesky pre-conditioners are implemented
 *
 * @author Wenzhi Gao, Shanghai University of Finance and Economics
 * @date Aug 29th, 2022
 *
 */

#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include "adpcg.h"
#include "vec.h"
#include "schurmat.h"

// #define CG_DEBUG

/* Wrapper for the internal methods */
static void cg_dummy_diag( schurMat *A, vec *d ) {
    
    return;
}

static DSDP_INT cg_cond( schurMat *A ) {
    
    return (A->stype == SCHUR_TYPE_SPARSE) ? TRUE : A->denseM->isillCond;
}

static DSDP_INT cg_diagpcd( vec *x, vec *y ) {
    
    vec_vdiv(x, y);
    return DSDP_RETCODE_OK;
}

static DSDP_INT cg_cholpcd( schurMat *A, vec *x, vec *aux ) {
    
    schurMatSolve(A, 1, x->x, aux->x);
    return DSDP_RETCODE_OK;
}

static DSDP_INT cg_Ainv( schurMat *A, vec *v, vec *aux ) {
    
    schurMatSolve(A, 1, v->x, aux->x);
    return DSDP_RETCODE_OK;
}

/** @brief Get time stamp for now
 *  @return Current time stamp
 */
static double my_clock( void ) {
    
    struct timeval t;
    gettimeofday(&t, NULL);
    
    return (1.0e-6 * t.tv_usec + t.tv_sec);
}

/** @brief Logging for conjugate gradient
 *  @param[in] cg CG solver
 *
 */
static void cg_logging( adpcg *cg ) {
    
    printf("| Avg.sv: %3.3f | Avg.fc: %3.3f | currentTime: %3.3f | "
           "latestTime: %3.3f | pcond %3d | nused: %3d | status: %3d | niter: %2d/%2d | rnrm: %3.3e |"
           "nfactors: %3d | nrounds: %2d | nsolverd: %2d | nsolves: %2d | indef: %d |\n",
           cg->avgsvtime, cg->avgfctime, cg->currenttime, cg->latesttime,
           cg->ptype, cg->nused, cg->status, cg->niter, cg->maxiter, cg->rnrm,
           cg->nfactors, cg->nrounds, cg->nsolverd, cg->nsolves, cg->A_cond(cg->A));
}

/** @brief Prepare pre-conditioner for the CG solver
 *  @param[in] cg CG solver
 *  @return CG_OK if the pre-conditioner is successfully collected
 *
 *  The method invokes internal preparation routine of pre-conditioner.
 *  The time computing pre-conditioner is counted into avgsvtime
 */
static DSDP_INT cg_prepare_preconditioner( adpcg *cg ) {
    
    DSDP_INT retcode = CG_OK; double t;
    
    if ( cg->ptype == CG_PRECOND_DIAG ) {
        cg->A_getdiag(cg->A, cg->diag);
    } else {
        t = my_clock();
        retcode = cg->A_chol(cg->A);
        cg->avgfctime = (cg->nfactors * cg->avgfctime + \
                         my_clock() - t) / (cg->nfactors + 1);
        cg->nfactors += 1;
    }
    
    cg->nused = 0; return retcode;
}

/** @brief Apply pre-conditioning
 *  @param[in] cg CG solver
 *  @param[in] v Overwritten by P \ v
 *  @return CG_OK if pre-conditioning is done successfully
 */
static DSDP_INT cg_precondition( adpcg *cg, vec *v ) {
    
    return ( cg->ptype == CG_PRECOND_CHOL ) ? \
             cg->cholpcd(cg->chol, v, cg->aux) : cg->diagpcd(v, cg->diag);
}

/** @brief Decide whether to update pre-conditioner
 *  @param[in] cg CG solver
 *
 * 1.If the system is classified as ill-conditioned or indefinite by some user-defined criterion
 * 2.If the diagonal pre-conditioner is used
 * 3.If latesttime > 1.5 average solution time
 * 4.If ADP-CG is asked to perform direct solve
 * 5.If average solution time > average factorization time
 * 6.If the "nused" property of the current pre-conditioner exceeds the user-defined threshold
 *
 */
static DSDP_INT cg_decision( adpcg *cg ) {
    
    if ( cg->A_cond(cg->A) ) {
        return CG_TRUE;
    }
    
    if ( cg->ptype == CG_PRECOND_DIAG ) {
        return CG_TRUE;
    }
    
    if ( cg->latesttime > 1.5 * cg->avgsvtime + 0.3 * cg->avgfctime ) {
#ifdef CG_DEBUG
        printf("| Current pre-conditioner is considered out-dated with %3.3f > %3.3f. \n",
               cg->latesttime, 1.5 * cg->avgsvtime + 0.3 * cg->avgfctime);
#endif
        return CG_TRUE;
    }
    
    if ( cg->status == CG_STATUS_DIRECT ) {
        return CG_TRUE;
    }

    if ( cg->avgsvtime > cg->avgfctime ) {
        return CG_TRUE;
    }
    
    if ( cg->reuse > 0 && cg->nused > cg->reuse ) {
        return CG_TRUE;
    }
    
    return CG_FALSE;
}

/** @brief Implement conjugate gradinet
 *  @param[in] cg CG solver
 *  @param[in] b RHS vector
 *  @param[in] warm If there is warm start?
 *
 *  Implement pre-conditioned CG with restart
 */
extern DSDP_INT cg_iteration( adpcg *cg, vec *b, DSDP_INT warm ) {
    
    DSDP_INT status = CG_STATUS_UNKNOWN, retcode = CG_OK, i, f;
    double tol = cg->tol, dTMd, rTPinvr, alpha, beta, rnrm;
    
    vec *r = cg->r, *rnew = cg->rnew, *x = cg->x, *d = cg->d, \
        *pinvr = cg->pinvr, *Ad = cg->Ad;
    schurMat *A = cg->A;
    
    cg->v_reset(r);
    
    /* Process restart */
    if ( cg->niter != 0 ) {
        f = (DSDP_INT) cg->niter / 3 + 1;
    } else {
        f = (DSDP_INT) cg->maxiter / 5;
    }
    
    f = MAX(f, 20);
    
    if ( warm ) {
        cg->v_copy(b, r); cg->Av(A, x, d);
        cg->v_axpy(-1.0, d, r); cg->v_norm(r, &alpha); /* alpha = ||Ax0 - b|| */
        if ( alpha < 1e-06 * tol ) {
            status = CG_STATUS_SOLVED; cg->v_copy(x, b);
            cg->status = status; return status;
        }
        cg->v_norm(b, &rnrm);
        if ( rnrm < alpha ) { /* If quality of warm start if worse than 0 */
            cg->v_reset(x); cg->v_copy(b, r); alpha = rnrm;
        }
        
    } else {
        cg->v_reset(x); cg->v_copy(b, r);
    }
    
    cg->v_norm(r, &alpha); rnrm = alpha;
    
    if ( alpha == 0.0 ) {
        status = CG_STATUS_SOLVED;
        cg->status = status; return status;
    }
 
    /* Relative error */
    alpha = MIN(100, alpha);
    tol = MAX(tol * alpha * 0.1, tol * 0.1);
    cg->v_copy(r, d);
    
    retcode = cg_precondition(cg, d);
    if ( retcode != CG_OK ) {
        status = CG_STATUS_FAILED;
        cg->status = status;
        return status;
    }
    
    cg->Av(A, d, Ad); cg->v_copy(d, pinvr);
    
    for ( i = 0; i < cg->maxiter; ++i ) {
        
        cg->v_dot(r, pinvr, &rTPinvr);
        cg->v_dot(d, Ad, &dTMd);
        alpha = rTPinvr / dTMd;
        cg->v_axpy(alpha, d, x);
        
        if ( ((cg->restart == -1 && i % f == 3) ||
             (cg->restart > 0 && i % cg->restart == 1)) ) {
            /* Restart CG solver */
            cg->v_copy(b, r); cg->Av(A, x, d);
            cg->v_axpy(-1.0, d, r);
            cg->v_copy(r, d);
            
            retcode = cg_precondition(cg, d);
            if ( retcode != CG_OK ) {
                status = CG_STATUS_FAILED;
                cg->status = status; break;
            }
            
            cg->Av(A, d, Ad);
            cg->v_copy(r, pinvr);
            
            retcode = cg_precondition(cg, pinvr);
            if ( retcode != CG_OK ) {
                status = CG_STATUS_FAILED;
                cg->status = status; break;
            }
            continue;
        }
        
        cg->v_zaxpby(rnew, 1.0, r, -alpha, Ad);
        cg->v_copy(rnew, pinvr);
        retcode = cg_precondition(cg, pinvr);
        
        if ( retcode != CG_OK ) {
            status = CG_STATUS_FAILED;
            cg->status = status; break;
        }
        
        cg->v_dot(rnew, pinvr, &beta);
        beta = beta / rTPinvr;
        cg->v_axpby(1.0, pinvr, beta, d);
        cg->Av(A, d, Ad);
        cg->v_copy(rnew, r);
        cg->v_norm(r, &rnrm);
        
        if ( rnrm != rnrm ) {
            status = CG_STATUS_FAILED;
            cg->status = status; break;
        }
        
        if ( rnrm < tol ) {
            status = CG_STATUS_SOLVED;
            cg->status = status; break;
        }
    }
    
    if ( status == CG_STATUS_FAILED ) {
        return status;
    }
    
    vec_copy(x, b);
    
    if ( i >= cg->maxiter ) {
        status = CG_STATUS_MAXITER;
    }
    
    cg->niter = i; cg->rnrm = rnrm;
    cg->status = status; return status;
}

/** @brief Initialize the conjugate gradient solver
 *  @param[in] cg  Adaptive CG solver
 *  Initlaize all the pointers to NULL and values to 0
 */
extern void cg_init( adpcg *cg ) {
    
    if ( !cg ) { return; }
    
    memset(cg, 0, sizeof(adpcg));
    cg->ptype = CG_PRECOND_DIAG;    cg->maxiter = -1;
    cg->status = CG_STATUS_UNKNOWN; cg->reuse = -1;
    cg->restart = -1; cg->tol = 1e-08;
    
    /* Link CG methods */
    /* Vector */
    cg->v_init = vec_init; cg->v_alloc = vec_alloc;
    cg->v_free = vec_free; cg->v_copy = vec_copy;
    cg->v_reset = vec_reset; cg->v_norm = vec_norm;
    cg->v_axpy = vec_axpy; cg->v_axpby = vec_axpby;
    cg->v_zaxpby = vec_zaxpby; cg->v_dot = vec_dot;
    
    /* Matrix  */
    cg->A_chol = schurMatFactorize;
    cg->A_cond = cg_cond;
    cg->A_getdiag = cg_dummy_diag;
    
    /* Matrix-vector */
    cg->diagpcd = cg_diagpcd; cg->cholpcd = cg_cholpcd;
    cg->Av = schurMatMx; cg->Ainv = cg_Ainv;
    
    /* Check after linking */
    if ( !cg->v_init || !cg->v_alloc || !cg->v_free ||
         !cg->v_copy || !cg->v_reset || !cg->v_norm ||
         !cg->v_axpy || !cg->v_axpby || !cg->v_zaxpby ||
         !cg->v_dot ) {
        cgerr("Vector operations are not fully initialized. \n");
    }
    
    if ( !cg->A_chol || !cg->A_cond || !cg->A_getdiag ) {
        cgerr("Matrix operations are not fully initialized. \n");
    }
    
    if ( !cg->diagpcd || !cg->cholpcd || !cg->Av || !cg->Ainv ) {
        cgerr("Matrix-vector operations are not fully initialized. \n");
    }
    
    return;
}

/** @brief Allocate internal memory for CG solver
 *  @param[in] cg Adaptive CG solver
 *  @param[in] n Dimension of the linear system
 *  @param[in] vsize Size of vector structure
 *  @return CG_OK if memory is successfully allocated

 Allocate the internal memory for the CG solver.
 */
extern DSDP_INT cg_alloc( adpcg *cg, DSDP_INT n, DSDP_INT vsize ) {
    
    DSDP_INT retcode = CG_OK;
    if ( !cg ) { retcode = CG_ERR; return retcode; }
    
    void (*init) (vec *) = cg->v_init;
    DSDP_INT (*alloc) (vec *, DSDP_INT) = cg->v_alloc;
    
    if ( n <= 0 || vsize <= 0 ) {
        cgerr("Invalid dimension.\n");
        retcode = CG_ERR; return retcode;
    }
    
    cg->n = n;
    
    cg->r = (vec *) calloc(1, vsize); cg->rnew = (vec *) calloc(1, vsize);
    cg->d = (vec *) calloc(1, vsize); cg->pinvr = (vec *) calloc(1, vsize);
    cg->Ad = (vec *) calloc(1, vsize); cg->x = (vec *) calloc(1, vsize);
    cg->aux = (vec *) calloc(1, vsize); cg->btmp = (vec *) calloc(1, vsize);
    
    if ( !cg->r || !cg->rnew || !cg->d || !cg->pinvr ||
         !cg->Ad || !cg->x || !cg->aux || !cg->btmp ) {
        retcode = CG_ERR; return retcode;
    }
    
    init(cg->r); init(cg->rnew); init(cg->d); init(cg->pinvr);
    init(cg->Ad); init(cg->x); init(cg->aux); init(cg->btmp);
    
    retcode = alloc(cg->r, n);     if ( retcode != CG_OK ) { return retcode; }
    retcode = alloc(cg->rnew, n);  if ( retcode != CG_OK ) { return retcode; }
    retcode = alloc(cg->d, n);     if ( retcode != CG_OK ) { return retcode; }
    retcode = alloc(cg->pinvr, n); if ( retcode != CG_OK ) { return retcode; }
    retcode = alloc(cg->Ad, n);    if ( retcode != CG_OK ) { return retcode; }
    retcode = alloc(cg->x, n);     if ( retcode != CG_OK ) { return retcode; }
    retcode = alloc(cg->aux, n);   if ( retcode != CG_OK ) { return retcode; }
    retcode = alloc(cg->btmp, n);  if ( retcode != CG_OK ) { return retcode; }
    
    return retcode;
}

/** @brief Link pointers to LHS matrix, diagonal, and Cholesky pre-conditioner
 *  @param[in] cg CG Solver
 *  @param[in] A Left hand side matrix
 *  @param[in] diag Diagonal matrix
 *  @param[in] chol Cholesky factor
 *
 *  Register pointers for coefficient data, diagonal matrix and Cholesky factor
 */
extern void cg_register( adpcg *cg, schurMat *A, vec *diag, schurMat *chol ) {
    
    if ( !cg ) { return; }
    if ( !A || !diag || !chol ) { return; }
    cg->A = A; cg->diag = diag; cg->chol = chol;
    
    return;
}

/** @brief Free the internal memory of CG solver
 *  @param[in] cg CG solver
 *
 * Free all the internal memory allocated by adaptive CG solver.
 * The solver pointer itself has to be freed by user.
 */
extern void cg_free( adpcg *cg ) {
    
    if ( !cg ) { return; }
    void (*v_free) (vec *) = cg->v_free;
    
    v_free(cg->r); v_free(cg->rnew); v_free(cg->d);
    v_free(cg->pinvr); v_free(cg->Ad); v_free(cg->x);
    v_free(cg->aux); v_free(cg->btmp);

    DSDP_FREE(cg->r); DSDP_FREE(cg->rnew); DSDP_FREE(cg->d);
    DSDP_FREE(cg->pinvr); DSDP_FREE(cg->Ad); DSDP_FREE(cg->x);
    DSDP_FREE(cg->aux); DSDP_FREE(cg->btmp);

    memset(cg, 0, sizeof(adpcg));
    
    return;
}

/** @brief Set parameters for the CG solver
 *  @param[in] cg CG solver
 *  @param[in] tol Relative solution tolerance
 *  @param[in] reuse The maximum reuse number
 *  @param[in] maxiter Maximum of iteration
 *  @param[in] restart Restart frequency of CG. -1 if automatically decided
 *
 * Set CG parameters
 */
extern void cg_setparam( adpcg *cg, double tol,
                         DSDP_INT reuse, DSDP_INT maxiter, DSDP_INT restart ) {
    
    if ( !cg ) { return; }
    if ( tol <= 0.0 || maxiter <= 0 ) {
        cgerr("Invalid parameter range for tol and maxiter. \n");
    }
    
    cg->tol = tol; cg->reuse = reuse;
    cg->maxiter = maxiter; cg->restart = restart;

    return;
}

/** @brief Extract CG statistics after some solve
 *  @param[in] cg CG solver
 *  @param[out] status CG solution status
 *  @param[out] niter Number of CG iterations
 *  @param[out] rnrm Residual norm
 *  @param[out] avgsvtime Current avarage solution time
 *  @param[out] avgfctime Current average factorization time
 *  @param[out] nused Current number of iterations the pre-conditioner is used
 *  @param[out] nmaixter  Number of iterations CG fails to solve the system
 *  @param[out] nfactors Number of factorizations
 *  @param[out] nrounds Number of rounds
 *  @param[out] nsolverd Number of solves in the current round
 *  @param[out] nsolves Number of solves
 *
 *  Collect different CG statistics. If some statistic is not needed, just let it be NULL.
 */
extern void cg_getstats( adpcg *cg, DSDP_INT *status, DSDP_INT *niter, double *rnrm,
                         double *avgsvtime, double *avgfctime, DSDP_INT *nused,
                         DSDP_INT *nmaixter, DSDP_INT *nfactors, DSDP_INT *nrounds,
                         DSDP_INT *nsolverd, DSDP_INT *nsolves ) {
    
    if ( !cg ) { return; }
    if ( status ) { *status = cg->status; }
    if ( niter ) { *niter = cg->niter; }
    if ( rnrm ) { *rnrm = cg->rnrm; }
    if ( avgsvtime ) { *avgsvtime = cg->avgsvtime; }
    if ( avgfctime ) { *avgfctime = cg->avgfctime; }
    if ( nused ) { *nused = cg->nused; }
    if ( nmaixter ) { *nmaixter = cg->nmaxiter; }
    if ( nsolves ) { *nsolves = cg->nsolves; }
    if ( nsolverd ) { *nsolverd = cg->nsolverd; }
    if ( nfactors ) { *nfactors = cg->nfactors; }

    return;
}

/** @brief Start a round of solves
 *  @param[in] cg CG solver
 *
 *  Several rules decide whether to update the current pre-conditioner
 */
extern DSDP_INT cg_start( adpcg *cg ) {
    
    DSDP_INT retcode = CG_OK, decision = CG_FALSE;
    /* Update pre-conditioner */
    if ( cg->A->stype == SCHUR_TYPE_SPARSE ) {
        cg->status = CG_STATUS_DIRECT;
        cg->ptype = CG_PRECOND_CHOL;
    }
    
    decision = cg_decision(cg);
    if ( decision ) {
        cg_prepare_preconditioner(cg);
#ifdef CG_DEBUG
        printf("| CG round %d decides to update pre-conditioner. \n", cg->nrounds);
#endif
    } else {
#ifdef CG_DEBUG
        printf("| CG round %d decides NOT to update pre-conditioner. \n", cg->nrounds);
#endif
    }
    
    return retcode;
}

/** @brief Finish a round of solve
 *  @param[in] cg CG solver
 *
 * At the end of each round,
 * 1. nused increases by 1
 * 2. currenttime overwrites latesttime
 * 3. nsolverd is reset
 */
extern void cg_finish( adpcg *cg ) {
    
    cg->nused += 1; cg->latesttime = ( cg->nsolverd ) ? cg->currenttime / cg->nsolverd : 0.0;
    cg->currenttime = 0.0; cg->nsolverd = 0; cg->nrounds += 1;
    
    if ( cg->nmaxiter > 20 ) {
        cg->status = CG_STATUS_DIRECT;
    }
    
    if ( cg->latesttime > 2.0 * cg->avgfctime && cg->nfactors > 0 ) {
        printf("| ADP-CG Warning: No longer re-using pre-conditioners. \n");
        cg->status = CG_STATUS_DIRECT;
    }
    
    return;
}

/** @brief Solve linear system using adaptive pre-conditioned conjugate gradient with restart
 *  @param[in] cg CG Solver
 *  @param[in] b  RHS vector. Overwritten when solved
 *  @param[in] x0 Initial point
 *
 *  Solve the linear system by adaptive pre-conditioning and restart.
 */
extern DSDP_INT cg_solve( adpcg *cg, vec *b, vec *x0 ) {
    
    DSDP_INT retcode = CG_OK, status, warm = CG_FALSE;
    double t;
    cg->niter = 0;
    
    if ( cg->status == CG_STATUS_DIRECT ) {
        return cg->Ainv(cg->A, b, cg->aux);
    }
    
    if ( cg->ptype == CG_PRECOND_CHOL &&
         cg->nused == 0 ) {
        return cg->Ainv(cg->A, b, cg->aux);
    }
    
    t = my_clock(); cg->v_copy(b, cg->btmp);
    
    if ( x0 ) {
        cg->v_copy(x0, cg->x); warm = CG_TRUE;
    }
    
    status = cg_iteration(cg, b, warm);
    
    if ( status == CG_STATUS_FAILED ) {
        retcode = CG_ERR; return retcode;
    }
      
    if ( status == CG_STATUS_MAXITER ) {
        /* Allow for regret if in the first solve */
        cg->nmaxiter += 1;
        if ( cg->nsolverd == 0 ) {
//            printf("| Regret happens at round %d. \n", cg->nrounds);
            cg->ptype = CG_PRECOND_CHOL;
            retcode = cg_prepare_preconditioner(cg);
            vec_copy(cg->btmp, b);
            return ( retcode == CG_OK ) ? cg_solve(cg, b, cg->x) : retcode;
        }
    }
    
    t = my_clock() - t; cg->currenttime += t; cg->nsolverd += 1;
    cg->avgsvtime = (cg->avgsvtime * cg->nsolves + t) / (cg->nsolves + 1);
    cg->nsolves += 1;
    
#ifdef CG_DEBUG
    cg_logging(cg);
#endif
    
    return retcode;
}
