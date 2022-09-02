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

#define adpcg_free(var) do {free((var)); (var) = NULL;} while (0)

/** @brief Get time stamp for now
 *  @return Current time stamp
 */

static double my_clock( void ) {
    
    struct timeval t;
    gettimeofday(&t, NULL);
    
    return (1.0e-6 * t.tv_usec + t.tv_sec);
}

/** @brief Prepare pre-conditioner for the CG solver
 *  @param[in] cg CG solver
 *  @return CG_OK if the pre-conditioner is successfully collected
 *
 *  The method invokes internal preparation routine of pre-conditioner.
 *  The time computing pre-conditioner is counted into avgsvtime
 */
static cgint cg_prepare_preconditioner( adpcg *cg ) {
    
    cgint retcode = CG_OK; double t;
    
    if ( cg->ptype == CG_PRECOND_DIAG ) {
        cg->A_getdiag(cg->diag, cg->A);
    } else {
        t = my_clock();
        retcode = cg->A_chol(cg->A);
        cg->avgfctime = (cg->nfactors * cg->avgfctime + \
                         t - my_clock()) / (cg->nfactors + 1);
        cg->nfactors += 1; cg->nused = 0;
    }
    
    return retcode;
}

/** @brief Apply pre-conditioning
 *  @param[in] cg CG solver
 *  @param[in] v Overwritten by P \ v
 *  @return CG_OK if pre-conditioning is done successfully
 */
static cgint cg_precondition( adpcg *cg, void *v ) {
    
    return ( cg->ptype == CG_PRECOND_CHOL ) ? \
             cg->cholpcd(cg->chol, v, cg->aux) : cg->diagpcd(cg->diag, v);
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
static cgint cg_decision( adpcg *cg ) {
    
    if ( cg->A_cond(cg->A) ) {
        return CG_TRUE;
    }
    
    if ( cg->ptype == CG_PRECOND_DIAG ) {
        return CG_TRUE;
    }
    
    if ( cg->latesttime > 1.5 * cg->avgsvtime ) {
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
extern cgint cg_iteration( adpcg *cg, void *b, cgint warm ) {
    
    cgint status = CG_STATUS_UNKNOWN, retcode = CG_OK, i;
    double tol = cg->tol, dTMd, rTPinvr, alpha, beta, rnrm;
    
    void *r = cg->r, *rnew = cg->rnew, *x = cg->x, *d = cg->d, \
         *pinvr = cg->pinvr, *Ad = cg->Ad, *A = cg->A;
    
    cg->v_reset(r);
    
    if ( warm ) {
        cg->v_copy(b, r); cg->Av(A, x, d);
        cg->v_axpy(-1.0, d, r); cg->v_norm(r, &alpha);
        if ( alpha < 1e-06 * tol ) {
            status = CG_STATUS_SOLVED; cg->v_copy(x, b);
            return status;
        }
    } else {
        cg->v_reset(x); cg->v_copy(b, r);
    }
    
    cg->v_norm(r, &alpha); cg->v_norm(b, &rnrm);
    
    if ( rnrm < alpha ) {
        cg->v_reset(x); cg->v_copy(b, r); alpha = rnrm;
    }
    
    if ( alpha == 0.0 ) {
        status = CG_STATUS_SOLVED;
        return status;
    }
 
    alpha = MIN(100, alpha);
    tol = MAX(tol * alpha * 0.1, tol * 0.1);
    cg->v_copy(r, d);
    
    retcode = cg_precondition(cg, d);
    if ( retcode != CG_OK ) {
        status = CG_STATUS_FAILED;
        return status;
    }
    
    cg->Av(A, d, Ad); cg->v_copy(d, pinvr);
    
    for ( i = 0; i < cg->maxiter; ++i ) {
        
        cg->v_dot(r, pinvr, &rTPinvr);
        cg->v_dot(d, Ad, &dTMd);
        alpha = rTPinvr / dTMd;
        cg->v_axpy(alpha, d, x);
        
        if ( (cg->restart == -1 && i % 50 == 1) ||
             (cg->restart > 0 && i % cg->restart == 1) ) {
            /* Restart CG solver */
            cg->v_copy(b, r); cg->Av(A, x, d);
            cg->v_axpy(-1.0, d, r);
            cg->v_copy(r, d);
            
            retcode = cg_precondition(cg, d);
            if ( retcode != CG_OK ) {
                status = CG_STATUS_FAILED; break;
            }
            
            cg->Av(A, d, Ad);
            cg->v_copy(r, pinvr);
            
            retcode = cg_precondition(cg, pinvr);
            if ( retcode != CG_OK ) {
                status = CG_STATUS_FAILED; break;
            }
            continue;
        }
        
        cg->v_zaxpby(rnew, 1.0, r, -alpha, Ad);
        cg->v_copy(rnew, pinvr);
        retcode = cg_precondition(cg, pinvr);
        
        if ( retcode != CG_OK ) {
            status = CG_STATUS_FAILED; break;
        }
        
        cg->v_dot(rnew, pinvr, &beta);
        beta = beta / rTPinvr;
        cg->v_axpby(1.0, pinvr, beta, d);
        cg->Av(A, d, Ad);
        cg->v_copy(rnew, r);
        cg->v_norm(r, &rnrm);
        
        if ( rnrm != rnrm ) {
            status = CG_STATUS_FAILED; break;
        }
        
        if ( rnrm < tol ) {
            status = CG_STATUS_SOLVED; break;
        }
    }
    
    if ( status == CG_STATUS_FAILED ) {
        return status;
    }
    
    if ( i >= cg->maxiter ) {
        status = CG_STATUS_MAXITER;
    }
    
    return status;
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
    cg->restart = -1;
    
    /* Link CG methods */
    
    /* Vector */
    cg->v_init = NULL; cg->v_alloc = NULL; cg->v_free = NULL;
    cg->v_copy = NULL; cg->v_reset = NULL; cg->v_norm = NULL;
    cg->v_axpy = NULL; cg->v_axpby = NULL; cg->v_zaxpby = NULL;
    cg->v_dot = NULL;
    
    /* Matrix  */
    cg->A_chol = NULL; cg->A_cond = NULL; cg->A_getdiag = NULL;
    
    /* Matrix-vector */
    cg->diagpcd = NULL; cg->cholpcd = NULL; cg->Av = NULL;
    cg->Ainv = NULL;
    
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
extern cgint cg_alloc( adpcg *cg, cgint n, cgint vsize ) {
    
    cgint retcode = CG_OK;
    if ( !cg ) { retcode = CG_ERR; return retcode; }
    
    void (*init) (void *) = cg->v_init;
    cgint (*alloc) (void *, cgint) = cg->v_alloc;
    
    if ( n <= 0 || vsize <= 0 ) {
        cgerr("Invalid dimension.\n");
        retcode = CG_ERR; return retcode;
    }
    
    cg->n = n;
    
    cg->r = (void *) calloc(1, vsize); cg->rnew = (void *) calloc(1, vsize);
    cg->d = (void *) calloc(1, vsize); cg->pinvr = (void *) calloc(1, vsize);
    cg->Ad = (void *) calloc(1, vsize); cg->x = (void *) calloc(1, vsize);
    cg->aux = (void *) calloc(1, vsize);
    
    if ( !cg->r || !cg->rnew || !cg->d ||
         !cg->pinvr || !cg->Ad || !cg->x || !cg->aux ) {
        retcode = CG_ERR; return retcode;
    }
    
    init(cg->r); init(cg->rnew); init(cg->d); init(cg->pinvr);
    init(cg->Ad); init(cg->x); init(cg->aux);
    
    retcode = alloc(cg->r, n);     if ( retcode != CG_OK ) { return retcode; }
    retcode = alloc(cg->rnew, n);  if ( retcode != CG_OK ) { return retcode; }
    retcode = alloc(cg->d, n);     if ( retcode != CG_OK ) { return retcode; }
    retcode = alloc(cg->pinvr, n); if ( retcode != CG_OK ) { return retcode; }
    retcode = alloc(cg->Ad, n);    if ( retcode != CG_OK ) { return retcode; }
    retcode = alloc(cg->x, n);     if ( retcode != CG_OK ) { return retcode; }
    retcode = alloc(cg->aux, n);   if ( retcode != CG_OK ) { return retcode; }
    
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
extern void cg_register( adpcg *cg, void *A, void *diag, void *chol ) {
    
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
    void (*v_free) (void *) = cg->v_free;
    
    v_free(cg->r); v_free(cg->rnew); v_free(cg->d);
    v_free(cg->pinvr); v_free(cg->Ad); v_free(cg->x);
    v_free(cg->aux);

    adpcg_free(cg->r); adpcg_free(cg->rnew); adpcg_free(cg->d);
    adpcg_free(cg->pinvr); adpcg_free(cg->Ad); adpcg_free(cg->x);
    adpcg_free(cg->aux);

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
                         cgint reuse, cgint maxiter, cgint restart ) {
    
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
extern void cg_getstats( adpcg *cg, cgint *status, cgint *niter, double *rnrm,
                         double *avgsvtime, double *avgfctime, cgint *nused,
                         cgint *nmaixter, cgint *nfactors, cgint *nrounds,
                         cgint *nsolverd, cgint *nsolves ) {
    
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
extern cgint cg_start( adpcg *cg ) {
    
    cgint retcode = CG_OK, decision = CG_FALSE;
    /* Update pre-conditioner */
    decision = cg_decision(cg);
    if ( decision ) {
        cg_prepare_preconditioner(cg);
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
    
    cg->nused += 1; cg->latesttime = cg->currenttime / cg->nsolverd;
    cg->currenttime = 0.0; cg->nsolverd = 0;
    
    if ( cg->nmaxiter > 20 ) {
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
extern cgint cg_solve( adpcg *cg, void *b, void *x0 ) {
    
    cgint retcode = CG_OK, status, warm = CG_FALSE;
    double t;
    
    if ( cg->status == CG_STATUS_DIRECT ) {
        return cg->Ainv(cg->A, b, cg->aux);
    }
    
    if ( cg->ptype == CG_PRECOND_CHOL &&
         cg->nused == 0 ) {
        return cg->Ainv(cg->A, b, cg->aux);
    }
    
    t = my_clock();
    cg->v_copy(b, cg->aux);
    
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
        if ( cg->ptype == CG_PRECOND_DIAG && cg->nsolverd == 0 ) {
            cg->ptype = CG_PRECOND_CHOL;
            retcode = cg_prepare_preconditioner(cg);
            return ( retcode == CG_OK ) ? cg_solve(cg, b, x0) : retcode;
        }
    }
    
    cg->currenttime += my_clock() - t;
    cg->nsolverd += 1; cg->nsolves += 1;
    return retcode;
}
