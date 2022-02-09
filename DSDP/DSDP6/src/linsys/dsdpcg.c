#include <assert.h>
#include "sparsemat.h"
#include "densemat.h"
#include "dsdpcg.h"

static void diagPrecond( vec *P, vec *x ) {
    // Diagonal preconditioning
    assert( P->dim == x->dim );
    for (DSDP_INT i = 0; i < x->dim; ++i) {
        x->x[i] /= P->x[i];
    }
    return;
}

static void cholPrecond( dsMat *P, vec *x ) {
    // Cholesky preconditioning
    assert( P->dim == x->dim );
    if (!P->isFactorized) {
        denseMatFactorize(P);
    }
    denseArrSolveInp(P, 1, x->x);
    return;
}

static void preCond( CGSolver *cgSolver, vec *x ) {
    
    if (cgSolver->pType == CG_PRECOND_DIAG) {
        diagPrecond(cgSolver->preCond, x);
    } else {
        cholPrecond(cgSolver->preCond, x);
    }
    return;
}

extern DSDP_INT dsdpCGInit( CGSolver *cgSolver ) {
    
    DSDP_INT retcode = DSDP_RETCODE_OK;
    
    cgSolver->M     = NULL;
    cgSolver->r     = NULL;
    cgSolver->rnew  = NULL;
    cgSolver->d     = NULL;
    cgSolver->Pinvr = NULL;
    cgSolver->x     = NULL;
    cgSolver->Md    = NULL;
    
    cgSolver->dim     = 0;
    cgSolver->tol     = 1e-04;
    cgSolver->resinrm = 0.0;
    cgSolver->niter   = 0;
    cgSolver->maxiter = 100;
    
    return retcode;
}

extern DSDP_INT dsdpCGAlloc( CGSolver *cgSolver, DSDP_INT m ) {
    
    DSDP_INT retcode = DSDP_RETCODE_OK;
    
    assert( m > 0 );
    
    cgSolver->dim   = m;
    cgSolver->r     = (vec *) calloc(1, sizeof(vec));
    cgSolver->rnew  = (vec *) calloc(1, sizeof(vec));
    cgSolver->d     = (vec *) calloc(1, sizeof(vec));
    cgSolver->Pinvr = (vec *) calloc(1, sizeof(vec));
    cgSolver->Md    = (vec *) calloc(1, sizeof(vec));
    cgSolver->x     = (vec *) calloc(1, sizeof(vec));
    
    vec_init(cgSolver->r);
    vec_init(cgSolver->rnew);
    vec_init(cgSolver->d);
    vec_init(cgSolver->Pinvr);
    vec_init(cgSolver->Md);
    vec_init(cgSolver->x);
    
    vec_alloc(cgSolver->r, m);
    vec_alloc(cgSolver->rnew, m);
    vec_alloc(cgSolver->d, m);
    vec_alloc(cgSolver->Pinvr, m);
    vec_alloc(cgSolver->Md, m);
    vec_alloc(cgSolver->x, m);

    return retcode;
}

extern DSDP_INT dsdpCGFree( CGSolver *cgSolver ) {
    
    DSDP_INT retcode = DSDP_RETCODE_OK;
    
    vec_free(cgSolver->r);
    vec_free(cgSolver->rnew);
    vec_free(cgSolver->d);
    vec_free(cgSolver->Pinvr);
    vec_free(cgSolver->Md);
    vec_free(cgSolver->x);
    
    DSDP_FREE(cgSolver->r);
    DSDP_FREE(cgSolver->rnew);
    DSDP_FREE(cgSolver->d);
    DSDP_FREE(cgSolver->Pinvr);
    DSDP_FREE(cgSolver->Md);
    DSDP_FREE(cgSolver->x);
    
    cgSolver->M     = NULL;
    cgSolver->tol   = 0.0;
    cgSolver->dim   = 0;
    cgSolver->niter = 0;
    
    return retcode;
}

extern DSDP_INT dsdpCGSetTol( CGSolver *cgSolver, double tol ) {
    assert( tol > 0 );
    cgSolver->tol = MAX(1e-20, tol);
    return DSDP_RETCODE_OK;
}

extern DSDP_INT dsdpCGSetMaxIter( CGSolver *cgSolver, DSDP_INT maxiter ) {
    assert( maxiter >= 1 );
    cgSolver->niter = MAX(maxiter, 1);
    return DSDP_RETCODE_OK;
}

extern DSDP_INT dsdpCGSetM( CGSolver *cgSolver, dsMat *M ) {
    assert( M->dim == cgSolver->dim );
    cgSolver->M = M;
    return DSDP_RETCODE_OK;
}

extern DSDP_INT dsdpCGSetPreCond( CGSolver *cgSolver, DSDP_INT pType, void *preCond ) {
    
    assert( pType == CG_PRECOND_DIAG || pType == CG_PRECOND_CHOL );
    cgSolver->pType = pType;
    cgSolver->preCond = preCond;
    
    if (pType == CG_PRECOND_CHOL) {
        if (((dsMat *) preCond)->isillCond) {
            return DSDP_RETCODE_FAILED;
        }
    }
    return DSDP_RETCODE_OK;
}

extern DSDP_INT dsdpGetCGSolStatistic( CGSolver *cgSolver, DSDP_INT *status, double *resinorm ) {
    *status = cgSolver->status;
    *resinorm = cgSolver->resinrm;
    return DSDP_RETCODE_OK;
}

/* Implement (pre-conditioned) conjugate gradient for Schur system */
extern DSDP_INT dsdpCGSolve( CGSolver *cgSolver, vec *b, vec *x0 ) {
    /* CG solver for schur system. Solve P^-1 M x = P^-1 b
       On exit b is over-written
     */
    DSDP_INT retcode = DSDP_RETCODE_OK;
    
    /* Initialize
     x = x0;
     r = b - M * x;
     d = P \ r;
     Md = M * d;
     Pinvr = P \ r;
    */
    
    dsMat *M = cgSolver->M;
    vec *r = cgSolver->r, *rnew = cgSolver->rnew, *x = cgSolver->x,
        *d = cgSolver->d, *Pinvr = cgSolver->Pinvr, *Md = cgSolver->Md;
    
    double tol = cgSolver->tol, dTMd = 0.0, rTPinvr = 0.0,
           alpha = 0.0, beta = 0.0, resinorm = DSDP_INFINITY;
    cgSolver->status = CG_STATUS_UNKNOWN;
    DSDP_INT iter = 0;
    
    vec_reset(r);
    
    if (x0) {
        vec_copy(x0, x);
        // r = b - M * x;
        vec_copy(b, r);
        denseMataAxpby(M, -1.0, x, 1.0, r);
    } else {
        vec_reset(x);
        // r = b;
        vec_copy(b, r);
    }
    
    // d = P \ r;
    vec_copy(r, d);
    preCond(cgSolver, d);
    // Md = M * d;
    denseMataAxpby(M, 1.0, d, 0.0, Md);
    // Pinvr = P \ r;
    vec_copy(d, Pinvr);
    
    // Start CG
    for (iter = 0; iter < cgSolver->maxiter; ++iter) {
        // alpha = r' * Pinvr / (d' * Md);s
        vec_dot(r, Pinvr, &rTPinvr);
        vec_dot(d, Md, &dTMd);
        alpha = rTPinvr / dTMd;
        // x = x + alpha * d;
        vec_axpy(alpha, d, x);
        // rnew = r - alpha * Md
        vec_zaxpby(rnew, 1.0, r, -alpha, Md);
        // Pinvr = P \ rnew;
        vec_copy(rnew, Pinvr);
        preCond(cgSolver, Pinvr);
        // beta = rnew' * Pinvr / rTPinvr;
        vec_dot(rnew, Pinvr, &beta);
        beta = beta / rTPinvr;
        // d = Pinvr + beta * d;
        vec_axpby(1.0, Pinvr, beta, d);
        // Md = M * d;
        denseMataAxpby(M, 1.0, d, 0.0, Md);
        // r = rnew
        vec_copy(rnew, r);
        // norm(r)
        vec_norm(r, &resinorm);
        
        printf("%3d %10.3e %20.10e \n", iter + 1, alpha, resinorm);
        
        if (resinorm <= tol) {
            cgSolver->status = CG_STATUS_SOLVED;
            break;
        }
    }
    
    if (iter == cgSolver->maxiter - 1 && resinorm > tol) {
        cgSolver->status = CG_STATUS_MAXITER;
    }
    
    vec_copy(x, b);

    return retcode;
}
