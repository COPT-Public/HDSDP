#include <assert.h>
#include "sparsemat.h"
#include "densemat.h"
#include "dsdpcg.h"

static void diagPrecond( vec *P, vec *x ) {
    // Diagonal preconditioning
    vec_vdiv(x, P); return;
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
        diagPrecond(cgSolver->vecPre, x);
    } else {
        cholPrecond(cgSolver->cholPre, x);
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
    
    cgSolver->vecPre  = NULL;
    cgSolver->cholPre = NULL;
    cgSolver->pType = CG_PRECOND_DIAG;
    cgSolver->reuse = 0;
    cgSolver->nused = 1024;
    
    cgSolver->dim     = 0;
    cgSolver->tol     = 1e-04;
    cgSolver->resinrm = 0.0;
    cgSolver->niter   = 0;
    cgSolver->maxiter = 20;
    cgSolver->nfailed = 0;
    
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
    
    cgSolver->M = NULL;
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
    
    cgSolver->pType = 0;
    cgSolver->cholPre = NULL;
    cgSolver->vecPre = NULL;
    
    cgSolver->tol     = 0.0;
    cgSolver->resinrm = 0.0;
    cgSolver->dim     = 0;
    cgSolver->niter   = 0;
    cgSolver->maxiter = 0;
    cgSolver->status  = 0;
    cgSolver->reuse   = 0;
    cgSolver->nused   = 0;
    cgSolver->nfailed = 0;
    
    return retcode;
}

extern DSDP_INT dsdpCGSetTol( CGSolver *cgSolver, double tol ) {
    assert( tol > 0 );
    cgSolver->tol = MAX(1e-20, tol);
    return DSDP_RETCODE_OK;
}

extern DSDP_INT dsdpCGSetMaxIter( CGSolver *cgSolver, DSDP_INT maxiter ) {
    assert( maxiter >= 1 );
    cgSolver->maxiter = MAX(maxiter, 1);
    return DSDP_RETCODE_OK;
}

extern DSDP_INT dsdpCGSetM( CGSolver *cgSolver, dsMat *M ) {
    assert( M->dim == cgSolver->dim );
    cgSolver->M = M;
    return DSDP_RETCODE_OK;
}

extern DSDP_INT dsdpCGSetDPre( CGSolver *cgSolver, vec *preCond ) {
    cgSolver->vecPre = preCond;
    return DSDP_RETCODE_OK;
}

extern DSDP_INT dsdpCGSetCholPre( CGSolver *cgSolver, dsMat *preCond ) {
    cgSolver->cholPre = preCond;
    return DSDP_RETCODE_OK;
}

extern DSDP_INT dsdpCGSetPType( CGSolver *cgSolver, DSDP_INT pType ) {
    assert( pType == CG_PRECOND_DIAG || pType == CG_PRECOND_CHOL);
    cgSolver->pType = pType;
    return DSDP_RETCODE_OK;
}

extern DSDP_INT dsdpCGGetStatus( CGSolver *cgSolver, DSDP_INT *status ) {
    *status = cgSolver->status;
    return DSDP_RETCODE_OK;
}

extern DSDP_INT dsdpCGGetSolStatistic( CGSolver *cgSolver, DSDP_INT *iter, double *resinorm ) {
    *iter = cgSolver->niter;
    *resinorm = cgSolver->resinrm;
    return DSDP_RETCODE_OK;
}

extern DSDP_INT dsdpCGSetPreReuse( CGSolver *cgSolver, DSDP_INT reuse ) {
    assert( reuse >= 0 );
    cgSolver->reuse = reuse;
    return DSDP_RETCODE_OK;
}

extern DSDP_INT dsdpCGResetPreReuse( CGSolver *cgSolver ) {
    cgSolver->nused = 0;
    return DSDP_RETCODE_OK;
}

extern DSDP_INT dsdpCGprepareP( CGSolver *cgSolver ) {
    
    /* TODO: Necessary to try diagonal scaling every iteration ?*/
    if (cgSolver->status == CG_STATUS_INDEFINITE) {
        denseMatResetFactor(cgSolver->M);
        denseMatFactorize(cgSolver->M);
        return DSDP_RETCODE_OK;
    }
    
    if (cgSolver->reuse == 0) {
        cgSolver->status = CG_STATUS_INDEFINITE;
    }
    
    // Prepare pre-conditioner in CG Solver
    if (cgSolver->pType == CG_PRECOND_DIAG) {
        return DSDP_RETCODE_OK;
    } else if (cgSolver->pType == CG_PRECOND_CHOL && cgSolver->nused > cgSolver->reuse){
        denseMatResetFactor(cgSolver->M);
        denseMatFactorize(cgSolver->M);
        dsdpCGResetPreReuse(cgSolver);
    } else {
        if (cgSolver->M->isillCond) {
            cgSolver->status = CG_STATUS_INDEFINITE;
            return DSDP_RETCODE_OK;
        }
        if (cgSolver->nused >= cgSolver->reuse) {
            dsdpCGSetPType(cgSolver, CG_PRECOND_DIAG);
        } else {
            cgSolver->nused += 1;
        }
        
    }
    return DSDP_RETCODE_OK;
}

/* Implement (pre-conditioned) conjugate gradient for Schur system */
extern DSDP_INT dsdpCGSolve( CGSolver *cgSolver, vec *b, vec *x0 ) {
    /* CG solver for schur system. Solve P^-1 M x = P^-1 b
       On exit b is over-written
     */
    DSDP_INT retcode = DSDP_RETCODE_OK;
    
    if (cgSolver->status == CG_STATUS_INDEFINITE) {
        denseArrSolveInp(cgSolver->M, 1, b->x);
        return retcode;
    }
    
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
        
        // Get initial residual
        vec_norm(r, &alpha);
        
        if (alpha < 1e-06 * tol) {
            cgSolver->status = CG_STATUS_SOLVED;
            vec_copy(x, b);
            return retcode;
        }
    } else {
        vec_reset(x);
        // r = b;
        vec_copy(b, r);
    }
    
    vec_norm(r, &alpha);
    
    if (alpha == 0.0) {
        cgSolver->status = CG_STATUS_SOLVED;
        return retcode;
    }
    
    alpha = MIN(100, alpha); tol = MAX(tol * alpha * 0.1, tol * 1e-01);
    
    // d = P \ r;
    vec_copy(r, d); preCond(cgSolver, d);
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
        
        // printf("%3d %10.3e %20.10e \n", iter + 1, alpha, resinorm);
        
        if (resinorm <= tol) {
            cgSolver->status = CG_STATUS_SOLVED;
            break;
        }
    }
    
    cgSolver->resinrm = resinorm;
    if (iter >= cgSolver->maxiter - 1 && resinorm > tol) {
        cgSolver->status = CG_STATUS_MAXITER;
    }
    
    vec_copy(x, b);

    return retcode;
}
