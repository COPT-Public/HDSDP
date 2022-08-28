#include <stdlib.h>
#include "cg.h"
#include "dsdplapack.h"
#include "schurmat.h"
#include "vec.h"


static void diagPrecond( vec *P, vec *x ) {
    // Diagonal preconditioning
    vec_vdiv(x, P); return;
}

static void cholPrecond( schurMat *P, vec *x, vec *aux ) {
    // Cholesky preconditioning
    assert( P->m == x->dim );
    if (!P->isFactorized) {
        schurMatFactorize(P);
    }   
    schurMatSolve(P, 1, x->x, aux->x);
    return;
}

static void preCond( CGSolver *cgSolver, vec *x ) {
    
    if (cgSolver->pType == CG_PRECOND_DIAG) {
        diagPrecond(cgSolver->vecPre, x);
    } else {
        cholPrecond(cgSolver->cholPre, x, cgSolver->aux);
    }
    return;
}

extern void dsdpCGinit( CGSolver *cgSolver ) {
        
    cgSolver->M     = NULL; cgSolver->r   = NULL;
    cgSolver->rnew  = NULL; cgSolver->d   = NULL;
    cgSolver->Pinvr = NULL; cgSolver->x   = NULL;
    cgSolver->Md    = NULL; cgSolver->aux = NULL;
    
    cgSolver->vecPre  = NULL;
    cgSolver->cholPre = NULL;
    cgSolver->pType = CG_PRECOND_DIAG;
    cgSolver->reuse = 0; cgSolver->nused = 1024;
    
    cgSolver->dim = 0; cgSolver->tol = 0.0;
    cgSolver->resinrm = 0.0; cgSolver->niter = 0;
    cgSolver->maxiter = 0; cgSolver->nfailed = 0;
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
    cgSolver->aux   = (vec *) calloc(1, sizeof(vec));
    
    vec_init(cgSolver->r);     vec_alloc(cgSolver->r, m);
    vec_init(cgSolver->rnew);  vec_alloc(cgSolver->rnew, m);
    vec_init(cgSolver->d);     vec_alloc(cgSolver->d, m);
    vec_init(cgSolver->Pinvr); vec_alloc(cgSolver->Pinvr, m);
    vec_init(cgSolver->Md);    vec_alloc(cgSolver->Md, m);
    vec_init(cgSolver->x);     vec_alloc(cgSolver->x, m);
    vec_init(cgSolver->aux);   vec_alloc(cgSolver->aux, m);

    return retcode;
}

extern void dsdpCGFree( CGSolver *cgSolver ) {
    
    cgSolver->M = NULL;
    vec_free(cgSolver->r);     DSDP_FREE(cgSolver->r);
    vec_free(cgSolver->rnew);  DSDP_FREE(cgSolver->rnew);
    vec_free(cgSolver->d);     DSDP_FREE(cgSolver->d);
    vec_free(cgSolver->Pinvr); DSDP_FREE(cgSolver->Pinvr);
    vec_free(cgSolver->Md);    DSDP_FREE(cgSolver->Md);
    vec_free(cgSolver->x);     DSDP_FREE(cgSolver->x);
    vec_free(cgSolver->aux);   DSDP_FREE(cgSolver->aux);
    
    cgSolver->pType = 0;
    cgSolver->cholPre = NULL;
    cgSolver->vecPre = NULL;
    
    cgSolver->tol     = 0.0; cgSolver->resinrm = 0.0;
    cgSolver->dim     = 0;   cgSolver->niter   = 0;
    cgSolver->maxiter = 0;   cgSolver->status  = 0;
    cgSolver->reuse   = 0;   cgSolver->nused   = 0;
    cgSolver->nfailed = 0;
}

extern void dsdpCGSetTol( CGSolver *cgSolver, double tol ) {
    
    cgSolver->tol = MAX(1e-20, tol);
}

extern void dsdpCGSetMaxIter( CGSolver *cgSolver, DSDP_INT maxiter ) {
    
    cgSolver->maxiter = MAX(maxiter, 1);
}

extern void dsdpCGSetM( CGSolver *cgSolver, schurMat *M ) {
    
    cgSolver->M = M;
}

extern void dsdpCGSetDPre( CGSolver *cgSolver, vec *preCond ) {
    
    cgSolver->vecPre = preCond;
}

extern void dsdpCGSetCholPre( CGSolver *cgSolver, schurMat *preCond ) {
    
    cgSolver->cholPre = preCond;
}

extern void dsdpCGSetPType( CGSolver *cgSolver, DSDP_INT pType ) {
    
    cgSolver->pType = pType;
}

extern void dsdpCGGetStatus( CGSolver *cgSolver, DSDP_INT *status ) {
    
    *status = cgSolver->status;
}

extern void dsdpCGGetSolStatistic( CGSolver *cgSolver, DSDP_INT *iter, double *resinorm ) {
    
    *iter = cgSolver->niter;
    *resinorm = cgSolver->resinrm;
}

extern void dsdpCGSetPreReuse( CGSolver *cgSolver, DSDP_INT reuse ) {
    
    cgSolver->reuse = reuse;
}

extern void dsdpCGResetPreReuse( CGSolver *cgSolver ) {
    
    cgSolver->nused = 0;
}

extern void dsdpCGprepareP( CGSolver *cgSolver ) {
    
    if (cgSolver->status == CG_STATUS_INDEFINITE) {
        schurMatReset(cgSolver->M, FALSE); // Clean factor
        schurMatFactorize(cgSolver->M);
        return;
    }
    
    if (cgSolver->reuse == 0) {
        cgSolver->status = CG_STATUS_INDEFINITE;
    }
    
    // Prepare pre-conditioner in CG Solver
    if (cgSolver->pType == CG_PRECOND_DIAG) {
        return;
    } else if (cgSolver->pType == CG_PRECOND_CHOL && cgSolver->nused > cgSolver->reuse){
        schurMatReset(cgSolver->M, FALSE);
        schurMatFactorize(cgSolver->M);
        if (cgSolver->M->denseM->isillCond) {
            cgSolver->M->isillCond = TRUE;
        }
        dsdpCGResetPreReuse(cgSolver);
    } else {
        if (cgSolver->M->isillCond || cgSolver->M->denseM->isillCond) {
            cgSolver->status = CG_STATUS_INDEFINITE;
            return;
        }
        if (cgSolver->nused >= cgSolver->reuse) {
            dsdpCGSetPType(cgSolver, CG_PRECOND_DIAG);
        } else {
            cgSolver->nused += 1;
        }
        
    }
}

extern void dsdpCGStoreRHS( CGSolver *cgSolver, vec *bin ) {
    
    vec_copy(bin, cgSolver->aux);
}

extern void dsdpCGRestoreRHS( CGSolver *cgSolver, vec *bout ) {
    
    vec_copy(cgSolver->aux, bout);
}

/* Implement (pre-conditioned) conjugate gradient for Schur system */
extern DSDP_INT dsdpCGSolve( CGSolver *cgSolver, vec *b, vec *x0 ) {
    // CG solver for schur system. Solve P^-1 M x = P^-1 b. On exit b is over-written
    DSDP_INT retcode = DSDP_RETCODE_OK;
    
    if (cgSolver->M->stype == SCHUR_TYPE_SPARSE) {
        schurMatFactorize(cgSolver->M);
        schurMatSolve(cgSolver->M, 1, b->x, cgSolver->r->x);
        cgSolver->status = CG_STATUS_SOLVED;
        return retcode;
    }
    
    
    if (cgSolver->status == CG_STATUS_INDEFINITE) {
        schurMatSolve(cgSolver->M, 1, b->x, cgSolver->r->x);
        return retcode;
    }

    /* Initialize
     x = x0;
     r = b - M * x;
     d = P \ r;
     Md = M * d;
     Pinvr = P \ r;
    */
    
    schurMat *M = cgSolver->M;
    vec *r = cgSolver->r, *rnew = cgSolver->rnew, *x = cgSolver->x,
        *d = cgSolver->d, *Pinvr = cgSolver->Pinvr, *Md = cgSolver->Md;
    
    double tol = cgSolver->tol, dTMd = 0.0, rTPinvr = 0.0,
           alpha = 0.0, beta = 0.0, resinorm = DSDP_INFINITY;
    
    DSDP_INT iter = 0; vec_reset(r); cgSolver->status = CG_STATUS_UNKNOWN;
    
    if (x0) {
        // Warm start CG solver
        vec_copy(x0, x); vec_copy(b, r); schurMatMx(M, x, cgSolver->d);
        vec_axpy(-1.0, cgSolver->d, r); vec_norm(r, &alpha); // Get initial residual
        if (alpha < 1e-06 * tol) {
            cgSolver->status = CG_STATUS_SOLVED; vec_copy(x, b);
            return retcode;
        }
    } else {
        vec_reset(x); vec_copy(b, r);
    }
    
    vec_norm(r, &alpha); vec_norm(b, &resinorm);
    
    if (resinorm < alpha) {
        vec_reset(x); vec_copy(b, r);
        alpha = resinorm;
    }
    
    if (alpha == 0.0) {
        cgSolver->status = CG_STATUS_SOLVED;
        return retcode;
    }
    
    alpha = MIN(100, alpha); tol = MAX(tol * alpha * 0.1, tol * 1e-01);
    
    vec_copy(r, d); preCond(cgSolver, d); // d = P \ r;
    
    schurMatMx(M, d, Md); // Md = M * d;
    vec_copy(d, Pinvr); // Pinvr = P \ r;
    
    // Start CG
    for (iter = 0; iter < cgSolver->maxiter; ++iter) {
        vec_dot(r, Pinvr, &rTPinvr); // alpha = r' * Pinvr / (d' * Md);s
        vec_dot(d, Md, &dTMd); alpha = rTPinvr / dTMd;
        vec_axpy(alpha, d, x); // x = x + alpha * d;
        vec_zaxpby(rnew, 1.0, r, -alpha, Md); // rnew = r - alpha * Md
        vec_copy(rnew, Pinvr); // Pinvr = P \ rnew;
        preCond(cgSolver, Pinvr);
        vec_dot(rnew, Pinvr, &beta); // beta = rnew' * Pinvr / rTPinvr;
        beta = beta / rTPinvr;
        vec_axpby(1.0, Pinvr, beta, d); // d = Pinvr + beta * d;
        schurMatMx(M, d, Md); // Md = M * d;
        vec_copy(rnew, r); // r = rnew
        vec_norm(r, &resinorm); // norm(r) // printf("%3d %10.3e %20.10e \n", iter + 1, alpha, resinorm);
        if (resinorm <= tol) {
            cgSolver->status = CG_STATUS_SOLVED; break;
        }
    }
    cgSolver->resinrm = resinorm;
    if (iter >= cgSolver->maxiter - 1 && resinorm > tol) {
        cgSolver->status = CG_STATUS_MAXITER;
    }
    
    vec_copy(x, b);
    return retcode;
}
