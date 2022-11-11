#include "linsys.h"
#include "qdldl.h"

typedef struct {
    
    int    n; ///< Dimension of the linear system
    
    int    *iWork; ///< Working array
    double *dWork; ///< Working array
    int    *P; ///< Symbolic tree
    
    int    *Lnz;
    int    *Lp; ///< CSC representation of L
    int    *Li;
    double *Lx;
    
} qdldl_linsys;

static int ldlCreate( void **pldl, int n ) {
    
    int retcode = RETCODE_OK;
    qdldl_linsys *qdldl = NULL;
    
    POTLP_INIT(qdldl, qdldl_linsys, 1);
    
    if ( !qdldl ) {
        retcode = RETCODE_FAILED;
        goto exit_cleanup;
    }
    
    qdldl->n = n;
    *pldl = qdldl;

exit_cleanup:
    return retcode;
}


static int ldlSymbolic( void *ldl, int *Ap, int *Ai ) {
    
    int retcode = RETCODE_OK;
    qdldl_linsys *qdldl = (qdldl_linsys *) ldl;
    
    POTLP_INIT(qdldl->P, int, qdldl->n);
    POTLP_INIT(qdldl->iWork, int, qdldl->n * 4);
    POTLP_INIT(qdldl->dWork, double, qdldl->n * 3);
    POTLP_INIT(qdldl->Lp, int, qdldl->n + 1);
    POTLP_INIT(qdldl->Lnz, int, qdldl->n);
    
    if ( !qdldl->P || !qdldl->Lp || !qdldl->iWork ||
         !qdldl->dWork || !qdldl->Lnz ) {
        retcode = RETCODE_FAILED;
        goto exit_cleanup;
    }
    
    int ldlret = QDLDL_etree(qdldl->n, Ap, Ai, qdldl->iWork, qdldl->Lnz, qdldl->P);
    
    if ( ldlret == -1 || ldlret == -2 ) {
        retcode = RETCODE_FAILED;
        goto exit_cleanup;
    }
    
    POTLP_INIT(qdldl->Li, int, ldlret);
    POTLP_INIT(qdldl->Lx, double, ldlret);
    
    if ( !qdldl->Li || !qdldl->Lx ) {
        retcode = RETCODE_FAILED;
        goto exit_cleanup;
    }
    
exit_cleanup:
    return retcode;
}

static int ldlNumeric( void *ldl, int *Ap, int *Ai, double *Ax ) {
    
    int retcode = RETCODE_OK;
    qdldl_linsys *qdldl = (qdldl_linsys *) ldl;
    
    int n = qdldl->n;
    
    double *D = qdldl->dWork;
    double *Dinv = qdldl->dWork + n;
    double *fWork = Dinv + n;
    
    int *bWork = qdldl->iWork;
    int *iWork = bWork + n;
    
    int nPos = QDLDL_factor(n, Ap, Ai, Ax, qdldl->Lp, qdldl->Li, qdldl->Lx,
                            D, Dinv, qdldl->Lnz, qdldl->P, bWork, iWork, fWork);
    
    if ( nPos == -1 ) {
        retcode = RETCODE_FAILED;
        goto exit_cleanup;
    }
    
exit_cleanup:
    return retcode;
}

static void ldlSolve( void *ldl, double *bx ) {
    
    qdldl_linsys *qdldl = (qdldl_linsys *) ldl;
    QDLDL_solve(qdldl->n, qdldl->Lp, qdldl->Li, qdldl->Lx,
                qdldl->dWork + qdldl->n, bx);
    
    return;
}

static void ldlDestroy( void **pldl ) {
    
    if ( !pldl ) {
        return;
    }
    
    qdldl_linsys *qdldl = (qdldl_linsys *) (*pldl);
    
    if ( qdldl ) {
        
        POTLP_FREE(qdldl->P);
        POTLP_FREE(qdldl->iWork);
        POTLP_FREE(qdldl->dWork);
        POTLP_FREE(qdldl->Lnz);
        POTLP_FREE(qdldl->Lp);
        POTLP_FREE(qdldl->Li);
        POTLP_FREE(qdldl->Lx);
        
        POTLP_ZERO(qdldl, qdldl_linsys, 1);
    }
    
    POTLP_FREE(*pldl);
    
    return;
}

extern pot_int potLinsysCreate( pot_linsys **ppotLinsys ) {
    
    pot_int retcode = RETCODE_OK;
    
    if ( !ppotLinsys ) {
        retcode = RETCODE_FAILED;
        goto exit_cleanup;
    }
    
    pot_linsys *potLinsys = NULL;
    POTLP_INIT(potLinsys, pot_linsys, 1);
    
    if ( !potLinsys ) {
        retcode = RETCODE_FAILED;
        goto exit_cleanup;
    }
    
    POTLP_ZERO(potLinsys, pot_linsys, 1);
    
    potLinsys->LCreate = ldlCreate;
    potLinsys->LDestroy = ldlDestroy;
    potLinsys->LSFac = ldlSymbolic;
    potLinsys->LNFac = ldlNumeric;
    potLinsys->LSolve = ldlSolve;
    
    *ppotLinsys = potLinsys;
    
exit_cleanup:
    return retcode;
}

extern pot_int potLinsysInit( pot_linsys *potLinsys, pot_int nCol ) {
    
    pot_int retcode = RETCODE_OK;
    
    potLinsys->nCol = nCol;
    retcode = potLinsys->LCreate(&potLinsys->solver, nCol);
    
    if ( retcode != RETCODE_OK ) {
        retcode = RETCODE_FAILED;
        goto exit_cleanup;
    }
    
exit_cleanup:
    return retcode;
}

extern pot_int potLinsysSymFactorize( pot_linsys *potLinsys, pot_int *colMatBeg, pot_int *colMatIdx ) {
    
    pot_int retcode = RETCODE_OK;
    retcode = potLinsys->LSFac(potLinsys->solver, colMatBeg, colMatIdx);
    
    if ( retcode != RETCODE_OK ) {
        retcode = RETCODE_FAILED;
        goto exit_cleanup;
    }
    
exit_cleanup:
    return retcode;
}

extern pot_int potLinsysNumFactorize( pot_linsys *potLinsys, int *colMatBeg, int *colMatIdx, double *colMatElem ) {
    
    pot_int retcode = RETCODE_OK;
    retcode = potLinsys->LNFac(potLinsys->solver, colMatBeg, colMatIdx, colMatElem);
    
    if ( retcode != RETCODE_OK ) {
        retcode = RETCODE_FAILED;
        goto exit_cleanup;
    }
    
exit_cleanup:
    return retcode;
}

extern void potLinsysSolve( pot_linsys *potLinsys, double *rhsVec, double *solVec ) {
    
    if ( solVec ) {
        POTLP_MEMCPY(solVec, rhsVec, double, potLinsys->nCol);
        potLinsys->LSolve(potLinsys->solver, solVec);
    } else {
        potLinsys->LSolve(potLinsys->solver, rhsVec);
    }
    
    return;
}

extern void potLinsysClear( pot_linsys *potLinsys ) {
    
    if ( !potLinsys ) {
        return;
    }
    
    potLinsys->LDestroy(&potLinsys->solver);
    POTLP_ZERO(potLinsys, pot_linsys, 1);
    return;
}

extern void potLinsysDestroy( pot_linsys **ppotLinsys ) {
    
    if ( !ppotLinsys ) {
        return;
    }
    
    potLinsysClear(*ppotLinsys);
    POTLP_FREE(*ppotLinsys);
    
    return;
}
