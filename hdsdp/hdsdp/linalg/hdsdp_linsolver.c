#include "interface/hdsdp_utils.h"

#include "linalg/def_hdsdp_linsolver.h"
#include "linalg/hdsdp_linsolver.h"

/* Sparse direct solver interface */
static hdsdp_retcode pardisoLinSolverCreate( void **pchol, int nCol ) {
    
    hdsdp_retcode retcode = HDSDP_RETCODE_OK;
    
    HDSDP_NULLCHECK(pchol);
    pardiso_linsys *pds = NULL;
    HDSDP_INIT(pds, pardiso_linsys, 1);
    HDSDP_MEMCHECK(pds);
    
    pds->nCol = nCol;
    HDSDP_ZERO(pds->pt, void *, 64);
    HDSDP_ZERO(pds->iparm, int, 64);
    
    int mType = PARDISO_SYM_POSDEFINITE;
    pardisoinit(pds->pt, &mType, pds->iparm);
    
    /* Non-default setting */
    set_pardiso_param(pds->iparm, PARDISO_PARAM_NONDEFAULT, 1);
    /* Symbolic ordering */
    set_pardiso_param(pds->iparm, PARDISO_PARAM_SYMBOLIC, PARDISO_PARAM_SYMBOLIC_MMD);
    /* Pivoting perturbation */
    set_pardiso_param(pds->iparm, PARDISO_PARAM_PERTURBATION, 8);
    /* Solve in-place*/
    set_pardiso_param(pds->iparm, PARDISO_PARAM_INPLACE, 1);
    /* Use 0-based index*/
    set_pardiso_param(pds->iparm, PARDISO_PARAM_INDEX, PARDISO_PARAM_INDEX_C);
    /* Get diagonal elements */
    set_pardiso_param(pds->iparm, PARDISO_PARAM_DIAGONAL, PARDISO_PARAM_DIAGONAL_ON);
    
    *pchol = (void *) pds;
    
exit_cleanup:
    
    return retcode;
}

/* TODO: Enable reproducibility using pardiso CNR mode */
static void pardisoLinSolverSetThreads( void *chol, void *pThreads ) {
    
    pardiso_linsys *pds = (pardiso_linsys *) chol;
    
    int nThreads = *((int *) pThreads);
    set_pardiso_param(pds->iparm, PARDISO_PARAM_THREADS, nThreads);
    
    return;
}

#define SHOW_ORDERING(ord) printf("Using ordering %s. \n", ord)
static hdsdp_retcode pardisoLinSolverSymbolic( void *chol, int *colMatBeg, int *colMatIdx ) {
    
    hdsdp_retcode retcode = HDSDP_RETCODE_OK;
    
    pardiso_linsys *pds = (pardiso_linsys *) chol;
    
    pds->colMatBeg = colMatBeg;
    pds->colMatIdx = colMatIdx;
    
    HDSDP_INIT(pds->dWork, double, pds->nCol);
    HDSDP_MEMCHECK(pds->dWork);
    
    int maxfct = 1, mnum = 1, mtype = PARDISO_SYM_POSDEFINITE, phase = PARDISO_PHASE_SYM;
    int idummy = 0, msg = 0, pdsret = PARDISO_RET_OK;
    
    /* Use amd first */
    int amdFactorNnz = 0;
    set_pardiso_param(pds->iparm, PARDISO_PARAM_SYMBOLIC, PARDISO_PARAM_SYMBOLIC_MMD);
    set_pardiso_param(pds->iparm, PARDISO_PARAM_FACNNZ, -1);
    
    pardiso(pds->pt, &maxfct, &mnum, &mtype, &phase,
            &pds->nCol, NULL, colMatBeg, colMatIdx, &idummy, &idummy,
            pds->iparm, &msg, NULL, NULL, &pdsret);
    
    if ( pdsret != PARDISO_RET_OK ) {
        retcode = HDSDP_RETCODE_FAILED;
        goto exit_cleanup;
    }
    
    amdFactorNnz = get_pardiso_output(pds->iparm, PARDISO_PARAM_FACNNZ);
    
    /* Try nested dissection */
    int ndFactorNnz = 0;
    set_pardiso_param(pds->iparm, PARDISO_PARAM_SYMBOLIC, PARDISO_PARAM_SYMBOLIC_ND);
    set_pardiso_param(pds->iparm, PARDISO_PARAM_FACNNZ, -1);
    
    pardiso(pds->pt, &maxfct, &mnum, &mtype, &phase,
            &pds->nCol, NULL, colMatBeg, colMatIdx, &idummy, &idummy,
            pds->iparm, &msg, NULL, NULL, &pdsret);
    
    if ( pdsret != PARDISO_RET_OK ) {
        retcode = HDSDP_RETCODE_FAILED;
        goto exit_cleanup;
    }
    
    ndFactorNnz = get_pardiso_output(pds->iparm, PARDISO_PARAM_FACNNZ);
    
    if ( ndFactorNnz > amdFactorNnz ) {
        /* Finally decide to use amd */
        set_pardiso_param(pds->iparm, PARDISO_PARAM_SYMBOLIC, PARDISO_PARAM_SYMBOLIC_MMD);
        set_pardiso_param(pds->iparm, PARDISO_PARAM_FACNNZ, -1);
        pardiso(pds->pt, &maxfct, &mnum, &mtype, &phase,
                &pds->nCol, NULL, colMatBeg, colMatIdx, &idummy, &idummy,
                pds->iparm, &msg, NULL, NULL, &pdsret);
#ifdef HDSDP_LINSYS_DEBUG
        printf("Using AMD ordering. nFactorNnz = %d \n", amdFactorNnz);
#endif
    } else {
#ifdef HDSDP_LINSYS_DEBUG
        printf("Using ND ordering. nFactorNnz = %d \n", ndFactorNnz);
#endif
    }
    
exit_cleanup:
    
    return retcode;
}

/* Stable Cholesky solve */
static hdsdp_retcode pardisoLinSolverStableNumeric( void *chol, int *colMatBeg, int *colMatIdx, double *colMatElem ) {
    
    hdsdp_retcode retcode = HDSDP_RETCODE_OK;
    
    pardiso_linsys *pds = (pardiso_linsys *) chol;
    
    pds->colMatBeg = colMatBeg;
    pds->colMatIdx = colMatIdx;
    pds->colMatElem = colMatElem;
    
    int maxfct = 1, mnum = 1, mtype = PARDISO_SYM_POSDEFINITE, phase = PARDISO_PHASE_FAC;
    int idummy = 0, msg = 0, pdsret = PARDISO_RET_OK;
    
    set_pardiso_param(pds->iparm, PARDISO_PARAM_SCALING, 1);
    set_pardiso_param(pds->iparm, PARDISO_PARAM_MATCHING, 1);
    
    pardiso(pds->pt, &maxfct, &mnum, &mtype, &phase,
            &pds->nCol, colMatElem, colMatBeg, colMatIdx, &idummy, &idummy,
            pds->iparm, &msg, NULL, NULL, &pdsret);
    
    if ( pdsret != PARDISO_RET_OK ) {
        retcode = HDSDP_RETCODE_FAILED;
        goto exit_cleanup;
    }

exit_cleanup:
    return retcode;
}

static hdsdp_retcode pardisoLinSolverNumeric( void *chol, int *colMatBeg, int *colMatIdx, double *colMatElem ) {
    
    hdsdp_retcode retcode = HDSDP_RETCODE_OK;
    
    pardiso_linsys *pds = (pardiso_linsys *) chol;
    
    pds->colMatBeg = colMatBeg;
    pds->colMatIdx = colMatIdx;
    pds->colMatElem = colMatElem;
    
    int maxfct = 1, mnum = 1, mtype = PARDISO_SYM_POSDEFINITE, phase = PARDISO_PHASE_FAC;
    int idummy = 0, msg = 0, pdsret = PARDISO_RET_OK;
    
    pardiso(pds->pt, &maxfct, &mnum, &mtype, &phase,
            &pds->nCol, colMatElem, colMatBeg, colMatIdx, &idummy, &idummy,
            pds->iparm, &msg, NULL, NULL, &pdsret);
    
    if ( pdsret != PARDISO_RET_OK ) {
        retcode = HDSDP_RETCODE_FAILED;
        goto exit_cleanup;
    }

exit_cleanup:
    return retcode;
}

static hdsdp_retcode pardisoLinSolverPsdCheck( void *chol, int *colMatBeg, int *colMatIdx, double *colMatElem, int *isPsd ) {
    
    hdsdp_retcode retcode = HDSDP_RETCODE_OK;
    
    pardiso_linsys *pds = (pardiso_linsys *) chol;
    
    pds->colMatBeg = colMatBeg;
    pds->colMatIdx = colMatIdx;
    pds->colMatElem = colMatElem;
    
    int maxfct = 1, mnum = 1, mtype = PARDISO_SYM_POSDEFINITE, phase = PARDISO_PHASE_FAC;
    int idummy = 0, msg = 0, pdsret = PARDISO_RET_OK;
    
    pardiso(pds->pt, &maxfct, &mnum, &mtype, &phase,
            &pds->nCol, colMatElem, colMatBeg, colMatIdx, &idummy, &idummy,
            pds->iparm, &msg, NULL, NULL, &pdsret);
    
    if ( pdsret == PARDISO_RET_OK ) {
        *isPsd = 1;
    } else if ( pdsret == PARDISO_RET_INDEFINITE ) {
        *isPsd = 0;
    } else {
        retcode = HDSDP_RETCODE_FAILED;
        goto exit_cleanup;
    }

exit_cleanup:
    return retcode;
}

/* If nRhs > 1, then the solution must be overwritten */
static void pardisoLinSolverForwardN( void *chol, int nRhs, double *rhsVec, double *solVec ) {
    
    pardiso_linsys *pds = (pardiso_linsys *) chol;
    
    int maxfct = 1, mnum = 1, mtype = PARDISO_SYM_POSDEFINITE, phase = PARDISO_PHASE_FORWARD;
    int idummy = 0, msg = 0, pdsret = PARDISO_RET_OK;
    
    /* Do not overwrite */
    if ( solVec ) {
        set_pardiso_param(pds->iparm, PARDISO_PARAM_INPLACE, 0);
        pardiso(pds->pt, &maxfct, &mnum, &mtype, &phase,
                &pds->nCol, pds->colMatElem, pds->colMatBeg, pds->colMatIdx, &idummy, &nRhs,
                pds->iparm, &msg, rhsVec, solVec, &pdsret);
    } else {
        
        assert( nRhs == 1 );
        set_pardiso_param(pds->iparm, PARDISO_PARAM_INPLACE, 1);
        pardiso(pds->pt, &maxfct, &mnum, &mtype, &phase,
                &pds->nCol, pds->colMatElem, pds->colMatBeg, pds->colMatIdx, &idummy, &nRhs,
                pds->iparm, &msg, rhsVec, pds->dWork, &pdsret);
    }
    
    return;
}

static void pardisoLinSolverBackwardN( void *chol, int nRhs, double *rhsVec, double *solVec ) {
    
    pardiso_linsys *pds = (pardiso_linsys *) chol;
    
    int maxfct = 1, mnum = 1, mtype = PARDISO_SYM_POSDEFINITE, phase = PARDISO_PHASE_BACKWARD;
    int idummy = 0, msg = 0, pdsret = PARDISO_RET_OK;
    
    /* Do not overwrite */
    if ( solVec ) {
        set_pardiso_param(pds->iparm, PARDISO_PARAM_INPLACE, 0);
        pardiso(pds->pt, &maxfct, &mnum, &mtype, &phase,
                &pds->nCol, pds->colMatElem, pds->colMatBeg, pds->colMatIdx, &idummy, &nRhs,
                pds->iparm, &msg, rhsVec, solVec, &pdsret);
    } else {
        
        assert( nRhs == 1 );
        set_pardiso_param(pds->iparm, PARDISO_PARAM_INPLACE, 1);
        pardiso(pds->pt, &maxfct, &mnum, &mtype, &phase,
                &pds->nCol, pds->colMatElem, pds->colMatBeg, pds->colMatIdx, &idummy, &nRhs,
                pds->iparm, &msg, rhsVec, pds->dWork, &pdsret);
    }
    
    return;
}

static void pardisoLinSolverSolveN( void *chol, int nRhs, double *rhsVec, double *solVec ) {
    
    pardiso_linsys *pds = (pardiso_linsys *) chol;
    
    int maxfct = 1, mnum = 1, mtype = PARDISO_SYM_POSDEFINITE, phase = PARDISO_PHASE_SOLVE;
    int idummy = 0, msg = 0, pdsret = PARDISO_RET_OK;
    
    /* Do not overwrite */
    if ( solVec ) {
        set_pardiso_param(pds->iparm, PARDISO_PARAM_INPLACE, 0);
        pardiso(pds->pt, &maxfct, &mnum, &mtype, &phase,
                &pds->nCol, pds->colMatElem, pds->colMatBeg, pds->colMatIdx, &idummy, &nRhs,
                pds->iparm, &msg, rhsVec, solVec, &pdsret);
    } else {
        
        assert( nRhs == 1 );
        set_pardiso_param(pds->iparm, PARDISO_PARAM_INPLACE, 1);
        pardiso(pds->pt, &maxfct, &mnum, &mtype, &phase,
                &pds->nCol, pds->colMatElem, pds->colMatBeg, pds->colMatIdx, &idummy, &nRhs,
                pds->iparm, &msg, rhsVec, pds->dWork, &pdsret);
    }
    
    return;
}

static hdsdp_retcode pardisoLinSolverGetDiag( void *chol, double *diagVec ) {
    
    hdsdp_retcode retcode = HDSDP_RETCODE_OK;
    pardiso_linsys *pds = (pardiso_linsys *) chol;
    
    int one = 1, pdsret = PARDISO_RET_OK;
    pardiso_getdiag((const void **) pds->pt, (void *) diagVec, (void *) pds->dWork, &one, &pdsret);
    
    if ( pdsret != PARDISO_RET_OK ) {
        retcode = HDSDP_RETCODE_FAILED;
        goto exit_cleanup;
    }
    
exit_cleanup:
    return retcode;
}

static void pardisoLinSolverInvert( void *chol, double *dFullMatrix, double *dAuxiMatrix ) {
    
    pardiso_linsys *pds = (pardiso_linsys *) chol;
    HDSDP_ZERO(dAuxiMatrix, double, pds->nCol * pds->nCol);
    
    double *pElem = dAuxiMatrix;
    for ( int iRow = 0; iRow < pds->nCol; ++iRow ) {
        *pElem = 1.0;
        pElem += iRow + pds->nCol;
    }
    
    pardisoLinSolverSolveN(chol, pds->nCol, dAuxiMatrix, dFullMatrix);
    
    return;
}

static void pardisoLinSolverClear( void *chol ) {
    
    if ( !chol ) {
        return;
    }
    
    pardiso_linsys *pds = (pardiso_linsys *) chol;
    
    int maxfct = 1, mnum = 1, mtype = PARDISO_SYM_POSDEFINITE, phase = PARDISO_PHASE_FREE;
    int idummy = 0, msg = 0, pdsret = PARDISO_RET_OK;
    
    pardiso(pds->pt, &maxfct, &mnum, &mtype, &phase,
            &pds->nCol, pds->colMatElem, pds->colMatBeg, pds->colMatIdx, &idummy, NULL,
            pds->iparm, &msg, NULL, pds->dWork, &pdsret);
    
    HDSDP_FREE(pds->dWork);
    HDSDP_ZERO(pds, pardiso_linsys, 1);
    
    return;
}

static void pardisoLinSolverDestroy( void **pchol ) {
    
    if ( !pchol ) {
        return;
    }
    
    pardisoLinSolverClear(*pchol);
    HDSDP_FREE(*pchol);
    
    return;
}

/* Dense direct solver interface */
static hdsdp_retcode lapackLinSolverCreate( void **pchol, int nCol ) {
    
    hdsdp_retcode retcode = HDSDP_RETCODE_OK;
    
    HDSDP_NULLCHECK(pchol);
    lapack_linsys *lap = NULL;
    HDSDP_INIT(lap, lapack_linsys, 1);
    HDSDP_MEMCHECK(lap);
    
    lap->nCol = nCol;
    HDSDP_ZERO(lap, lapack_linsys, 1);
    
    HDSDP_INIT(lap->dFullMatElem, double, nCol * nCol);
    HDSDP_MEMCHECK(lap->dFullMatElem);
    
exit_cleanup:
    return retcode;
}

static void lapackLinSolverSetParamsDummy( void *chol, void *param ) {
    
    (void) param;
    
    return;
}

static hdsdp_retcode lapackLinSolverSymbolic( void *chol, int *dummy1, int *dummy2 ) {
    
    (void) dummy1;
    (void) dummy2;
    
    return HDSDP_RETCODE_OK;
}


static hdsdp_retcode lapackLinSolverNumeric( void *chol, int *dummy1, int *dummy2, double *dFullMatrix ) {
    
    hdsdp_retcode retcode = HDSDP_RETCODE_OK;
    
    int info = 0;
    
    
    return retcode;
}
