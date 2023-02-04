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

static hdsdp_retcode pardisoLinSolverSolveN( void *chol, int nRhs, double *rhsVec, double *solVec ) {
    
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
    
    return HDSDP_RETCODE_OK;
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
    lapack_linsys *lap = (lapack_linsys *) chol;
    
    HDSDP_MEMCPY(lap->dFullMatElem, dFullMatrix, double, lap->nCol * lap->nCol);
    
    int info = LAPACK_RET_OK;
    char uplolow = LAPACK_UPLOW_LOW;
    dpotrf(&uplolow, &lap->nCol, lap->dFullMatElem, &lap->nCol, &info);
    
    if ( info != LAPACK_RET_OK ) {
        retcode = HDSDP_RETCODE_FAILED;
        goto exit_cleanup;
    }
    
exit_cleanup:
    return retcode;
}

static hdsdp_retcode lapackLinSolverPsdCheck( void *chol, int *dummy1, int *dummy2, double *dFullMatrix, int *isPsd ) {
    
    hdsdp_retcode retcode = HDSDP_RETCODE_OK;
    lapack_linsys *lap = (lapack_linsys *) chol;
    
    HDSDP_MEMCPY(lap->dFullMatElem, dFullMatrix, double, lap->nCol * lap->nCol);
    
    int info = LAPACK_RET_OK;
    char uplolow = LAPACK_UPLOW_LOW;
    dpotrf(&uplolow, &lap->nCol, lap->dFullMatElem, &lap->nCol, &info);
    
    if ( info == LAPACK_RET_OK ) {
        *isPsd = 1;
    } else if ( info > 0 ) {
        *isPsd = 0;
    } else {
        retcode = HDSDP_RETCODE_FAILED;
        goto exit_cleanup;
    }
 
exit_cleanup:
    return retcode;
}

static void lapackLinSolverForwardN( void *chol, int nRhs, double *rhsVec, double *solVec ) {
    
    lapack_linsys *lap = (lapack_linsys *) chol;
    
    int info = LAPACK_RET_OK;
    char uplolow = LAPACK_UPLOW_LOW, sideleft = LAPACK_SIDE_LEFT, notrans = LAPACK_NOTRANS, nonunit = LAPACK_DIAG_NONUNIT;
    double done = 1.0;
    
    if ( solVec ) {
        HDSDP_MEMCPY(solVec, rhsVec, double, nRhs * lap->nCol);
        dtrsm(&sideleft, &uplolow, &notrans, &nonunit, &lap->nCol, &nRhs,
              &done, lap->dFullMatElem, &lap->nCol, solVec, &lap->nCol);
    } else {
        dtrsm(&sideleft, &uplolow, &notrans, &nonunit, &lap->nCol, &nRhs,
              &done, lap->dFullMatElem, &lap->nCol, rhsVec, &lap->nCol);
    }
 
    return;
}

static void lapackLinSolverBackwardN( void *chol, int nRhs, double *rhsVec, double *solVec ) {
    
    lapack_linsys *lap = (lapack_linsys *) chol;
    
    int info = LAPACK_RET_OK;
    char uplolow = LAPACK_UPLOW_LOW, sideleft = LAPACK_SIDE_LEFT, trans = LAPACK_TRANS, nonunit = LAPACK_DIAG_NONUNIT;
    double done = 1.0;
    
    if ( solVec ) {
        HDSDP_MEMCPY(solVec, rhsVec, double, nRhs * lap->nCol);
        dtrsm(&sideleft, &uplolow, &trans, &nonunit, &lap->nCol, &nRhs,
              &done, lap->dFullMatElem, &lap->nCol, solVec, &lap->nCol);
    } else {
        dtrsm(&sideleft, &uplolow, &trans, &nonunit, &lap->nCol, &nRhs,
              &done, lap->dFullMatElem, &lap->nCol, rhsVec, &lap->nCol);
    }
 
    return;
}

static hdsdp_retcode lapackLinSolverSolveN( void *chol, int nRhs, double *rhsVec, double *solVec ) {
    
    lapack_linsys *lap = (lapack_linsys *) chol;
    
    int info = LAPACK_RET_OK;
    char uplolow = LAPACK_UPLOW_LOW;
    
    if ( solVec ) {
        HDSDP_MEMCPY(solVec, rhsVec, double, nRhs * lap->nCol);
        dpotrs(&uplolow, &lap->nCol, &nRhs, lap->dFullMatElem,
               &lap->nCol, solVec, &lap->nCol, &info);
    } else {
        dpotrs(&uplolow, &lap->nCol, &nRhs, lap->dFullMatElem,
               &lap->nCol, rhsVec, &lap->nCol, &info);
    }
    
    assert( info == LAPACK_RET_OK );
    
    return HDSDP_RETCODE_OK;
}

static hdsdp_retcode lapackLinSolverGetDiag( void *chol, double *diagVec ) {
    
    lapack_linsys *lap = (lapack_linsys *) chol;
    
    for ( int iRow = 0; iRow < lap->nCol; ++iRow ) {
        diagVec[iRow] = lap->dFullMatElem[iRow * lap->nCol + iRow];
    }
    
    return HDSDP_RETCODE_OK;
}

static void lapackLinSolverInvert( void *chol, double *dFullMatrix, double *dAuxiMatrix ) {
    
    (void) dAuxiMatrix;
    lapack_linsys *lap = (lapack_linsys *) chol;
    HDSDP_MEMCPY(dFullMatrix, lap->dFullMatElem, double, lap->nCol * lap->nCol);
    
    int info = LAPACK_RET_OK;
    char uplolow = LAPACK_UPLOW_LOW;
    
    dpotri(&uplolow, &lap->nCol, dFullMatrix, &lap->nCol, &info);
    HUtilMatSymmetrize(lap->nCol, dFullMatrix);
    
    assert( info == LAPACK_RET_OK );
    
    return;
}

static void lapackLinSolverClear( void *chol ) {
    
    if ( !chol ) {
        return;
    }
    
    lapack_linsys *lap = (lapack_linsys *) chol;
    HDSDP_FREE(lap->dFullMatElem);
    
    HDSDP_ZERO(lap, lapack_linsys, 1);
    
    return;
}

static void lapackLinSolverDestroy( void **pchol ) {
    
    if ( !pchol ) {
        return;
    }
    
    lapackLinSolverClear(*pchol);
    HDSDP_FREE(*pchol);
    
    return;
}

/* Dense iterative solver interface */
static hdsdp_retcode conjGradLinSolverCreate( void **pchol, int nCol ) {
    
    hdsdp_retcode retcode = HDSDP_RETCODE_OK;
    HDSDP_NULLCHECK(pchol);
    
    iterative_linsys *cg = NULL;
    HDSDP_INIT(cg, iterative_linsys, 1);
    HDSDP_MEMCHECK(cg);
    HDSDP_ZERO(cg, iterative_linsys, 1);
    
    cg->nCol = nCol;
    
    HDSDP_INIT(cg->iterResi, double, nCol);
    HDSDP_INIT(cg->iterResiNew, double, nCol);
    HDSDP_INIT(cg->iterDirection, double, nCol);
    HDSDP_INIT(cg->preInvResi, double, nCol);
    HDSDP_INIT(cg->MTimesDirection, double, nCol);
    HDSDP_INIT(cg->iterVec, double, nCol);
    HDSDP_INIT(cg->iterVecAuxi, double, nCol);
    HDSDP_INIT(cg->rhsBuffer, double, nCol);
    
    cg->useJacobi = 1;
    HDSDP_INIT(cg->JacobiPrecond, double, nCol);
    HDSDP_CALL(lapackLinSolverCreate((void **)&cg->lap, nCol));
    
    /* Set parameters */
    cg->params.absTol = 1e-06;
    cg->params.relTol = 1e-06;
    cg->params.maxIter = HDSDP_MAX(50, nCol / 20);
    cg->params.nRestartFreq = 20;
        
exit_cleanup:
    return retcode;
}

static void conjGradLinSolverSetParam( void *chol, void *iterParams ) {
    
    iterative_linsys *cg = (iterative_linsys *) chol;
    iterative_params *params = (iterative_params *) iterParams;
    
    cg->params.absTol = params->absTol;
    cg->params.relTol = params->relTol;
    cg->params.maxIter = params->maxIter;
    cg->params.nRestartFreq = params->nRestartFreq;
    
    return;
}

static hdsdp_retcode conjGradNumeric( void *chol, int *dummy1, int *dummy2, double *dFullMatrix ) {
    
    (void) chol;
    (void) dummy1;
    (void) dummy2;
    
    iterative_linsys *cg = (iterative_linsys *) chol;
    cg->pFullMatElem = dFullMatrix;
    
    return HDSDP_RETCODE_OK;
}

static hdsdp_retcode conjGradSolve( void *chol, double *rhsVec, double *solVec ) {
    
    hdsdp_retcode retcode = HDSDP_RETCODE_OK;
    /* TODO: Implement the conjugate gradient solver */
    
    return retcode;
}

static hdsdp_retcode conjGradSolveN( void *chol, int nRhs, double *rhsVec, double *solVec ) {
    
    hdsdp_retcode retcode = HDSDP_RETCODE_OK;
    iterative_linsys *cg = (iterative_linsys *) chol;
    
    for ( int iCol = 0; iCol < nRhs; ++iCol ) {
        HDSDP_CALL(conjGradSolve(chol, rhsVec + iCol * cg->nCol, solVec + iCol * cg->nCol))
    }
    
exit_cleanup:
    
    return retcode;
}

static void conjGradClear( void *chol ) {
    
    if ( !chol ) {
        return;
    }
    
    iterative_linsys *cg = (iterative_linsys *) chol;
    
    HDSDP_FREE(cg->iterResi);
    HDSDP_FREE(cg->iterResiNew);
    HDSDP_FREE(cg->iterDirection);
    HDSDP_FREE(cg->preInvResi);
    HDSDP_FREE(cg->MTimesDirection);
    HDSDP_FREE(cg->iterVec);
    HDSDP_FREE(cg->iterVecAuxi);
    HDSDP_FREE(cg->rhsBuffer);
    HDSDP_FREE(cg->JacobiPrecond);
    
    lapackLinSolverDestroy((void **) &cg->lap);
    
    return;
}

static void conjGradDestroy( void **pchol ) {
    
    if ( !pchol ) {
        return;
    }
    
    conjGradClear(*pchol);
    HDSDP_FREE(*pchol);
    
    return;
}
