#ifdef HEADERPATH
#include "interface/hdsdp_utils.h"
#include "linalg/def_hdsdp_linsolver.h"
#include "linalg/hdsdp_linsolver.h"
#include "linalg/dense_opts.h"
#include "linalg/vec_opts.h"
#else
#include "hdsdp_utils.h"
#include "def_hdsdp_linsolver.h"
#include "hdsdp_linsolver.h"
#include "dense_opts.h"
#include "vec_opts.h"
#endif

#include <math.h>

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
    
    /* Get diagonal of Cholesky factor */
    hdsdp_retcode retcode = HDSDP_RETCODE_OK;
    pardiso_linsys *pds = (pardiso_linsys *) chol;
    
    int one = 1, pdsret = PARDISO_RET_OK;
    pardiso_getdiag((const void **) pds->pt, (void *) diagVec, (void *) pds->dWork, &one, &pdsret);
    
    if ( pdsret != PARDISO_RET_OK ) {
        retcode = HDSDP_RETCODE_FAILED;
        goto exit_cleanup;
    }
    
    /* For pardiso we need to take squareroot */
    for ( int i = 0; i < pds->nCol; ++i ) {
        diagVec[i] = sqrt(diagVec[i]);
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
        pElem += pds->nCol + 1;
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
    lapack_flinsys *lap = NULL;
    HDSDP_INIT(lap, lapack_flinsys, 1);
    HDSDP_MEMCHECK(lap);
    
    HDSDP_ZERO(lap, lapack_flinsys, 1);
    lap->nCol = nCol;
    
    HDSDP_INIT(lap->dFullMatElem, double, nCol * nCol);
    HDSDP_MEMCHECK(lap->dFullMatElem);
    
    *pchol = lap;
    
exit_cleanup:
    return retcode;
}

static void lapackLinSolverSetParamsDummy( void *chol, void *param ) {
    
    (void) chol;
    (void) param;
    
    return;
}

static hdsdp_retcode lapackLinSolverSymbolic( void *chol, int *dummy1, int *dummy2 ) {
    
    (void) chol;
    (void) dummy1;
    (void) dummy2;
    
    return HDSDP_RETCODE_OK;
}

static hdsdp_retcode lapackLinSolverNumeric( void *chol, int *dummy1, int *dummy2, double *dFullMatrix ) {
    
    hdsdp_retcode retcode = HDSDP_RETCODE_OK;
    lapack_flinsys *lap = (lapack_flinsys *) chol;
    
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
    lapack_flinsys *lap = (lapack_flinsys *) chol;
    
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
    
    lapack_flinsys *lap = (lapack_flinsys *) chol;
    
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
    
    lapack_flinsys *lap = (lapack_flinsys *) chol;
    
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
    
    lapack_flinsys *lap = (lapack_flinsys *) chol;
    
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
    
    lapack_flinsys *lap = (lapack_flinsys *) chol;
    
    for ( int iRow = 0; iRow < lap->nCol; ++iRow ) {
        diagVec[iRow] = lap->dFullMatElem[iRow * lap->nCol + iRow];
    }
    
    return HDSDP_RETCODE_OK;
}

static void lapackLinSolverInvert( void *chol, double *dFullMatrix, double *dAuxiMatrix ) {
    
    (void) dAuxiMatrix;
    lapack_flinsys *lap = (lapack_flinsys *) chol;
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
    
    lapack_flinsys *lap = (lapack_flinsys *) chol;
    HDSDP_FREE(lap->dFullMatElem);
    
    HDSDP_ZERO(lap, lapack_flinsys, 1);
    
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
    HDSDP_INIT(cg->rhsBuffer, double, nCol);
    
    cg->useJacobi = 1;
    HDSDP_INIT(cg->JacobiPrecond, double, nCol);
    HDSDP_CALL(lapackLinSolverCreate((void **)&cg->lap, nCol));
    
    /* Set parameters */
    cg->params.absTol = 1e-06;
    cg->params.relTol = 1e-06;
    cg->params.maxIter = HDSDP_MAX(50, nCol / 20);
    cg->params.nRestartFreq = 20;
    
    *pchol = cg;
        
exit_cleanup:
    return retcode;
}

static void conjGradLinSolverSetParam( void *chol, void *iterParams ) {
    
    iterative_linsys *cg = (iterative_linsys *) chol;
    iterative_params *params = (iterative_params *) iterParams;
    
    if ( params->absTol > 0 ) {
        cg->params.absTol = params->absTol;
    }
    
    if ( params->relTol > 0 ) {
        cg->params.relTol = params->relTol;
    }
    
    if ( params->maxIter > 0 ) {
        cg->params.maxIter = params->maxIter;
    }
    
    if ( params->nRestartFreq > 0 )
    cg->params.nRestartFreq = params->nRestartFreq;
    
    return;
}

static hdsdp_retcode conjGradLinSolverSymbolic( void *chol, int *dummy1, int *dummy2 ) {
    
    (void) chol;
    (void) dummy1;
    (void) dummy2;
    
    return HDSDP_RETCODE_OK;
}

static hdsdp_retcode conjGradLinSolverNumeric( void *chol, int *dummy1, int *dummy2, double *dFullMatrix ) {
    
    (void) chol;
    (void) dummy1;
    (void) dummy2;
    
    iterative_linsys *cg = (iterative_linsys *) chol;
    cg->fullMatElem = dFullMatrix;
    
    return HDSDP_RETCODE_OK;
}

static hdsdp_retcode conjGradLinSolverPsdCheck( void *chol, int *dummy1, int *dummy2, double *dummy3, int *dummy4 ) {
    
    (void) chol;
    (void) dummy1;
    (void) dummy2;
    (void) dummy3;
    (void) dummy4;
    
    return HDSDP_RETCODE_FAILED;
}

static void conjGradLinSolverForwardN( void *chol, int dummy1, double *dummy2, double *dummy3 ) {
    
    (void) chol;
    (void) dummy1;
    (void) dummy2;
    (void) dummy3;
    
    return;
}

static void conjGradLinSolverBackwardN( void *chol, int dummy1, double *dummy2, double *dummy3 ) {
    
    (void) chol;
    (void) dummy1;
    (void) dummy2;
    (void) dummy3;
    
    return;
}

static hdsdp_retcode conjGradBuildPreconditioner( iterative_linsys *cg ) {
    
    hdsdp_retcode retcode = HDSDP_RETCODE_OK;
    
    /*
     Jacobi preconditioner just uses the diagonal elements of a matrix
     https://netlib.org/linalg/html_templates/node55.html
     TODO: consider sparsity-approximate inverse method mentioned in
     https://math.stackexchange.com/questions/978051/find-diagonal-of-inverse-matrix
     */
    if ( cg->useJacobi ) {
        for ( int iCol = 0; iCol < cg->nCol; ++iCol ) {
            cg->JacobiPrecond[iCol] = cg->fullMatElem[iCol + iCol * cg->nCol];
        }
    } else {
        /* Use Cholesky as the preconditioner */
        HDSDP_CALL(lapackLinSolverNumeric((void *) cg->lap, NULL, NULL, cg->fullMatElem));
    }
    
exit_cleanup:
    
    return retcode;
}

static void conjGradApplyPreconditioner( iterative_linsys *cg, double *rhsVec ) {
    
    /* Apply preconditioner rhsVec <- P \ rhsVec */
    
    if ( cg->useJacobi ) {
        vvrscl(&cg->nCol, cg->JacobiPrecond, rhsVec);
    } else {
        lapackLinSolverSolveN((void *) cg, 1, rhsVec, NULL);
    }
    
    return;
}

static hdsdp_retcode conjGradSolve( iterative_linsys *cg, double *rhsVec, double *solVec ) {
    
    hdsdp_retcode retcode = HDSDP_RETCODE_OK;
    
    int nCGMaxIter = cg->params.maxIter;
    cg->useJacobi = cg->params.useJacobi;
    
    /* CG working scalers */
    double dDotMTimesd = 0.0;
    double resiDotPreInvResi = 0.0;
    double alpha = 0.0;
    double beta = 0.0;
    double resiNorm = 0.0;
    double rhsNorm = 0.0;
    
    double *iterResi = cg->iterResi;
    double *iterResiNew = cg->iterResiNew;
    double *iterDirection = cg->iterDirection;
    double *MTimesDirection = cg->MTimesDirection;
    double *preInvResi = cg->preInvResi;
    double *iterVec = cg->iterVec;
    
    double cgStartTime = HUtilGetTimeStamp();
    double cgDuration = 0.0;
    
    int iter = 0;
    
    /* Get restart frequency */
    int nRestartFreq = HDSDP_MAX(cg->params.nRestartFreq, 20);
    
    int incx = 1;
    double one = 1.0;
    double minusOne = -1.0;
    
    HDSDP_ZERO(iterVec, double, cg->nCol);
    HDSDP_MEMCPY(iterResi, rhsVec, double, cg->nCol);
    
    resiNorm = nrm2(&cg->nCol, iterResi, &incx);
    rhsNorm = nrm2(&cg->nCol, rhsVec, &incx);
    
    double cgTol = HDSDP_MIN(cg->params.absTol, rhsNorm * cg->params.relTol);
    
    /* Initial value is already a qualified solution */
    if ( resiNorm < cgTol ) {
        cg->solStatus = ITERATIVE_STATUS_OK;
        goto exit_cleanup;
    }
    
    /* Build preconditioner */
    HDSDP_CALL(conjGradBuildPreconditioner(cg));
    HDSDP_MEMCPY(iterDirection, iterResi, double, cg->nCol);
    
    /* Initial iteration preparation */
    conjGradApplyPreconditioner(cg, iterDirection);
    HDSDP_MEMCPY(preInvResi, iterDirection, double, cg->nCol);
    fds_symv(cg->nCol, 1.0, cg->fullMatElem, iterVec, 0.0, MTimesDirection);
    
    for ( iter = 0; iter < nCGMaxIter; ++iter ) {
        
        resiDotPreInvResi = dot(&cg->nCol, preInvResi, &incx, iterResi, &incx);
        dDotMTimesd = dot(&cg->nCol, iterDirection, &incx, MTimesDirection, &incx);
        alpha = resiDotPreInvResi / dDotMTimesd;
        axpy(&cg->nCol, &alpha, iterDirection, &incx, iterVec, &incx);
        
        if ( iter % nRestartFreq == 5 ) {
            /* Restart CG solver */
#ifdef HDSDP_CONJGRAD_DEBUG
            printf("Conjugate gradient restarts at iteration %d. \n", iter);
#endif
            HDSDP_MEMCPY(iterResi, rhsVec, double, cg->nCol);
            fds_symv(cg->nCol, 1.0, cg->fullMatElem, iterVec, 0.0, iterDirection);
            axpy(&cg->nCol, &minusOne, iterDirection, &incx, iterResi, &incx);
            HDSDP_MEMCPY(iterDirection, iterResi, double, cg->nCol);
            conjGradApplyPreconditioner(cg, iterDirection);
            fds_symv(cg->nCol, 1.0, cg->fullMatElem, iterDirection, 0.0, MTimesDirection);
            HDSDP_MEMCPY(preInvResi, iterResi, double, cg->nCol);
            conjGradApplyPreconditioner(cg, preInvResi);
        }
        
        HDSDP_MEMCPY(iterResiNew, iterResi, double, cg->nCol);
        double minusAlpha = -alpha;
        axpy(&cg->nCol, &minusAlpha, MTimesDirection, &incx, iterResiNew, &incx);
        HDSDP_MEMCPY(preInvResi, iterResiNew, double, cg->nCol);
        conjGradApplyPreconditioner(cg, preInvResi);
        beta = dot(&cg->nCol, iterResiNew, &incx, preInvResi, &incx);
        beta = beta / resiDotPreInvResi;
        axpby(&cg->nCol, &one, preInvResi, &incx, &beta, iterDirection, &incx);
        fds_symv(cg->nCol, 1.0, cg->fullMatElem, iterDirection, 0.0, MTimesDirection);
        HDSDP_MEMCPY(iterResi, iterResiNew, double, cg->nCol);
        resiNorm = nrm2(&cg->nCol, iterResi, &incx);
        
        if ( resiNorm != resiNorm ) {
            cg->solStatus = ITERATIVE_STATUS_NUMERICAL;
            retcode = HDSDP_RETCODE_FAILED;
            goto exit_cleanup;
        }
        
#ifdef HDSDP_CONJGRAD_DEBUG
        printf("Conjugate iteration %d. Resi: %10.3e", iter, resiNorm);
#endif
        
        if ( resiNorm < cgTol ) {
            cg->solStatus = ITERATIVE_STATUS_OK;
            retcode = HDSDP_RETCODE_FAILED;
            goto exit_cleanup;
        }
    }
    
    if ( iter >= nCGMaxIter ) {
        cg->solStatus = ITERATIVE_STATUS_MAXITER;
    }
    
    HDSDP_MEMCPY(solVec, iterVec, double, cg->nCol);

exit_cleanup:
    
    /* Collect solution time statistics */
    cgDuration = HUtilGetTimeStamp() - cgStartTime;
    cg->nSolves += 1;
    cg->solveTime += cgDuration;
    cg->nIters = iter;
    
    return retcode;
}

static hdsdp_retcode conjGradLinSolverSolveN( void *chol, int nRhs, double *rhsVec, double *solVec ) {
    
    hdsdp_retcode retcode = HDSDP_RETCODE_OK;
    iterative_linsys *cg = (iterative_linsys *) chol;
    
    if ( solVec ) {
        for ( int iCol = 0; iCol < nRhs; ++iCol ) {
            HDSDP_CALL(conjGradSolve(chol, rhsVec + iCol * cg->nCol, solVec + iCol * cg->nCol))
        }
    } else {
        /* Overwrite solution */
        for ( int iCol = 0; iCol < nRhs; ++iCol ) {
            HDSDP_MEMCPY(cg->rhsBuffer, rhsVec + iCol * cg->nCol, double, cg->nCol);
            HDSDP_CALL(conjGradSolve(chol, cg->rhsBuffer, rhsVec + iCol * cg->nCol));
        }
    }
    
exit_cleanup:
    
    return retcode;
}

static hdsdp_retcode conjGradLinSolverGetDiag( void *chol, double *dummy ) {
    
    (void) chol;
    (void) dummy;
    
    return HDSDP_RETCODE_FAILED;
}

static void conjGradLinSolverInvert( void *chol, double *dummy1, double *dummy2 ) {
    
    (void) chol;
    (void) dummy1;
    (void) dummy2;
    
    return;
}

static void conjGradLinSolverClear( void *chol ) {
    
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
    HDSDP_FREE(cg->rhsBuffer);
    HDSDP_FREE(cg->JacobiPrecond);
    
    lapackLinSolverDestroy((void **) &cg->lap);
    
    return;
}

static void conjGradLinSolverDestroy( void **pchol ) {
    
    if ( !pchol ) {
        return;
    }
    
    conjGradLinSolverClear(*pchol);
    HDSDP_FREE(*pchol);
    
    return;
}

extern hdsdp_retcode HFpLinsysCreate( hdsdp_linsys_fp **pHLin, int nCol, linsys_type Ltype ) {
    
    hdsdp_retcode retcode = HDSDP_RETCODE_OK;
    
    HDSDP_NULLCHECK(pHLin);
    
    hdsdp_linsys_fp *HLinsys = NULL;
    HDSDP_INIT(HLinsys, hdsdp_linsys_fp, 1);
    HDSDP_MEMCHECK(HLinsys);
    
    HLinsys->LinType = Ltype;
    HLinsys->nCol = nCol;
    
    switch ( Ltype ) {
        case HDSDP_LINSYS_DENSE_DIRECT:
            HLinsys->cholCreate = lapackLinSolverCreate;
            HLinsys->cholSetParam = lapackLinSolverSetParamsDummy;
            HLinsys->cholSymbolic = lapackLinSolverSymbolic;
            HLinsys->cholNumeric = lapackLinSolverNumeric;
            HLinsys->cholPsdCheck = lapackLinSolverPsdCheck;
            HLinsys->cholFSolve = lapackLinSolverForwardN;
            HLinsys->cholBSolve = lapackLinSolverBackwardN;
            HLinsys->cholSolve = lapackLinSolverSolveN;
            HLinsys->cholGetDiag = lapackLinSolverGetDiag;
            HLinsys->cholInvert = lapackLinSolverInvert;
            HLinsys->cholDestroy = lapackLinSolverDestroy;
            break;
        case HDSDP_LINSYS_SPARSE_DIRECT:
            HLinsys->cholCreate = pardisoLinSolverCreate;
            HLinsys->cholSetParam = pardisoLinSolverSetThreads;
            HLinsys->cholSymbolic = pardisoLinSolverSymbolic;
            HLinsys->cholNumeric = pardisoLinSolverNumeric;
            HLinsys->cholPsdCheck = pardisoLinSolverPsdCheck;
            HLinsys->cholFSolve = pardisoLinSolverForwardN;
            HLinsys->cholBSolve = pardisoLinSolverBackwardN;
            HLinsys->cholSolve = pardisoLinSolverSolveN;
            HLinsys->cholGetDiag = pardisoLinSolverGetDiag;
            HLinsys->cholInvert = pardisoLinSolverInvert;
            HLinsys->cholDestroy = pardisoLinSolverDestroy;
            break;
        case HDSDP_LINSYS_DENSE_ITERATIVE:
            HLinsys->cholCreate = conjGradLinSolverCreate;
            HLinsys->cholSetParam = conjGradLinSolverSetParam;
            HLinsys->cholSymbolic = conjGradLinSolverSymbolic;
            HLinsys->cholNumeric = conjGradLinSolverNumeric;
            HLinsys->cholPsdCheck = conjGradLinSolverPsdCheck;
            HLinsys->cholFSolve = conjGradLinSolverForwardN;
            HLinsys->cholBSolve = conjGradLinSolverBackwardN;
            HLinsys->cholSolve = conjGradLinSolverSolveN;
            HLinsys->cholGetDiag = conjGradLinSolverGetDiag;
            HLinsys->cholInvert = conjGradLinSolverInvert;
            HLinsys->cholDestroy = conjGradLinSolverDestroy;
            break;
        /* Not implemented */
        case HDSDP_LINSYS_SMALL_DIRECT:
            retcode = HDSDP_RETCODE_FAILED;
            goto exit_cleanup;
        default:
            retcode = HDSDP_RETCODE_FAILED;
            goto exit_cleanup;
    }
    
    *pHLin = HLinsys;
    
    HDSDP_CALL(HLinsys->cholCreate(&HLinsys->chol, nCol));
    
exit_cleanup:
    
    return retcode;
}

extern void HFpLinsysSetParam( hdsdp_linsys_fp *HLin, double relTol, double absTol, int nThreads, int maxIter, int nRestartFreq ) {
    
    if ( HLin->LinType == HDSDP_LINSYS_DENSE_ITERATIVE ) {
        
        iterative_params params;
        params.absTol = absTol;
        params.relTol = relTol;
        params.maxIter = relTol;
        params.nRestartFreq = nRestartFreq;
        HLin->cholSetParam(HLin->chol, &params);
        
    } else if ( HLin->LinType == HDSDP_LINSYS_SPARSE_DIRECT ) {
        
        void *pnThreads = (void *) (&nThreads);
        HLin->cholSetParam(HLin->chol, pnThreads);
        
    }
    
    return;
}

extern hdsdp_retcode HFpLinsysSymbolic( hdsdp_linsys_fp *HLin, int *colMatBeg, int *colMatIdx ) {
    
    hdsdp_retcode retcode = HDSDP_RETCODE_OK;
    HDSDP_CALL(HLin->cholSymbolic(HLin->chol, colMatBeg, colMatIdx));
    
exit_cleanup:
    
    return retcode;
}

extern hdsdp_retcode HFpLinsysNumeric( hdsdp_linsys_fp *HLin, int *colMatBeg, int *colMatIdx, double *colMatElem ) {
    
    hdsdp_retcode retcode = HDSDP_RETCODE_OK;
    HDSDP_CALL(HLin->cholNumeric(HLin->chol, colMatBeg, colMatIdx, colMatElem));
    
exit_cleanup:
    
    return retcode;
}

extern hdsdp_retcode HFpLinsysPsdCheck( hdsdp_linsys_fp *HLin, int *colMatBeg, int *colMatIdx, double *colMatElem, int *isPsd ) {
    
    hdsdp_retcode retcode = HDSDP_RETCODE_OK;
    HDSDP_CALL(HLin->cholPsdCheck(HLin->chol, colMatBeg, colMatIdx, colMatElem, isPsd));
    
exit_cleanup:
    
    return retcode;
}


extern void HFpLinsysFSolve( hdsdp_linsys_fp *HLin, int nRhs, double *rhsVec, double *solVec ) {
    
    HLin->cholFSolve(HLin->chol, nRhs, rhsVec, solVec);
    
    return;
}

extern void HFpLinsysBSolve( hdsdp_linsys_fp *HLin, int nRhs, double *rhsVec, double *solVec ) {
    
    HLin->cholBSolve(HLin->chol, nRhs, rhsVec, solVec);
    
    return;
}


extern hdsdp_retcode HFpLinsysSolve( hdsdp_linsys_fp *HLin, int nRhs, double *rhsVec, double *solVec ) {
    
    hdsdp_retcode retcode = HDSDP_RETCODE_OK;
    HDSDP_CALL(HLin->cholSolve(HLin->chol, nRhs, rhsVec, solVec));
    
exit_cleanup:
    
    return retcode;
}

extern hdsdp_retcode HFpLinsysGetDiag( hdsdp_linsys_fp *HLin, double *diagElem ) {
    
    hdsdp_retcode retcode = HDSDP_RETCODE_OK;
    HLin->cholGetDiag(HLin->chol, diagElem);
    
exit_cleanup:
    
    return retcode;
}


extern void HFpLinsysInvert( hdsdp_linsys_fp *HLin, double *dFullMatrix, double *dAuxiMatrix ) {
    
    HLin->cholInvert(HLin->chol, dFullMatrix, dAuxiMatrix);
    
    return;
}

extern void HFpLinsysClear( hdsdp_linsys_fp *HLin ) {
    
    if ( !HLin ) {
        return;
    }
    
    HLin->nCol = 0;
    HLin->cholDestroy(&HLin->chol);
    
    return;
}

extern void HFpLinsysDestroy( hdsdp_linsys_fp **HLin ) {
    
    if ( !HLin ) {
        return;
    }
    
    HFpLinsysClear(*HLin);
    HDSDP_FREE(*HLin);
    
    return;
}
