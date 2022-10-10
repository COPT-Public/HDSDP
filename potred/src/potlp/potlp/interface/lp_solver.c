/** @file lpdata.c
 *  @brief Implement the LP data interface for potential reduction
 *
 * @TODO: Add more detailed comments
 */

#include "lp_solver.h"
#include "lp_solver.h"
#include "pot_solver.h"
#include "pot_structs.h"
#include "pot_utils.h"
#include "pot_vector.h"
#include "pot_constr_mat.h"
#include "pot_objfunc.h"
#include "vec_mat.h"

#include <sys/time.h>
#include <math.h>

static double my_clock( void ) {
    
    struct timeval t;
    gettimeofday(&t, NULL);
    return (1.0e-6 * t.tv_usec + t.tv_sec);
}

static void potLpRuizGetColMaxElem( potlp_solver *potlp ) {
    
    
    
    return;
}

/* Abstract Ruiz scaling */
static pot_int potLpRuizScale( potlp_solver *potlp ) {
    
    pot_int retcode = RETCODE_OK;
    
    pot_int nCol = potlp->nCol;
    pot_int nRow = potlp->nConstr;
    pot_int nColMatElem = potlp->colMatBeg[nCol];
    
    double *ruizWorkDiagRow = NULL;
    double *ruizWorkDiagCol = NULL;
    double *ruizWorkColMatElem = NULL;
    double *ruizWorkColMatElemTrans = NULL;
    double *ruizWorkColVal = NULL;
    double *ruizWorkColValTrans = NULL;
    double *ruizWorkRhs = NULL;
    double *ruizWorkRhsTrans = NULL;
    
    double *ruizScalDiagRow = potlp->ruizRow;
    double *ruizScalDiagCol = potlp->ruizCol;
    
    /* Allocate memory for ruiz */
    POTLP_INIT(ruizWorkColMatElem, double, nColMatElem);
    POTLP_INIT(ruizWorkColMatElemTrans, double, nColMatElem);
    POTLP_INIT(ruizWorkColVal, double, nCol);
    POTLP_INIT(ruizWorkColValTrans, double, nCol);
    POTLP_INIT(ruizWorkRhs, double, nRow);
    POTLP_INIT(ruizWorkRhsTrans, double, nRow);
    
    if ( !ruizScalDiagRow || !ruizScalDiagCol ||
         !ruizWorkColMatElem || !ruizWorkColMatElemTrans || !ruizWorkColVal || !ruizWorkColValTrans ||
         !ruizWorkRhs || !ruizWorkRhsTrans || !ruizWorkDiagCol || !ruizWorkDiagRow ) {
        retcode = RETCODE_FAILED;
        goto exit_cleanup;
    }
    
    POTLP_MEMCPY(ruizWorkColMatElem, potlp->colMatElem, double, nColMatElem);
    POTLP_MEMCPY(ruizWorkColMatElemTrans, potlp->colMatElem, double, nColMatElem);
    POTLP_MEMCPY(ruizWorkColVal, potlp->lpObj, double, nCol);
    POTLP_MEMCPY(ruizWorkColValTrans, potlp->lpObj, double, nCol);
    POTLP_MEMCPY(ruizWorkRhs, potlp->lpRHS, double, nRow);
    POTLP_MEMCPY(ruizWorkRhsTrans, potlp->lpRHS, double, nRow);
    
    /* Start ruiz scaling */
    int maxRuizIter = potlp->intParams[INT_PARAM_MAXRUIZITER];
    
    for ( int i = 0; i < maxRuizIter; ++i ) {
        
        
        
        
        
        
    }
    
exit_cleanup:
    
    POTLP_FREE(ruizWorkColMatElem);
    POTLP_FREE(ruizWorkColMatElemTrans);
    POTLP_FREE(ruizWorkColVal);
    POTLP_FREE(ruizWorkColValTrans);
    POTLP_FREE(ruizWorkRhs);
    POTLP_FREE(ruizWorkRhsTrans);
    
    return retcode;
}

/* Constraint methods */
static void potLpConstrMatImplPrepareX( void *AMatData, pot_vec *xInit ) {
    
    potVecReset(xInit);
    
    for ( int i = xInit->n - xInit->ncone; i < xInit->n; ++i ) {
        xInit->x[i] = 1.0;
    }
    
    return;
}

static void potLpConstrMatImplProject( void *AMatData, pot_vec *xVec ) {
    
    double eTx = potVecSumCone(xVec);
    potVecConeAddConstant(xVec, -eTx / xVec->ncone);
    
    return;
}

static void potLpConstrMatImplScalProject( void *AMatData, pot_vec *xVec, pot_vec *yVec ) {
    
    potVecNormalize(xVec);
    double xTy = potVecSumScalCone(xVec, yVec);
    potVecConeAxpy(-xTy, xVec, yVec);
    
    return;
}

static void potLpConstrMatImplMonitor( void *AMatData, void *info ) {
    
    return;
}

static void potLpObjFISetupRes( potlp_solver *potlp, pot_vec *xVec ) {
    
    pot_int nRow = potlp->nConstr;
    pot_int nCol = potlp->nCol;
    
    double *rowDual = xVec->x;
    double *colVal = rowDual + nRow;
    double *colSlack = colVal + nCol;
    double *kappaVar = colSlack + nCol;
    double *tauVar = kappaVar + 1;
    
    double *lpObj = potlp->lpObj;
    double *lpRHS = potlp->lpRHS;
    
    double *pRes = potlp->pRes;
    double *dRes = potlp->dRes;
    double *compl = potlp->cplRes;
    
    pot_int *Ap = potlp->colMatBeg;
    pot_int *Ai = potlp->colMatIdx;
    double  *Ax = potlp->colMatElem;
    
    /* Primal and dual objective */
    potlp->pObjVal = dot(&nCol, colVal, &potIntConstantOne, potlp->lpObj, &potIntConstantOne);
    potlp->dObjVal = dot(&nRow, rowDual, &potIntConstantOne, potlp->lpRHS, &potIntConstantOne);
    
    for ( int i = 0; i < nRow; ++i ) {
        pRes[i] = -lpRHS[i] * tauVar[0];
    }
    
    for ( int i = 0; i < nCol; ++i ) {
        dRes[i] = lpObj[i] * tauVar[0];
    }
    
    /* Primal infeasibility */
    spMatAxpy(nCol, Ap, Ai, Ax, 1.0, colVal, pRes);
    
    /* Dual infeasibility */
    spMatATxpy(nCol, Ap, Ai, Ax, -1.0, rowDual, dRes);
    axpy(&nCol, &potDblConstantMinusOne, colSlack, &potIntConstantOne, dRes, &potIntConstantOne);
    
    /* Complementarity */
    potlp->kappa = *kappaVar;
    potlp->tau = *tauVar;
    *compl = potlp->dObjVal - potlp->pObjVal - *kappaVar;
    
    /* Get statistics */
    potlp->pInfeas = nrm2(&nRow, pRes, &potIntConstantOne);
    potlp->dInfeas = nrm2(&nCol, dRes, &potIntConstantOne);
    potlp->complInfeas = *compl;
    potlp->complGap = potlp->pObjVal - potlp->dObjVal;
    
    return;
}

static void potLpObjFISetupGrad( potlp_solver *potlp, pot_vec *gVec ) {
    
    pot_int nRow = potlp->nConstr;
    pot_int nCol = potlp->nCol;
        
    double *gRowDual = gVec->x;
    double *gColVal = gRowDual + nRow;
    double *gColSlack = gColVal + nCol;
    double *gKappaVar = gColSlack + nCol;
    double *gTauVar = gKappaVar + 1;
    
    double *lpObj = potlp->lpObj;
    double *lpRHS = potlp->lpRHS;
    
    double *pRes = potlp->pRes;
    double *dRes = potlp->dRes;
    double compl = *potlp->cplRes;
    
    pot_int *Ap = potlp->colMatBeg;
    pot_int *Ai = potlp->colMatIdx;
    double  *Ax = potlp->colMatElem;
    
    /* Gradient of row dual */
    axpy(&nRow, &compl, lpRHS, &potIntConstantOne, gRowDual, &potIntConstantOne);
    spMatAxpy(nCol, Ap, Ai, Ax, -1.0, dRes, gRowDual);
    
    /* Gradient of column value */
    double ncompl = -compl;
    axpy(&nCol, &ncompl, lpObj, &potIntConstantOne, gColVal, &potIntConstantOne);
    spMatATxpy(nCol, Ap, Ai, Ax, 1.0, pRes, gColVal);
    
    /* Gradient of column slack */
    for ( int i = 0; i < nCol; ++i ) {
        gColSlack[i] = - dRes[i];
    }
    
    /* Gradient of kappa */
    gKappaVar[0] = - compl;
    
    /* Gradient of tau */
    gTauVar[0] = dot(&nCol, lpObj, &potIntConstantOne, dRes, &potIntConstantOne);
    gTauVar[0] -= dot(&nRow, lpRHS, &potIntConstantOne, pRes, &potIntConstantOne);
    
    return;
}

static void potLpObjFISetupHVec( potlp_solver *potlp, pot_vec *xVec, pot_vec *fHvec ) {
    
    /* TODO: Unify the three methods */
    pot_int nRow = potlp->nConstr;
    pot_int nCol = potlp->nCol;
    
    double *rowDual = xVec->x;
    double *colVal = rowDual + nRow;
    double *colSlack = colVal + nCol;
    double *kappaVar = colSlack + nCol;
    double *tauVar = kappaVar + 1;
    
    double *lpObj = potlp->lpObj;
    double *lpRHS = potlp->lpRHS;
    
    double *pRes = potlp->pRes;
    double *dRes = potlp->dRes;
    double *compl = potlp->cplRes;
    
    pot_int *Ap = potlp->colMatBeg;
    pot_int *Ai = potlp->colMatIdx;
    double  *Ax = potlp->colMatElem;
    
    double pObjVal = dot(&nCol, colVal, &potIntConstantOne, potlp->lpObj, &potIntConstantOne);
    double dObjVal = dot(&nRow, rowDual, &potIntConstantOne, potlp->lpRHS, &potIntConstantOne);
    
    for ( int i = 0; i < nRow; ++i ) {
        pRes[i] = -lpRHS[i] * tauVar[0];
    }
    
    for ( int i = 0; i < nCol; ++i ) {
        dRes[i] = lpObj[i] * tauVar[0];
    }
    
    /* Primal infeasibility */
    spMatAxpy(nCol, Ap, Ai, Ax, 1.0, colVal, pRes);
    
    /* Dual infeasibility */
    spMatATxpy(nCol, Ap, Ai, Ax, -1.0, rowDual, dRes);
    axpy(&nCol, &potDblConstantMinusOne, colSlack, &potIntConstantOne, dRes, &potIntConstantOne);
    
    *compl = dObjVal - pObjVal - *kappaVar;
    potLpObjFISetupGrad(potlp, fHvec);

    return;
}

/* Objective methods */
static double potLpObjFImplVal( void *objFData, pot_vec *xVec ) {
    
    potlp_solver *potlp = (potlp_solver *) objFData;
    potLpObjFISetupRes(potlp, xVec);
    
    double pInfeas = potlp->pInfeas;
    double dInfeas = potlp->dInfeas;
    double complInfeas = potlp->complInfeas;
    
    double objFVal = pInfeas * pInfeas + dInfeas * dInfeas + complInfeas * complInfeas;
    
    return 0.5 * objFVal;
}

static void potLpObjFImplGrad( void *objFData, pot_vec *xVec, pot_vec *fGrad ) {
    
    potVecReset(fGrad);
    potlp_solver *potlp = (potlp_solver *) objFData;
    potLpObjFISetupGrad(potlp, fGrad);
    
    return;
}

static void potLpObjFImplHess( void *objFData, pot_vec *xVec, double *fHess ) {
    
    assert( 0 );
    return;
}

static void potLpObjFImplHVec( void *objFData, pot_vec *xVec, pot_vec *fHvec ) {

    potVecReset(fHvec);
    potlp_solver *potlp = (potlp_solver *) objFData;
    potLpObjFISetupHVec(potlp, xVec, fHvec);
    
    return;
}

static void potLpObjFImplMonitor( void *objFData, void *info ) {
    
    potlp_solver *potlp = (potlp_solver *) objFData;
    potlp->nIter += 1;
    
    if ( potlp->nIter == 1 ) {
        printf("Iteration log. \n");
        printf("%8s  %10s  %10s  %10s  %10s  %10s  %10s \n", "nIter", "pObj", "dObj", "rGap", "pInf", "dInf", "k/t");
    }
    
    double absGap = fabs(potlp->pObjVal - potlp->dObjVal);
    double relGap = absGap / (1.0 + fabs(potlp->pObjVal) + fabs(potlp->dObjVal));
    double pInfeas = potlp->pInfeas / potlp->tau;
    double dInfeas = potlp->dInfeas / potlp->tau;
    
    if ( potlp->nIter % 500 == 0 || potlp->nIter == 1 ) {
        printf("%8lld  %10.3e  %10.3e  %10.3e  %10.3e  %10.3e  %10.3e |%5.2f [s] \n",
               potlp->nIter, potlp->pObjVal, potlp->dObjVal, relGap, pInfeas,
               dInfeas, potlp->kappa / potlp->tau,
               my_clock() - potlp->startT);
    }
    
    double absFeasTol = 1e-04;
//    double relFeasTol = 1e-03;
    double absOptTol = 1e-04;
    double relOptTol = 1e-03;
    
    if ( absGap < absOptTol && relGap < relOptTol && pInfeas < absFeasTol && dInfeas < absFeasTol ) {
        int *intInfo = (int *) info;
        *intInfo = 1;
    }
    
    if ( potlp->nIter >= potlp->intParams[INT_PARAM_MAXITER]) {
        int *intInfo = (int *) info;
        *intInfo = 1;
    }
    
    return;
}

/** @brief Scale the problem data
 *
 */
static void LPSolverIScale( potlp_solver *potlp ) {
    
    double rhsOneNorm = 0.0;
    double objOneNorm = 0.0;
    
    double *lpObj = potlp->lpObj;
    double *lpRHS = potlp->lpRHS;
    
    for ( int i = 0; i < potlp->nCol; ++i ) {
        objOneNorm += fabs(lpObj[i]);
    }
    
    for ( int i = 0; i < potlp->nConstr; ++i ) {
        rhsOneNorm += fabs(lpRHS[i]);
    }
    
    potlp->lpObjNorm = objOneNorm;
    potlp->lpRHSNorm = rhsOneNorm;
    
    objOneNorm += 1.0;
    rhsOneNorm += 1.0;
    
    for ( int i = 0; i < potlp->nCol; ++i ) {
        lpObj[i] = lpObj[i] / objOneNorm;
    }
    
    for ( int i = 0; i < potlp->nConstr; ++i ) {
        lpRHS[i] = lpRHS[i] / rhsOneNorm;
    }
    
    return;
}

extern pot_int LPSolverCreate( potlp_solver **ppotlp ) {
    
    pot_int retcode = RETCODE_OK;
    
    if ( !ppotlp ) {
        retcode = RETCODE_FAILED;
        goto exit_cleanup;
    }
    
    potlp_solver *potlp = NULL;
    
    POTLP_INIT(potlp, potlp_solver, 1);
    
    if ( !potlp ) {
        retcode = RETCODE_FAILED;
        goto exit_cleanup;
    }
    
    memset(potlp, 0, sizeof(potlp_solver));
    
    memcpy(potlp->dblParams, defaultDblParam, sizeof(double) * NUM_DBL_PARAM);
    memcpy(potlp->intParams, defaultIntParam, sizeof(int) * NUM_INT_PARAM);
    
    POT_CALL(potLPCreate(&potlp->potIterator));
    POT_CALL(potConstrMatCreate(&potlp->potConstrMat));
    POT_CALL(potObjFCreate(&potlp->potObjF));
    
    potlp->nIter = 0;
    potlp->startT = my_clock();
    
    *ppotlp = potlp;
    
exit_cleanup:
    return retcode;
}

extern pot_int LPSolverInit( potlp_solver *potlp, pot_int nCol, pot_int nRow ) {
    
    pot_int retcode = RETCODE_OK;
    assert( nCol > 0 && nRow > 0 );
    
    potlp->nCol = nCol;
    potlp->nConstr = nRow;
    
    /* x; s; kappa; tau*/
    pot_int coneDim = 2 * nCol + 2;
    
    /* y */
    pot_int varDim = coneDim + nRow;
    
    /* Construct internal solver */
    POT_CALL(potLPInit(potlp->potIterator, varDim, coneDim));
    POT_CALL(potConstrMatInit(potlp->potConstrMat, 1, varDim));
    POT_CALL(potObjFInit(potlp->potObjF, varDim));
    
    /* Link objective, constraint and data to the iterator */
    POT_CALL(potLPSetLinearConstrs(potlp->potIterator, potlp->potConstrMat));
    POT_CALL(potLPSetObj(potlp->potIterator, potlp->potObjF));
    
    /* Line data and methods to objective */
    potlp->potObjF->objFData = potlp;
    potlp->potObjF->objFVal = potLpObjFImplVal;
    potlp->potObjF->objFGrad = potLpObjFImplGrad;
    potlp->potObjF->objFHess = potLpObjFImplHess;
    potlp->potObjF->objFHVec = potLpObjFImplHVec;
    potlp->potObjF->objFMonitor = potLpObjFImplMonitor;
    
    /* Line data and methods to constraint  */
    potlp->potConstrMat->AMatData = potlp;
    potlp->potConstrMat->AMatPrepareX = potLpConstrMatImplPrepareX;
    potlp->potConstrMat->AMatProject = potLpConstrMatImplProject;
    potlp->potConstrMat->AMatScalProject = potLpConstrMatImplScalProject;
    potlp->potConstrMat->AMatMonitor = potLpConstrMatImplMonitor;
    
exit_cleanup:
    return retcode;
}

extern pot_int LPSolverSetData( potlp_solver *potlp, pot_int *Ap, pot_int *Ai, double *Ax, double *lpObj, double *lpRHS ) {
    
    pot_int retcode = RETCODE_OK;
    
    pot_int nCol = potlp->nCol;
    pot_int nRow = potlp->nConstr;
    pot_int nNz = Ap[nCol];
    
    POTLP_INIT(potlp->colMatBeg, pot_int, nCol + 1);
    POTLP_INIT(potlp->colMatIdx, pot_int, nNz);
    POTLP_INIT(potlp->colMatElem, double, nNz);
    
    POTLP_INIT(potlp->lpObj, double, nCol);
    POTLP_INIT(potlp->lpRHS, double, nRow);
    POTLP_INIT(potlp->pdcRes, double, nRow + nCol + 1);
    
    POTLP_INIT(potlp->ruizCol, double, nRow + nCol + 1);
    POTLP_INIT(potlp->ruizRow, double, nRow + nCol + 1);
    
    if ( !potlp->colMatBeg || !potlp->colMatIdx || !potlp->colMatElem ||
         !potlp->lpObj || !potlp->lpRHS || !potlp->pdcRes ||
         !potlp->ruizCol || !potlp->ruizRow ) {
        retcode = RETCODE_FAILED;
        goto exit_cleanup;
    }
    
    potlp->pRes = potlp->pdcRes;
    potlp->dRes = potlp->pRes + nRow;
    potlp->cplRes = potlp->dRes + nCol;
    
    /* Copy data */
    memcpy(potlp->colMatBeg, Ap, sizeof(pot_int) * (nCol + 1));
    memcpy(potlp->colMatIdx, Ai, sizeof(pot_int) * nNz);
    memcpy(potlp->colMatElem, Ax, sizeof(double) * nNz);
    memcpy(potlp->lpObj, lpObj, sizeof(double) * nCol);
    memcpy(potlp->lpRHS, lpRHS, sizeof(double) * nRow);
    
exit_cleanup:
    return retcode;
}

extern pot_int LPSolverOptimize( potlp_solver *potlp ) {
    
    pot_int retcode = RETCODE_OK;
    
    if ((0)) {
        LPSolverIScale(potlp);
    }
    
    retcode = potReductionSolve(potlp->potIterator);
    
exit_cleanup:
    return retcode;
}

extern void LPSolverClear( potlp_solver *potlp ) {
    
    if ( !potlp ) {
        return;
    }
    
    POTLP_FREE(potlp->colMatBeg);
    POTLP_FREE(potlp->colMatIdx);
    POTLP_FREE(potlp->colMatElem);
    
    POTLP_FREE(potlp->lpRHS);
    POTLP_FREE(potlp->lpObj);
    
    POTLP_FREE(potlp->ruizCol);
    POTLP_FREE(potlp->ruizRow);
    
    POTLP_FREE(potlp->pdcRes);
    
    potConstrMatDestroy(&potlp->potConstrMat);
    potObjFDestroy(&potlp->potObjF);
    potLPDestroy(&potlp->potIterator);
    
    memset(potlp, 0, sizeof(potlp_solver));
    
    return;
}

extern void LPSolverDestroy( potlp_solver **ppotlp ) {
    
    if ( !ppotlp ) {
        return;
    }
    
    LPSolverClear(*ppotlp);
    POTLP_FREE(*ppotlp);
    
    return;
}
