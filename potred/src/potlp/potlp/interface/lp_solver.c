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
    
    pot_int nRow = potlp->nRow;
    pot_int nCol = potlp->nCol;
    
    /* Post-multiply Ruiz */
    double *scalVals = potlp->scalVals;
    potVecExport(xVec, scalVals);
    vvscl(&xVec->n, potlp->ruizCol, scalVals);
    
    double *rowDual = scalVals;
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
    potlp->pInfeas = potlp->pInfeas / *tauVar;
    potlp->dInfeas = nrm2(&nCol, dRes, &potIntConstantOne);
    potlp->dInfeas = potlp->dInfeas / *tauVar;
    potlp->complInfeas = *compl;
    potlp->complGap = potlp->pObjVal - potlp->dObjVal;
    
    /* Pre-multiply Ruiz */
    int nRuizRow = nRow + nCol + 1;
    vvscl(&nRuizRow, potlp->ruizRow, potlp->pdcRes);
    
    return;
}

static void potLpObjFISetupGrad( potlp_solver *potlp, pot_vec *gVec ) {
    
    pot_int nRow = potlp->nRow;
    pot_int nCol = potlp->nCol;
    
    /* Pre-multiply Ruiz */
    int nRuizRow = nRow + nCol + 1;
    vvscl(&nRuizRow, potlp->ruizRow, potlp->pdcRes);
        
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
    
    /* Post-multiply Ruiz */
    vvscl(&gVec->n, potlp->ruizCol, gVec->x);
    
    return;
}

static void potLpObjFISetupHVec( potlp_solver *potlp, pot_vec *xVec, pot_vec *fHvec ) {
    
    pot_int nRow = potlp->nRow;
    pot_int nCol = potlp->nCol;
    
    /* Post-multiply Ruiz */
    double *scalVals = potlp->scalVals;
    potVecExport(xVec, scalVals);
    vvscl(&xVec->n, potlp->ruizCol, scalVals);
    
    double *rowDual = scalVals;
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
        pRes[i] = -lpRHS[i] * (*tauVar);
    }
    
    for ( int i = 0; i < nCol; ++i ) {
        dRes[i] = lpObj[i] * (*tauVar);
    }
    
    /* Primal infeasibility */
    spMatAxpy(nCol, Ap, Ai, Ax, 1.0, colVal, pRes);
    
    /* Dual infeasibility */
    spMatATxpy(nCol, Ap, Ai, Ax, -1.0, rowDual, dRes);
    axpy(&nCol, &potDblConstantMinusOne, colSlack, &potIntConstantOne, dRes, &potIntConstantOne);
    
    *compl = dObjVal - pObjVal - *kappaVar;
    
    /* Pre-multiply Ruiz */
    int nRuizRow = nRow + nCol + 1;
    vvscl(&nRuizRow, potlp->ruizRow, potlp->pdcRes);
    
    potLpObjFISetupGrad(potlp, fHvec);

    return;
}

/* Objective methods */
static double potLpObjFImplVal( void *objFData, pot_vec *xVec ) {
    
    potlp_solver *potlp = (potlp_solver *) objFData;
    potLpObjFISetupRes(potlp, xVec);
    int nRes = potlp->nRow + potlp->nCol + 1;
    double objFVal = nrm2(&nRes, potlp->pdcRes, &potIntConstantOne);
    
    return 0.5 * objFVal * objFVal;
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
        printf("Potential reduction log \n\n");
        printf("%8s  %10s  %10s  %10s  %10s  %10s  %10s \n", "nIter", "pObj", "dObj", "rGap", "pInf", "dInf", "k/t");
    }
    
    double objScalCoeff = (1.0 + potlp->lpObjNorm) * (1.0 + potlp->lpRHSNorm) / potlp->tau;
    double pObjVal = potlp->pObjVal * objScalCoeff;
    double dObjVal = potlp->dObjVal * objScalCoeff;
    double absGap = fabs(pObjVal - dObjVal);
    double relGap = absGap / (1.0 + fabs(pObjVal) + fabs(dObjVal));
    double pInfeas = potlp->pInfeas;
    double dInfeas = potlp->dInfeas;
    
    int logFreq = 0;
    
    if ( potlp->nIter < 10000 ) {
        logFreq = 1000;
    } else if ( potlp->nIter < 50000 ) {
        logFreq = 5000;
    } else if ( potlp->nIter < 200000 ) {
        logFreq = 10000;
    } else {
        logFreq = 50000;
    }
    
    double relFeasTol = 1e-05;
    double relOptTol = 1e-04;
    
    int *intInfo = NULL;
    if ( relGap < relOptTol && pInfeas < relFeasTol && dInfeas < relFeasTol ) {
        intInfo = (int *) info;
        *intInfo = 1;
    }
    
    if ( potlp->nIter >= potlp->intParams[INT_PARAM_MAXITER]) {
        intInfo = (int *) info;
        *intInfo = 1;
    }
    
    if ( potlp->nIter % logFreq == 0 || potlp->nIter == 1 || intInfo ) {
        
        double elapsedTime = my_clock() - potlp->startT;
        if ( elapsedTime < 100.0 ) {
            printf("%8lld  %10.3e  %10.3e  %10.3e  %10.3e  %10.3e  %10.3e |%5.1f [s] \n",
                   potlp->nIter, pObjVal, dObjVal, relGap, pInfeas,
                   dInfeas, potlp->kappa / potlp->tau, elapsedTime);
        } else if ( elapsedTime < 3600  ){
            printf("%8lld  %10.3e  %10.3e  %10.3e  %10.3e  %10.3e  %10.3e |%5.1f [m] \n",
                   potlp->nIter, pObjVal, dObjVal, relGap, pInfeas,
                   dInfeas, potlp->kappa / potlp->tau,
                   elapsedTime / 60.0 );
        } else {
            printf("%8lld  %10.3e  %10.3e  %10.3e  %10.3e  %10.3e  %10.3e |%5.1f [h] \n",
                   potlp->nIter, pObjVal, dObjVal, relGap, pInfeas,
                   dInfeas, potlp->kappa / potlp->tau,
                   elapsedTime / 3600.0 );
        }
    }
    
    return;
}

#define RUIZ_DEBUG(format, info) // printf(format, info);
static pot_int LPSolverIRuizScale( potlp_solver *potlp ) {
    
    pot_int retcode = RETCODE_OK;
    
    pot_int nCol = potlp->nCol;
    pot_int nRow = potlp->nRow;
    pot_int nColMatElem = potlp->colMatBeg[nCol];
    pot_int *colMatBeg = potlp->colMatBeg;
    pot_int *colMatIdx = potlp->colMatIdx;
    
    pot_int nRuizRow = nRow + nCol + 1;
    pot_int nRuizCol = nRow + 2 * nCol + 2;
    
    double *ruizWorkDiagRow = NULL;
    double *ruizWorkDiagCol = NULL;
    double *ruizWorkColMatElem = NULL;
    double *ruizWorkColMatElemTrans = NULL;
    double *ruizWorkObjVal = NULL;
    double *ruizWorkObjValTrans = NULL;
    double *ruizWorkRhs = NULL;
    double *ruizWorkRhsTrans = NULL;
    double *ruizWorkEye = NULL;
    double ruizWorkOne = 1.0;
    
    double *ruizScalDiagRow = potlp->ruizRow;
    double *ruizScalDiagCol = potlp->ruizCol;
    
    /* Allocate memory for ruiz */
    POTLP_INIT(ruizWorkDiagRow, double, nRuizRow);
    POTLP_INIT(ruizWorkDiagCol, double, nRuizCol);
    
    POTLP_INIT(ruizWorkColMatElem, double, nColMatElem);
    POTLP_INIT(ruizWorkColMatElemTrans, double, nColMatElem);
    POTLP_INIT(ruizWorkObjVal, double, nCol);
    POTLP_INIT(ruizWorkObjValTrans, double, nCol);
    POTLP_INIT(ruizWorkRhs, double, nRow);
    POTLP_INIT(ruizWorkRhsTrans, double, nRow);
    POTLP_INIT(ruizWorkEye, double, nCol);
    
    if ( !ruizScalDiagRow || !ruizScalDiagCol || !ruizWorkColMatElem || !ruizWorkColMatElemTrans ||
         !ruizWorkObjVal || !ruizWorkObjValTrans || !ruizWorkRhs || !ruizWorkRhsTrans ||
         !ruizWorkDiagCol || !ruizWorkDiagRow || !ruizWorkEye ) {
        retcode = RETCODE_FAILED;
        goto exit_cleanup;
    }
    
    /* Initialize working data */
    POTLP_MEMCPY(ruizWorkColMatElem, potlp->colMatElem, double, nColMatElem);
    POTLP_MEMCPY(ruizWorkColMatElemTrans, potlp->colMatElem, double, nColMatElem);
    POTLP_MEMCPY(ruizWorkObjVal, potlp->lpObj, double, nCol);
    POTLP_MEMCPY(ruizWorkObjValTrans, potlp->lpObj, double, nCol);
    POTLP_MEMCPY(ruizWorkRhs, potlp->lpRHS, double, nRow);
    POTLP_MEMCPY(ruizWorkRhsTrans, potlp->lpRHS, double, nRow);
    
    for ( int i = 0; i < nCol; ++i ) {
        ruizWorkEye[i] = 1.0;
    }
    
    /* Initialize scalers */
    for ( int i = 0; i < nRuizRow; ++i ) {
        ruizScalDiagRow[i] = 1.0;
    }
    
    for ( int i = 0; i < nRuizCol; ++i ) {
        ruizScalDiagCol[i] = 1.0;
    }
    
    int maxRuizIter = potlp->intParams[INT_PARAM_MAXRUIZITER];
    
    /* Start ruiz scaling on the matrix
        [ 0   A  0  0 -b
          A'  0 -I  0  c
          b' -c' 0 -1  0 ]
    */
    
    RUIZ_DEBUG("Start Ruiz-scaling %s\n", "");
    
    for ( int i = 0; i < maxRuizIter; ++i ) {
        
        POTLP_ZERO(ruizWorkDiagRow, double, nRuizRow);
        
        /* row[0] <- max( [ A  ] )
           row[1] <- max( [ A' ] )
           row[2] <- max( [ b' -c' ] )
         */
        spMatMaxRowAbs(nCol, colMatBeg, colMatIdx, ruizWorkColMatElem, ruizWorkDiagRow);
        spMatMaxColAbs(nCol, colMatBeg, colMatIdx, ruizWorkColMatElemTrans, ruizWorkDiagRow + nRow);
        
        int maxRHSIdx = idamax(&nRow, ruizWorkRhsTrans, &potIntConstantOne);
        int maxObjIdx = idamax(&nCol, ruizWorkObjValTrans, &potIntConstantOne);
        
        double maxRHS = fabs(ruizWorkRhsTrans[maxRHSIdx]);
        double maxObj = fabs(ruizWorkObjValTrans[maxObjIdx]);
        double compElem = POTLP_MAX(maxRHS, maxObj);
        
        /* row[0] <- max( [ row[0], -b    ] )
           row[1] <- max( [ row[1], -I, c ] )
           row[2] <- max( [ row[2], -1    ] )
         */
        
        double *ruizWorkDiagRowRHS = ruizWorkDiagRow;
        for ( int j = 0; j < nRow; ++j ) {
            double rhsElem = fabs(ruizWorkRhs[j]);
            ruizWorkDiagRowRHS[j] = POTLP_MAX(rhsElem, ruizWorkDiagRowRHS[j]);
        }
        
        double *ruizWorkDiagRowObj = ruizWorkDiagRow + nRow;
        for ( int j = 0; j < nCol; ++j ) {
            double objElem = fabs(ruizWorkObjVal[j]);
            ruizWorkDiagRowObj[j] = POTLP_MAX(objElem, ruizWorkDiagRowObj[j]);
            ruizWorkDiagRowObj[j] = POTLP_MAX(ruizWorkEye[j], ruizWorkDiagRowObj[j]);
        }
        
        compElem = POTLP_MAX(compElem, ruizWorkOne);
        ruizWorkDiagRow[nRow + nCol] = compElem;
        
        /* col[0] <- max( [  A' ] )
           col[1] <- max( [  A  ] )
           col[2] <- max( [ -I  ] )
           col[3] <- max( [ -1  ] )
           col[4] <- max( [ -b; c] )
         */
        
        POTLP_ZERO(ruizWorkDiagCol, double, nRuizCol);
        spMatMaxRowAbs(nCol, colMatBeg, colMatIdx, ruizWorkColMatElemTrans, ruizWorkDiagCol);
        spMatMaxColAbs(nCol, colMatBeg, colMatIdx, ruizWorkColMatElem, ruizWorkDiagCol + nRow);
        
        double *ruizWorkDiagColEye = ruizWorkDiagCol + nRow + nCol;
        for ( int j = 0; j < nCol; ++j ) {
            ruizWorkDiagColEye[j] = ruizWorkEye[j];
        }
        
        ruizWorkDiagCol[nRow + nCol + nCol] = ruizWorkOne;
        
        maxRHSIdx = idamax(&nRow, ruizWorkRhs, &potIntConstantOne);
        maxObjIdx = idamax(&nCol, ruizWorkObjVal, &potIntConstantOne);
        maxRHS = fabs(ruizWorkRhs[maxRHSIdx]);
        maxObj = fabs(ruizWorkObjVal[maxObjIdx]);
        ruizWorkDiagCol[nRow + 2 * nCol + 1] = POTLP_MAX(maxRHS, maxObj);
        
        /* col[0] <- max( [ col[0],  b ] )
           col[1] <- max( [ col[1], -c ] )
           col[2] <- col[2]
           col[3] <- col[3],
           col[4] <- col[4] */
        
        double *ruizWorkDiagColRHS = ruizWorkDiagCol;
        for ( int j = 0; j < nRow; ++j ) {
            double rhsElem = fabs(ruizWorkRhsTrans[j]);
            ruizWorkDiagColRHS[j] = POTLP_MAX(rhsElem, ruizWorkDiagColRHS[j]);
        }
        
        double *ruizWorkDiagColObj = ruizWorkDiagCol + nRow;
        for ( int j = 0; j < nCol; ++j ) {
            double objElem = fabs(ruizWorkObjValTrans[j]);
            ruizWorkDiagColRHS[j] = POTLP_MAX(objElem, ruizWorkDiagColRHS[j]);
        }
        
        /* Get sqrt operation */
        double maxRuizDiagDeviate = 0.0;
        double ruizDiagDeviate = 0.0;
        for ( int j = 0; j < nRuizRow; ++j ) {
            ruizWorkDiagRow[j] = sqrtl(ruizWorkDiagRow[j]);
            ruizScalDiagRow[j] = ruizScalDiagRow[j] * ruizWorkDiagRow[j];
            ruizDiagDeviate = fabs(ruizWorkDiagRow[j] - 1.0);
            maxRuizDiagDeviate = POTLP_MAX(maxRuizDiagDeviate, ruizDiagDeviate);
        }
        
        for ( int j = 0; j < nRuizCol; ++j ) {
            ruizWorkDiagCol[j] = sqrtl(ruizWorkDiagCol[j]);
            ruizScalDiagCol[j] = ruizScalDiagCol[j] * ruizWorkDiagCol[j];
            ruizDiagDeviate = fabs(ruizWorkDiagCol[j] - 1.0);
            maxRuizDiagDeviate = POTLP_MAX(maxRuizDiagDeviate, ruizDiagDeviate);
        }
        
        RUIZ_DEBUG("Ruiz Deviation %e \n", maxRuizDiagDeviate);
        
        if ( maxRuizDiagDeviate < 1e-08 ) {
            RUIZ_DEBUG("Ruiz Successfully Ends in %d iterations \n", i);
            break;
        }
        
        /* Row scaling */
        spMatRowScal(nCol, colMatBeg, colMatIdx, ruizWorkColMatElem, ruizWorkDiagRowRHS);
        vvrscl(&nRow, ruizWorkDiagRowRHS, ruizWorkRhs);
        
        spMatColScal(nCol, colMatBeg, colMatIdx, ruizWorkColMatElemTrans, ruizWorkDiagRowObj);
        vvrscl(&nCol, ruizWorkDiagRowObj, ruizWorkEye);
        vvrscl(&nCol, ruizWorkDiagRowObj, ruizWorkObjVal);
        
        double *ruizWorkDiagRowCompl = ruizWorkDiagRow + nRow + nCol;
        rscl(&nRow, ruizWorkDiagRowCompl, ruizWorkRhsTrans, &potIntConstantOne);
        rscl(&nCol, ruizWorkDiagRowCompl, ruizWorkObjValTrans, &potIntConstantOne);
        ruizWorkOne = ruizWorkOne / (*ruizWorkDiagRowCompl);
    
        /* Column scaling */
        spMatRowScal(nCol, colMatBeg, colMatIdx, ruizWorkColMatElemTrans, ruizWorkDiagColRHS);
        vvrscl(&nRow, ruizWorkDiagColRHS, ruizWorkRhsTrans);
        
        spMatColScal(nCol, colMatBeg, colMatIdx, ruizWorkColMatElem, ruizWorkDiagColObj);
        vvrscl(&nRow, ruizWorkDiagColObj, ruizWorkObjValTrans);
        
        vvrscl(&nCol, ruizWorkDiagColEye, ruizWorkEye);
        
        ruizWorkOne = ruizWorkOne / ruizWorkDiagCol[nRow + nCol + nCol];
        
        double *ruizWorkDiagColCompl = ruizWorkDiagCol + nRow + nCol + nCol + 1;
        rscl(&nRow, ruizWorkDiagColCompl, ruizWorkRhs, &potIntConstantOne);
        rscl(&nCol, ruizWorkDiagColCompl, ruizWorkObjVal, &potIntConstantOne);
    }
    
    RUIZ_DEBUG("Ruiz-scaling Ends %s\n", "");
    
    for ( int i = 0; i < nRuizRow; ++i ) {
        ruizScalDiagRow[i] = 1.0 / ruizScalDiagRow[i];
    }
    
    for ( int i = 0; i < nRuizCol; ++i ) {
        ruizScalDiagCol[i] = 1.0 / ruizScalDiagCol[i];
    }
    
exit_cleanup:
    
    POTLP_FREE(ruizWorkDiagRow);
    POTLP_FREE(ruizWorkDiagCol);
    POTLP_FREE(ruizWorkColMatElem);
    POTLP_FREE(ruizWorkColMatElemTrans);
    POTLP_FREE(ruizWorkObjVal);
    POTLP_FREE(ruizWorkObjValTrans);
    POTLP_FREE(ruizWorkRhs);
    POTLP_FREE(ruizWorkRhsTrans);
    POTLP_FREE(ruizWorkEye);
    
    return retcode;
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
    
    for ( int i = 0; i < potlp->nRow; ++i ) {
        rhsOneNorm += fabs(lpRHS[i]);
    }
    
    potlp->lpObjNorm = objOneNorm;
    potlp->lpRHSNorm = rhsOneNorm;
    
    double objScaler = objOneNorm + 1.0;
    double rhsScaler = rhsOneNorm + 1.0;
    
    for ( int i = 0; i < potlp->nCol; ++i ) {
        lpObj[i] = lpObj[i] / objScaler;
    }
    
    for ( int i = 0; i < potlp->nRow; ++i ) {
        lpRHS[i] = lpRHS[i] / rhsScaler;
    }
    
    return;
}

static void LPSovlerIPrintLPStats( potlp_solver *potlp ) {
    
    printf("\nOptimizing an LP with %d variables and %d constraints. \n", potlp->nCol, potlp->nRow);
    printf("LP Data norm |b| = %5.2e |c| = %5.2e. \n", potlp->lpRHSNorm, potlp->lpObjNorm);
    
    return;
}

static void LPSolverIRetrieveSolution( potlp_solver *potlp ) {
    
    pot_solver *pot = potlp->potIterator;
    
    int nCol = potlp->nCol;
    int nRow = potlp->nRow;
    
    double *scalVals = potlp->scalVals;
    potVecExport(pot->xVec, scalVals);
    vvscl(&pot->xVec->n, potlp->ruizCol, scalVals);
    
    double *lpObj = potlp->lpObj;
    double *lpRHS = potlp->lpRHS;
    
    double *rowDual = scalVals;
    double *colVal = rowDual + nRow;
    double *colDual = colVal + nCol;
    double tau = *( colDual + nCol + 1 );
    
    double primalScaler = tau / (potlp->lpRHSNorm + 1.0);
    double dualScaler = tau / (potlp->lpObjNorm + 1.0);
    
    /* Scale back */
    rscl(&nCol, &primalScaler, colVal, &potIntConstantOne);
    rscl(&nCol, &dualScaler, colDual, &potIntConstantOne);
    rscl(&nRow, &dualScaler, rowDual, &potIntConstantOne);
    
    /* Verify errors */
    double pObjVal = dot(&nCol, colVal, &potIntConstantOne, lpObj, &potIntConstantOne);
    double dObjVal = dot(&nRow, rowDual, &potIntConstantOne, lpRHS, &potIntConstantOne);
    pObjVal = pObjVal * (potlp->lpObjNorm + 1.0);
    dObjVal = dObjVal * (potlp->lpRHSNorm + 1.0);
    
    double *pRes = potlp->pRes;
    double *dRes = potlp->dRes;
    
    for ( int i = 0; i < nRow; ++i ) {
        pRes[i] = - lpRHS[i] * (potlp->lpRHSNorm + 1.0);
    }
    
    for ( int i = 0; i < nCol; ++i ) {
        dRes[i] = colDual[i] - lpObj[i] * (potlp->lpObjNorm + 1.0);
    }
    
    /* Primal infeasibility b - A * x  */
    spMatAxpy(nCol, potlp->colMatBeg, potlp->colMatIdx, potlp->colMatElem, 1.0, colVal, pRes);
    
    /* Dual infeasibility s - c + A' * y */
    spMatATxpy(nCol, potlp->colMatBeg, potlp->colMatIdx, potlp->colMatElem, 1.0, rowDual, dRes);
    
    double pInfeas = nrm2(&nRow, pRes, &potIntConstantOne);
    double relpInfeas = pInfeas / (potlp->lpRHSNorm + 1.0);
    double dInfeas = nrm2(&nCol, dRes, &potIntConstantOne);
    double reldInfeas = dInfeas / (potlp->lpObjNorm + 1.0);
    double compGap = fabs(pObjVal - dObjVal);
    double rGap = compGap / ( 1 + fabs(pObjVal) + fabs(dObjVal) );
    
    potlp->pInfeas = pInfeas;
    potlp->dInfeas = dInfeas;
    potlp->complGap = compGap;
    
    printf("\nLP Solution statistic \n");
    printf("    pObj: %+8.2e   Obj: %+8.2e \n", pObjVal, dObjVal);
    printf("    pInf:  %8.2e   Rel:  %8.2e \n", pInfeas, relpInfeas);
    printf("    dInf:  %8.2e   Rel:  %8.2e \n", dInfeas, reldInfeas);
    printf("    Gap :  %8.2e   Rel:  %8.2e \n", compGap, rGap);
    
    POTLP_MEMCPY(potlp->colVal, colVal, double, nCol);
    POTLP_MEMCPY(potlp->colDual, colDual, double, nCol);
    POTLP_MEMCPY(potlp->rowDual, rowDual, double, nRow);
    
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
    
    POTLP_ZERO(potlp, potlp_solver, 1);
    
    POTLP_MEMCPY(potlp->dblParams, defaultDblParam, double, NUM_DBL_PARAM);
    POTLP_MEMCPY(potlp->intParams, defaultIntParam, int, NUM_INT_PARAM);
    
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
    potlp->nRow = nRow;
    
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
    pot_int nRow = potlp->nRow;
    pot_int nNz = Ap[nCol];
    
    POTLP_INIT(potlp->colMatBeg, pot_int, nCol + 1);
    POTLP_INIT(potlp->colMatIdx, pot_int, nNz);
    POTLP_INIT(potlp->colMatElem, double, nNz);
    
    POTLP_INIT(potlp->lpObj, double, nCol);
    POTLP_INIT(potlp->lpRHS, double, nRow);
    POTLP_INIT(potlp->pdcRes, double, nRow + nCol + 1);
    
    POTLP_INIT(potlp->ruizRow, double, nRow + nCol + 1);
    POTLP_INIT(potlp->ruizCol, double, nRow + 2 * nCol + 2);
    POTLP_INIT(potlp->scalVals, double, nRow + 2 * nCol + 2);
    
    POTLP_INIT(potlp->colVal, double, nCol);
    POTLP_INIT(potlp->colDual, double, nCol);
    POTLP_INIT(potlp->rowDual, double, nRow);
    
    if ( !potlp->colMatBeg || !potlp->colMatIdx || !potlp->colMatElem ||
         !potlp->lpObj || !potlp->lpRHS || !potlp->pdcRes ||
         !potlp->ruizCol || !potlp->ruizRow  || !potlp->scalVals ||
         !potlp->colVal || !potlp->colDual || !potlp->rowDual ) {
        retcode = RETCODE_FAILED;
        goto exit_cleanup;
    }
    
    potlp->pRes = potlp->pdcRes;
    potlp->dRes = potlp->pRes + nRow;
    potlp->cplRes = potlp->dRes + nCol;
    
    /* Copy data */
    POTLP_MEMCPY(potlp->colMatBeg, Ap, pot_int, nCol + 1);
    POTLP_MEMCPY(potlp->colMatIdx, Ai, pot_int, nNz);
    POTLP_MEMCPY(potlp->colMatElem, Ax, double, nNz);
    POTLP_MEMCPY(potlp->lpObj, lpObj, double, nCol);
    POTLP_MEMCPY(potlp->lpRHS, lpRHS, double, nRow);
    
exit_cleanup:
    return retcode;
}

extern pot_int LPSolverOptimize( potlp_solver *potlp ) {
    
    pot_int retcode = RETCODE_OK;
    
    LPSolverIScale(potlp);
    LPSovlerIPrintLPStats(potlp);
    
    POT_CALL(LPSolverIRuizScale(potlp));
    POT_CALL(potReductionSolve(potlp->potIterator));
    LPSolverIRetrieveSolution(potlp);
    
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
    
    POTLP_FREE(potlp->ruizRow);
    POTLP_FREE(potlp->ruizCol);
    POTLP_FREE(potlp->scalVals);
    POTLP_FREE(potlp->pdcRes);
    
    POTLP_FREE(potlp->colVal);
    POTLP_FREE(potlp->colDual);
    POTLP_FREE(potlp->rowDual);
    
    potConstrMatDestroy(&potlp->potConstrMat);
    potObjFDestroy(&potlp->potObjF);
    potLPDestroy(&potlp->potIterator);
    
    POTLP_ZERO(potlp, potlp_solver, 1);
    
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
