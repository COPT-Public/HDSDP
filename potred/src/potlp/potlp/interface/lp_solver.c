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

#include <math.h>


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
    double  *Ax = potlp->colMatVal;
    
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
    for ( int i = 0, j; i < nCol; ++i ) {
        for ( j = Ap[i]; j < Ap[i + 1]; ++j ) {
            pRes[Ai[j]] += colVal[i] * Ax[j];
        }
    }
    
    /* Dual infeasibility */
    for ( int i = 0, j; i < nCol; ++i ) {
        double aTy = 0.0;
        for ( j = Ap[i]; j < Ap[i + 1]; ++j ) {
            aTy += rowDual[Ai[j]] * Ax[j];
        }
        dRes[i] -= aTy;
        dRes[i] -= colSlack[i];
    }
    
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
    double  *Ax = potlp->colMatVal;
    
    /* Gradient of row dual */
    axpy(&nRow, &compl, lpRHS, &potIntConstantOne, gRowDual, &potIntConstantOne);
    
    for ( int i = 0, j; i < nCol; ++i ) {
        for ( j = Ap[i]; j < Ap[i + 1]; ++j ) {
            gRowDual[Ai[j]] -= dRes[i] * Ax[j];
        }
    }
    
    /* Gradient of column value */
    double ncompl = -compl;
    axpy(&nCol, &ncompl, lpObj, &potIntConstantOne, gColVal, &potIntConstantOne);
    for ( int i = 0, j; i < nCol; ++i ) {
        double aTp = 0.0;
        for ( j = Ap[i]; j < Ap[i + 1]; ++j ) {
            aTp += pRes[Ai[j]] * Ax[j];
        }
        gColVal[i] += aTp;
    }
    
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
    double  *Ax = potlp->colMatVal;
    
    double pObjVal = dot(&nCol, colVal, &potIntConstantOne, potlp->lpObj, &potIntConstantOne);
    double dObjVal = dot(&nRow, rowDual, &potIntConstantOne, potlp->lpRHS, &potIntConstantOne);
    
    for ( int i = 0; i < nRow; ++i ) {
        pRes[i] = -lpRHS[i] * tauVar[0];
    }
    
    for ( int i = 0; i < nCol; ++i ) {
        dRes[i] = lpObj[i] * tauVar[0];
    }
    
    /* Primal infeasibility */
    for ( int i = 0, j; i < nCol; ++i ) {
        for ( j = Ap[i]; j < Ap[i + 1]; ++j ) {
            pRes[Ai[j]] += colVal[i] * Ax[j];
        }
    }
    
    /* Dual infeasibility */
    for ( int i = 0, j; i < nCol; ++i ) {
        double aTy = 0.0;
        for ( j = Ap[i]; j < Ap[i + 1]; ++j ) {
            aTy += rowDual[Ai[j]] * Ax[j];
        }
        dRes[i] -= aTy;
        dRes[i] -= colSlack[i];
    }
    
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
    printf("%10.3e  %10.3e  %10.3e  %10.3e  %10.3e \n",
           potlp->pObjVal, potlp->dObjVal, potlp->pInfeas / potlp->tau,
           potlp->dInfeas / potlp->tau, potlp->kappa / potlp->tau);
    
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
    
    POT_CALL(potLPCreate(&potlp->potIterator));
    POT_CALL(potConstrMatCreate(&potlp->potConstrMat));
    POT_CALL(potObjFCreate(&potlp->potObjF));
    
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
    POTLP_INIT(potlp->colMatVal, double, nNz);
    
    POTLP_INIT(potlp->lpObj, double, nCol);
    POTLP_INIT(potlp->lpRHS, double, nRow);
    POTLP_INIT(potlp->pdcRes, double, nCol + nRow + 1);
    
    POTLP_INIT(potlp->auxArray, double, nCol + 2 * nRow + 2);
    
    if ( !potlp->colMatBeg || !potlp->colMatIdx || !potlp->colMatVal ||
         !potlp->lpObj || !potlp->lpRHS || !potlp->pdcRes || !potlp->auxArray ) {
        retcode = RETCODE_FAILED;
        goto exit_cleanup;
    }
    
    potlp->pRes = potlp->pdcRes;
    potlp->dRes = potlp->pRes + nRow;
    potlp->cplRes = potlp->dRes + nCol;
    
    /* Copy data */
    memcpy(potlp->colMatBeg, Ap, sizeof(pot_int) * (nCol + 1));
    memcpy(potlp->colMatIdx, Ai, sizeof(pot_int) * nNz);
    memcpy(potlp->colMatVal, Ax, sizeof(double) * nNz);
    memcpy(potlp->lpObj, lpObj, sizeof(double) * nCol);
    memcpy(potlp->lpRHS, lpRHS, sizeof(double) * nRow);
    
exit_cleanup:
    return retcode;
}

extern pot_int LPSolverOptimize( potlp_solver *potlp ) {
    
    pot_int retcode = RETCODE_OK;
    
    if (0) {
        LPSolverIScale(potlp);
    }
    
    retcode = potReductionSolve(potlp->potIterator);
    
exit_cleanup:
    return retcode;
}

extern void LPSolverDestroy( potlp_solver **ppotlp ) {
    
    if ( !ppotlp ) {
        return;
    }
    
    potlp_solver *potlp = *ppotlp;
    
    POTLP_FREE(potlp->colMatBeg);
    POTLP_FREE(potlp->colMatIdx);
    POTLP_FREE(potlp->colMatVal);
    
    POTLP_FREE(potlp->lpRHS);
    POTLP_FREE(potlp->lpObj);
    
    POTLP_FREE(potlp->ruizCol);
    POTLP_FREE(potlp->ruizRow);
    
    POTLP_FREE(potlp->pdcRes);
    POTLP_FREE(potlp->auxArray);
    
    potConstrMatDestroy(&potlp->potConstrMat);
    potObjFDestroy(&potlp->potObjF);
    potLPDestroy(&potlp->potIterator);
    
    memset(potlp, 0, sizeof(potlp_solver));
    return;
}
