/** @file lpdata.c
 *  @brief Implement the LP data interface for potential reduction
 *
 * @TODO: Add more detailed comments
 */

#include "lp_solver.h"
#include "pot_solver.h"
#include "pot_structs.h"
#include "pot_utils.h"
#include "pot_vector.h"
#include "pot_constr_mat.h"
#include "pot_objfunc.h"

#include <math.h>

typedef struct {
    
    pot_int nCol; ///< Number of LP variables
    pot_int nConstr; ///< Number of constraints
    
    pot_int *colMatBeg;
    pot_int *colMatIdx;
    double  *colMatVal;
    
    double *lpRHS;
    double *lpObj;
    
    double lpRHSNorm;
    double lpObjNorm;
    
    double *ruizCol;
    double *ruizRow;
    
    double *pdcRes;
    double *pRes;
    double *dRes;
    double *cplRes;
    
    double pInfeas;
    double dInfeas;
    double complGap;
    
    pot_solver *potIterator;
    pot_constr_mat *potConstrMat;
    pot_fx *potObjF;
    
} potlp_solver;


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
    potVecConeAddConstant(xVec, -eTx);
    
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

/* Objective methods */
static double potLpObjFImplVal( void *ObjFData, pot_vec *xVec ) {
    
    
    return 0.0;
}

static void potLpObjFImplGrad( void *ObjFData, pot_vec *xVec, pot_vec *fGrad ) {
    
    return;
}

static void potLpObjFImplHess( void *objFData, pot_vec *xVec, double *fHess ) {
    
    return;
}

static void potLpObjFImplHVec( void *objFData, pot_vec *xVec, pot_vec *fHvec ) {
    
    return;
}

static void potLpObjFImplMonitor( void *objFData, void *info ) {
    
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
    
    if ( !*ppotlp ) {
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
    
    /* x;  s; kappa; tau*/
    pot_int coneDim = nCol + nRow + 2;
    
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
    
    potlp->pRes = potlp->pdcRes;
    potlp->dRes = potlp->pRes + nRow;
    potlp->cplRes = potlp->pRes + nCol;
    
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
    
    potConstrMatDestroy(&potlp->potConstrMat);
    potObjFDestroy(&potlp->potObjF);
    potLPDestroy(&potlp->potIterator);
    
    memset(potlp, 0, sizeof(potlp_solver));
    return;
    
}
