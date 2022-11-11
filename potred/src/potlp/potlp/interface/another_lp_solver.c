#include "another_lp_solver.h"
#include "pot_solver.h"
#include "pot_structs.h"
#include "pot_utils.h"
#include "pot_vector.h"
#include "pot_constr_mat.h"
#include "pot_objfunc.h"
#include "vec_mat.h"

#include <sys/time.h>
#include <math.h>

#ifdef MATLAB_MEX_FILE
#include <stdbool.h>
#ifdef __cplusplus
extern "C" {
    bool utIsInterruptPending();
}
#else
extern bool utIsInterruptPending();
#endif
#endif

static int isCtrlC = 0;

static void monitorCtrlC( int sigNum ) {
    
    isCtrlC = 1;
    return;
}

static void POT_FNAME(potLpResidualScal) ( potlp_solver *potlp ) {
    
    double *pRes = potlp->pRes;
    double *dRes = potlp->dRes;
    
    scal(&potlp->nRow, &potlp->pResOmega, pRes, &potIntConstantOne);
    scal(&potlp->nCol, &potlp->dResOmega, dRes, &potIntConstantOne);
    *potlp->cplRes = (*potlp->cplRes) * potlp->cplResOmega;
    
    return;
}

#define REWEIGHT_DEBUG(s, p, d, c) // printf("%s: pInf = %3.3f dInf = %3.3f Gap = %3.3f \n", s, p, d, c)
static void POT_FNAME(potLpReWeight) ( potlp_solver *potlp ) {
    

    /* Implement the weighted objective Heuristic */
    double rRate = potlp->dblParams[DBL_IPARAM_RESTARTRATE];
    double rMax = potlp->dblParams[DBL_IPARAM_RESTARTMAX];
    
    double pOmega = potlp->pResOmega;
    double dOmega = potlp->dResOmega;
    double cOmega = potlp->cplResOmega;
    
    double pInfeas = potlp->pInfeasRel;
    double dInfeas = potlp->dInfeasRel;
    double cInfeas = potlp->complGapRel;
    
//    double pInfeas = potlp->dInfeas;
//    double dInfeas = potlp->pInfeas;
//    double cInfeas = potlp->complGap;
    
    REWEIGHT_DEBUG("Before", pOmega, dOmega, cOmega);
    
    /* Are the residuals stuck? */
    int pStuck = 0, pFast = 0;
    int dStuck = 0, dFast = 0;
    int cStuck = 0, cFast = 0;
    
    double minInfeas = 0.0;
    minInfeas = POTLP_MIN(pInfeas, dInfeas);
    minInfeas = POTLP_MIN(minInfeas, cInfeas);
    
    if ( rRate > 0 ) {
        /* Heuristic update */
        if ( pInfeas > 5 * minInfeas ) { pStuck = 1; }
        else if ( pInfeas < 1.1 * minInfeas ) { pFast = 1; }
        if ( dInfeas > 5 * minInfeas ) { dStuck = 1; }
        else if ( dInfeas < 1.1 * minInfeas ) { dFast = 1; }
        if ( cInfeas > 5 * minInfeas ) { cStuck = 1; }
        else if ( cInfeas < 1.1 * minInfeas ) { cFast = 1; }
        
        /* Do if stuck */
        if ( pStuck ) { pOmega = POTLP_MIN(pOmega * rRate, rMax); }
        if ( dStuck ) { dOmega = POTLP_MIN(dOmega * rRate, rMax); }
        if ( cStuck ) { cOmega = POTLP_MIN(cOmega * rRate, rMax); }
        
        /* Do if fast */
        if ( pFast ) { };
        if ( dFast ) { };
        if ( cFast ) { };
        
        /* Normalize */
        double minOmega = POTLP_MIN(pOmega, dOmega);
        minOmega = POTLP_MIN(minOmega, cOmega);

        pOmega = pOmega / minOmega;
        dOmega = dOmega / minOmega;
        cOmega = cOmega / minOmega;
        
    } else {
        /* Make residuals the same */
        pOmega = pInfeas / minInfeas;
        dOmega = dInfeas / minInfeas;
        cOmega = cInfeas / minInfeas;
        
        pOmega = pOmega * pOmega;
        dOmega = dOmega * dOmega;
        cOmega = cOmega * cOmega;
        
        pOmega = POTLP_MIN(pOmega, rMax);
        dOmega = POTLP_MIN(dOmega, rMax);
        cOmega = POTLP_MIN(cOmega, rMax);
    }
    
    REWEIGHT_DEBUG("After", pOmega, dOmega, cOmega);
    
    potlp->pResOmega = pOmega;
    potlp->dResOmega = dOmega;
    potlp->cplResOmega = cOmega;
    
    potReductionRestart(potlp->potIterator);
    
    return;
}

static pot_int POT_FNAME(potLpNewtonStep) ( potlp_solver *potlp ) {
    
    pot_int retcode = RETCODE_OK;
    double simp = 2 * potlp->nCol + 2;
    POT_CALL(LpNewtonOneStep(potlp->ipm, potlp->lpObj, potlp->lpRHS, potlp->colMatBeg,
                             potlp->colMatIdx, potlp->colMatElem, potlp->colVal, potlp->rowDual,
                             potlp->colDual, &potlp->kappa, &potlp->tau, potlp->pRes, potlp->dRes,
                             potlp->pObjVal * potlp->tau, potlp->dObjVal * potlp->tau, simp));
    
    pot_solver *pot = potlp->potIterator;
    
    /* TODO: Write solution back */
    potVecCopy(pot->xVec, pot->xVecOld);
    double *xNew = pot->xVec->x;
    
    int nCol = potlp->nCol;
    int nRow = potlp->nRow;
    
    POTLP_MEMCPY(xNew, potlp->rowDual, double, nRow);
    xNew += nRow;
    
    POTLP_MEMCPY(xNew, potlp->colVal, double, nCol);
    xNew += nCol;
    
    POTLP_MEMCPY(xNew, potlp->colDual, double, nCol);
    xNew += nCol;
    
    xNew[0] = potlp->kappa;
    xNew[1] = potlp->tau;
    
    LPQMatScal(potlp->potQMatrix, pot->xVec->x);
    potReductionRestart(pot);
    
exit_cleanup:
    return retcode;
}

/* Constraint methods */
static void POT_FNAME(potLpConstrMatImplPrepareX)( void *AMatData, pot_vec *xInit ) {
    
    potVecReset(xInit);

    for ( int i = xInit->n - xInit->ncone; i < xInit->n; ++i ) {
        xInit->x[i] = 1.0;
    }
    
    return;
}

static void POT_FNAME(potLpConstrMatImplProject)( void *AMatData, pot_vec *xVec ) {
    
    double eTx = potVecSumCone(xVec);
    potVecConeAddConstant(xVec, -eTx / xVec->ncone);
    
    return;
}

static void POT_FNAME(potLpConstrMatImplScalProject)( void *AMatData, pot_vec *xVec, pot_vec *yVec ) {
    
    double xTy = potVecSumScalCone(xVec, yVec);
    potVecConeAxpy(-xTy, xVec, yVec);
    
    return;
}

static void POT_FNAME(potLpConstrMatImplMonitor)( void *AMatData, void *info ) {
    
    return;
}

static void POT_FNAME(potLpObjFISetupRes)( potlp_solver *potlp, pot_vec *xVec ) {
    
    LPQMatMultiply(potlp->potQMatrix, potlp->isColBasic, xVec->x, potlp->pdcRes);
    POT_FNAME(potLpResidualScal)(potlp);
    
    return;
}

static void POT_FNAME(potLpObjFISetupGrad)( potlp_solver *potlp, pot_vec *gVec ) {
    
    POT_FNAME(potLpResidualScal)(potlp);
    LPQMatTransMultiply(potlp->potQMatrix, potlp->isColBasic, potlp->pdcRes, gVec->x);
    
    return;
}

static void POT_FNAME(potLpObjFISetupHVec)( potlp_solver *potlp, pot_vec *xVec, pot_vec *fHvec ) {
    
    LPQMatMultiply(potlp->potQMatrix, potlp->isColBasic, xVec->x, potlp->pdcRes);
    POT_FNAME(potLpResidualScal)(potlp);
    POT_FNAME(potLpResidualScal)(potlp);
    LPQMatTransMultiply(potlp->potQMatrix, potlp->isColBasic, potlp->pdcRes, fHvec->x);
    
    return;
}

/* Objective methods */
static double POT_FNAME(potLpObjFImplVal)( void *objFData, pot_vec *xVec ) {
    
    potlp_solver *potlp = (potlp_solver *) objFData;
    POT_FNAME(potLpObjFISetupRes)(potlp, xVec);
    int nRes = potlp->nRow + potlp->nCol + 1;
    double objFVal = nrm2(&nRes, potlp->pdcRes, &potIntConstantOne);
    
    return 0.5 * objFVal * objFVal;
}

static void POT_FNAME(potLpObjFImplGrad)( void *objFData, pot_vec *xVec, pot_vec *fGrad ) {
    
    potVecReset(fGrad);
    potlp_solver *potlp = (potlp_solver *) objFData;
    POT_FNAME(potLpObjFISetupGrad)(potlp, fGrad);
    
    return;
}

static void POT_FNAME(potLpObjFImplHess)( void *objFData, pot_vec *xVec, double *fHess ) {
    
    assert( 0 );
    return;
}


static void POT_FNAME(potLpObjFImplHVec)( void *objFData, pot_vec *xVec, pot_vec *fHvec ) {
    
    potVecReset(fHvec);
    potlp_solver *potlp = (potlp_solver *) objFData;
    POT_FNAME(potLpObjFISetupHVec)(potlp, xVec, fHvec);
    
    return;
}

static void POT_FNAME(LPSolverIRetrieveSolution)( potlp_solver *potlp );
static void POT_FNAME(potLpObjFImplMonitor)( void *objFData, void *info ) {
    
    potlp_solver *potlp = (potlp_solver *) objFData;
    potlp->nIter += 1;
    
    if ( potlp->nIter == 1 ) {
        printf("Potential reduction log \n\n");
        printf("%8s  %10s  %10s    %12s  %10s  %10s  %10s |   "
               "T  [u]\n", "nIter", "pObj", "dObj", "Rel/Best Gap", "pInf", "dInf", "k/t");
    }
    
    int logFreq = 0;
    
    if ( potlp->nIter < 1000 ) {
        logFreq = 50;
    } else if ( potlp->nIter < 5000 ) {
        logFreq = 500;
    } else if ( potlp->nIter < 20000 ) {
        logFreq = 1000;
    } else {
        logFreq = 5000;
    }
    
    int *intInfo = NULL;
    
#ifdef MATLAB_MEX_FILE
    if (utIsInterruptPending()) {
        monitorCtrlC(0);
    }
#endif
    
    if ( isCtrlC ) {
        potlp->Lpstatus = POTLP_USER_INTERRUPT;
        intInfo = (int *) info;
        *intInfo = 1;
        isCtrlC = 0; /* Reset the ctrl c status */
    }
    
    if ( potlp->nIter % logFreq == 0 || potlp->nIter == 1 || intInfo || potlp->potIterator->useCurvature ) {
        
        POT_FNAME(LPSolverIRetrieveSolution)(potlp);
        
        /* Interior point solver */
        int useIPM = 1;
        if ( useIPM ) {
            int icode = POT_FNAME(potLpNewtonStep)(potlp);
            if ( icode != RETCODE_OK ) {
                potlp->Lpstatus = POTLP_INTERNAL_ERROR;
                intInfo = (int *) info;
                *intInfo = 1;
                return;
            }
        }
        
        double relOptTol = potlp->dblParams[DBL_PARAM_RELOPTTOL];
        double relFeasTol = potlp->dblParams[DBL_PARAM_RELFEASTOL];
        double elapsedTime = potUtilGetTimeStamp() - potlp->startT;
        
        double pObjVal = potlp->pObjVal;
        double dObjVal = potlp->dObjVal;
        double pInfeas = potlp->pInfeasRel;
        double dInfeas = potlp->dInfeasRel;
        double relGap = potlp->complGapRel;
        
        double pObjBestNew = potlp->pObjBest;
        double dObjBestNew = potlp->dObjBest;
        
#define Gap(x, y) fabs((x) - (y)) / (fabs(x) + fabs(y) + 1.0);
        
        if ( pInfeas < relFeasTol ) {
            pObjBestNew = pObjVal;
        }
        
        if ( dInfeas < relFeasTol ) {
            dObjBestNew = dObjVal;
        }
        
        double relBestGap = Gap(potlp->pObjBest, potlp->dObjBest);
        double oldPNewDGap = Gap(potlp->pObjBest, dObjBestNew);
        double oldDNewPGap = Gap(pObjBestNew, potlp->dObjBest);
        double newPNewDGap = Gap(pObjBestNew, dObjBestNew);
        
        double minGap = POTLP_MIN(relBestGap, oldPNewDGap);
        minGap = POTLP_MIN(minGap, oldDNewPGap);
        minGap = POTLP_MIN(minGap, newPNewDGap);
        
        if ( minGap == oldPNewDGap ) {
            potlp->dObjBest = dObjBestNew;
        } else if ( minGap == oldDNewPGap ) {
            potlp->pObjBest = pObjBestNew;
        } else if ( minGap == newPNewDGap ) {
            potlp->dObjBest = dObjBestNew;
            potlp->pObjBest = pObjBestNew;
        } else {
            /* No change in gap */
        }
                
        if ( minGap < relOptTol && pInfeas < relFeasTol && dInfeas < relFeasTol ) {
            potlp->Lpstatus = POTLP_OPTIMAL;
            intInfo = (int *) info;
            *intInfo = 1;
        }
        
        if ( elapsedTime >= potlp->dblParams[DBL_PARAM_TIMELIMIT] ) {
            potlp->Lpstatus = POTLP_TIMELIMIT;
            intInfo = (int *) info;
            *intInfo = 1;
        }
        
        if ( potlp->nIter >= potlp->intParams[INT_PARAM_MAXITER] ) {
            potlp->Lpstatus = POTLP_MAXITER;
            intInfo = (int *) info;
            *intInfo = 1;
        }
        
        if ( potlp->kappa > 1e+08 * potlp->tau && potlp->potIterator->fVal < 1e-06 ) {
            potlp->Lpstatus = POTLP_INFEAS_OR_UNBOUNDED;
            intInfo = (int *) info;
            *intInfo = 1;
        }
        
        if ( elapsedTime < 100.0 ) {
            printf("%8lld  %10.3e  %10.3e  %5.1e %5.1e %10.3e  %10.3e  %10.3e |%5.1f [s] \n",
                   potlp->nIter, pObjVal, dObjVal, relGap, minGap, pInfeas,
                   dInfeas, potlp->kappa / potlp->tau, elapsedTime);
        } else if ( elapsedTime < 3600.0  ){
            printf("%8lld  %10.3e  %10.3e  %5.1e %5.1e %10.3e  %10.3e  %10.3e |%5.1f [m] \n",
                   potlp->nIter, pObjVal, dObjVal, relGap, minGap, pInfeas,
                   dInfeas, potlp->kappa / potlp->tau,
                   elapsedTime / 60.0 );
        } else {
            printf("%8lld  %10.3e  %10.3e  %5.1e %5.1e %10.3e  %10.3e  %10.3e |%5.1f [h] \n",
                   potlp->nIter, pObjVal, dObjVal, relGap, minGap, pInfeas,
                   dInfeas, potlp->kappa / potlp->tau,
                   elapsedTime / 3600.0 );
        }
#ifdef MATLAB_MEX_FILE
        /* Let MATLAB show log immediately */
        mexEvalString("drawnow;");
#endif
    }
    
    int resiFreq = 500;
    if ( potlp->intParams[INT_PARAM_RSCALFREQ] > 0 ) {
        resiFreq = potlp->intParams[INT_PARAM_RSCALFREQ];
    }
    
    if ( potlp->nIter % resiFreq == 0 && !potlp->potIterator->useCurvature ) {
        POT_FNAME(potLpReWeight)(potlp);
    }
    
    return;
}

static void POT_FNAME(LPSolverIScale)( potlp_solver *potlp ) {
    
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
    potlp->objScaler = 1.0;
    potlp->rhsScaler = 1.0;
    
    if ( potlp->intParams[INT_PARAM_COEFSCALE] ) {
        
        int iMaxAbsb = idamax(&potlp->nRow, lpRHS, &potIntConstantOne);
        int iMaxAbsc = idamax(&potlp->nCol, lpObj, &potIntConstantOne);
        double maxAbsb = fabs(lpRHS[iMaxAbsb]);
        double maxAbsc = fabs(lpObj[iMaxAbsc]);
        
        double scal = POTLP_MAX(maxAbsb, maxAbsc);
        
        if ( scal > 1e+04 ) {
            scal = 1e+04;
        }
        
        if ( scal < 1e-04 ) {
            scal = 1.0;
        }
        
        double objScaler = scal; // objOneNorm + 1.0;
        double rhsScaler = scal; // rhsOneNorm + 1.0;
        
        potlp->objScaler = objScaler;
        potlp->rhsScaler = rhsScaler;
        
        for ( int i = 0; i < potlp->nCol; ++i ) {
            lpObj[i] = lpObj[i] / objScaler;
        }

        for ( int i = 0; i < potlp->nRow; ++i ) {
            lpRHS[i] = lpRHS[i] / rhsScaler;
        }
    }
    
    return;
}

#if 0
/* Test scaling from SCS. The original LP data is destroyed */
static pot_int POT_FNAME(LPSolverIScalInplace)( potlp_solver *potlp ) {
    
    pot_int retcode = RETCODE_OK;
    
    printf("\n[Warning!] Inplace scaling is on. The original LP is destroyed. \n");
    potlp->intParams[INT_PARAM_MAXRUIZITER] = 0;
    potlp->intParams[INT_PARAM_COEFSCALE] = 0;
    
    pot_int nRow = potlp->nRow;
    pot_int nCol = potlp->nCol;
    
    pot_int *colMatBeg = potlp->colMatBeg;
    pot_int *colMatIdx = potlp->colMatIdx;
    double *colMatElem = potlp->colMatElem;
    double *lpObj = potlp->lpObj;
    double *lpRHS = potlp->lpRHS;
    
    double *inpScalWorkRow = NULL;
    double *inpScalWorkCol = NULL;
    
    POTLP_INIT(inpScalWorkRow, double, nRow);
    POTLP_INIT(inpScalWorkCol, double, nCol);
    
    if ( !inpScalWorkRow || !inpScalWorkCol ) {
        retcode = RETCODE_FAILED;
        goto exit_cleanup;
    }
    
    /* Do 0 Ruiz */
    POT_CALL(spMatRuizScal(nRow, nCol, colMatBeg, colMatIdx, colMatElem,
                           inpScalWorkRow, inpScalWorkCol, 0));
    /* Do 1 L2 */
//    POT_CALL(spMatL2Scal(nRow, nCol, colMatBeg, colMatIdx, colMatElem,
//                         inpScalWorkRow, inpScalWorkCol));
    
    /* Scale b and c */
    vvscl(&nRow, inpScalWorkRow, lpRHS);
    vvscl(&nCol, inpScalWorkCol, lpObj);
    
    int iMaxAbsb = idamax(&nRow, lpRHS, &potIntConstantOne);
    double maxAbsb = fabs(lpRHS[iMaxAbsb]);
    int iMaxAbsc = idamax(&nCol, lpObj, &potIntConstantOne);
    double maxAbsc = fabs(lpObj[iMaxAbsc]);
    
    double scal = POTLP_MAX(maxAbsb, maxAbsc);
    
    if ( scal > 1e+04 ) {
        scal = 1e+04;
    }
    
    if ( scal < 1e-04 ) {
        scal = 1.0;
    }
    
    rscl(&nRow, &scal, lpRHS, &potIntConstantOne);
    rscl(&nCol, &scal, lpObj, &potIntConstantOne);
    
exit_cleanup:
    POTLP_FREE(inpScalWorkRow);
    POTLP_FREE(inpScalWorkCol);
    return retcode;
}

static pot_int POT_FNAME(LPSolverISetupQMatrix)( potlp_solver *potlp ) {
    
    pot_int retcode = RETCODE_OK;
    POT_CALL(LPQMatSetup(potlp->potQMatrix, potlp->nCol, potlp->nRow,
                         potlp->colMatBeg, potlp->colMatIdx, potlp->colMatElem,
                         potlp->lpObj, potlp->lpRHS));
    
exit_cleanup:
    return retcode;
}
#endif

static pot_int POT_FNAME(LPSolverIRuizScale)( potlp_solver *potlp ) {
    
    pot_int retcode = RETCODE_OK;
    POT_CALL(LPQMatRuizScal(potlp->potQMatrix, potlp->intParams[INT_PARAM_MAXRUIZITER]));
    
exit_cleanup:
    return retcode;
}

static pot_int POT_FNAME(LPSolverIL2Scale)( potlp_solver *potlp ) {
    
    pot_int retcode = RETCODE_OK;
    POT_CALL(LPQMatL2Scal(potlp->potQMatrix));
    
exit_cleanup:
    return retcode;
}

static void POT_FNAME(LPSovlerIPrintLPStats)( potlp_solver *potlp ) {
    
    printf("\nOptimizing an LP with %d variables and %d constraints. \n", potlp->nCol, potlp->nRow);
    printf("LP Data norm |b| = %5.2e |c| = %5.2e. \n", potlp->lpRHSNorm, potlp->lpObjNorm);
    
    return;
}

static void POT_FNAME(LPSolverIParamAdjust)( potlp_solver *potlp ) {
    
    pot_solver *pot = potlp->potIterator;
    
    potlp->pObjBest = POTLP_INFINITY;
    potlp->dObjBest = -POTLP_INFINITY;
    
    if ( potlp->intParams[INT_PARAM_CURVATURE] == 0 ) {
        pot->allowCurvature = 0;
    } else {
        pot->allowCurvature = 1;
    }
    
    pot->curvLimit = potlp->intParams[INT_PARAM_CURVATURE] - 1;
    pot->curvMinInterval = potlp->intParams[INT_PARAM_CURVINTERVAL];
    
    return;
}

static void POT_FNAME(LPSolverIHeurInitialize)( potlp_solver *potlp ) {
    
    int nColQ = potlp->potQMatrix->nColQ;
    int *isColBasic = potlp->isColBasic;
    
    /* Initialize basis status */
    if ( isColBasic ) {
        for ( int i = 0; i < nColQ; ++i ) {
            isColBasic[i] = 1;
        }
    }
    
    /* Initialize residual weight */
    potlp->pResOmega = 1.0;
    potlp->dResOmega = 1.0;
    potlp->cplResOmega = 1.0;
    
    return;
}

static void POT_FNAME(LPSolverIRetrieveSolution)( potlp_solver *potlp ) {
    
    pot_solver *pot = potlp->potIterator;
    
    int nCol = potlp->nCol;
    int nRow = potlp->nRow;
    
    potVecExport(pot->xVec, potlp->scalVals);
    LPQMatScalBack(potlp->potQMatrix, potlp->scalVals);
    
    double *lpObj = potlp->lpObj;
    double *lpRHS = potlp->lpRHS;
    
    double *rowDual = potlp->scalVals;
    double *colVal = rowDual + nRow;
    double *colDual = colVal + nCol;
    double kappa = *( colDual + nCol );
    double tau = *( colDual + nCol + 1 );
    
    double objScaler = potlp->objScaler;
    double rhsScaler = potlp->rhsScaler;
    
    double primalScaler = 1.0 / rhsScaler;
    double dualScaler = 1.0 / objScaler;
    
    /* Scale back */
    rscl(&nCol, &primalScaler, colVal, &potIntConstantOne);
    rscl(&nCol, &dualScaler, colDual, &potIntConstantOne);
    rscl(&nRow, &dualScaler, rowDual, &potIntConstantOne);
    
    /* Verify errors */
    double pObjVal = dot(&nCol, colVal, &potIntConstantOne, lpObj, &potIntConstantOne);
    double dObjVal = dot(&nRow, rowDual, &potIntConstantOne, lpRHS, &potIntConstantOne);
    pObjVal = pObjVal * objScaler / tau;
    dObjVal = dObjVal * rhsScaler / tau;
    
    double *pRes = potlp->pRes;
    double *dRes = potlp->dRes;
    
    double rhsScalerTau = rhsScaler * tau;
    for ( int i = 0; i < nRow; ++i ) {
        pRes[i] = - lpRHS[i] * rhsScalerTau;
    }
    
    double objScalerTau = objScaler * tau;
    for ( int i = 0; i < nCol; ++i ) {
        dRes[i] = - colDual[i] + lpObj[i] * objScalerTau;
    }
    
    /* Scaler when evaluating errors */
    objScaler = potlp->lpObjNorm + 1.0;
    rhsScaler = potlp->lpRHSNorm + 1.0;
    
    /* Primal infeasibility A * x  - b * tau */
    spMatAxpy(nCol, potlp->colMatBeg, potlp->colMatIdx, potlp->colMatElem, 1.0, colVal, pRes);
    
    /* Dual infeasibility - s + c * tau - A' * y */
    spMatATxpy(nCol, potlp->colMatBeg, potlp->colMatIdx, potlp->colMatElem, -1.0, rowDual, dRes);
    
    double pInfeas = nrm2(&nRow, pRes, &potIntConstantOne) / tau;
    double relpInfeas = pInfeas / rhsScaler;
    double dInfeas = nrm2(&nCol, dRes, &potIntConstantOne) / tau;
    double reldInfeas = dInfeas / objScaler;
    double compGap = fabs(pObjVal - dObjVal);
    double rGap = compGap / ( 1 + fabs(pObjVal) + fabs(dObjVal) );
    
    potlp->kappa = kappa;
    potlp->tau = tau;
    potlp->pObjVal = pObjVal;
    potlp->dObjVal = dObjVal;
    potlp->pInfeas = pInfeas;
    potlp->pInfeasRel = relpInfeas;
    potlp->dInfeas = dInfeas;
    potlp->dInfeasRel = reldInfeas;
    potlp->complGap = compGap;
    potlp->complGapRel = rGap;
    
    POTLP_MEMCPY(potlp->colVal, colVal, double, nCol);
    POTLP_MEMCPY(potlp->colDual, colDual, double, nCol);
    POTLP_MEMCPY(potlp->rowDual, rowDual, double, nRow);
    
    return;
}

static void POT_FNAME(LPSolverIPrintSolStatistics)( potlp_solver *potlp ) {
    
    double solTime = potUtilGetTimeStamp() - potlp->startT;
    
    /* Potential statistics */
    int nCurvs = 0;
    double curvT = 0.0;
    potReductionGetStatistics(potlp->potIterator, &nCurvs, &curvT);
    double avgTcurv = (nCurvs == 0) ? 0.0 : curvT / nCurvs;
    printf("\nPotential Reduction statistic \n");
    printf("In all [%d] curvs take [%5.3f] seconds ([%5.3f]s/curv) \n", nCurvs, curvT, avgTcurv);
    
    /* LP Statistics */
    if ( potlp->Lpstatus == POTLP_OPTIMAL ) {
        printf("\nLP Status: %s \n", "Optimal");
    } else if ( potlp->Lpstatus == POTLP_MAXITER ) {
        printf("\nLP Status: %s \n", "Maxiter");
    } else if ( potlp->Lpstatus == POTLP_INFEAS_OR_UNBOUNDED ) {
        printf("\nLP Status: %s \n", "Infeasible or Unbounded");
    } else if ( potlp->Lpstatus == POTLP_TIMELIMIT ) {
        printf("\nLP Status: %s \n", "Timelimit");
    } else if ( potlp->Lpstatus == POTLP_USER_INTERRUPT ) {
        printf("\nLP Status: %s \n", "User Interrupt");
    } else if ( potlp->Lpstatus == POTLP_INTERNAL_ERROR ) {
        printf("\nLP Status: %s \n", "Internal Error");
    } else if ( potlp->Lpstatus == POTLP_UNKNOWN ){
        printf("\nLP Status: %s \n", "Unknown");
    } else {
        assert( 0 );
    }
    
    printf("LP Solution statistic \n");
    printf("pObj: %+8.3e   dObj: %+8.3e \n", potlp->pObjVal, potlp->dObjVal);
    printf("pInf:  %8.3e    Rel:  %8.3e \n", potlp->pInfeas, potlp->pInfeasRel);
    printf("dInf:  %8.3e    Rel:  %8.3e \n", potlp->dInfeas, potlp->dInfeasRel);
    printf("Gap :  %8.3e    Rel:  %8.3e \n", potlp->complGap, potlp->complGapRel);
    printf("Time:  %5.3f seconds \n\n", solTime);
    
    return;
}

extern pot_int POT_FNAME(LPSolverCreate)( potlp_solver **ppotlp ) {
    
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
    
    potUtilGetDefaultParams(potlp->dblParams, potlp->intParams);

    POT_CALL(LPQMatCreate(&potlp->potQMatrix));
    POT_CALL(potLPCreate(&potlp->potIterator));
    POT_CALL(potConstrMatCreate(&potlp->potConstrMat));
    POT_CALL(potObjFCreate(&potlp->potObjF));
    POT_CALL(LpNewtonCreate(&potlp->ipm));
    
    potlp->nIter = 0;
    potlp->startT = potUtilGetTimeStamp();
    
    potlp->act.sa_handler = monitorCtrlC;
    sigaction(SIGINT, &potlp->act, NULL);
    
    *ppotlp = potlp;
    
exit_cleanup:
    return retcode;
}

extern pot_int POT_FNAME(LPSolverInit)( potlp_solver *potlp, pot_int nCol, pot_int nRow ) {
    
    pot_int retcode = RETCODE_OK;
    assert( nCol > 0 && nRow > 0 );
    
    potlp->nCol = nCol;
    potlp->nRow = nRow;
    
    /* x; s; kappa; tau */
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
    
    /* Link data and methods to objective */
    potlp->potObjF->objFData = potlp;
    potlp->potObjF->objFVal = POT_FNAME(potLpObjFImplVal);
    potlp->potObjF->objFGrad = POT_FNAME(potLpObjFImplGrad);
    potlp->potObjF->objFHess = POT_FNAME(potLpObjFImplHess);
    potlp->potObjF->objFHVec = POT_FNAME(potLpObjFImplHVec);
    potlp->potObjF->objFMonitor = POT_FNAME(potLpObjFImplMonitor);
    
    /* Line data and methods to constraint  */
    potlp->potConstrMat->AMatData = potlp;
    potlp->potConstrMat->AMatPrepareX = POT_FNAME(potLpConstrMatImplPrepareX);
    potlp->potConstrMat->AMatProject = POT_FNAME(potLpConstrMatImplProject);
    potlp->potConstrMat->AMatScalProject = POT_FNAME(potLpConstrMatImplScalProject);
    potlp->potConstrMat->AMatMonitor = POT_FNAME(potLpConstrMatImplMonitor);
    
    potlp->Lpstatus = POTLP_UNKNOWN;
    
exit_cleanup:
    return retcode;
}

extern pot_int POT_FNAME(LPSolverSetData)( potlp_solver *potlp, pot_int *Ap, pot_int *Ai, double *Ax, double *lpObj, double *lpRHS ) {
    
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
    
    POTLP_INIT(potlp->colVal, double, nCol);
    POTLP_INIT(potlp->colDual, double, nCol);
    POTLP_INIT(potlp->rowDual, double, nRow);
    
    POTLP_INIT(potlp->scalVals, double, 2 * nCol + 2 * nRow + 2);
//    POTLP_INIT(potlp->isColBasic, int, 2 * nCol + 2 * nRow + 2);
    
    if ( !potlp->colMatBeg || !potlp->colMatIdx || !potlp->colMatElem ||
         !potlp->lpObj || !potlp->lpRHS || !potlp->pdcRes ||
         !potlp->colVal || !potlp->colDual || !potlp->rowDual ||
         !potlp->scalVals ) {
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
    
    /* Prepare Q Matrix */
    POT_CALL(LPQMatInit(potlp->potQMatrix, nCol, nRow, potlp->colMatBeg));
    
    /* Prepare the interior point solver */
    POT_CALL(LpNewtonInit(potlp->ipm, nCol, nRow, Ap, Ai, Ax));
    
exit_cleanup:
    return retcode;
}

extern void POT_FNAME(LPSolverParamsPrint)( potlp_solver *potlp ) {
    
    potUtilPrintParams(potlp->dblParams, potlp->intParams);
    potUtilPrintIParams(potlp->dblParams, potlp->intParams);
    
    return;
}

extern pot_int POT_FNAME(LPSolverOptimize)( potlp_solver *potlp ) {
    
    pot_int retcode = RETCODE_OK;
    
    /* Pre-solve */
    POT_FNAME(LPSolverIParamAdjust)(potlp);

    /* Test SCS strategy but destroying the original LP */
//    POT_FNAME(LPSolverIScalInplace(potlp));
    POT_FNAME(LPSolverIScale)(potlp);
    POT_FNAME(LPSovlerIPrintLPStats)(potlp);
    POT_CALL(POT_FNAME(LPSolverISetupQMatrix(potlp)));
    POT_CALL(POT_FNAME(LPSolverIRuizScale(potlp)));
//    POT_CALL(POT_FNAME(LPSolverIL2Scale(potlp)));
    POT_FNAME(LPSolverIHeurInitialize)(potlp);
    
    /* Solve */
    POT_CALL(potReductionSolve(potlp->potIterator));
    
    /* Post-solve */
    POT_FNAME(LPSolverIRetrieveSolution(potlp));
    POT_FNAME(LPSolverIPrintSolStatistics)(potlp);
    
exit_cleanup:
    return retcode;
}

extern void POT_FNAME(LPSolverClear)( potlp_solver *potlp ) {
    
    if ( !potlp ) {
        return;
    }
    
    POTLP_FREE(potlp->colMatBeg);
    POTLP_FREE(potlp->colMatIdx);
    POTLP_FREE(potlp->colMatElem);
    
    POTLP_FREE(potlp->lpRHS);
    POTLP_FREE(potlp->lpObj);
    
    POTLP_FREE(potlp->pdcRes);
    
    POTLP_FREE(potlp->colVal);
    POTLP_FREE(potlp->colDual);
    POTLP_FREE(potlp->rowDual);
    
    POTLP_FREE(potlp->scalVals);
    POTLP_FREE(potlp->isColBasic);
    
    potConstrMatDestroy(&potlp->potConstrMat);
    potObjFDestroy(&potlp->potObjF);
    potLPDestroy(&potlp->potIterator);
    LpNewtonDestroy(&potlp->ipm);
    
    POTLP_ZERO(potlp, potlp_solver, 1);
    
    return;
}

extern void POT_FNAME(LPSolverGetSolution)( potlp_solver *potlp, double *colVal, double *rowDual, double *colDual ) {
        
    if ( !potlp ) {
        return;
    }
    
    pot_int nCol = potlp->nCol;
    pot_int nRow = potlp->nRow;
    
    if ( colVal ) {
        rscl(&nCol, &potlp->tau, potlp->colVal, &potIntConstantOne);
        POTLP_MEMCPY(colVal, potlp->colVal, double, nCol);
    }
    
    if ( rowDual ) {
        rscl(&nRow, &potlp->tau, potlp->rowDual, &potIntConstantOne);
        POTLP_MEMCPY(rowDual, potlp->rowDual, double, nRow);
    }
    
    if ( colDual ) {
        rscl(&nCol, &potlp->tau, potlp->colDual, &potIntConstantOne);
        POTLP_MEMCPY(colDual, potlp->colDual, double, nCol);
    }
    
    return;
}

extern void POT_FNAME(LPSolverDestroy)( potlp_solver **ppotlp ) {
    
    if ( !ppotlp ) {
        return;
    }
    
    POT_FNAME(LPSolverClear)(*ppotlp);
    POTLP_FREE(*ppotlp);
    
    return;
}
