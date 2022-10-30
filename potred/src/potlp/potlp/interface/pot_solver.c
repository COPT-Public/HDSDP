#include "pot_def.h"
#include "pot_structs.h"
#include "pot_vector.h"
#include "pot_solver.h"
#include "pot_objfunc.h"
#include "pot_constr_mat.h"
#include "pot_lanczos.h"

#include <math.h>

static void potAssembleGrad( double rho, double f, pot_vec *gkVec, pot_vec *gVec, pot_vec *xVec ) {
    
    potVecCopy(gVec, gkVec);
    potVecAxinvpBy(-f, rho, xVec, 1.0, gkVec);
    
    return;
}

static double quadform2by2( double projH[4], double projh[2], double alpha[2] ) {
    
    double linObj = 0.0, quadObj = 0.0;
    
    if ( projh ) {
        linObj = projh[0] * alpha[0] + projh[1] * alpha[1];
    }
    if ( projH ) {
        quadObj += 0.5 * projH[0] * alpha[0] * alpha[0] + projH[1] * alpha[0] * alpha[1];
        quadObj += 0.5 * projH[3] * alpha[1] * alpha[1];
    }
    
    return quadObj + linObj;
}

static void slv2by2( double HplusM[4], const double rhsVec[2], double solVec[2] ) {
    
    double a = HplusM[0], b = HplusM[2], c = HplusM[1], d = HplusM[3];
    double det = a * d - b * c;
    solVec[0] =  d * rhsVec[0] - b * rhsVec[1];
    solVec[1] = -c * rhsVec[0] + a * rhsVec[1];
    solVec[0] = solVec[0] / det;
    solVec[1] = solVec[1] / det;
    
    return;
}

/** @brief Compute the minimum eigenvalue of a 2 by 2 matrix
 *
 */
static double eigs2by2( double GinvL[4] ) {
    
    double a = GinvL[0], b = GinvL[1], c = GinvL[2], d = GinvL[3];
    double trace = (a + d) / 2;
    double det = a * d - b * c;
    return trace - sqrtl(trace * trace - det);
}

static double potReductionTrustRegionSolve( double alpha[2], double projH[4], double projh[2], double projG[4],
                                            double trustRadius, double tol ) {
    
    if ( projH[3] == 0.0 ) {
        
        double G11 = projG[0];
        double alpha11 = sqrt(2.0 * trustRadius / G11);
        
        double modelVal = 0.5 * alpha11 * alpha11 * projH[0] + projh[0] * alpha11;
        
        if ( modelVal > 0.0 ) {
            alpha11 = -alpha11;
            modelVal = 0.5 * alpha11 * alpha11 * projH[0] + projh[0] * alpha11;
        }
        
        alpha[0] = alpha11;
        alpha[1] = 0.0;
        
        return modelVal;
    }
    
    /* 2 by 2 eigen solver */
    double GinvL[4] = {0.0};
    slv2by2(projG, &projH[0], &GinvL[0]);
    slv2by2(projG, &projH[2], &GinvL[2]);
    double lamLowBound = -eigs2by2(GinvL);
    
    if ( lamLowBound < 0.0 ) {
        lamLowBound = 0.0;
    }
    
    double lamUpBound = (lamLowBound > 1.0) ? lamLowBound : 1.0;
    double HplusLamG[4] = {0.0};
    
    /* Reach out for an upper-bound */
    while (1) {
        
        lamUpBound = lamUpBound * 2.0;
        
        HplusLamG[0] = projH[0] + lamUpBound * projG[0];
        HplusLamG[1] = projH[1] + lamUpBound * projG[1];
        HplusLamG[2] = projH[2] + lamUpBound * projG[2];
        HplusLamG[3] = projH[3] + lamUpBound * projG[3];
        
        slv2by2(HplusLamG, projh, alpha);
        alpha[0] = -alpha[0]; alpha[1] = -alpha[1];
        
        double alphaMalpha = quadform2by2(projG, NULL, alpha);
        
        if ( alphaMalpha - trustRadius <= 0 ) {
            lamLowBound = lamUpBound * 0.5;
            break;
        }
    }
    
    /* Trace back by Bisection */
    double diff = lamUpBound - lamLowBound;
    
    while ( diff > 1e-08 ) {
        
        double lamCurrent = lamLowBound + diff * 0.5;
        HplusLamG[0] = projH[0] + lamCurrent * projG[0];
        HplusLamG[1] = projH[1] + lamCurrent * projG[1];
        HplusLamG[2] = projH[2] + lamCurrent * projG[2];
        HplusLamG[3] = projH[3] + lamCurrent * projG[3];
        
        slv2by2(HplusLamG, projh, alpha);
        alpha[0] = -alpha[0]; alpha[1] = -alpha[1];
        
        double alphaMalpha = quadform2by2(projG, NULL, alpha);
        double normDiff = alphaMalpha - trustRadius;
        
        if ( normDiff > tol ) {
            lamLowBound = lamCurrent;
        } else if ( normDiff < -tol ) {
            lamUpBound = lamCurrent;
        } else {
            break;
        }
        
        diff = lamUpBound - lamLowBound;
    }
    
    double modelVal = quadform2by2(projH, projh, alpha);
    return modelVal > 0 ? -POTLP_INFINITY : modelVal;
}

static double potReductionComputePotValue( double rhoVal, double fVal, double zVal, pot_vec *xVec ) {
    
    return rhoVal * logl(fVal - zVal) - potVecLogDet(xVec);
}

static pot_int potReductionOneStep( pot_solver *pot ) {
    
    pot_int retcode = RETCODE_OK;
    
    pot_fx *objFunc = pot->objFunc;
    pot_constr_mat *AMat = pot->AMat;
        
    pot_vec *xPrev = pot->xVecOld;
    pot_vec *xPres = pot->xVec;
    pot_vec *xNorm = pot->xVecNorm;
    pot_vec *dXStep = pot->xStepVec;
    pot_vec *mkVec = pot->mkVec;
    pot_vec *gVec = pot->gVec;
    pot_vec *gkVec = pot->gkVec;
    pot_vec *auxVec1 = pot->auxVec1;
    pot_vec *auxVec2 = pot->auxVec2;
    
    /* Scale xPres */
    potVecCopy(xPres, xNorm);
    potVecConeNormalize(xNorm);
    
    /* [f, g] = fpot(A, ATA, x_pres); */
    pot->fVal = potObjFVal(objFunc, xPres);
    potObjFGrad(objFunc, xPres, gVec);
    
    double fVal = pot->fVal;
    double rhoVal = pot->rhoVal;
    int lczCode = RETCODE_OK;
    
    /* Build the gradient of the potential function */
    potAssembleGrad(rhoVal, fVal, gkVec, gVec, xPres);
    
    /* Prepare the second direction */
    if ( pot->useCurvature && pot->allowCurvature ) {
        /* Use negative curvature */
        lczCode = potLanczosSolve(pot->lczTool, gVec, mkVec);
        potVecConeScal(xPres, mkVec);
    } else {
        /* Use momentum mk = x_pres - x_prev; */
        potVecDiff(mkVec, xPrev, xPres);
    }
    
    /* Gradient projection */
    potConstrMatProj(AMat, gkVec, NULL);
    potConstrMatProj(AMat, mkVec, NULL);
    
    /* Normalize */
    double gkNorm = potVecNormalize(gkVec);
    gkNorm = gkNorm * rhoVal / fVal;
    potVecNormalize(mkVec);
    
    /* Scaled norm */
    double gkXXgk = potVecScaledDot(xPres, gkVec, gkVec);
    double mkXXmk = potVecScaledDot(xPres, mkVec, mkVec);
    double gkXXmk = potVecScaledDot(xPres, gkVec, mkVec);
    
    /* Set up fHgk and fHmk */
    potObjFHVec(objFunc, gkVec, auxVec1);
    potObjFHVec(objFunc, mkVec, auxVec2);
    
    /* Prepare Matrix H */
    double gkfHgk = potVecDot(gkVec, auxVec1);
    double gkfHmk = potVecDot(gkVec, auxVec2);
    double mkfHmk = potVecDot(mkVec, auxVec2);
    
    /* Prepare Matrix G */
    double gTgk = potVecDot(gVec, gkVec);
    double gTmk = potVecDot(gVec, mkVec);
    
    /* Prepare vector h */
    double gkTgk = gkNorm;
    double gkTmk = gkNorm * potVecDot(gkVec, mkVec);
    
    /* Assemble H */
    pot->projHessMat[0] =  (-rhoVal * (gTgk / fVal) * (gTgk / fVal));
    pot->projHessMat[0] += (gkXXgk) + (rhoVal * gkfHgk / fVal);
    pot->projHessMat[1] =  (-rhoVal * (gTmk * gTgk) / (fVal * fVal));
    pot->projHessMat[1] += (gkXXmk) + (rhoVal * gkfHmk / fVal);
    pot->projHessMat[2] =  pot->projHessMat[1];
    pot->projHessMat[3] =  (-rhoVal * (gTmk / fVal) * (gTmk / fVal));
    pot->projHessMat[3] += (mkXXmk) + (rhoVal * mkfHmk / fVal);
    
    /* Assemble h */
    pot->projgVec[0] = gkTgk;
    pot->projgVec[1] = gkTmk;
    
    /* Assemble G */
    pot->projGMat[0] = gkXXgk;
    pot->projGMat[1] = gkXXmk;
    pot->projGMat[2] = gkXXmk;
    pot->projGMat[3] = mkXXmk;
    
    /* Solve the trust region subproblem */
    double potVal = pot->potVal;
    double potReduce = 0.0;
    
    while (1) {
        
        double alphaStep[2] = {0.0};
        double modelVal = potReductionTrustRegionSolve(alphaStep, pot->projHessMat, pot->projgVec,
                                                       pot->projGMat, pot->betaRadius * pot->betaRadius / 4,
                                                       1e-08);
        
        if ( modelVal > 0.0 || pot->betaRadius < 1e-05 ) {
            retcode = RETCODE_FAILED;
            goto exit_cleanup;
        }
        
        /* Assemble the direction */
        potVeczAxpby(dXStep, alphaStep[0], gkVec, alphaStep[1], mkVec);
        potVecAxpy(1.0, xPres, dXStep);
        
        double fValTmp = potObjFVal(objFunc, dXStep);
        double zVal = pot->zVal;
        double potValTmp = potReductionComputePotValue(rhoVal, fValTmp, zVal, dXStep);
        
        potReduce = potValTmp - potVal;
        
        double potReduceRatio = potReduce / modelVal;
        
        if ( potReduceRatio < 0.2 ) {
            pot->betaRadius *= 0.25;
            printf("\n");
        } else if ( potReduceRatio > 0.75 ) {
            pot->betaRadius *= 2.0;
            pot->betaRadius = POTLP_MIN(pot->betaRadius, 1.0);
            pot->potVal = potValTmp;
            break;
        } else {
            pot->betaRadius *= 1.0;
            pot->potVal = potValTmp;
            break;
        }
    }
    
    potVecCopy(xPres, xPrev);
    potVecCopy(dXStep, xPres);
    
    if ( (potReduce > -1e-03 && pot->useCurvature) || lczCode != RETCODE_OK ) {
        pot->allowCurvature = 0;
        
        if ( lczCode != RETCODE_OK ) {
            printf("Curvature is shut down due to failed Lanczos \n");
        } else {
            printf("Curvature is shut down due to insufficient descent \n");
        }
    }
    
    pot->useCurvature = 0;
        
    if ( potReduce > -1e-03 && pot->allowCurvature ) {
        pot->useCurvature = 1;
    }
    
exit_cleanup:
    
    return retcode;
}


static void potLPPotentialHVec( void *pot, pot_vec *vVec, pot_vec *vVecP ) {

    pot_solver *p = pot;
    /* Assemble potential function Hessian-vector product
      
      v    <- Proj_x ( v )
      v    <- [ I  0 ]  [ x1 ]
              [ 0  X ]  [ x2 ]
     
      u[0] <- (g' * v) * g
      u[1] <-  Hess * v
       
      u[0] + u[1] <- S ( u[0] + u[1] )
     
      u[2] <-  [  0 ]
               [ v2 ]
     
      Pv  <- -rho / fVal * u[0] + rho * u[1] + f * u[2]
      Pv  <- Proj_x ( v )
     */
    
    pot_constr_mat *AMat = p->AMat;
    pot_fx *objFunc = p->objFunc;
    
    double rhoVal = p->rhoVal;
    double fVal = p->fVal;
    
    double u0Coeff = rhoVal / fVal;
    double u1Coeff = -rhoVal;
    double u2Coeff = -fVal;
    
    /* Do not modify their contents */
    pot_vec *gVec = p->gVec;
    pot_vec *xVec = p->xVec;
    
    /* xVec contains normalized part of cone */
    pot_vec *xVecNorm = p->xVecNorm;
    
    /* Auxiliary */
    pot_vec *auxVec1 = p->auxVec1;
    pot_vec *auxVec2 = p->auxVec2;
    pot_vec *auxVec3 = p->xStepVec;
    
    /* Reset v and make a backup */
    potVecCopy(vVec, auxVec3);
    potVecReset(vVecP);
    potVecReset(auxVec2);
    
    /* x <- Proj_x ( v ) */
    potConstrMatScalProj(AMat, xVecNorm, vVec, NULL);
    
    /* Store u[2] */
    potVecConeAxpy(u2Coeff, vVec, auxVec2);
    
    /* Scale by x */
    potVecConeScal(xVec, vVec);
    
    /* Add u[0] */
    double gTv = potVecDot(gVec, vVec);
    potVecAxpy(gTv * u0Coeff, gVec, vVecP);
    
    /* Add u[1] */
    potObjFHVec(objFunc, vVec, auxVec1);
    potVecAxpy(u1Coeff, auxVec1, vVecP);
    
    /* Scale by x */
    potVecConeScal(xVec, vVecP);
    
    /* Add u[2] */
    potVecAxpy(1.0, auxVec2, vVecP);
    
    /* Pv <- Proj_x ( v ) */
    potConstrMatScalProj(AMat, xVecNorm, vVecP, NULL);
    
    /* Copy back */
    potVecCopy(auxVec3, vVec);
    
    return;
}

extern pot_int potLPCreate( pot_solver **ppot ) {
    
    pot_int retcode = RETCODE_OK;
    
    pot_solver *pot = NULL;
    POTLP_INIT(pot, pot_solver, 1);
    
    if (!pot) {
        retcode = RETCODE_FAILED;
        goto exit_cleanup;
    }
    
    POTLP_ZERO(pot, pot_solver, 1);
    
    pot->potVal = POTLP_INFINITY;
    pot->zVal = 0.0;
    
    POT_CALL(potLanczosCreate(&pot->lczTool));
    
    *ppot = pot;
    
exit_cleanup:
    return retcode;
}

extern pot_int potLPInit( pot_solver *pot, pot_int vDim, pot_int vConeDim ) {
    
    pot_int retcode = RETCODE_OK;
    
    if ( pot->potVal != POTLP_INFINITY || pot->xVec || pot->xVecOld ||
         pot->gVec || pot->mkVec || pot->xStepVec || pot->HessMat ||
         pot->auxVec1 || pot->auxVec2 ) {
        goto exit_cleanup;
    }

    POT_CALL(potVecCreate(&pot->xVec));
    POT_CALL(potVecInit(pot->xVec, vDim, vConeDim));
    
    POT_CALL(potVecCreate(&pot->xVecOld));
    POT_CALL(potVecInit(pot->xVecOld, vDim, vConeDim));
    
    POT_CALL(potVecCreate(&pot->xVecNorm));
    POT_CALL(potVecInit(pot->xVecNorm, vDim, vConeDim));
    
    POT_CALL(potVecCreate(&pot->gVec));
    POT_CALL(potVecInit(pot->gVec, vDim, vConeDim));
    
    POT_CALL(potVecCreate(&pot->gkVec));
    POT_CALL(potVecInit(pot->gkVec, vDim, vConeDim));
    
    POT_CALL(potVecCreate(&pot->mkVec));
    POT_CALL(potVecInit(pot->mkVec, vDim, vConeDim));
    
    POT_CALL(potVecCreate(&pot->xStepVec));
    POT_CALL(potVecInit(pot->xStepVec, vDim, vConeDim));
    
    POT_CALL(potVecCreate(&pot->auxVec1));
    POT_CALL(potVecInit(pot->auxVec1, vDim, vConeDim));
    
    POT_CALL(potVecCreate(&pot->auxVec2));
    POT_CALL(potVecInit(pot->auxVec2, vDim, vConeDim));
    
    pot->rhoVal = vConeDim + sqrt(vConeDim);
    POT_CALL(potLanczosInit(pot->lczTool, vDim, vConeDim));
    pot->lczTool->MMat = pot;
    pot->lczTool->lczMatVec = potLPPotentialHVec;
    
#ifdef POT_DEBUG
    POTLP_INIT(pot->HessMat, double, vDim * vDim);
    if (!pot->HessMat) {
        retcode = RETCODE_FAILED;
    }
#endif
    
exit_cleanup:
    return retcode;
}

extern pot_int potLPSetObj( pot_solver *pot, pot_fx *objFunc ) {
    
    pot_int retcode = RETCODE_OK;
    
    if (pot->objFunc) {
        retcode = RETCODE_FAILED;
        goto exit_cleanup;
    }
    
    pot->objFunc = objFunc;

exit_cleanup:
    return retcode;
}

extern pot_int potLPSetLinearConstrs( pot_solver *pot, pot_constr_mat *AMat ) {
    
    pot_int retcode = RETCODE_OK;
    
    if (pot->AMat) {
        retcode = RETCODE_FAILED;
        goto exit_cleanup;
    }
    
    pot->AMat = AMat;
    
exit_cleanup:
    return retcode;
}

extern pot_int potReductionSolve( pot_solver *pot ) {
    
    pot_int retcode = RETCODE_OK;
    
    /* Reset information */
    pot->betaRadius = 1.0;
    
    potVecReset(pot->xVec);
    potVecReset(pot->xVecOld);
    potVecReset(pot->gVec);
    potVecReset(pot->gkVec);
    potVecReset(pot->mkVec);
    potVecReset(pot->auxVec1);
    potVecReset(pot->auxVec2);
    
    /* Use curvature in the initial iteration */
    pot->useCurvature = 1;
    pot_int info = 0;
    
    /* Initialize */
    potConstrMatPrepareX(pot->AMat, pot->xVec);
    potConstrMatPrepareX(pot->AMat, pot->xVecOld);
    
    for ( int i = 0; ; ++i ) {
        
        if ( i % 1000 == 0 && i < 10000 && pot->allowCurvature ) {
            pot->useCurvature = 1;
        }
        
        retcode = potReductionOneStep(pot);
        potObjFMonitor(pot->objFunc, &info);
        
        if ( retcode != RETCODE_OK ) {
            break;
        }
        
        if ( info ) {
            break;
        }
    }
    
    return retcode;
}

extern void potLPClear( pot_solver *pot ) {
    
    if (!pot) {
        return;;
    }
    
    pot->objFunc = NULL;
    pot->AMat = NULL;
    pot->fVal = POTLP_INFINITY;
    
    potVecDestroy(&pot->xVec);
    potVecDestroy(&pot->xVecOld);
    potVecDestroy(&pot->xVecNorm);
    potVecDestroy(&pot->gVec);
    potVecDestroy(&pot->mkVec);
    potVecDestroy(&pot->xStepVec);
    potVecDestroy(&pot->auxVec1);
    potVecDestroy(&pot->auxVec2);
    
    pot->objFunc = NULL;
    pot->AMat = NULL;
    
    potLanczosDestroy(&pot->lczTool);
    
    POTLP_FREE(pot->HessMat);
    POTLP_ZERO(pot, pot_solver, 1);
    
    return;
}

extern void potLPDestroy( pot_solver **ppot ) {
    
    if (!ppot) {
        return;
    }
    
    potLPClear(*ppot);
    POTLP_FREE(*ppot);
    
    return;
}
