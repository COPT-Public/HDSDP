#include "pot_utils.h"
#include "pot_vector.h"
#include "pot_constr_mat.h"
#include "pot_objfunc.h"
#include "pot_lanczos.h"
#include "pot_param.h"

#include "math.h"

static void potAssembleGrad( double rho, double f, pot_vec *gkVec, pot_vec *gVec, pot_vec *xVec ) {
    
    potVecCopy(gVec, gkVec);
    potVecAxinvpBy(-f, rho, xVec, 1.0, gkVec);
    
    return;
}

#if 0
static int chol2by2( double psdG[4], double cholL[4] ) {
    /*
       [a, b]
       [b, c]
     */
    double G11 = psdG[0], G12 = psdG[1], G22 = psdG[3];
    
    if ( G11 <= 0.0 ) {
        return 0;
    }
    
    cholL[0] = sqrtl(G11);
    cholL[1] = G12 / cholL[0];
    cholL[2] = 0.0;
    cholL[3] = G22 - (cholL[1] * cholL[1]);
    
    if ( cholL[3] <= 0.0 ) {
        return 0;
    } else {
        cholL[3] = sqrtl(cholL[3]);
    }
    
    return 1;
}

static void trsm2by2( double cholL[4], const double rhsMat[4], double solMat[4] ) {
    
    solMat[0] = rhsMat[0] / cholL[0];
    solMat[1] = (rhsMat[1] - cholL[1] * solMat[0]) / cholL[3];
    
    solMat[2] = rhsMat[2] / cholL[0];
    solMat[3] = (rhsMat[3] - cholL[1] * solMat[2]) / cholL[3];
    
    return;
}

static void transpose2by2( double LinvQ[4] ) {
    
    double tmp = LinvQ[2];
    LinvQ[2] = LinvQ[1];
    LinvQ[1] = tmp;
}
#endif

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
    
    while ( diff > 1e-06 ) {
        
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
    pot_vec *dXStep = pot->xStepVec;
    pot_vec *mkVec = pot->mkVec;
    pot_vec *gVec = pot->gVec;
    pot_vec *gkVec = pot->gkVec;
    pot_vec *auxVec1 = pot->auxVec1;
    pot_vec *auxVec2 = pot->auxVec2;
    
    /* [f, g] = fpot(A, ATA, x_pres); */
    pot->fVal = potObjFVal(objFunc, xPres);
    potObjFGrad(objFunc, xPres, gVec);
    
    double fVal = pot->fVal;
    double rhoVal = pot->rhoVal;
    
    /* Build the gradient of the potential function */
    potAssembleGrad(rhoVal, fVal, gkVec, gVec, xPres);
    
    /* Prepare the second direction */
    if ( 0 ) {
        /* Use negative curvature */
        POT_CALL(potLanczosSolve(pot->lczTool, mkVec));
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
        
        if ( modelVal > 0.0 ) {
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
        } else if ( potReduceRatio > 0.75 ) {
            pot->betaRadius *= 2.0;
            pot->betaRadius = (pot->betaRadius > 0.99995) ? 0.99995 : pot->betaRadius;
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
    
    if ( potReduce > -0.01 ) {
        pot->useCurvature = 1;
    }
    
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
    
    int maxIter = pot->intParams[INT_PARAM_MAXITER];
    pot_int info = 0;
    
    /* Initialize */
    potConstrMatPrepareX(pot->AMat, pot->xVec);
    potConstrMatPrepareX(pot->AMat, pot->xVecOld);
    
    printf("Iteration log. \n");
    printf("%10s  %10s  %10s  %10s  %10s \n", "pObj", "dObj", "pInf", "dInf", "k/t");
    
    for ( int i = 0; i < maxIter; ++i ) {
        
        retcode = potReductionOneStep(pot);
        
        if ( i % 1000 == 1) {
            potObjFMonitor(pot->objFunc, &info);
        }
        
        if ( retcode != RETCODE_OK ) {
            break;
        }
        
        if ( info ) {
            break;
        }
    }
    
    return retcode;
}
