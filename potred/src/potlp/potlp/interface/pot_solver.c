#include "pot_def.h"
#include "pot_param.h"
#include "pot_structs.h"
#include "pot_vector.h"
#include "pot_utils.h"
#include "pot_objfunc.h"
#include "pot_constr_mat.h"
#include "pot_lanczos.h"

#include <math.h>

/* TODO: Merge pot_utils into pot_solver */
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
