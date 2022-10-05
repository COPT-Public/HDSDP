#include "pot_utils.h"
#include "pot_vector.h"
#include "pot_constr_mat.h"
#include "pot_objfunc.h"
#include "pot_param.h"

static pot_int potReductionFindNegativeCurvature( pot_solver *pot, pot_vec *negCurv ) {
    
    pot_int retcode = RETCODE_OK;
    
    
    return retcode;
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
    
    /* [f, g] = fpot(A, ATA, x_pres); */
    pot->fVal = potObjFVal(objFunc, xPres);
    potObjFGrad(objFunc, xPres, gkVec);
    
    /* Prepare momentum mk = x_pres - x_prev; */
    potVecDiff(mkVec, xPres, xPrev);
    
    /* Gradient projection */
    potConstrMatProj(AMat, gkVec, NULL);
    
    /* Decide whether to update curvature */
    
    
    
    
    
    
    
exit_cleanup:
    
    return retcode;
}

