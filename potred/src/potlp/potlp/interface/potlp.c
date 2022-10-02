#include "potlp.h"
#include "pot_param.h"
#include "pot_structs.h"
#include "pot_vector.h"
#include "pot_utils.h"

extern pot_int potLPCreate( pot_solver **ppot ) {
    
    pot_int retcode = RETCODE_OK;
    
    pot_solver *pot = NULL;
    POTLP_INIT(pot, pot_solver, 1);
    
    if (!pot) {
        retcode = RETCODE_FAILED;
        goto exit_cleanup;
    }
    
    memset(pot, 0, sizeof(pot_int));
    memcpy(pot->dblParams, defaultDblParam, sizeof(double) * NUM_DBL_PARAM);
    memcpy(pot->intParams, defaultIntParam, sizeof(int) * NUM_INT_PARAM);
    
    pot->potVal = POTLP_INFINITY;
    pot->zVal = -POTLP_INFINITY;
    *ppot = pot;
    
exit_cleanup:
    return retcode;
}

extern pot_int potLPSetDim( pot_solver *pot, pot_int vDim, pot_int vConeDim ) {
    
    pot_int retcode = RETCODE_OK;
    
    if ( pot->fVal != POTLP_INFINITY || pot->xVec || pot->xVecOld ||
         pot->gVec || pot->mVec || pot->xStepVec || pot->HessMat ||
         pot->auxVec1 || pot->auxVec2 ) {
        goto exit_cleanup;
    }

    POT_CALL(potVecInit(pot->xVec, vDim, vConeDim));
    POT_CALL(potVecInit(pot->xVecOld, vDim, vConeDim));
    POT_CALL(potVecInit(pot->gVec, vDim, vConeDim));
    POT_CALL(potVecInit(pot->mVec, vDim, vConeDim));
    POT_CALL(potVecInit(pot->xStepVec, vDim, vConeDim));
    POT_CALL(potVecInit(pot->auxVec1, vDim, vConeDim));
    POT_CALL(potVecInit(pot->auxVec2, vDim, vConeDim));
    
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
    potVecDestroy(&pot->gVec);
    potVecDestroy(&pot->mVec);
    potVecDestroy(&pot->xStepVec);
    potVecDestroy(&pot->auxVec1);
    potVecDestroy(&pot->auxVec2);
    
    pot->objFunc = NULL;
    pot->AMat = NULL;
    
    POTLP_FREE(pot->HessMat);
    memset(pot, 0, sizeof(pot_solver));
    
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
