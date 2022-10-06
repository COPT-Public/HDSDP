#include "pot_lanczos.h"
#include "pot_vector.h"
#include "pot_utils.h"

#define SYEV_LWORK 30
#define SYEV_IWORK 12

extern pot_int potLanczosCreate( pot_lanczos **ppotLanczos ) {
    
    pot_int retcode = RETCODE_OK;
    
    if ( !ppotLanczos ) {
        retcode = RETCODE_FAILED;
        goto exit_cleanup;
    }
    
    pot_lanczos *potLanczos = NULL;
    POTLP_INIT(potLanczos, pot_lanczos, 1);
    
    if ( !potLanczos ) {
        retcode = RETCODE_FAILED;
        goto exit_cleanup;
    }
    
    memset(potLanczos, 0, sizeof(pot_lanczos));
    *ppotLanczos = potLanczos;
    
exit_cleanup:
    return retcode;
}

extern pot_int potLanczosInit( pot_lanczos *potLanczos, pot_int nCols ) {
    
    pot_int retcode = RETCODE_OK;
    potLanczos->n = nCols;
    
    if ( potLanczos->vVec   || potLanczos->wVec     || potLanczos->z1Vec ||
         potLanczos->z2Vec  || potLanczos->vaVec    || potLanczos->VMat  ||
         potLanczos->HMat   || potLanczos->YMat     || potLanczos->UMat  ||
         potLanczos->dArray || potLanczos->eiDblMat || potLanczos->eiIntMat ) {
        
        retcode = RETCODE_FAILED;
        goto exit_cleanup;
    }
    
    potLanczos->maxIter = 1000;
    
    POT_CALL(potVecCreate(&potLanczos->vVec));
    POT_CALL(potVecInit(potLanczos->vVec, nCols, 0));
    
    POT_CALL(potVecCreate(&potLanczos->wVec));
    POT_CALL(potVecInit(potLanczos->wVec, nCols, 0));
    
    POT_CALL(potVecCreate(&potLanczos->z1Vec));
    POT_CALL(potVecInit(potLanczos->z1Vec, nCols, 0));
    
    POT_CALL(potVecCreate(&potLanczos->z2Vec));
    POT_CALL(potVecInit(potLanczos->z2Vec, nCols, 0));
    
    POT_CALL(potVecCreate(&potLanczos->vaVec));
    POT_CALL(potVecInit(potLanczos->vaVec, nCols, 0));
    
    POTLP_INIT(potLanczos->VMat, double, nCols * (potLanczos->maxIter + 1));
    POTLP_INIT(potLanczos->HMat, double, (potLanczos->maxIter + 1) * (potLanczos->maxIter + 1));
    POTLP_INIT(potLanczos->YMat, double, potLanczos->maxIter * 2);
    POTLP_INIT(potLanczos->dArray, double, potLanczos->maxIter);
    POTLP_INIT(potLanczos->UMat, double, potLanczos->maxIter * potLanczos->maxIter);
    POTLP_INIT(potLanczos->eiDblMat, double, potLanczos->maxIter * SYEV_LWORK);
    POTLP_INIT(potLanczos->eiIntMat, pot_int, potLanczos->maxIter * SYEV_IWORK);
    
    if ( !potLanczos->YMat || !potLanczos->HMat || !potLanczos->YMat ||
         !potLanczos->dArray || !potLanczos->UMat || !potLanczos->eiDblMat ||
         !potLanczos->eiIntMat ) {
        retcode = RETCODE_FAILED;
        goto exit_cleanup;
    }
    
exit_cleanup:
    
    return retcode;
}

extern void potLanczosInitData( pot_lanczos *potLanczos, void *MMat, void (*lczMatVec) (void *, pot_vec *, pot_vec *) ) {
    
    if ( potLanczos->MMat || potLanczos->lczMatVec ) {
        return;
    }
    
    potLanczos->MMat = MMat;
    potLanczos->lczMatVec = lczMatVec;
    
    return;
}

extern pot_int potLanczosSolve( pot_lanczos *potLanczos, pot_vec *nCurv ) {
    
    pot_int retcode = RETCODE_OK;
    
    /* TODO: Implement the Lanczos procedure */
    
    return retcode;
}

extern void potLanczosClear( pot_lanczos *potLanczos ) {
    
    if ( !potLanczos ) {
        return;
    }
    
    potLanczos->n = 0;
    potLanczos->maxIter = 0;
    
    potLanczos->MMat = NULL;
    
    potVecDestroy(&potLanczos->vVec);
    potVecDestroy(&potLanczos->wVec);
    potVecDestroy(&potLanczos->z1Vec);
    potVecDestroy(&potLanczos->z2Vec);
    potVecDestroy(&potLanczos->vaVec);
    
    POTLP_FREE(potLanczos->VMat);
    POTLP_FREE(potLanczos->HMat);
    POTLP_FREE(potLanczos->YMat);
    POTLP_FREE(potLanczos->dArray);
    POTLP_FREE(potLanczos->UMat);
    POTLP_FREE(potLanczos->eiDblMat);
    POTLP_FREE(potLanczos->eiIntMat);
    
    potLanczos->lczMatVec = NULL;
    
    return;
}

extern void potLanczosDestroy( pot_lanczos **ppotLanczos ) {
    
    if ( !ppotLanczos ) {
        return;
    }
    
    potLanczosClear(*ppotLanczos);
    POTLP_FREE(*ppotLanczos);
    
    return;
}
