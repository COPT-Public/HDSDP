/** @file pot\_vector.c
 *  @brief Implement the vector methods
 *
 * @TODO: Add more detailed comments
 */

#include "pot_vector.h"

extern pot_int potVecInit( pot_vec *pVec, pot_int vDim, pot_int vConeDim ) {
    
    pot_int retcode = RETCODE_OK;
    
    if ( !pVec ) {
        retcode = RETCODE_FAILED;
        goto exit_cleanup;
    }
    
    if ( pVec->n > 0 || pVec->x ) {
        retcode = RETCODE_FAILED;
        goto exit_cleanup;
    }
    
    if ( vConeDim > vDim || vConeDim < 0 ) {
        retcode = RETCODE_FAILED;
        goto exit_cleanup;
    }
    
    POTLP_INIT(pVec->x, double, vDim);
    
    if ( !pVec->x ) {
        retcode = RETCODE_FAILED;
        goto exit_cleanup;
    }
    
    pVec->n = vDim; pVec->ncone = vConeDim;
    
exit_cleanup:
    return retcode;
}

extern void potVecCopy( pot_vec *srcVec, pot_vec *dstVec ) {
    
    memcpy(dstVec->x, srcVec->x, sizeof(double) * srcVec->n);
    
    return;
}

extern void potVecDestroy( pot_vec *pVec ) {
    
    if (!pVec) {
        return;
    }
    
    pVec->n = 0; pVec->ncone = 0; pVec->nrm = -1.0;
    POTLP_FREE(pVec->x);
    
    return;
}
