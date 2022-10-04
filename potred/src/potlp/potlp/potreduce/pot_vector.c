/** @file pot\_vector.c
 *  @brief Implement the vector methods
 *
 * @TODO: Add more detailed comments
 */

#include "pot_vector.h"
#include "vec_mat.h"

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

extern void potVecDiff( pot_vec *pVecOut, pot_vec *pVecInPrev, pot_vec *pVecInPres ) {
    
    assert( pVecOut->n == pVecInPres->n );
    assert( pVecInPres->n == pVecInPrev->n );
    
    for ( int i = 0; i < pVecOut->n; ++i ) {
        pVecOut->x[i] = pVecInPres->x[i] - pVecInPrev->x[i];
    }
    
    pVecOut->nrm = -1.0;
}

extern void potVecNormalize( pot_vec *pVec ) {
    
    assert( pVec->nrm == -1.0 );
    pVec->nrm = nrm2(&pVec->n, pVec->x, &potIntConstantOne);
    rscl(&pVec->n, pVec->x, &pVec->nrm, &potIntConstantOne);
    
}

extern void potVecCopy( pot_vec *srcVec, pot_vec *dstVec ) {
    
    assert( srcVec->n == dstVec->n );
    memcpy(dstVec->x, srcVec->x, sizeof(double) * srcVec->n);
    
    dstVec->nrm = srcVec->nrm;
    
    return;
}

extern void potVecClear( pot_vec *pVec ) {
    
    if (!pVec) {
        return;
    }
    
    pVec->n = 0; pVec->ncone = 0; pVec->nrm = -1.0;
    POTLP_FREE(pVec->x);
    
    return;
}

extern void potVecDestroy( pot_vec **ppVec ) {
    
    if (!ppVec) {
        return;
    }
    
    potVecClear(*ppVec);
    POTLP_FREE(*ppVec);
    
    return;
}
