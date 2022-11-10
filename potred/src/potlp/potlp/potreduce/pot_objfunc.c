#include "pot_objfunc.h"

extern pot_int potObjFCreate( pot_fx **ppotObjF ) {
    
    pot_int retcode = RETCODE_OK;
    
    if ( !ppotObjF ) {
        retcode = RETCODE_FAILED;
        goto exit_cleanup;
    }
    
    pot_fx *potObjF = NULL;
    
    POTLP_INIT(potObjF, pot_fx, 1);
    
    if ( !potObjF ) {
        retcode = RETCODE_FAILED;
        goto exit_cleanup;
    }
    
    POTLP_ZERO(potObjF, pot_fx, 1);
    *ppotObjF = potObjF;
    
exit_cleanup:
    return retcode;
}

extern pot_int potObjFInit( pot_fx *potObjF, pot_int nCols ) {
    
    pot_int retcode = RETCODE_OK;
    
    if ( !potObjF ) {
        retcode = RETCODE_FAILED;
        goto exit_cleanup;
    }
    
    if ( potObjF->n > 0 || potObjF->objFData ) {
        retcode = RETCODE_FAILED;
        goto exit_cleanup;
    }
    
    potObjF->n = nCols;
    potObjF->objFData = NULL;
    
    /* Initialize methods */
    potObjF->objFGrad = NULL;
    potObjF->objFHess = NULL;
    potObjF->objFHVec = NULL;
    potObjF->objFMonitor = NULL;
    
exit_cleanup:
    return retcode;
}

extern double potObjFVal( pot_fx *potObjF, pot_vec *xVec ) {
    
    return potObjF->objFVal(potObjF->objFData, xVec);
}

extern void potObjFGrad( pot_fx *potObjF, pot_vec *xVec, pot_vec *fGrad ) {
    
    potObjF->objFGrad(potObjF->objFData, xVec, fGrad);
    
    return;
}

extern void potObjFHess( pot_fx *potObjF, pot_vec *xVec, double *fHess ) {
    
    potObjF->objFHess(potObjF->objFData, xVec, fHess);
    
    return;
}

extern void potObjFHVec( pot_fx *potObjF, pot_vec *vVec, pot_vec *fHvVec ) {
    
    potObjF->objFHVec(potObjF->objFData, vVec, fHvVec);
    
    return;
}

extern void potObjFMonitor( pot_fx *potObjF, void *info ) {
    
    potObjF->objFMonitor(potObjF->objFData, info);
    
    return;
}

extern void potObjFClear( pot_fx *potObjF ) {
    
    if ( !potObjF ) {
        return;
    }
    
    POTLP_ZERO(potObjF, pot_fx, 1);
    return;
}

extern void potObjFDestroy( pot_fx **ppotObjF ) {
    
    if ( !ppotObjF ) {
        return;
    }
    
    potObjFClear(*ppotObjF);
    POTLP_FREE(*ppotObjF);
    
    return;
}
