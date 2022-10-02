#include "pot_objfunc.h"

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
    
    potObjF->objFData = NULL;
    potObjF->xVec = NULL;
    
    /* Initialize methods */
    potObjF->objFInit = NULL;
    potObjF->objFGrad = NULL;
    potObjF->objFHess = NULL;
    potObjF->objFHVec = NULL;
    potObjF->objFMonitor = NULL;
    potObjF->objFDestroy = NULL;
    
    if ( !potObjF->objFInit || !potObjF->objFGrad ||
         !potObjF->objFHess || !potObjF->objFHVec ||
         !potObjF->objFMonitor || !potObjF->objFDestroy ) {
        retcode = RETCODE_FAILED;
        goto exit_cleanup;
    }
    
exit_cleanup:
    return retcode;
}

extern void potObjFSetX( pot_fx *potObjF, pot_vec *xVec ) {
    
    if ( potObjF->xVec ) {
        return;
    }
    
    potObjF->xVec = xVec;
    
    return;
}

extern pot_int potObjFInitData( pot_fx *potObjF, void *inputData ) {
    
    return potObjF->objFInit(&potObjF->objFData, inputData);
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

extern void potObjFHVec( pot_fx *potObjF, pot_vec *xVec, pot_vec *fHVec ) {
    
    potObjF->objFGrad(potObjF->objFData, xVec, fHVec);
    
    return;
}

extern void potObjFMonitor( pot_fx *potObjF, pot_int *info ) {
    
    potObjF->objFMonitor(potObjF->objFData, info);
}

extern void potObjFDestroy( pot_fx *potObjF ) {
    
    if ( !potObjF ) {
        return;
    }
    
    potObjF->objFDestroy(&potObjF->objFData);
    memset(potObjF, 0, sizeof(pot_fx));
    
    return;
}
