#include "pot_constr_mat.h"
#include "pot_vector.h"

extern pot_int potConstrMatInit( pot_constr_mat *potConstrMat, pot_int nRows, pot_int nCols ) {
    
    pot_int retcode = RETCODE_OK;
    
    if ( !potConstrMat ) {
        retcode = RETCODE_FAILED;
        goto exit_cleanup;
    }
    
    if ( potConstrMat->m > 0 || potConstrMat->n > 0 || potConstrMat->AMatData ) {
        retcode = RETCODE_FAILED;
        goto exit_cleanup;
    }
    
    /* Initialize methods */
    potConstrMat->AMatData = NULL;
    potConstrMat->AMatInit = NULL;
    potConstrMat->AMatPrepareX = NULL;
    potConstrMat->AMatProject = NULL;
    potConstrMat->AMatMonitor = NULL;
    potConstrMat->AMatDestroy = NULL;
    
    if ( !potConstrMat->AMatData || !potConstrMat->AMatInit ||
         !potConstrMat->AMatPrepareX || !potConstrMat->AMatProject ||
         !potConstrMat->AMatMonitor || !potConstrMat->AMatDestroy ) {
        retcode = RETCODE_FAILED;
        goto exit_cleanup;
    }
    
exit_cleanup:
    return retcode;
}

extern pot_int potConstrMatInitData( pot_constr_mat *potConstrMat, void *inputData ) {
    
    return potConstrMat->AMatInit(&potConstrMat->AMatData, inputData);
}

extern void potConstrMatPrepareX( pot_constr_mat *potConstrMat, pot_vec *xVec ) {
    
    potConstrMat->AMatPrepareX(potConstrMat->AMatData, xVec);
}

extern void potConstrMatProj( pot_constr_mat *potConstrMat, pot_vec *xVec, pot_vec *yVec ) {
    
    if ( !yVec ) {
        potConstrMat->AMatProject(potConstrMat->AMatData, xVec);
    } else {
        potVecCopy(xVec, yVec);
        potConstrMat->AMatProject(potConstrMat->AMatData, yVec);
    }
    
    return;
}

extern void potConstrMatMonitor( pot_constr_mat *potConstrMat ) {
    
    return;
}

extern void potConstrMatDestroy( pot_constr_mat *potConstrMat ) {
    
    if ( !potConstrMat ) {
        return;
    }
    potConstrMat->AMatDestroy(&potConstrMat->AMatData);
    memset(potConstrMat, 0, sizeof(pot_constr_mat));
    
    return;
}
