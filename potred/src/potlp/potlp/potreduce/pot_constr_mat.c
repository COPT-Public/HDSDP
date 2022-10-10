#include "pot_constr_mat.h"
#include "pot_vector.h"

extern pot_int potConstrMatCreate( pot_constr_mat **ppotConstrMat ) {
    
    pot_int retcode = RETCODE_OK;
    
    if ( !ppotConstrMat ) {
        retcode = RETCODE_FAILED;
        goto exit_cleanup;
    }
    
    pot_constr_mat *potConstrMat = NULL;
    POTLP_INIT(potConstrMat, pot_constr_mat, 1);
    
    if ( !potConstrMat ) {
        retcode = RETCODE_FAILED;
        goto exit_cleanup;
    }
    
    memset(potConstrMat, 0, sizeof(pot_constr_mat));
    *ppotConstrMat = potConstrMat;
    
exit_cleanup:
    return retcode;
}

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
    
    potConstrMat->AMatData = NULL;
    
    /* Initialize methods */
    potConstrMat->AMatPrepareX = NULL;
    potConstrMat->AMatProject = NULL;
    potConstrMat->AMatScalProject = NULL;
    potConstrMat->AMatMonitor = NULL;
    
exit_cleanup:
    return retcode;
}

extern void potConstrMatPrepareX( pot_constr_mat *potConstrMat, pot_vec *xVec ) {
    
    potConstrMat->AMatPrepareX(potConstrMat->AMatData, xVec);
}

extern void potConstrMatProj( pot_constr_mat *potConstrMat, pot_vec *xVec, pot_vec *xVecP ) {
    
    if ( !xVecP ) {
        potConstrMat->AMatProject(potConstrMat->AMatData, xVec);
    } else {
        potVecCopy(xVec, xVecP);
        potConstrMat->AMatProject(potConstrMat->AMatData, xVecP);
    }
    
    return;
}

extern void potConstrMatScalProj( pot_constr_mat *potConstrMat, pot_vec *xVec, pot_vec *yVec, pot_vec *yVecP ) {
    
    if ( !yVecP ) {
        potConstrMat->AMatScalProject(potConstrMat->AMatData, xVec, yVec);
    } else {
        potVecCopy(yVec, yVecP);
        potConstrMat->AMatScalProject(potConstrMat->AMatData, xVec, yVecP);
    }
    
    return;
}

extern void potConstrMatMonitor( pot_constr_mat *potConstrMat, void *info ) {
    
    potConstrMat->AMatMonitor(potConstrMat->AMatData, info);
}

extern void potConstrMatClear( pot_constr_mat *potConstrMat ) {
    
    memset(potConstrMat, 0, sizeof(pot_constr_mat));
    
    return;
}

extern void potConstrMatDestroy( pot_constr_mat **ppotConstrMat ) {
    
    if ( !ppotConstrMat ) {
        return;
    }
    
    potConstrMatClear(*ppotConstrMat);
    POTLP_FREE(*ppotConstrMat);
    
    return;
}
