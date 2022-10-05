/** @file lpdata.c
 *  @brief Implement the LP data interface for potential reduction
 *
 * @TODO: Add more detailed comments
 */

#include "lp_solver.h"
#include "pot_structs.h"

typedef struct {
    
    pot_int nCol; ///< Number of LP variables
    pot_int nConstr; ///< Number of constraints
    
    pot_int *colMatBeg;
    pot_int *colMatIdx;
    double  *colMatVal;
    
    double *lpRHS;
    double *lpObj;
    
    double lpRHSNorm;
    double lpObjNorm;
    
    double *ruizCol;
    double *ruizRow;
    
    double *pdcRes;
    double *pRes;
    double *dRes;
    double *cplRes;
    
    double pInfeas;
    double dInfeas;
    double complGap;
    
} potlp_solver;


/* Constraint methods */
static void potLpConstrMatImplPrepareX( void *AMatData, pot_vec *xInit ) {
    
    
    return;
}

static void potLpConstrMatImplProject( void *AMatData, pot_vec *xVec ) {
    
    
    
    return;
}

static void potLpConstrMatImplScalProject( void *AMatData, pot_vec *xVec, pot_vec *yVec ) {
    
    
    return;
}

static void potLpConstrMatImplMonitor( void *AMatData, void *info ) {
    
    
    return;
}

/* Objective methods */
static double potLpObjFImplVal( void *ObjFData, pot_vec *xVec ) {
    
    
    return 0.0;
}

static void potLpObjFImplGrad( void *ObjFData, pot_vec *xVec, pot_vec *fGrad ) {
    
    return;
}

static void potLpObjFImplHess( void *objFData, pot_vec *xVec, double *fHess ) {
    
    return;
}

static void potLpObjFImplHVec( void *objFData, pot_vec *xVec, pot_vec *fHvec ) {
    
    return;
}

static void potLpObjFImplMonitor( void *objFData, void *info ) {
    
    return;
}
