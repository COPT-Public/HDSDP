/** @file another\_lp\_solver.h
 *  @brief Another implementation of the potential reduction solver
 * This LP solver treats scaling differently by constructing the large HSD matrix
 *
 *  H = [ 0    A   0  0   -b ]
 *     [ -A'   0  -I   0   c ]
 *     [  b'  -c'  0  -1   0 ]
 */

#ifndef another_lp_solver_h
#define another_lp_solver_h

#include "pot_structs.h"
#include "lp_qmatrix.h"

typedef struct {
    
    pot_int nCol; ///< Number of LP variables
    pot_int nRow; ///< Number of constraints
    
    pot_int *colMatBeg;
    pot_int *colMatIdx;
    double  *colMatElem;
    
    double *lpRHS;
    double *lpObj;
    
    double lpRHSNorm;
    double lpObjNorm;
    
    double objScaler;
    double rhsScaler;
    
    double *scalVals;
    
    double *pdcRes;
    double *pRes;
    double *dRes;
    double *cplRes;
    
    double pObjVal;
    double dObjVal;
    
    double pInfeas;
    double pInfeasRel;
    double dInfeas;
    double dInfeasRel;
    double complGap;
    double complGapRel;
    
    double kappa;
    double tau;
    
    double *colVal;
    double *colDual;
    double *rowDual;
    
    /* Heuristics */
    double pResOmega;
    double dResOmega;
    double cplResOmega;
    int *isColBasic;
    
    lp_qmatrix *potQMatrix;
    pot_solver *potIterator;
    pot_constr_mat *potConstrMat;
    pot_fx *potObjF;
    
    int intParams[NUM_INT_PARAM];
    double dblParams[NUM_DBL_PARAM];
    
    int64_t nIter; ///< Number of iterations
    double  startT; ///< Start time
    
} potlp_anosolver;

#endif /* another_lp_solver_h */
