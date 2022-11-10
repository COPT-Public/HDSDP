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
#include "lp_newton.h"

/* To allow multiple definitions of the same interface */
#define POT_FNAME(x) Ano##x
#ifdef potlp_solver
#undef potlp_solver
#endif
#define potlp_solver potlp_anosolver

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
    double pObjBest;
    double dObjVal;
    double dObjBest;
    
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
    
    /* TODO: Implement the best obj heuristic
     double *colValBest;
     double *colDualBest;
     double *rowDualBest;
     */
    
    /* Heuristics */
    double pResOmega;
    double dResOmega;
    double cplResOmega;
    int *isColBasic;
    
    lp_qmatrix *potQMatrix;
    pot_solver *potIterator;
    pot_constr_mat *potConstrMat;
    pot_fx *potObjF;
    
    /* Interior point solver */
    lp_newton *ipm;
    
    int intParams[NUM_INT_PARAM];
    double dblParams[NUM_DBL_PARAM];
    
    int Lpstatus;
    int64_t nIter; ///< Number of iterations
    double  startT; ///< Start time
    
    struct sigaction act; /// Control C detector 
    
} potlp_anosolver;

extern pot_int POT_FNAME(LPSolverCreate)( potlp_solver **ppotlp );
extern pot_int POT_FNAME(LPSolverInit)( potlp_solver *potlp, pot_int nCol, pot_int nRow );
extern pot_int POT_FNAME(LPSolverSetData)( potlp_solver *potlp, pot_int *Ap, pot_int *Ai,
                                double *Ax, double *lpObj, double *lpRHS );
extern void POT_FNAME(LPSolverParamsPrint)( potlp_solver *potlp );
extern pot_int POT_FNAME(LPSolverOptimize)( potlp_solver *potlp );
extern void POT_FNAME(LPSolverGetSolution)( potlp_solver *potlp, double *colVal, double *rowDual, double *colDual );
extern void POT_FNAME(LPSolverDestroy)( potlp_solver **ppotlp );

#endif /* another_lp_solver_h */
