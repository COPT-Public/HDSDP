#ifndef lpdata_h
#define lpdata_h

#include "pot_structs.h"

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
    
    double *ruizCol;
    double *ruizRow;
    double *scalVals;
    
    double *pdcRes;
    double *pRes;
    double *dRes;
    double *cplRes;
    
    double pObjVal;
    double dObjVal;
    
    double pInfeas;
    double dInfeas;
    double complInfeas;
    double complGap;
    
    double kappa;
    double tau;
    
    double *colVal;
    double *colDual;
    double *rowDual;
            
    pot_solver *potIterator;
    pot_constr_mat *potConstrMat;
    pot_fx *potObjF;
    
    int intParams[NUM_INT_PARAM];
    double dblParams[NUM_DBL_PARAM];
    
    int64_t nIter; ///< Number of iterations
    double  startT; ///< Start time
    
} potlp_solver;

extern pot_int LPSolverCreate( potlp_solver **ppotlp );
extern pot_int LPSolverInit( potlp_solver *potlp, pot_int nCol, pot_int nRow );
extern pot_int LPSolverSetData( potlp_solver *potlp, pot_int *Ap, pot_int *Ai,
                                double *Ax, double *lpObj, double *lpRHS );
extern pot_int LPSolverOptimize( potlp_solver *potlp );
extern void LPSolverDestroy( potlp_solver **ppotlp );

#endif /* lpdata_h */
