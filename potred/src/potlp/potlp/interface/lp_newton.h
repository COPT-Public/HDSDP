#ifndef lp_newton_h
#define lp_newton_h

#include "pot_solver.h"
#include "linsys.h"

typedef struct {
    
    int nCol;
    int nRow;
    
    int    *AugBeg;
    int    *AugIdx;
    double *AugElem;
    
    double *colBackup;
    
    pot_linsys *lpLinsys; ///< Indefinite augmented system
    
    /* Algorithm parameters */
    double alpha;
    double beta;
    double gamma;
    double mu;
    
    /* Intermediate arrays */
    double *dd; ///< Size n
    double *xse; ///< Size n
    double *d1; ///< Size m + n
    double *d2; ///< Size m + n
    double *daux; ///< Size m + n
    
    /* Consecutive memory for [dx; dy; ds] */
    double *dx;
    double *dy;
    double *ds;
    
    /* Signal ill-conditioning Newton*/
    int badNewton;

} lp_newton;

extern pot_int LpNewtonCreate( lp_newton **pnewton );
extern pot_int LpNewtonInit( lp_newton *newton, pot_int nCol, pot_int nRow, pot_int *colMatBeg, pot_int *colMatIdx, double *colMatElem );
extern pot_int LpNewtonOneStep( lp_newton *newton, double *lpObj, double *lpRHS, int *colMatBeg, int *colMatIdx, double *colMatElem,
                                double *colVal, double *rowDual, double *colDual, double *kappa, double *tau,
                                double *pRes, double *dRes, double pObjVal, double dObjVal, double spxSize );
extern void LpNewtonClear( lp_newton *newton );
extern void LpNewtonDestroy( lp_newton **pnewton );


#endif /* lp_newton_h */
