#ifndef lp_qmatrix_h
#define lp_qmatrix_h

#include "pot_structs.h"
#include "pot_vector.h"

typedef struct {
    
    pot_int nColQ;
    pot_int nRowQ;
    
    pot_int *QMatBeg;
    pot_int *QMatIdx;
    double  *QMatElem;
    
    double  *sclCol;
    double  *sclRow;

} lp_qmatrix;

extern pot_int LPQMatCreate( lp_qmatrix **pQMat );
extern pot_int LPQMatInit( lp_qmatrix *QMat, pot_int nCol, pot_int nRow, pot_int *colMatBeg );
extern pot_int LPQMatSetup( lp_qmatrix *QMat, pot_int nCol, pot_int nRow, pot_int *colMatBeg,
                            pot_int *colMatIdx, double *colMatElem, double *lpObj, double *lpRHS );
extern pot_int LPQMatRuizScal( lp_qmatrix *QMat, int maxIter );
extern pot_int LPQMatPCScal( lp_qmatrix *QMat, int maxIter );
extern pot_int LPQMatL2Scal( lp_qmatrix *QMat );
extern void LPQMatMultiply( lp_qmatrix *QMat, int *isColBasic, double *xVal, double *qxVal );
extern void LPQMatTransMultiply( lp_qmatrix *QMat, int *isColBasic, double *yVal, double *qtyVal );
extern void LPQMatScal( lp_qmatrix *QMat, double *xVal );
extern void LPQMatScalBack( lp_qmatrix *QMat, double *xVal );
extern void LPQMatClear( lp_qmatrix *QMat );
extern void LPQMatDestroy( lp_qmatrix **pQMat );

#endif /* lp_qmatrix_h */
