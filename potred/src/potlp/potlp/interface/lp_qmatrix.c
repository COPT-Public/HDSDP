/** @brief lp\_qmatrix.c
 * The Q matrix implements the following functionalities and heuristics
 *
 * 1. HSD matrx setup and scaling
 * 2. Matrix operations with residual balancing
 * 3. Matrix operations with reduced support
 *
 */

#include "lp_qmatrix.h"
#include "vec_mat.h"

#include <math.h>

extern pot_int LPQMatCreate( lp_qmatrix **pQMat ) {
    
    pot_int retcode = RETCODE_OK;
    
    if ( !pQMat ) {
        retcode = RETCODE_FAILED;
        goto exit_cleanup;
    }
    
    lp_qmatrix *QMat = NULL;
    POTLP_INIT(QMat, lp_qmatrix, 1);
    
    if ( !QMat ) {
        retcode = RETCODE_FAILED;
        goto exit_cleanup;
    }
    
    POTLP_ZERO(QMat, lp_qmatrix, 1);
    *pQMat = QMat;
    
exit_cleanup:
    return retcode;
}

extern pot_int LPQMatInit( lp_qmatrix *QMat, pot_int nCol, pot_int nRow, pot_int *colMatBeg ) {
    
    pot_int retcode = RETCODE_OK;
    
    if ( !QMat ) {
        retcode = RETCODE_FAILED;
        goto exit_cleanup;
    }
    
    if ( QMat->nColQ > 0 || QMat->nRowQ > 0 ) {
        retcode = RETCODE_FAILED;
        goto exit_cleanup;
    }
    
    /* [y, x, s, kappa, tau ] */
    QMat->nColQ = 2 * nCol + nRow + 2;
    /* [ pRes, dRes, comp ] */
    QMat->nRowQ = nCol + nRow + 1;
    
    POTLP_INIT(QMat->sclCol, double, QMat->nColQ);
    POTLP_INIT(QMat->sclRow, double, QMat->nRowQ);
    
    if ( !QMat->sclCol || !QMat->sclRow ) {
        retcode = RETCODE_FAILED;
        goto exit_cleanup;
    }
    /* Allocate memory for HSD Q matrix
          [  0    A   0   0  -b ]
          [ -A'   0  -I   0   c ]
          [  b'  -c'  0  -1   0 ]
    */
    
    int nColMatElem = colMatBeg[nCol];
    /* 2 * nnz(A) + 2 * m + 3 * n [+ 1] kappa would give an extra entry  */
    int nQMatElem = 2 * nColMatElem + 2 * nRow + 3 * nCol + 1;
    
    POTLP_INIT(QMat->QMatBeg, pot_int, QMat->nColQ + 1);
    POTLP_INIT(QMat->QMatIdx, pot_int, nQMatElem);
    POTLP_INIT(QMat->QMatElem, double, nQMatElem);
    
    if ( !QMat->QMatBeg || !QMat->QMatIdx || !QMat->QMatElem ) {
        retcode = RETCODE_FAILED;
        goto exit_cleanup;
    }
    
exit_cleanup:
    return retcode;
}

extern pot_int LPQMatSetup( lp_qmatrix *QMat, pot_int nCol, pot_int nRow, pot_int *colMatBeg,
                            pot_int *colMatIdx, double *colMatElem, double *lpObj, double *lpRHS ) {
    
    pot_int retcode = RETCODE_OK;
    
    if ( QMat->nRowQ != nRow + nCol + 1 || QMat->nColQ != 2 * nCol + nRow + 2 ) {
        retcode = RETCODE_FAILED;
        goto exit_cleanup;
    }
    
    /* Start concatenating */
    retcode = spMatBuildQMat(QMat->nRowQ, QMat->nColQ, QMat->QMatBeg,
                             QMat->QMatIdx, QMat->QMatElem, nRow, nCol,
                             colMatBeg, colMatIdx, colMatElem, lpRHS, lpObj);
    
    if ( retcode != RETCODE_OK ) {
        goto exit_cleanup;
    }
    
exit_cleanup:
    return retcode;
}

#define RUIZ_DEBUG(format, info) // printf(format, info);
extern pot_int LPQMatRuizScal( lp_qmatrix *QMat, int maxIter ) {
    
    pot_int retcode = RETCODE_OK;
    
    pot_int nRowQ = QMat->nRowQ;
    pot_int nColQ = QMat->nColQ;
    
    pot_int *QMatBeg = QMat->QMatBeg;
    pot_int *QMatIdx = QMat->QMatIdx;
    double  *QMatElem = QMat->QMatElem;
    
    double *ruizScalDiagRow = QMat->sclRow;
    double *ruizScalDiagCol = QMat->sclCol;
    double *ruizWorkDiagRow = NULL;
    double *ruizWorkDiagCol = NULL;
    
    /* Allocate workspace */
    POTLP_INIT(ruizWorkDiagRow, double, nRowQ);
    POTLP_INIT(ruizWorkDiagCol, double, nColQ);
    
    if ( !ruizWorkDiagRow || !ruizWorkDiagCol ) {
        retcode = RETCODE_FAILED;
        goto exit_cleanup;
    }
    
    /* Initialize scalers */
    for ( int i = 0; i < nRowQ; ++i ) {
        ruizScalDiagRow[i] = 1.0;
    }
    
    for ( int i = 0; i < nColQ; ++i ) {
        ruizScalDiagCol[i] = 1.0;
    }
    
    RUIZ_DEBUG("Start Ruiz-scaling %s\n", "");
    
    for ( int i = 0; i < maxIter; ++i ) {
        
        POTLP_ZERO(ruizWorkDiagRow, double, nRowQ);
        spMatMaxRowAbs(nColQ, QMatBeg, QMatIdx, QMatElem, ruizWorkDiagRow);
        POTLP_ZERO(ruizWorkDiagCol, double, nColQ);
        spMatMaxColAbs(nColQ, QMatBeg, QMatIdx, QMatElem, ruizWorkDiagCol);
        
        /* sqrt operation */
        double maxRuizDiagDeviate = 0.0;
        double ruizDiagDeviate = 0.0;
        
        for ( int j = 0; j < nRowQ; ++j ) {
            ruizWorkDiagRow[j] = sqrtl(ruizWorkDiagRow[j]);
            ruizScalDiagRow[j] = ruizScalDiagRow[j] * ruizWorkDiagRow[j];
            ruizDiagDeviate = fabs(ruizWorkDiagRow[j] - 1.0);
            maxRuizDiagDeviate = POTLP_MAX(maxRuizDiagDeviate, ruizDiagDeviate);
        }
        
        for ( int j = 0; j < nColQ; ++j ) {
            ruizWorkDiagCol[j] = sqrtl(ruizWorkDiagCol[j]);
            ruizScalDiagCol[j] = ruizScalDiagCol[j] * ruizWorkDiagCol[j];
            ruizDiagDeviate = fabs(ruizWorkDiagCol[j] - 1.0);
            maxRuizDiagDeviate = POTLP_MAX(maxRuizDiagDeviate, ruizDiagDeviate);
        }
        
        RUIZ_DEBUG("Ruiz Deviation %e \n", maxRuizDiagDeviate);
        
        if ( maxRuizDiagDeviate < 1e-08 ) {
            RUIZ_DEBUG("Ruiz Successfully Ends in %d iterations \n", i);
            break;
        }
        
        /* Scaling */
        spMatRowScal(nColQ, QMatBeg, QMatIdx, QMatElem, ruizWorkDiagRow);
        spMatColScal(nColQ, QMatBeg, QMatIdx, QMatElem, ruizWorkDiagCol);
    }
    
    RUIZ_DEBUG("Ruiz-scaling Ends %s\n", "");
    
    for ( int i = 0; i < nRowQ; ++i ) {
        ruizScalDiagRow[i] = 1.0 / ruizScalDiagRow[i];
    }
    
    for ( int i = 0; i < nColQ; ++i ) {
        ruizScalDiagCol[i] = 1.0 / ruizScalDiagCol[i];
    }
    
exit_cleanup:
    
    POTLP_FREE(ruizWorkDiagRow);
    POTLP_FREE(ruizWorkDiagCol);
    
    return retcode;
}

/** @brief Q matrix multiplication Q x
 *
 *The internal isColBasic allows skipping some of the columns that are considered nonbasic
 *
 * @param[in]  QMat The Q Matrix
 * @param[in]  isColBasic Basis status of the columns
 * @param[in]  xVal x in Qx
 * @param[out] qxVal Qx. Overwritten on exit
 */
extern void LPQMatMultiply( lp_qmatrix *QMat, int *isColBasic, double *xVal, double *qxVal ) {
    
    pot_int *QMatBeg = QMat->QMatBeg;
    pot_int *QMatIdx = QMat->QMatIdx;
    double  *QMatElem = QMat->QMatElem;
    
    POTLP_ZERO(qxVal, double, QMat->nRowQ);
    
    if ( isColBasic ) {

        for ( int i = 0, j; i < QMat->nColQ; ++i ) {
            if ( !isColBasic[i] ) {
                continue;
            }
            
            for ( j = QMatBeg[i]; j < QMatBeg[i + 1]; ++j ) {
                qxVal[QMatIdx[j]] += QMatElem[j] * xVal[i];
            }
        }
        
    } else {
        spMatAxpy(QMat->nColQ, QMatBeg, QMatIdx, QMatElem, 1.0, xVal, qxVal);
    }
    
    return;
}

/** @brief Q matrix transpose multiplication Q'  y
 *
 *The internal isColBasic allows skipping some of the columns that are considered nonbasic
 *
 * @param[in]  QMat The Q Matrix
 * @param[in]  isColBasic Basis status of the columns
 * @param[in]  yVal y in Q' y
 * @param[out] qtyVal Q' y. Overwritten on exit
 */
extern void LPQMatTransMultiply( lp_qmatrix *QMat, int *isColBasic, double *yVal, double *qtyVal ) {
    
    pot_int *QMatBeg = QMat->QMatBeg;
    pot_int *QMatIdx = QMat->QMatIdx;
    double  *QMatElem = QMat->QMatElem;
    
    POTLP_ZERO(qtyVal, double, QMat->nColQ);
    
    if ( isColBasic ) {
        for ( int i = 0, j; i < QMat->nColQ; ++i ) {
            double qty = 0.0;
            if ( isColBasic[i] ) {
                for ( j = QMatBeg[i]; j < QMatBeg[i + 1]; ++j ) {
                    qty += QMatElem[j] * yVal[QMatIdx[j]];
                }
            }
            qtyVal[i] = qty;
        }
    } else {
        spMatATxpy(QMat->nColQ, QMatBeg, QMatIdx, QMatElem, 1.0, yVal, qtyVal);
    }
    
    return;
}

extern void LPQMatScalBack( lp_qmatrix *QMat, double *xVal ) {
    
    vvscl(&QMat->nColQ, QMat->sclCol, xVal);
    
    return;
}

extern void LPQMatClear( lp_qmatrix *QMat ) {
    
    if ( !QMat ) {
        return;
    }
    
    POTLP_FREE(QMat->QMatBeg);
    POTLP_FREE(QMat->QMatIdx);
    POTLP_FREE(QMat->QMatElem);
    
    POTLP_FREE(QMat->sclCol);
    POTLP_FREE(QMat->sclRow);
    
    POTLP_ZERO(QMat, lp_qmatrix, 1);
    
    return;
}

extern void LPQMatDestroy( lp_qmatrix **pQMat ) {
    
    if ( !pQMat ) {
        return;
    }
    
    LPQMatClear(*pQMat);
    POTLP_FREE(*pQMat);
    
    return;
}
