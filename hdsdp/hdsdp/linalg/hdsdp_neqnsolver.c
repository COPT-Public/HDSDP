#ifdef HEADERPATH
#include "linalg/def_hdsdp_neqnsolver.h"
#include "linalg/hdsdp_neqnsolver.h"
#include "linalg/vec_opts.h"
#include "interface/hdsdp_utils.h"
#else
#include "def_hdsdp_neqnsolver.h"
#include "hdsdp_neqnsolver.h"
#include "vec_opts.h"
#include "hdsdp_utils.h"
#endif

static double HNEquationIGetDenseColThreshold( int nCol ) {
    
    if ( nCol > 100000 ) {
        return 0.001;
    }
    
    if ( nCol > 50000 ) {
        return 0.002;
    }
    
    if ( nCol > 10000 ) {
        return 0.01;
    }
    
    if ( nCol > 2000 ) {
        return 0.05;
    }
    
    if ( nCol > 500 ) {
        return 0.1;
    }
    
    return 1.0;
}

static hdsdp_retcode HNEquationIAnalyze( hdsdp_normal_linsys *HNeq ) {
    
    hdsdp_retcode retcode = HDSDP_RETCODE_OK;
    
    /* Get sparsity of each column */
    HDSDP_INIT(HNeq->colSparsity, int, HNeq->nCol);
    for ( int iCol = 0; iCol < HNeq->nCol; ++iCol ) {
        HNeq->colSparsity[iCol] = HNeq->colMatBeg[iCol + 1] - HNeq->colMatBeg[iCol];
    }
    
    /* Sort sparsity of columns */
    HUtilDescendSortIntByInt(HNeq->colSparsity, HNeq->colSparsitySortedId, 0, HNeq->nCol);
    
    /* Threshold columns */
    double dColThresh = HNEquationIGetDenseColThreshold(HNeq->nCol);
    int iDenseColThresh = (int) dColThresh * HNeq->nCol;
    
    /* Detect dense columns up to threshold */
    int nDenseCol = 0;
    for ( int iCol = 0; iCol < iDenseColThresh; ++iCol ) {
        if ( HNeq->colSparsity[HNeq->nCol - iCol - 1] < iDenseColThresh ) {
            break;
        }
        nDenseCol += 1;
    }
    
    hdsdp_printf("Num. Dense cols: %d \n", nDenseCol);
    int nSparseCol = HNeq->nCol - nDenseCol;
    
    /* Allocate memory for dense and sparse columns respectively */
    
    
exit_cleanup:
    return retcode;
}

extern hdsdp_retcode HNEquationCreate( hdsdp_normal_linsys **pHNeq ) {
    
    hdsdp_retcode retcode = HDSDP_RETCODE_OK;
    
    if ( !pHNeq ) {
        retcode = HDSDP_RETCODE_FAILED;
        goto exit_cleanup;
    }
    
    hdsdp_normal_linsys *HNeq = NULL;
    HDSDP_INIT(HNeq, hdsdp_normal_linsys, 1);
    HDSDP_MEMCHECK(HNeq);
    HDSDP_ZERO(HNeq, hdsdp_normal_linsys, 1);
    *pHNeq = HNeq;
    
exit_cleanup:
    return retcode;
}

extern hdsdp_retcode HNEquationInit( hdsdp_normal_linsys *HNeq, int nRow, int nCol ) {
    
    hdsdp_retcode retcode = HDSDP_RETCODE_OK;
    
    HNeq->nRow = nRow;
    HNeq->nCol = nCol;
    
    HDSDP_INIT(HNeq->colSparsity, int, nCol);
    HDSDP_MEMCHECK(HNeq->colSparsity);
    
    HDSDP_INIT(HNeq->colSparsitySortedId, int, nCol);
    HDSDP_MEMCHECK(HNeq->colSparsitySortedId);
    
    HDSDP_INIT(HNeq->dColBuffer, double, nCol);
    HDSDP_MEMCHECK(HNeq->dColBuffer);
    
    HDSDP_INIT(HNeq->dRowBuffer, double, nRow);
    HDSDP_MEMCHECK(HNeq->dRowBuffer);
    
    HDSDP_CALL(HFpLinsysCreate(&HNeq->spchol, nCol, HDSDP_LINSYS_SPARSE_DIRECT));
    
exit_cleanup:
    return retcode;
}

extern hdsdp_retcode HNEquationSetData( hdsdp_normal_linsys *HNeq, int *colMatBeg, int *colMatIdx, double *colMatElem ) {
    
    hdsdp_retcode retcode = HDSDP_RETCODE_OK;
    
    HNeq->colMatBeg = colMatBeg;
    HNeq->colMatIdx = colMatIdx;
    HNeq->colMatElem = colMatElem;
    
    hdsdp_printf("Analyzing normal equation. \n");
    HDSDP_CALL(HNEquationIAnalyze(HNeq));
 
exit_cleanup:
    return retcode;
}

extern void HNEquationClear( hdsdp_normal_linsys *HNeq ) {
    
    if ( !HNeq ) {
        return;
    }
    
    HDSDP_FREE(HNeq->colSparsity);
    HDSDP_FREE(HNeq->colSparsitySortedId);
    
    HDSDP_FREE(HNeq->spColIdx);
    HDSDP_FREE(HNeq->colMatSparseBeg);
    HDSDP_FREE(HNeq->colMatSparseIdx);
    HDSDP_FREE(HNeq->colMatSparseElem);
    
    HDSDP_FREE(HNeq->dsColIdx);
    HDSDP_FREE(HNeq->colMatDenseElem);
    
    HDSDP_FREE(HNeq->schurMatBeg);
    HDSDP_FREE(HNeq->schurMatIdx);
    HDSDP_FREE(HNeq->schurMatElem);
    
    HDSDP_FREE(HNeq->dColBuffer);
    HDSDP_FREE(HNeq->dRowBuffer);
    
    HFpLinsysDestroy(&HNeq->spchol);
    
    HDSDP_ZERO(HNeq, hdsdp_normal_linsys, 1);
    
    return;
}

extern void HNEquationDestroy( hdsdp_normal_linsys **pHNeq ) {
    
    if ( !pHNeq ) {
        return;
    }
    
    HNEquationClear(*pHNeq);
    HDSDP_FREE(*pHNeq);
    
    return;
}

