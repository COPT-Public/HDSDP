#include <stdio.h>
#include "dsdputils.h"
#include "dsdpsolver.h"

/* Recording prototype methods that are no longer in use */

/* Check if a dense matrix is rank-one
   Depreciated due to inaccuracy
 */
static DSDP_INT isDenseRank1InAcc( dsMat *dataMat, DSDP_INT *isRank1 ) {
    // Check if a dense matrix is rank-one
    // This version is fast but not accurate due to potential numerical error
    DSDP_INT retcode = DSDP_RETCODE_OK;
    double *A    = dataMat->array;
    DSDP_INT n   = dataMat->dim;
    DSDP_INT r1  = TRUE;
    DSDP_INT col = 0;
    
    double benchCol  = 0.0;
    double scaleCol  = 0.0;
    double benchCol2 = 0.0;
    double scaleCol2 = 0.0;
    double diff      = 0.0;
    
    // Get the first column that contains non-zero elements
    for (DSDP_INT i = 0; i < n; ++i) {
        col = i;
        if (packIdx(A, n, i, i) != 0) {
            break;
        }
    }
    
    assert( col != n - 1 ); // or it is a zero matrix
    
    // Check the scaling coefficient
    for (DSDP_INT i = col + 1; i < n; ++i) {
        scaleCol = packIdx(A, n, i, col);
        for (DSDP_INT j = col; j < n; ++j) {
            benchCol2 = packIdx(A, n, j, col);
            if (i <= j) {
                scaleCol2 = packIdx(A, n, j, i);
            } else {
                scaleCol2 = packIdx(A, n, i, j);
            }
            diff = benchCol * scaleCol2 - benchCol2 * scaleCol;
            if (fabs(diff) > 1e-04 * MAX(1, benchCol)) {
                r1 = FALSE;
                break;
            }
        }
    }
    
    *isRank1 = r1;
    
    return retcode;
}
