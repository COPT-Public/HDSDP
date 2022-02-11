#include <stdio.h>
#include "dsdphsd.h"
#include "dsdputils.h"
#include "dsdpsolver.h"

static char etype[] = "Prototype methods";
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

/* Set up Schur matrix by simple loop
   Depreciated due to low speed
 */

static DSDP_INT setupSDPSchurBlock( HSDSolver *dsdpSolver, DSDP_INT blockid ) {
    /* Set up the schur matrix for a given block. This routine updates
     
     1. Schur matrix Msdp
     2. asinv
     3. u
     4. d4 ( = ASinvRyASinv by the event indicator)
    */
    
    DSDP_INT retcode = DSDP_RETCODE_OK;
    assert( blockid < dsdpSolver->nBlock );
    
    DSDP_INT m       = dsdpSolver->m;
    DSDP_INT mattype = MAT_TYPE_UNKNOWN;
    sdpMat *sdpData = dsdpSolver->sdpData[blockid];
    
    // Temporary storage array for SinvASinv
    rkMat *rkdata = dsdpSolver->rkaux[blockid];
    dsMat *dsdata = dsdpSolver->dsaux[blockid];
    
    DSDP_INT dim = 0;
    
    if (dsdpSolver->eventMonitor[EVENT_IN_PHASE_B]) {
        dim = m;
    } else {
        dim = m + 1;
    }
    
    void *data = NULL;
    DSDP_INT i, j;
    for (i = 0; i < dim; ++i) {
        // Compute SinvASinv
        mattype = sdpData->types[i];
        
        if (mattype == MAT_TYPE_ZERO) {
            continue;
        } else if (mattype == MAT_TYPE_RANKK) {
            retcode = getSinvASinv(dsdpSolver, blockid, i, rkdata);
            data = (void *) rkdata;
            mattype = MAT_TYPE_RANKK;
            checkCode;
        } else {
            retcode = getSinvASinv(dsdpSolver, blockid, i, dsdata);
            data = (void *) dsdata;
            mattype = MAT_TYPE_DENSE;
            checkCode;
        }
        
        for (j = 0; j <= i; ++j) {
            getTraceASinvASinv(dsdpSolver, blockid, j, mattype, i, data);
        }
    }
    
    if (dsdpSolver->eventMonitor[EVENT_IN_PHASE_B] &&
        !dsdpSolver->eventMonitor[EVENT_INVALID_GAP]) {
        
        double maxdiag = 0.0;
        for (DSDP_INT i = 0; i < m; ++i) {
            maxdiag = MAX(packIdx(dsdpSolver->Msdp->array, m, i, i), maxdiag);
        }
        
        dsdpSolver->Mscaler = maxdiag;
        
        for (DSDP_INT i = 0; i < m; ++i) {
            packIdx(dsdpSolver->Msdp->array, m, i, i) += \
            MIN(maxdiag * 1e-06, 1e-08);
        }
    }
    
    // Perturb in Phase A
//    if (dsdpSolver->eventMonitor[EVENT_IN_PHASE_A]) {
//        for (DSDP_INT i = 0; i < m; ++i) {
//            packIdx(dsdpSolver->Msdp->array, m, i, i) += \
//            MIN(maxdiag * 1e-05, 1e-05);
//        }
//    }
    
    assert( retcode == DSDP_RETCODE_OK );
    return retcode;
}

/* Evaluate a single element of the Schur matrix
   Depreciated due to low speed
 */
extern DSDP_INT getTraceASinvASinvSlow( HSDSolver *dsdpSolver, DSDP_INT blockid, DSDP_INT constrid,
                                    DSDP_INT mattype, DSDP_INT constrid2, void *SinvASinv ) {
    
    // Compute trace between SinvASinv and some A
    // constrid is the position of (A in A * (.)) and constrid2 is the position of (A in (.) * SinvASinv)
    DSDP_INT retcode = DSDP_RETCODE_OK;
    DSDP_INT Atype = dsdpSolver->sdpData[blockid]->types[constrid];
    DSDP_INT m = dsdpSolver->m;
    void *A = dsdpSolver->sdpData[blockid]->sdpData[constrid];
    double trace = 0.0;
    double *M = dsdpSolver->Msdp->array;
    
    // assert( constrid <= constrid2 );
    // assert( constrid <= m && constrid2 <= m);
    
    if (mattype == MAT_TYPE_ZERO || Atype == MAT_TYPE_ZERO) {
        return retcode;
    }
    
    retcode = rkMatrkTrace((rkMat *) SinvASinv, (rkMat *) A, &trace);
    
    if (mattype == MAT_TYPE_RANKK) {
        switch (Atype) {
            case MAT_TYPE_DENSE:
                retcode = rkMatdenseTrace((rkMat *) SinvASinv, (dsMat *) A, &trace);
                // retcode = r1MatdenseTrace((r1Mat *) SinvASinv, (dsMat *) A, &trace);
                break;
            case MAT_TYPE_SPARSE:
                retcode = rkMatspsTrace((rkMat *) SinvASinv, (spsMat *) A, &trace);
                // retcode = r1MatspsTrace((r1Mat *) SinvASinv, (spsMat *) A, &trace);
                break;
            case MAT_TYPE_RANKK:
                retcode = rkMatrkTrace((rkMat *) SinvASinv, (rkMat *) A, &trace);
                // retcode = r1Matr1Trace((r1Mat *) SinvASinv, (r1Mat *) A, &trace);
                break;
            default:
                error(etype, "Invalid matrix type. \n");
                break;
        }
    } else if (mattype == MAT_TYPE_DENSE) {
        switch (Atype) {
            case MAT_TYPE_DENSE:
                retcode = denseDsTrace((dsMat *) SinvASinv, (dsMat *) A, &trace);
                break;
            case MAT_TYPE_SPARSE:
                retcode = denseSpsTrace((dsMat *) SinvASinv, (spsMat *) A, &trace);
                break;
            case MAT_TYPE_RANKK:
                retcode = rkMatdenseTrace((rkMat *) A, (dsMat *) SinvASinv, &trace);
                // retcode = r1MatdenseTrace((r1Mat *) A, (dsMat *) SinvASinv, &trace);
                break;
            default:
                error(etype, "Invalid matrix type. \n");
                break;
        }
    } else {
        error(etype, "Invalid matrix type. \n");
    }
    
    checkCode;
    
    // Perturbation of M diagonal element
    if (constrid == constrid2 &&
        dsdpSolver->eventMonitor[EVENT_IN_PHASE_B]) {
        if (dsdpSolver->mu < 1) {
            dsdpSolver->Msdp->isillCond = TRUE;
        }
        
        if (dsdpSolver->mu < 1e-05) {
            trace += 1e-08;
        }
    }
    
    // Update Schur/auxiliary vectors
    if (constrid2 == m) {
        // The first A is C
        if (constrid == m) {
            dsdpSolver->csinvcsinv += trace;
        } else {
            dsdpSolver->u->x[constrid] += trace;
        }
    } else {
        packIdx(M, m, constrid2, constrid) += trace;
    }
    
    return retcode;
}


extern DSDP_INT getSinvASinvSlow( HSDSolver *dsdpSolver, DSDP_INT blockid, DSDP_INT constrid,
                              void *SinvASinv ) {
    
    // Given S and A, the routine computes A, asinv and trace(S, Sinv A Sinv)
    DSDP_INT retcode, typeA = dsdpSolver->sdpData[blockid]->types[constrid];
    spsMat *S = dsdpSolver->S[blockid];
    void *A = dsdpSolver->sdpData[blockid]->sdpData[constrid];
    double *asinv = NULL, *asinvrysinv = NULL, tracediag = 0.0;
    
    if (constrid < dsdpSolver->m) {
        asinv = &dsdpSolver->asinv->x[constrid];
        asinvrysinv = &dsdpSolver->d4->x[constrid];
    } else {
        asinv = &dsdpSolver->csinv;
        asinvrysinv = &dsdpSolver->csinvrysinv;
    }
    
    switch (typeA) {
        case MAT_TYPE_RANKK:
            retcode = spsSinvRkSinvSolve(S, A, SinvASinv, &tracediag);
            break;
        case MAT_TYPE_SPARSE:
            retcode = spsSinvSpSinvSolve(S, A, SinvASinv, &tracediag);
            break;
        case MAT_TYPE_DENSE:
            retcode = spsSinvDsSinvSolve(S, A, SinvASinv, &tracediag);
            break;
        default:
            error(etype, "Invalid matrix type. \n");
            break;
    }
    
    *asinv += tracediag;
    
    if (dsdpSolver->eventMonitor[EVENT_IN_PHASE_A]) {
        switch (typeA) {
            case MAT_TYPE_RANKK:
                retcode = rkMatdiagTrace(SinvASinv, dsdpSolver->Ry, &tracediag);
                break;
            default:
                retcode = denseDiagTrace(SinvASinv, dsdpSolver->Ry, &tracediag);
                break;
        }
        *asinvrysinv += tracediag;
    }
    return retcode;
}
