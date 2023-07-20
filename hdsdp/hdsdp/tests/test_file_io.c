#include <stdio.h>

#ifdef HEADERPATH
#include "interface/hdsdp.h"
#include "interface/hdsdp_utils.h"
#include "interface/hdsdp_user_data.h"
#include "interface/hdsdp_file_io.h"
#include "interface/hdsdp_conic.h"
#include "interface/hdsdp_schur.h"
#else
#include "hdsdp.h"
#include "hdsdp_utils.h"
#include "hdsdp_user_data.h"
#include "hdsdp_file_io.h"
#include "hdsdp_conic.h"
#include "hdsdp_schur.h"
#endif

#include <math.h>

hdsdp_retcode test_schur_consistency( hdsdp_kkt *kkt ) {
    
    hdsdp_retcode retcode = HDSDP_RETCODE_OK;
        
    double *kktBuffer2345 = NULL;
    double *kktBuffer3 = NULL;
    double *kktBuffer4 = NULL;
    
    int nKKTNzs = 0;
    
    if ( kkt->isKKTSparse ) {
        nKKTNzs = kkt->kktMatBeg[kkt->nRow];
    } else {
        nKKTNzs = kkt->nRow * kkt->nRow;
    }
    
    HDSDP_INIT(kktBuffer2345, double, nKKTNzs);
    HDSDP_MEMCHECK(kktBuffer2345);
    
    HDSDP_INIT(kktBuffer3, double, nKKTNzs);
    HDSDP_MEMCHECK(kktBuffer3);
    
    HDSDP_INIT(kktBuffer4, double, nKKTNzs);
    HDSDP_MEMCHECK(kktBuffer4);
    
    /* Use hybrid strategy as a benchmark */
    HDSDP_CALL(HKKTBuildUpFixed(kkt, KKT_TYPE_HOMOGENEOUS, KKT_M3));
    HDSDP_MEMCPY(kktBuffer3, kkt->kktMatElem, double, nKKTNzs);
    
    HDSDP_CALL(HKKTBuildUpFixed(kkt, KKT_TYPE_HOMOGENEOUS, KKT_M4));
    HDSDP_MEMCPY(kktBuffer4, kkt->kktMatElem, double, nKKTNzs);
    
    HDSDP_CALL(HKKTBuildUp(kkt, KKT_TYPE_INFEASIBLE));
    HDSDP_MEMCPY(kktBuffer2345, kkt->kktMatElem, double, nKKTNzs);
    
    double err3_4 = 0.0;
    double err3_2345 = 0.0;
    double err4_2345 = 0.0;
    
    for ( int iElem = 0; iElem < nKKTNzs; ++iElem ) {
        err3_4 += fabs(kktBuffer3[iElem] - kktBuffer4[iElem]);
        err3_2345 += fabs(kktBuffer3[iElem] - kktBuffer2345[iElem]);
        err4_2345 += fabs(kktBuffer4[iElem] - kktBuffer2345[iElem]);
    }
    
    printf("KKT consistency check: | 3-4: %6.3e | 3-2345: %6.3e | 4-2345: %6.3e \n",
           err3_4, err3_2345, err4_2345);
    
exit_cleanup:
    
    HDSDP_FREE(kktBuffer2345);
    HDSDP_FREE(kktBuffer3);
    HDSDP_FREE(kktBuffer4);
    
    return retcode;
}

int test_file_io( char *fname ) {
    
    
    hdsdp_retcode retcode = HDSDP_RETCODE_OK;
    
    int nConstrs = 0;
    int nBlks = 0;
    int *BlkDims = NULL;
    double *rowRHS = NULL;
    int **coneMatBeg = NULL;
    int **coneMatIdx = NULL;
    double **coneMatElem = NULL;
    int nLpCols = 0;
    int *LpMatBeg = NULL;
    int *LpMatIdx = NULL;
    double *LpMatElem = NULL;
    int nCols = 0;
    int nElem = 0;
    double *rowDual = NULL;
    double logdet = 0.0;
    int nCones = 1;
    
    double *kktLhsBuffer = NULL;
    
    hdsdp_cone **SDPCones = NULL;
    user_data *SDPData = NULL;
    hdsdp_cone *SDPCone = NULL;
    hdsdp_kkt *kkt = NULL;
    
    double timeStart = HUtilGetTimeStamp();
    
    HDSDP_CALL(HReadSDPA(fname, &nConstrs, &nBlks, &BlkDims, &rowRHS, &coneMatBeg,
                         &coneMatIdx, &coneMatElem, &nCols, &nLpCols, &LpMatBeg,
                         &LpMatIdx, &LpMatElem, &nElem));
    
    printf("Reading SDPA file in %f seconds. \n", HUtilGetTimeStamp() - timeStart);
    
    HDSDP_CALL(HUserDataCreate(&SDPData));
    HDSDP_INIT(rowDual, double, nConstrs);
    HDSDP_MEMCHECK(rowDual);
    
    HDSDP_INIT(kktLhsBuffer, double, nConstrs);
    HDSDP_MEMCHECK(kktLhsBuffer);
    
    HDSDP_INIT(SDPCones, hdsdp_cone *, nBlks);
    HDSDP_MEMCHECK(SDPCones);
    
    for ( int iBlk = 0; iBlk < nBlks; ++iBlk ) {
        HUserDataSetConeData(SDPData, HDSDP_CONETYPE_DENSE_SDP, nConstrs, BlkDims[iBlk],
                             coneMatBeg[iBlk], coneMatIdx[iBlk], coneMatElem[iBlk]);
        cone_type cone = HUserDataChooseCone(SDPData);
        HDSDP_CALL(HConeCreate(&SDPCone));
        HDSDP_CALL(HConeSetData(SDPCone, SDPData));
        HDSDP_CALL(HConeProcData(SDPCone));
//        HConeView(SDPCone);
        HDSDP_CALL(HConePresolveData(SDPCone));
//        HConeView(SDPCone);
        
        for ( int i = 0; i < nConstrs; ++i ) {
            rowDual[i] = (double) (i + 1) / nConstrs;
        }
        
        HConeSetStart(SDPCone, -1e+05);
        HConeUpdate(SDPCone, 1.5, rowDual);
//        HConeView(SDPCone);
        
        HDSDP_CALL(HConeGetLogBarrier(SDPCone, 1.5, rowDual, &logdet));
        printf("- Conic log det (S) = %e. \n", logdet);
        
        SDPCones[iBlk] = SDPCone;
        HUserDataClear(SDPData);
    }
    
    /* KKT setup */
    HDSDP_CALL(HKKTCreate(&kkt));
    HDSDP_CALL(HKKTInit(kkt, nConstrs, nBlks, SDPCones));
    
    HDSDP_CALL(HKKTBuildUp(kkt, KKT_TYPE_INFEASIBLE));
    
    /* KKT solve */
    double dCSinv = 0.0;
    double dCSinvRdCSinv = 0.0;
    double dCSinvCSinv = 0.0;
    
    HKKTExport(kkt, kktLhsBuffer, NULL, NULL, &dCSinvCSinv, &dCSinv, &dCSinvRdCSinv);
    HDSDP_CALL(HKKTFactorize(kkt));
    HDSDP_CALL(HKKTSolve(kkt, kktLhsBuffer, NULL));
    
    HKKTExport(kkt, NULL, kktLhsBuffer, NULL, &dCSinvCSinv, &dCSinv, &dCSinvRdCSinv);
    HDSDP_CALL(HKKTFactorize(kkt));
    HDSDP_CALL(HKKTSolve(kkt, kktLhsBuffer, NULL));
    
    HKKTExport(kkt, NULL, NULL, kktLhsBuffer, &dCSinvCSinv, &dCSinv, &dCSinvRdCSinv);
    HDSDP_CALL(HKKTFactorize(kkt));
    HDSDP_CALL(HKKTSolve(kkt, kktLhsBuffer, NULL));
    
    /* KKT consistency */
    // HDSDP_CALL(test_schur_consistency(kkt));
    
exit_cleanup:
    
    HDSDP_FREE(kktLhsBuffer);
    
    HUserDataDestroy(&SDPData);
    HKKTDestroy(&kkt);
    
    for ( int iBlk = 0; iBlk < nBlks; ++iBlk ) {
        SDPCone = SDPCones[iBlk];
        HConeDestroy(&SDPCone);
    }
    
    HDSDP_FREE(BlkDims);
    HDSDP_FREE(rowRHS);
    HDSDP_FREE(rowDual);
    HDSDP_FREE(SDPCones);
    
    for ( int iBlk = 0; iBlk < nBlks; ++iBlk ) {
        HDSDP_FREE(coneMatBeg[iBlk]);
        HDSDP_FREE(coneMatIdx[iBlk]);
        HDSDP_FREE(coneMatElem[iBlk]);
    }
    
    HDSDP_FREE(coneMatBeg);
    HDSDP_FREE(coneMatIdx);
    HDSDP_FREE(coneMatElem);
    
    if ( nLpCols > 0 ) {
        HDSDP_FREE(LpMatBeg);
        HDSDP_FREE(LpMatIdx);
        HDSDP_FREE(LpMatElem);
    }
    
    return retcode;
}

