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
    double *rowDualStep = NULL;
    double logdet = 0.0;
    
    double *kktLhsBuffer = NULL;
    
    hdsdp_cone **SDPCones = NULL;
    user_data *SDPData = NULL;
    hdsdp_cone *SDPCone = NULL;
    hdsdp_kkt *kkt = NULL;
    
    double timeStart = HUtilGetTimeStamp();
    printf("Filename: %s\n", fname);
    
    HDSDP_CALL(HReadSDPA(fname, &nConstrs, &nBlks, &BlkDims, &rowRHS, &coneMatBeg,
                         &coneMatIdx, &coneMatElem, &nCols, &nLpCols, &LpMatBeg,
                         &LpMatIdx, &LpMatElem, &nElem));
    
    printf("Reading SDPA file in %f seconds \n", HUtilGetTimeStamp() - timeStart);
    
    HDSDP_CALL(HUserDataCreate(&SDPData));
    
    HDSDP_INIT(rowDual, double, nConstrs);
    HDSDP_MEMCHECK(rowDual);
    
    HDSDP_INIT(rowDualStep, double, nConstrs);
    HDSDP_MEMCHECK(rowDualStep);
    
    HDSDP_INIT(kktLhsBuffer, double, nConstrs);
    HDSDP_MEMCHECK(kktLhsBuffer);
    
    HDSDP_INIT(SDPCones, hdsdp_cone *, nBlks);
    HDSDP_MEMCHECK(SDPCones);
    
    for ( int iBlk = 0; iBlk < nBlks; ++iBlk ) {
        HUserDataSetConeData(SDPData, HDSDP_CONETYPE_DENSE_SDP, nConstrs, BlkDims[iBlk],
                             coneMatBeg[iBlk], coneMatIdx[iBlk], coneMatElem[iBlk]);
        cone_type cone = HUserDataChooseCone(SDPData);
        HDSDP_CALL(HConeCreate(&SDPCone));
        
        SDPCones[iBlk] = SDPCone;
        
        HDSDP_CALL(HConeSetData(SDPCone, SDPData));
        HDSDP_CALL(HConeProcData(SDPCone));
//        HConeView(SDPCone);
        HDSDP_CALL(HConePresolveData(SDPCone));
        HConeView(SDPCone);
        
        for ( int i = 0; i < nConstrs; ++i ) {
            rowDual[i] = 0.0 * (double) (i + 1) / nConstrs;
            rowDualStep[i] = (double) (i + 1) / nConstrs;
        }
        
        HConeSetStart(SDPCone, -1e+03);
        HConeUpdate(SDPCone, 1.0, rowDual);
//        HConeView(SDPCone);
        
        HDSDP_CALL(HConeGetLogBarrier(SDPCone, 1.5, rowDual, BUFFER_DUALVAR, &logdet));
//        printf("- Conic log det (S) = %e. \n", logdet);
        
        double ratio = 0.0;
        HDSDP_CALL(HConeRatioTest(SDPCone, 1.0, rowDualStep, 0.9, BUFFER_DUALVAR, &ratio));
        printf("- Conic ratio test: %e. \n", ratio);
        HDSDP_CALL(HConeRatioTest(SDPCone, 1.0, rowDualStep, 0.9, BUFFER_DUALVAR, &ratio));
        printf("- Conic ratio test again: %e. \n", ratio);
        
        HUserDataClear(SDPData);
    }
    
    /* KKT setup */
    HDSDP_CALL(HKKTCreate(&kkt));
    HDSDP_CALL(HKKTInit(kkt, nConstrs, nBlks, SDPCones));
    
    HDSDP_CALL(HKKTBuildUp(kkt, KKT_TYPE_INFEASIBLE));
//    HDSDP_PROFILER(HKKTBuildUp(kkt, KKT_TYPE_INFEASIBLE), 100);
    
    /* KKT solve */
    double dCSinv = 0.0;
    double dCSinvRdCSinv = 0.0;
    double dCSinvCSinv = 0.0;
    
    HKKTExport(kkt, kktLhsBuffer, NULL, NULL, &dCSinvCSinv, &dCSinv, &dCSinvRdCSinv, NULL);
    HDSDP_CALL(HKKTFactorize(kkt));
    HDSDP_CALL(HKKTSolve(kkt, kktLhsBuffer, NULL));
    
    HKKTExport(kkt, NULL, kktLhsBuffer, NULL, &dCSinvCSinv, &dCSinv, &dCSinvRdCSinv, NULL);
    HDSDP_CALL(HKKTFactorize(kkt));
    HDSDP_CALL(HKKTSolve(kkt, kktLhsBuffer, NULL));
    
    HKKTExport(kkt, NULL, NULL, kktLhsBuffer, &dCSinvCSinv, &dCSinv, &dCSinvRdCSinv, NULL);
    HDSDP_CALL(HKKTFactorize(kkt));
    HDSDP_CALL(HKKTSolve(kkt, kktLhsBuffer, NULL));
    
    /* KKT consistency */
    HDSDP_CALL(HUtilKKTCheck(kkt));
    
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
    HDSDP_FREE(rowDualStep);
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

int test_solver( char *fname ) {
    
    
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
    
    double dErrors[6] = {0.0};
    
    hdsdp *hsolve = NULL;
    user_data **SDPDatas = NULL;
    user_data *SDPData = NULL;
    
    double timeStart = HUtilGetTimeStamp();
    printf("Filename: %s\n", fname);
    
    HDSDP_CALL(HReadSDPA(fname, &nConstrs, &nBlks, &BlkDims, &rowRHS, &coneMatBeg,
                         &coneMatIdx, &coneMatElem, &nCols, &nLpCols, &LpMatBeg,
                         &LpMatIdx, &LpMatElem, &nElem));
    
    if ( nLpCols > 0 ) {
        printf("LP Cone is being developed. \n");
        goto exit_cleanup;
    }
    
    printf("Reading SDPA file in %f seconds \n", HUtilGetTimeStamp() - timeStart);
    
    HDSDP_INIT(SDPDatas, user_data *, nBlks);
    HDSDP_MEMCHECK(SDPDatas);
    
    HDSDP_CALL(HDSDPCreate(&hsolve));
    HDSDP_CALL(HDSDPInit(hsolve, nConstrs, nBlks));
    
    for ( int iBlk = 0; iBlk < nBlks; ++iBlk ) {
        HDSDP_CALL(HUserDataCreate(&SDPData));
        HUserDataSetConeData(SDPData, HDSDP_CONETYPE_DENSE_SDP, nConstrs, BlkDims[iBlk],
                             coneMatBeg[iBlk], coneMatIdx[iBlk], coneMatElem[iBlk]);
        HDSDP_CALL(HDSDPSetCone(hsolve, iBlk, SDPData));
        SDPDatas[iBlk] = SDPData;
        SDPData = NULL;
    }
    
    HDSDPSetDualObjective(hsolve, rowRHS);
    HDSDP_CALL(HDSDPOptimize(hsolve, 1));
    
exit_cleanup:
    
    for ( int iBlk = 0; iBlk < nBlks; ++iBlk ) {
        HUserDataDestroy(&SDPDatas[iBlk]);
    }
    HDSDP_FREE(SDPDatas);
    
    HUserDataDestroy(&SDPData);
    HDSDPDestroy(&hsolve);
    
    HDSDP_FREE(BlkDims);
    HDSDP_FREE(rowRHS);
    
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

