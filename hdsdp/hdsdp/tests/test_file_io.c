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
            rowDual[i] = i + 1;
        }
        
        HConeSetStart(SDPCone, -1e+06);
        HConeUpdate(SDPCone, 1.5, rowDual);
//        HConeView(SDPCone);
        
        HDSDP_CALL(HConeGetLogBarrier(SDPCone, 1.5, rowDual, &logdet));
//        printf("- Conic log det (S) = %e. \n", logdet);
        
        SDPCones[iBlk] = SDPCone;
        HUserDataClear(SDPData);
    }
    
    /* KKT setup */
    HDSDP_CALL(HKKTCreate(&kkt));
    HDSDP_CALL(HKKTInit(kkt, nConstrs, nBlks, SDPCones));
    HDSDP_CALL(HKKTBuildUp(kkt, 1));
    
exit_cleanup:
    
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

