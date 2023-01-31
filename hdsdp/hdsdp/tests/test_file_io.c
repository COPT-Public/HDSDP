#include <stdio.h>

#include "interface/hdsdp.h"
#include "interface/hdsdp_utils.h"
#include "interface/hdsdp_file_io.h"

int test_file_io( void ) {
    
    
    hdsdp_retcode retcode = HDSDP_RETCODE_OK;
    
    char *fname = "/Users/gaowenzhi/Desktop/dsdp6/hsd/benchmark/sdplib/hinf4.dat-s";
    
    int nConstrs = 0;
    int nBlks = 0;
    int *nBlkDims = NULL;
    double *rowRHS = NULL;
    int **coneMatBeg = NULL;
    int **coneMatIdx = NULL;
    double **coneMatElem = NULL;
    int nLpCols = 0;
    int *LpMatBeg = NULL;
    int *LpMatIdx = NULL;
    double *LpMatElem = NULL;
    int nElem = 0;
    
    HDSDP_CALL(HReadSDPA(fname, &nConstrs, &nBlks, &nBlkDims, &rowRHS, &coneMatBeg,
                         &coneMatIdx, &coneMatElem, &nLpCols, &LpMatBeg, &LpMatIdx, &LpMatElem, &nElem));
    
    HDSDP_FREE(nBlkDims);
    HDSDP_FREE(rowRHS);
    
    for ( int iBlk = 0; iBlk < nBlks; ++iBlk ) {
        HDSDP_FREE(coneMatBeg[iBlk]);
        HDSDP_FREE(coneMatIdx[iBlk]);
        HDSDP_FREE(coneMatElem[iBlk]);
    }
    
    HDSDP_FREE(coneMatBeg);
    HDSDP_FREE(coneMatIdx);
    HDSDP_FREE(coneMatElem);
    
    HDSDP_FREE(LpMatBeg);
    HDSDP_FREE(LpMatIdx);
    HDSDP_FREE(LpMatElem);
    
exit_cleanup:
    
    return retcode;
}

