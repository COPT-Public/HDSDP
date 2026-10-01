#ifdef HEADERPATH
#include "interface/hdsdp.h"
#include "interface/hdsdp_utils.h"
#include "interface/hdsdp_user_data.h"
#include "interface/hdsdp_file_io.h"
#include "interface/hdsdp_conic.h"
#include "interface/hdsdp_schur.h"
#include "interface/hdsdp_lpsolve.h"
#include "linalg/sparse_opts.h"
#include "linalg/hdsdp_mineig.h"
#include "external/lp_mps.h"
#include "external/hdsdp_cs.h"
#else
#include "hdsdp.h"
#include "hdsdp_utils.h"
#include "hdsdp_user_data.h"
#include "hdsdp_file_io.h"
#include "hdsdp_conic.h"
#include "hdsdp_schur.h"
#include "hdsdp_mineig.h"
#include "sparse_opts.h"
#include "lp_mps.h"
#include "hdsdp_cs.h"
#include "hdsdp_lpsolve.h"
#endif

#include <stdio.h>
#include <string.h>
#include <math.h>

int test_file_io( char *fname );
int test_sdpa_io( char *fname );
int test_mat( char *path );
int test_mineig( char *path );
int test_primal_primal_dual_bench( char *fname );

#define FILE_TYPE_UNKOWN (0)
#define FILE_TYPE_MPS    (1)
#define FILE_TYPE_SDPA   (2)

static int get_file_type( char *fname ) {
    
    
    char *mpssuf = ".mps";
    char *sdpasuf = ".dat-s";
    
    char *ext = strrchr(fname, '.');
    if ( ext ) {
        if ( strcmp(mpssuf, ext) == 0 ) {
            return FILE_TYPE_MPS;
        } else if ( strcmp(sdpasuf, ext) == 0 ) {
            return FILE_TYPE_SDPA;
        }
    }
    
    return FILE_TYPE_UNKOWN;
}

static void get_sol_status( hdsdp_status status, char *sStatus ) {
    
    char *sOpt = "Optimal";
    char *sNumerical = "Numerical";
    char *sTimeLimit = "Timelimit";
    char *sError = "Error";
    char *sMaxIter = "MaxIter";
    char *sInvalid = "FatalErr";
    
    switch (status) {
        case HDSDP_PRIMAL_DUAL_OPTIMAL:
            strcpy(sStatus, sOpt);
            break;
        case HDSDP_NUMERICAL:
            strcpy(sStatus, sNumerical);
            break;
        case HDSDP_MAXITER:
            strcpy(sStatus, sMaxIter);
            break;
        case HDSDP_TIMELIMIT:
            strcpy(sStatus, sTimeLimit);
            break;
        case HDSDP_INTERNAL_ERROR:
            strcpy(sStatus, sError);
            break;
        default:
            strcpy(sStatus, sInvalid);
            break;
    }
    
    return;
}

static int is_psd( void *chol, int *Ap, int *Ai, double *Ax, int *ispd ) {
    
    hdsdp_retcode retcode = HDSDP_RETCODE_OK;
    hdsdp_linsys *ichol = (hdsdp_linsys *) chol;
    HDSDP_CALL(HFpLinsysPsdCheck(ichol, Ap, Ai, Ax, ispd));
    
exit_cleanup:
    return retcode;
}

static int file_io_mps( char *fname ) {
    
    hdsdp_retcode retcode = HDSDP_RETCODE_OK;
    
    char prob[128] = "?";
    int *Aeqp = NULL;
    int *Aeqi = NULL;
    double *Aeqx = NULL;
    
    int *AeqTransp = NULL;
    int *AeqTransi = NULL;
    double *AeqTransx = NULL;

    int *Aineqp = NULL;
    int *Aineqi = NULL;
    double *Aineqx = NULL;
    
    int *colUbIdx = NULL;
    double *colUbElem = NULL;
    
    int nCol;
    int nRow;
    int nEqRow;
    int nIneqRow;
    int nColUb;
    
    int nElem = 0;
    double *rowRhs = NULL;
    double *colObj = NULL;
    
    int *iTransBuffer = NULL;
    
    /* Reading the standard mps file */
    retcode = (hdsdp_retcode) potLpMpsRead(fname, prob, &nRow, &nEqRow, &nIneqRow, &nCol, &nElem,
                                           &Aeqp, &Aeqi, &Aeqx, &Aineqp, &Aineqi, &Aineqx, &rowRhs,
                                           &colObj, &nColUb, &colUbIdx, &colUbElem);
    
    if ( retcode != HDSDP_RETCODE_OK ) {
        goto exit_cleanup;
    }
    
    assert( nIneqRow == 0 && nColUb == 0 );
    
    HDSDP_INIT(AeqTransp, int, nRow + 1);
    HDSDP_INIT(AeqTransi, int, Aeqp[nCol]);
    HDSDP_INIT(AeqTransx, double, Aeqp[nCol]);
    HDSDP_INIT(iTransBuffer, int, nRow);
    /* Compute transpose of Aeq */
    for ( int iElem = 0; iElem < Aeqp[nCol]; ++iElem ) {
        iTransBuffer[Aeqi[iElem]] += 1;
    }
    
    dcs_cumsum(AeqTransp, iTransBuffer, nRow);
    
    int iPos = 0;
    
    for ( int iCol = 0; iCol < nCol; ++iCol ) {
        for ( int iElem = Aeqp[iCol]; iElem < Aeqp[iCol + 1]; ++iElem ) {
            AeqTransi[iPos = iTransBuffer[Aeqi[iElem]]++] = iCol;
            AeqTransx[iPos] = Aeqx[iElem];
        }
    }
    
    /* Call solver */
    hdsdp_lpsolver *lpsolve = NULL;
    HDSDP_CALL(HLpSolverCreate(&lpsolve));
    HDSDP_CALL(HLpSolverInit(lpsolve, nRow, nCol));
    HDSDP_CALL(HLpSolverSetData(lpsolve, Aeqp, Aeqi, Aeqx, AeqTransp, AeqTransi, AeqTransx, rowRhs, colObj));
    HDSDP_CALL(HLpSolverOptimize(lpsolve));
    
exit_cleanup:
    
    HLpSolverDestroy(&lpsolve);
    
    HDSDP_FREE(Aineqp);
    HDSDP_FREE(Aineqi);
    HDSDP_FREE(Aineqx);
    
    HDSDP_FREE(Aeqp);
    HDSDP_FREE(Aeqi);
    HDSDP_FREE(Aeqx);
    
    HDSDP_FREE(AeqTransp);
    HDSDP_FREE(AeqTransi);
    HDSDP_FREE(AeqTransx);
    HDSDP_FREE(iTransBuffer);
    
    HDSDP_FREE(colUbIdx);
    HDSDP_FREE(colUbElem);
    
    HDSDP_FREE(colObj);
    HDSDP_FREE(rowRhs);
    
    return (int) retcode;
}

static int file_io_sdpa( char *fname ) {
    
    hdsdp_retcode retcode = HDSDP_RETCODE_OK;
    
    int nConstrs = 0;
    int nCones = 0;
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
    
    hdsdp *hsolve = NULL;
    user_data **SDPDatas = NULL;
    user_data *SDPData = NULL;
    user_data *LPData = NULL;
    
    double timeStart = HUtilGetTimeStamp();
    printf("Filename: %s\n", fname);
    
    HDSDP_CALL(HReadSDPA(fname, &nConstrs, &nCones, &BlkDims, &rowRHS, &coneMatBeg,
                         &coneMatIdx, &coneMatElem, &nCols, &nLpCols, &LpMatBeg,
                         &LpMatIdx, &LpMatElem, &nElem));
    
    printf("Reading SDPA file in %f seconds \n", HUtilGetTimeStamp() - timeStart);
    
    HDSDP_INIT(SDPDatas, user_data *, nCones);
    HDSDP_MEMCHECK(SDPDatas);
    
    HDSDP_CALL(HDSDPCreate(&hsolve));
    
    if ( nLpCols > 0 ) {
        HDSDP_CALL(HDSDPInit(hsolve, nConstrs, nCones + 1));
    } else {
        HDSDP_CALL(HDSDPInit(hsolve, nConstrs, nCones));
    }
    
    for ( int iCone = 0; iCone < nCones; ++iCone ) {
        HDSDP_CALL(HUserDataCreate(&SDPData));
        HUserDataSetConeData(SDPData, HDSDP_CONETYPE_DENSE_SDP, nConstrs, BlkDims[iCone],
                             coneMatBeg[iCone], coneMatIdx[iCone], coneMatElem[iCone]);
        HDSDP_CALL(HDSDPSetCone(hsolve, iCone, SDPData));
        SDPDatas[iCone] = SDPData;
        SDPData = NULL;
    }
    
    if ( nLpCols > 0 ) {
        HDSDP_CALL(HUserDataCreate(&LPData));
        HUserDataSetConeData(LPData, HDSDP_CONETYPE_LP, nConstrs, nLpCols, LpMatBeg, LpMatIdx, LpMatElem);
        HDSDP_CALL(HDSDPSetCone(hsolve, nCones, LPData));
    }
    
    HDSDPSetDualObjective(hsolve, rowRHS);
    HDSDP_CALL(HDSDPOptimize(hsolve, 1));
    
exit_cleanup:
    
    for ( int iBlk = 0; iBlk < nCones; ++iBlk ) {
        HUserDataDestroy(&SDPDatas[iBlk]);
    }
    HDSDP_FREE(SDPDatas);
    
    HUserDataDestroy(&SDPData);
    HUserDataDestroy(&LPData);
    
    HDSDPDestroy(&hsolve);
    
    HDSDP_FREE(BlkDims);
    HDSDP_FREE(rowRHS);
    
    for ( int iBlk = 0; iBlk < nCones; ++iBlk ) {
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
    
    return (int) retcode;
}

int test_mineig( char *path ) {
    
    hdsdp_retcode retcode = HDSDP_RETCODE_OK;
    
    int nRow = 0;
    int nCol = 0;
    int *Ap = NULL;
    int *Ai = NULL;
    double *Ax = NULL;
    double *dRhs = NULL;
    double *dLhs = NULL;
    
    double dPerturb = 0.0;
    
    HDSDP_CALL(HUtilGetSparseMatrix(path, &nRow, &nCol, &Ap, &Ai, &Ax, &dRhs));
    HDSDP_INIT(dLhs, double, nCol);
    
    hdsdp_mineig *eig = NULL;
    hdsdp_linsys_fp *chol = NULL;
    
    HDSDP_CALL(HFpLinsysCreate(&chol, nCol, HDSDP_LINSYS_SPARSE_DIRECT));
    HDSDP_CALL(HFpLinsysSymbolic(chol, Ap, Ai));
    
    HDSDP_CALL(HMinEigCreate(&eig));
    HDSDP_CALL(HMinEigInit(eig, nCol, 30));
    
    HDSDP_CALL(HMinEigComputeMinEig(eig, Ap, Ai, Ax, (void *) chol, is_psd,
                                    MIN_EIG_LANCZOS, &dPerturb));
    
exit_cleanup:
    
    if ( retcode != HDSDP_RETCODE_OK ) {
        printf("Failed \n");
    }
    
    HDSDP_FREE(Ap);
    HDSDP_FREE(Ai);
    HDSDP_FREE(Ax);
    HDSDP_FREE(dLhs);
    HDSDP_FREE(dRhs);
    HMinEigDestroy(&eig);
    
    return (int) retcode;
}

int test_mat( char *path ) {
    
    hdsdp_retcode retcode = HDSDP_RETCODE_OK;
    
    int nRow = 0;
    int nCol = 0;
    int *Ap = NULL;
    int *Ai = NULL;
    double *Ax = NULL;
    double *dRhs = NULL;
    double *dLhs = NULL;
    
    HDSDP_CALL(HUtilGetSparseMatrix(path, &nRow, &nCol, &Ap, &Ai, &Ax, &dRhs));
    HDSDP_INIT(dLhs, double, nCol);
    
    double dResidual = 0.0;
    
    hdsdp_linsys_fp *chol = NULL;
    
    HDSDP_CALL(HFpLinsysCreate(&chol, nCol, HDSDP_LINSYS_SPARSE_DIRECT));
    HDSDP_CALL(HFpLinsysSymbolic(chol, Ap, Ai));
    
    double dTimeStart = HUtilGetTimeStamp();
    HDSDP_CALL(HFpLinsysNumeric(chol, Ap, Ai, Ax));
    HDSDP_CALL(HFpLinsysSolve(chol, 1, dRhs, dLhs));
    
    for ( int i = 0, j; i < nCol; ++i ) {
        for ( j = Ap[i]; j < Ap[i + 1]; ++j ) {
            dRhs[Ai[j]] -= dLhs[i] * Ax[j];
            if ( Ai[j] != i ) {
                dRhs[i] -= dLhs[Ai[j]] * Ax[j];
            }
        }
    }
    
    for ( int iCol = 0; iCol < nCol; ++iCol ) {
        dResidual += dRhs[iCol] * dRhs[iCol];
    }
    
    printf("Residual:     %e \n", sqrt(dResidual));
    printf("Elapsed Time: %f seconds \n ", HUtilGetTimeStamp() - dTimeStart);
    
exit_cleanup:
    
    if ( retcode != HDSDP_RETCODE_OK ) {
        printf("Failed \n");
    }
    
    HDSDP_FREE(Ap);
    HDSDP_FREE(Ai);
    HDSDP_FREE(Ax);
    HDSDP_FREE(dLhs);
    HDSDP_FREE(dRhs);
    HFpLinsysDestroy(&chol);
    
    return (int) retcode;
}


int test_file_io( char *fname ) {
    
    int retcode = HDSDP_RETCODE_OK;
    
    int iFileType = get_file_type(fname);
    
    if ( iFileType == FILE_TYPE_MPS ) {
        retcode = file_io_mps(fname);
    } else if ( iFileType == FILE_TYPE_SDPA ) {
        retcode = file_io_sdpa(fname);
    } else {
        hdsdp_printf("Unsupported file type. Only 'mps' and 'dat-s' are supported \n");
    }
    
    return retcode;
}

int test_sdpa_io( char *fname ) {
    
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
//        cone_type cone = HUserDataChooseCone(SDPData);
        HDSDP_CALL(HConeCreate(&SDPCone, iBlk));
        
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
    
    return (int) retcode;
}

int test_primal_primal_dual_bench( char *fname ) {
    
    /* Implement benchmark routine for primal-dual and primal IPM.
     For each problem, we solve the problem 5 times with primal-dual and primal solver. */
    
    hdsdp_retcode retcode = HDSDP_RETCODE_OK;
        
    char prob[128] = "?";
    int *Aeqp = NULL;
    int *Aeqi = NULL;
    double *Aeqx = NULL;
    
    int *AeqTransp = NULL;
    int *AeqTransi = NULL;
    double *AeqTransx = NULL;

    int *Aineqp = NULL;
    int *Aineqi = NULL;
    double *Aineqx = NULL;
    
    int *colUbIdx = NULL;
    double *colUbElem = NULL;
    
    int nCol;
    int nRow;
    int nEqRow;
    int nIneqRow;
    int nColUb;
    
    int nElem = 0;
    double *rowRhs = NULL;
    double *colObj = NULL;
    
    int *iTransBuffer = NULL;
    
    /* Reading the standard mps file */
    retcode = (hdsdp_retcode) potLpMpsRead(fname, prob, &nRow, &nEqRow, &nIneqRow, &nCol, &nElem,
                                           &Aeqp, &Aeqi, &Aeqx, &Aineqp, &Aineqi, &Aineqx, &rowRhs,
                                           &colObj, &nColUb, &colUbIdx, &colUbElem);
    
    if ( retcode != HDSDP_RETCODE_OK ) {
        goto exit_cleanup;
    }
    
    if ( nElem > 1e+07 ) {
        goto exit_cleanup;
    }
    
    assert( nIneqRow == 0 && nColUb == 0 );
    
    HDSDP_INIT(AeqTransp, int, nRow + 1);
    HDSDP_INIT(AeqTransi, int, Aeqp[nCol]);
    HDSDP_INIT(AeqTransx, double, Aeqp[nCol]);
    HDSDP_INIT(iTransBuffer, int, nRow);
    /* Compute transpose of Aeq */
    for ( int iElem = 0; iElem < Aeqp[nCol]; ++iElem ) {
        iTransBuffer[Aeqi[iElem]] += 1;
    }
    
    dcs_cumsum(AeqTransp, iTransBuffer, nRow);
    
    int iPos = 0;
    
    for ( int iCol = 0; iCol < nCol; ++iCol ) {
        for ( int iElem = Aeqp[iCol]; iElem < Aeqp[iCol + 1]; ++iElem ) {
            AeqTransi[iPos = iTransBuffer[Aeqi[iElem]]++] = iCol;
            AeqTransx[iPos] = Aeqx[iElem];
        }
    }
    
    /* Start benchmark */
    int nTest = 10;
    double dPrimalTime = 0.0;
    double dPrimalDualTime = 0.0;
    double dTimeStart = 0.0;
    
    char sPrimalStatus[20] = "";
    char sPrimalDualStatus[20] = "";
    
    int iSuccess = 1;
    int isPrimalUsed = 0;
    
    double dTWarm = HUtilGetTimeStamp();
    
    hdsdp_status iPrimalStatus = HDSDP_UNKNOWN;
    hdsdp_status iPrimalDualStatus = HDSDP_UNKNOWN;
    hdsdp_lpsolver *lpsolve = NULL;
    
    /* Warm-up */
    HDSDP_CALL(HLpSolverCreate(&lpsolve));
    HDSDP_CALL(HLpSolverInit(lpsolve, nRow, nCol));
    HDSDP_CALL(HLpSolverSetData(lpsolve, Aeqp, Aeqi, Aeqx, AeqTransp, AeqTransi, AeqTransx, rowRhs, colObj));
    lpsolve->params.iPrimalMethod = 0;
    retcode = HLpSolverOptimize(lpsolve);
    
    if ( HUtilGetTimeStamp() - dTWarm > 15.0 ) {
        nTest = 1;
    }
    
    if ( HLpSolverGetStatus(lpsolve) == HDSDP_UNKNOWN ) {
        goto exit_cleanup;
    }
    
    HLpSolverDestroy(&lpsolve);
    
    for ( int iTest = 0; iTest < nTest; ++iTest ) {
        
        /* Primal-dual solve */
        if ( iPrimalDualStatus == HDSDP_UNKNOWN || iPrimalDualStatus == HDSDP_PRIMAL_DUAL_OPTIMAL ) {
            
            HDSDP_CALL(HLpSolverCreate(&lpsolve));
            HDSDP_CALL(HLpSolverInit(lpsolve, nRow, nCol));
            HDSDP_CALL(HLpSolverSetData(lpsolve, Aeqp, Aeqi, Aeqx, AeqTransp, AeqTransi, AeqTransx, rowRhs, colObj));
            dTimeStart = HUtilGetTimeStamp();
            lpsolve->params.iPrimalMethod = 0;
            retcode = HLpSolverOptimize(lpsolve);
            dPrimalDualTime += HUtilGetTimeStamp() - dTimeStart;
            iPrimalDualStatus = HLpSolverGetStatus(lpsolve);
            
            if ( iPrimalDualStatus == HDSDP_UNKNOWN ) {
                break;
            }
            
            HLpSolverDestroy(&lpsolve);
        }
        
        /* Primal-solve */
        if ( iPrimalStatus == HDSDP_UNKNOWN || iPrimalStatus == HDSDP_PRIMAL_DUAL_OPTIMAL ) {
            
            HDSDP_CALL(HLpSolverCreate(&lpsolve));
            HDSDP_CALL(HLpSolverInit(lpsolve, nRow, nCol));
            HDSDP_CALL(HLpSolverSetData(lpsolve, Aeqp, Aeqi, Aeqx, AeqTransp, AeqTransi, AeqTransx, rowRhs, colObj));
            dTimeStart = HUtilGetTimeStamp();
            lpsolve->params.iPrimalMethod = 1;
            retcode = HLpSolverOptimize(lpsolve);
            dPrimalTime += HUtilGetTimeStamp() - dTimeStart;
            iPrimalStatus = HLpSolverGetStatus(lpsolve);
            
            if ( lpsolve->params.LpMethod == LP_ITER_PRIMAL ) {
                isPrimalUsed = 1;
            }
            HLpSolverDestroy(&lpsolve);
            if ( iPrimalStatus == HDSDP_UNKNOWN ) {
                break;
            }
            
        }
        
        if ( iPrimalStatus != HDSDP_PRIMAL_DUAL_OPTIMAL || iPrimalDualStatus != HDSDP_PRIMAL_DUAL_OPTIMAL ) {
            iSuccess = 0;
        }
        
        if ( iPrimalStatus != HDSDP_PRIMAL_DUAL_OPTIMAL && iPrimalDualStatus != HDSDP_PRIMAL_DUAL_OPTIMAL ) {
            break;
        }
        
        if ( !isPrimalUsed ) {
            break;
        }
        
    }
    
    get_sol_status(iPrimalStatus, sPrimalStatus);
    get_sol_status(iPrimalDualStatus, sPrimalDualStatus);
    
    /* Get summary */
    printf("\nTest summary\n");
    if ( iSuccess ) {
        printf("> Test result: Success \n");
    } else {
        printf("> Test result: Failed \n");
    }
    
    /* Print Primal-dual status */
    printf("> Primal dual status: %s \n", sPrimalDualStatus);
    
    if ( isPrimalUsed ) {
        printf("> Primal      status: %s \n", sPrimalStatus);
        /* Print speed ratio */
        if ( iSuccess ) {
            printf("> Primal speedup %f \n", (dPrimalDualTime - dPrimalTime) / dPrimalDualTime);
        }
    } else {
        printf("> Primal      status: %s \n", "Unused");
    }
    
exit_cleanup:
    
    HLpSolverDestroy(&lpsolve);
    
    HDSDP_FREE(Aineqp);
    HDSDP_FREE(Aineqi);
    HDSDP_FREE(Aineqx);
    
    HDSDP_FREE(Aeqp);
    HDSDP_FREE(Aeqi);
    HDSDP_FREE(Aeqx);
    
    HDSDP_FREE(AeqTransp);
    HDSDP_FREE(AeqTransi);
    HDSDP_FREE(AeqTransx);
    HDSDP_FREE(iTransBuffer);
    
    HDSDP_FREE(colUbIdx);
    HDSDP_FREE(colUbElem);
    
    HDSDP_FREE(colObj);
    HDSDP_FREE(rowRhs);
    
    return (int) retcode;
    
}
