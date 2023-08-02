#ifndef hdsdp_linsolver_h
#define hdsdp_linsolver_h

#ifdef HEADERPATH
#include "linalg/def_hdsdp_linsolver.h"
#else
#include "def_hdsdp_linsolver.h"
#endif

typedef hdsdp_linsys_fp hdsdp_linsys;

#ifdef __cplusplus
extern "C" {
#endif

extern hdsdp_retcode HFpLinsysCreate( hdsdp_linsys_fp **pHLin, int nCol, linsys_type Ltype );
extern void HFpLinsysSetParam( hdsdp_linsys_fp *HLin, double relTol, double absTol, int nThreads, int maxIter, int nRestartFreq );
extern hdsdp_retcode HFpLinsysSymbolic( hdsdp_linsys_fp *HLin, int *colMatBeg, int *colMatIdx );
extern hdsdp_retcode HFpLinsysNumeric( hdsdp_linsys_fp *HLin, int *colMatBeg, int *colMatIdx, double *colMatElem );
extern hdsdp_retcode HFpLinsysPsdCheck( hdsdp_linsys_fp *HLin, int *colMatBeg, int *colMatIdx, double *colMatElem, int *isPsd );
extern void HFpLinsysFSolve( hdsdp_linsys_fp *HLin, int nRhs, double *rhsVec, double *solVec );
extern void HFpLinsysBSolve( hdsdp_linsys_fp *HLin, int nRhs, double *rhsVec, double *solVec );
extern hdsdp_retcode HFpLinsysSolve( hdsdp_linsys_fp *HLin, int nRhs, double *rhsVec, double *solVec );
extern hdsdp_retcode HFpLinsysGetDiag( hdsdp_linsys_fp *HLin, double *diagElem );
extern void HFpLinsysInvert( hdsdp_linsys_fp *HLin, double *dFullMatrix, double *dAuxiMatrix );
extern void HFpLinsysClear( hdsdp_linsys_fp *HLin );
extern void HFpLinsysDestroy( hdsdp_linsys_fp **HLin );

#ifdef __cplusplus
}
#endif

#endif /* hdsdp_linsolver_h */
