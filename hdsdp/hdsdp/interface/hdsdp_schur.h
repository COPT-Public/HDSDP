#ifndef hdsdp_schur_h
#define hdsdp_schur_h

#ifdef HEADERPATH
#include "interface/def_hdsdp_schur.h"
#else
#include "def_hdsdp_schur.h"
#endif

extern hdsdp_retcode HKKTCreate( hdsdp_kkt **pHKKT );
extern hdsdp_retcode HKKTInit( hdsdp_kkt *HKKT, int nRow, int nCones, hdsdp_cone **cones );
extern hdsdp_retcode HKKTBuildUp( hdsdp_kkt *HKKT, int typeKKT );
extern hdsdp_retcode HKKTBuildUpExtraCone( hdsdp_kkt *HKKT, hdsdp_cone *cone, int typeKKT );
extern hdsdp_retcode HKKTBuildUpFixed( hdsdp_kkt *HKKT, int typeKKT, int kktStrategy );
extern void HKKTExport( hdsdp_kkt *HKKT, double *dKKTASinvVec, double *dKKTASinvRdSinvVec, double *dKKTASinvCSinvVec,
                       double *dCSinvCSinv, double *dCSinv, double *dCSinvRdCSinv, double *dTraceSinv );
extern hdsdp_retcode HKKTFactorize( hdsdp_kkt *HKKT );
extern hdsdp_retcode HKKTSolve( hdsdp_kkt *HKKT, double *dRhsVec, double *dLhsVec );
extern void HKKTRegularize( hdsdp_kkt *HKKT, double dKKTReg );
extern void HKKTRegisterPSDP( hdsdp_kkt *HKKT, double **dPrimalScalX );
extern void HKKTClear( hdsdp_kkt *HKKT );
extern void HKKTDestroy( hdsdp_kkt **pHKKT );

#endif /* hdsdp_schur_h */
