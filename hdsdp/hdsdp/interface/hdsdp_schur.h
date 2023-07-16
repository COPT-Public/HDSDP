#ifndef hdsdp_schur_h
#define hdsdp_schur_h

#ifdef HEADERPATH
#include "interface/def_hdsdp_schur.h"
#else
#include "def_hdsdp_schur.h"
#endif

extern hdsdp_retcode HKKTCreate( hdsdp_kkt **pHKKT );
extern hdsdp_retcode HKKTInit( hdsdp_kkt *HKKT, int nRow, int nCones, hdsdp_cone **cones );
extern hdsdp_retcode HKKTBuildUp( hdsdp_kkt *HKKT, int newM );
extern void HKKTRegularize( hdsdp_kkt *HKKT, double dKKTReg );
extern void HKKTClear( hdsdp_kkt *HKKT );
extern void HKKTDestroy( hdsdp_kkt **pHKKT );

#endif /* hdsdp_schur_h */
