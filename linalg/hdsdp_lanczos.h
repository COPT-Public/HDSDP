#ifndef hdsdp_lanczos_h
#define hdsdp_lanczos_h

#ifdef HEADERPATH
#include "linalg/def_hdsdp_lanczos.h"
#else
#include "def_hdsdp_lanczos.h"
#endif

#ifdef __cplusplus
extern "C" {
#endif

extern hdsdp_retcode HLanczosCreate( hdsdp_lanczos **pHLanczos );
extern hdsdp_retcode HLanczosInit( hdsdp_lanczos *HLanczos, int nCol, int nSpaceDim );
extern void HLanczosSetData( hdsdp_lanczos *HLanczos, void *MMat, void (*Mvec) (void *, double *, double *) );
extern hdsdp_retcode HLanczosSolve( hdsdp_lanczos *HLanczos, double *LanczosStart, double *dMaxStep );
extern void HLanczosClear( hdsdp_lanczos *HLanczos );
extern void HLanczosDestroy( hdsdp_lanczos **pHLanczos );

#ifdef __cplusplus
}
#endif

#endif /* hdsdp_lanczos_h */
