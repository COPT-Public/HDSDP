#ifndef hdsdp_psdp_h
#define hdsdp_psdp_h

#ifdef HEADERPATH
#include "interface/def_hdsdp_psdp.h"
#else
#include "def_hdsdp_psdp.h"
#endif

#ifdef __cplusplus
extern "C" {
#endif
extern hdsdp_retcode HPSDPCreate( hdsdp_psdp **pHpsdp );
extern hdsdp_retcode HPSDPInit( hdsdp_psdp *Hpsdp, hdsdp *HSolver );
extern hdsdp_retcode HPSDPOptimize( hdsdp_psdp *Hpsdp );
extern void HPSDPClear( hdsdp_psdp *Hpsdp );
extern void HPSDPDestroy( hdsdp_psdp **pHpsdp );

#ifdef __cplusplus
}
#endif


#endif /* hdsdp_psdp_h */
