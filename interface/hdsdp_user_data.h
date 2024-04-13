#ifndef hdsdp_user_data_h
#define hdsdp_user_data_h

#ifdef HEADERPATH
#include "interface/hdsdp.h"
#include "interface/def_hdsdp_conic.h"
#else
#include "hdsdp.h"
#include "def_hdsdp_conic.h"
#endif

typedef struct hdsdp_user_data user_data;

#ifdef __cplusplus
extern "C" {
#endif

extern hdsdp_retcode HUserDataCreate( user_data **pHdata );
extern void HUserDataSetConeData( user_data *Hdata, cone_type cone, int nRow, int nCol,
                                  int *coneMatBeg, int *coneMatIdx, double *coneMatElem );
extern cone_type HUserDataChooseCone( user_data *Hdata );
extern void HUserDataClear( user_data *Hdata );
extern void HUserDataDestroy( user_data **pHdata );

#ifdef __cplusplus
}
#endif

#endif /* hdsdp_user_data_h */
