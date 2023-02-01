#ifndef hdsdp_user_data_h
#define hdsdp_user_data_h

#include "interface/hdsdp.h"
#include "interface/def_hdsdp_conic.h"

typedef struct hdsdp_user_data user_data;

#ifdef __cplusplus
extern "C" {
#endif

extern hdsdp_retcode HUserDataCreate( user_data **pHdata );
extern void HUserDataSetConeInfo( user_data *Hdata, cone_type cone, int nRow, int nCol );
extern void HUserDataSetConeData( user_data *Hdata, int *coneMatBeg, int *coneMatIdx, double *coneMatElem );
extern cone_type HUserDataChooseCone( user_data *Hdata );
extern void HUserDataClear( user_data *Hdata );
extern void HUserDataDestroy( user_data **pHdata );

#ifdef __cplusplus
}
#endif

#endif /* hdsdp_user_data_h */
