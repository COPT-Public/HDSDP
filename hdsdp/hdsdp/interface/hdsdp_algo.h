#ifndef hdsdp_algo_h
#define hdsdp_algo_h

#ifdef HEADERPATH
#include "interface/hdsdp.h"
#else
#include "hdsdp.h"
#endif

#ifdef __cplusplus
extern "C" {
#endif

extern hdsdp_retcode HDSDPSolve( hdsdp *HSolver, int dOnly );

#ifdef __cplusplus
}
#endif

#endif /* hdsdp_algo_h */
