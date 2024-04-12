#ifndef hdsdp_neqnsolver_h
#define hdsdp_neqnsolver_h

#ifdef HEADERPATH
#include "linalg/def_hdsdp_neqnsolver.h"
#else
#include "def_hdsdp_neqnsolver.h"
#endif

#ifdef __cplusplus
extern "C" {
#endif

extern hdsdp_retcode HNEquationCreate( hdsdp_normal_linsys **pHNeq );


#ifdef __cplusplus
}
#endif

#endif /* hdsdp_neqnsolver_h */
