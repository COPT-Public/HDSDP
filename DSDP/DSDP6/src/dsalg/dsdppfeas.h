#ifndef dsdppfeas_h
#define dsdppfeas_h

#include "dsdphsd.h"
#include "dsdpdata.h"
#include "structs.h"
#include "dsdpsolver.h"
#include "dsdplog.h"

#ifdef __cplusplus
extern "C" {
#endif

extern DSDP_INT DSDPPFeasPhase( HSDSolver *dsdpSolver );

#ifdef __cplusplus
}
#endif

#endif /* dsdppfeas_h */
