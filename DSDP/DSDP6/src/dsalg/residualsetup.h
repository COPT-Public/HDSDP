#ifndef residualsetup_h
#define residualsetup_h

#include "dsdphsd.h"
#include "dsdpdata.h"
#include "structs.h"
#include "dsdpsolver.h"
#include "dsdplog.h"

/* Residual Ry is setup as proposed by Dr. Huangfu */
#ifdef __cplusplus
extern "C" {
#endif

extern DSDP_INT setupRes( HSDSolver *dsdpSolver );

#ifdef __cplusplus
}
#endif

#endif /* residualsetup_h */
