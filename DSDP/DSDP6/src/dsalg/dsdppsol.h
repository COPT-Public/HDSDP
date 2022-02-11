#ifndef dsdppsol_h
#define dsdppsol_h
/* Implement the primal solution extractor of DSDP */
#include "dsdphsd.h"
#include "dsdpdata.h"
#include "structs.h"
#include "dsdpsolver.h"
#include "hsd.h"

#ifdef __cplusplus
extern "C" {
#endif

extern DSDP_INT computePrimalX( HSDSolver *dsdpSolver );
extern DSDP_INT computeDIMACS( HSDSolver *dsdpSolver,
                               double *err1, double *err2,
                               double *err3, double *err4,
                               double *err5, double *err6 );

#ifdef __cplusplus
}
#endif

#endif /* dsdppsol_h */
