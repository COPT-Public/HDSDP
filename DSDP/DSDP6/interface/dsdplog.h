#ifndef dsdplog_h
#define dsdplog_h
#include "dsdphsd.h"
#include "dsdpsolver.h"
/* Implement problem statistic collection and logging */

#ifdef __cplusplus
extern "C" {
#endif

extern void dsdpshowdash          ( void );
extern void dsdpprintPhaseAheader ( void );
extern void dsdpprintPhaseBheader ( void );
extern void dsdpCheckNan          ( HSDSolver *dsdpSolver );

/* Phase A operations */
extern void     DSDPResetPhaseAMonitor     ( HSDSolver *dsdpSolver );
extern DSDP_INT DSDPCheckPhaseAConvergence ( HSDSolver *dsdpSolver, DSDP_INT *isOK );
extern DSDP_INT DSDPPhaseALogging          ( HSDSolver *dsdpSolver );
extern DSDP_INT printPhaseASummary         ( HSDSolver *dsdpSolver, double time );
extern DSDP_INT printPhaseABConvert        ( HSDSolver *dsdpSolver, DSDP_INT *goPb );

/* Phase B operations */
extern void     DSDPResetPhaseBMonitor     ( HSDSolver *dsdpSolver );
extern DSDP_INT DSDPCheckPhaseBConvergence ( HSDSolver *dsdpSolver, DSDP_INT *isOK );
extern DSDP_INT DSDPPhaseBLogging          ( HSDSolver *dsdpSolver );

#ifdef __cplusplus
}
#endif

#endif /* dsdplog_h */
