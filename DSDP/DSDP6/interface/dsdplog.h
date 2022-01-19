#ifndef dsdplog_h
#define dsdplog_h
#include "dsdphsd.h"
#include "dsdpsolver.h"
/* Implement problem statistic collection and logging */

extern void dsdpshowdash          ( void );
extern void dsdpprintPhaseAheader ( void );

/* Phase A operations */
extern DSDP_INT DSDPCheckPhaseAConvergence ( HSDSolver *dsdpSolver, DSDP_INT *isOK );
extern DSDP_INT DSDPPhaseALogging          ( HSDSolver *dsdpSolver );
extern DSDP_INT printPhaseASummary         ( HSDSolver *dsdpSolver, double time );

#endif /* dsdplog_h */
