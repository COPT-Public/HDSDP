#ifndef heurpool_h
#define heurpool_h

#include "dsdpsolver.h"

#ifndef DSDP_HEURS
#define DSDP_HEURS(x) DSDP_HEURISTIC_##x
#endif


extern void DSDP_HEURS( adjustSolverParams ) ( HSDSolver *dsdpSolver, double largeblock );
extern void DSDP_HEURS( adjustCScaler      ) ( HSDSolver *dsdpSolver );
extern void DSDP_HEURS( adjPRelaxPenalty   ) ( HSDSolver *dsdpSolver );
extern void DSDP_HEURS( adjDualPerturb     ) ( HSDSolver *dsdpSolver );

#endif /* heurpool_h */
