#ifndef dsdppresolve_h
#define dsdppresolve_h

/* Implement the presolving interface for DSDP */
#include "dsdphsd.h"
#include "dsdpsolver.h"

#ifdef __cplusplus
extern "C" {
#endif

extern DSDP_INT preSDPPrimal          ( HSDSolver *dsdpSolver );
extern DSDP_INT preSDPMatCScale       ( HSDSolver *dsdpSolver );
extern DSDP_INT preSDPDual            ( HSDSolver *dsdpSolver );
extern DSDP_INT preRank1Rdc           ( HSDSolver *dsdpSolver );
extern DSDP_INT preRankkRdc           ( HSDSolver *dsdpSolver );
extern DSDP_INT getMatIdx             ( HSDSolver *dsdpSolver );
extern DSDP_INT preSymbolic           ( HSDSolver *dsdpSolver );
extern DSDP_INT preStructureDetect    ( HSDSolver *dsdpSolver );
extern DSDP_INT DSDPPrepareMAssembler ( HSDSolver *dsdpSolver );
#ifdef __cplusplus
}
#endif

#endif /* dsdppresolve_h */
