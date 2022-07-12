#ifndef dsdpoutput_h
#define dsdpoutput_h

/* Implement the output interface for HDSDP */
#include "dsdpsolver.h"

#ifdef __cplusplus
extern "C" {
#endif

extern DSDP_INT dumpDualSymbolic( HSDSolver *dsdpSolver, char *fname );
extern DSDP_INT dumpDualSol     ( HSDSolver *dsdpSolver, char *fname );

#ifdef __cplusplus
}
#endif


#endif /* dsdpoutput_h */
