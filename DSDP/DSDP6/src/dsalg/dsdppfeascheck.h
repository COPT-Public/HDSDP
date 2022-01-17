#ifndef dsdppfeascheck_h
#define dsdppfeascheck_h

#include "dsdpdata.h"
#include "structs.h"
#include "dsdpsolver.h"
#include "hsd.h"

/* Implement the primal feasibility check */
extern DSDP_INT dsdpCheckPhaseAPfeas( HSDSolver *dsdpSolver, double dtaudelta, vec *dydelta, DSDP_INT *isPfeas );

#endif /* dsdppfeascheck_h */
