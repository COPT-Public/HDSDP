#ifndef dsdpproxmeasure_h
#define dsdpproxmeasure_h
/* Get the measure of proximity in DSDP */
#include "dsdphsd.h"
#include "residualsetup.h"
#include "dsdpdata.h"
#include "structs.h"
#include "dsdpsolver.h"
#include "hsd.h"


extern DSDP_INT dsdpgetPhaseAProxMeasure( HSDSolver *dsdpSolver, double newmu );



#endif /* dsdpproxmeasure_h */
