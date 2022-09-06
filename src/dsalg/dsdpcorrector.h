#ifndef dsdpcorrector_h
#define dsdpcorrector_h

/* Implement Dual corrector step in DSDP */
#include "dsdpsolver.h"

#ifdef __cplusplus
extern "C" {
#endif

extern DSDP_INT dInfeasCorrectorStep( HSDSolver *dsdpSolver, DSDP_INT isfinal );
extern DSDP_INT dualCorrectorStep( HSDSolver *dsdpSolver );

#ifdef __cplusplus
}
#endif



#endif /* dsdpCorrector_h */
