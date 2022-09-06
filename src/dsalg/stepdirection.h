#ifndef stepdirection_h
#define stepdirection_h

/* Stepsize direction computation in HSD system
 The following step directions will be recovered from the d1 and d2 (solutions to the Schur system)
 
 1. dtau = rM - (b2' * d2) / T
 2. dy   = d1 * dtau + d2
 3. dS   = Rd - ATdy + C * dtau
 4. dkpa = - kappa + mu / tau - (kappa / tau) * dtau

*/

#include "dsdpsolver.h"

#ifdef __cplusplus
extern "C" {
#endif

extern void getStepDirs( HSDSolver *dsdpSolver );

#ifdef __cplusplus
}
#endif

#endif /* stepdirection_h */
