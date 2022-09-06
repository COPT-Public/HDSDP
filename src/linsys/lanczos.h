/** @file lanczos.h
 *  @brief Headers and roitine list for Lanczos iteration
 *
 * Implement the Lanczos routine to compute the maximum \f$ \alpha \f$ such that
 * \f$ S + \alpha \Delta S \f$ is still in the cone.
 *
 *  @author Wenzhi Gao, Shanghai University of Finance and Economics
 *  @date Aug, 28th, 2022
 *
 */

#ifndef lanczos_h
#define lanczos_h

#include "structs.h"

#ifdef __cplusplus
extern "C" {
#endif

extern void     lczInit ( lczstep *lczSolver );
extern DSDP_INT lczAlloc ( lczstep *lczSolver, DSDP_INT n );
extern void     lczStepLength ( lczstep *lczSolver, spsMat *S, spsMat *dS, double *lbd, double *delta );
extern void     lczFree ( lczstep *lczSolver );

#ifdef __cplusplus
}
#endif

#endif /* lanczos_h */
