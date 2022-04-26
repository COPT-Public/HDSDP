#ifndef schurmatops_h
#define schurmatops_h

/* Schuar matrix assembly routine in DSDP
 This part completes the setup of the auxiliary arrays in SDP solution
 as well as the schur matrix
 
 Generally we need to setup
 
 Asinv for LP and ASinv for SDP
 Asinv^2ry for LP and Sinv Ry Sinv for SDP
 
 and then we can setup the Schur matrix by
 
 1. Computing SinvASinv for some A
 2. Compute tr(A, SinvASinv) for different A 
  
 */
#include "dsdphsd.h"
#include "residualsetup.h"
#include "dsdpdata.h"
#include "structs.h"
#include "dsdpsolver.h"
#include "hsd.h"

#ifdef __cplusplus
extern "C" {
#endif

extern DSDP_INT setupFactorize     ( HSDSolver *dsdpSolver );
extern DSDP_INT setupSchur   ( HSDSolver *dsdpSolver );

#ifdef __cplusplus
}
#endif

#endif /* schurmat_h */


