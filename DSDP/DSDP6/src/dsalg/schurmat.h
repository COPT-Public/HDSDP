#ifndef schurmat_h
#define schurmat_h

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


#endif /* schurmat_h */

