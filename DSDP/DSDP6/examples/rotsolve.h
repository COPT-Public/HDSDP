/** @file rotsolve.h
 * An example for HDSDP solving an SDP from rotation recovery
 *
 *
 *
 */

#ifndef rotsolve_h
#define rotsolve_h

#include "dsdphsd.h"
#include "dsdpparam.h"

#define MATRIX_DIM      (10)
#define NUM_CONSTR      (35)
#define NUM_LPVAR       (1)


extern DSDP_INT solve_rot( double *b );


#endif /* rotsolve_h */
