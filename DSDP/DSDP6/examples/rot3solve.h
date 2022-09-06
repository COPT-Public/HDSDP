/** @file rotsolve.h
 * An example for HDSDP solving an SDP from rotation recovery
 *
 */

#ifndef rotsolve_h
#define rotsolve_h

#ifndef ROTSOLVE
#define ROTSOLVE
#endif


#ifndef ROT_35
//#define ROT_35
#endif

#include "dsdphsd.h"
#include "dsdpparam.h"

#ifdef ROT_35
#define MATRIX_DIM      (10)
#define NUM_CONSTR      (35)
#define NUM_LPVAR       (1)
#else
#define MATRIX_DIM      (6)
#define NUM_CONSTR      (15)
#define NUM_LPVAR       (1)
#endif

extern DSDP_INT solveRot( double *b, char *fname );
extern DSDP_INT solveRotfromFile( int argc, char **argv );

#endif /* rotsolve_h */
