#ifndef rankonemat_h
#define rankonemat_h

/* Implement the rank one matrix */
#include <stdio.h>
#include "dsdphsd.h"
#include "vec.h"


#ifdef __cplusplus
extern "C" {
#endif

extern DSDP_INT r1MatInit     ( r1Mat *x                   );
extern DSDP_INT r1MatAlloc    ( r1Mat *x, const DSDP_INT n );
extern DSDP_INT r1MatCountNnz ( r1Mat *x                   );
extern DSDP_INT r1MatFree     ( r1Mat *x                   );
extern DSDP_INT r1MatFnorm    ( r1Mat *x, double *fnrm     );
extern DSDP_INT r1MatRscale   ( r1Mat *x, double r         );

#ifdef __cplusplus
}
#endif

#endif /* rankonemat_h */

