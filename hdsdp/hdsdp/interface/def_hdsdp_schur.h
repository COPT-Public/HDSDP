#ifndef def_hdsdp_schur_h
#define def_hdsdp_schur_h

#ifdef HEADERPATH
#include "interface/hdsdp.h"
#include "interface/hdsdp_utils.h"
#include "linalg/hdsdp_linsolver.h"
#else
#include "hdsdp.h"
#include "hdsdp_utils.h"
#include "hdsdp_linsolver.h"
#endif

/* Define the KKT solver for HDSDP
 
 In HDSDP, the most important part is how to set up the (reduced) KKT system, or the Schur complement
 
 M_{i, j} = Trace( A_i * S^-1 A_j * S^-1 )
 
 or M = A * S^-2 A^*
 
 */

typedef struct {
    
    int nRow;
    int nCones;
    
    
    
    
} hdsdp_kkt;

#endif /* def_hdsdp_schur_h */
