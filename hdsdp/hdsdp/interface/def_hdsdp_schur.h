#ifndef def_hdsdp_schur_h
#define def_hdsdp_schur_h

#ifdef HEADERPATH
#include "interface/hdsdp.h"
#include "interface/hdsdp_utils.h"
#include "interface/hdsdp_conic.h"
#include "linalg/hdsdp_linsolver.h"
#else
#include "hdsdp.h"
#include "hdsdp_utils.h"
#include "hdsdp_conic.h"
#include "hdsdp_linsolver.h"
#endif

/* Define the KKT solver for HDSDP
 
 In HDSDP, the most important part is how to set up the (reduced) KKT system, or the Schur complement
 
 M_{i, j} = Trace( A_i * S^-1 A_j * S^-1 )
 
 or M = A * S^-2 A^*
 
 */

#define KKT_M1 (0)
#define KKT_M2 (1)
#define KKT_M3 (2)
#define KKT_M4 (3)
#define KKT_M5 (4)

typedef struct {
    
    int nRow;
    int nCones;
    int maxConeDim;
    
    hdsdp_cone **cones;
    
    int isKKTSparse;
    hdsdp_linsys_fp *kktM;
    
    /* Buffer for the largest-size dual matrix */
    double *invBuffer;
    double *kktBuffer;
    double *kktBuffer2;
    
    /* KKT data */
    int *kktMatBeg;
    int *kktMatIdx;
    double *kktMatElem;
    double **kktDiag;
    
    /* KKT LHS and RHS components */
    double *dASinvVec;
    double *dASinvCSinvVec;
    double *dASinvRdSinvVec;
    
    double dCSinvCSinv;
    double dCSinvRdSinv;
    double dCSinv;
    double dTraceSinv;
    
    /* Primal IPM. No memory is allocated to inside the pointers */
    /* Primal matrix */
    double **dPrimalX;
    
} hdsdp_kkt;

#endif /* def_hdsdp_schur_h */
