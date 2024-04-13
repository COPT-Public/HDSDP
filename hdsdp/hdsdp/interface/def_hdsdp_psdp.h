/** @file def\_hdsdp\_psdp
 *  @brief Primal refinement algorithm
 *
 */
#ifndef def_hdsdp_psdp_h
#define def_hdsdp_psdp_h

#ifdef HEADERPATH
#include "interface/hdsdp.h"
#include "interface/hdsdp_utils.h"
#include "interface/hdsdp_conic.h"
#include "interface/hdsdp_schur.h"
#include "linalg/hdsdp_linsolver.h"
#else
#include "hdsdp.h"
#include "hdsdp_utils.h"
#include "hdsdp_conic.h"
#include "hdsdp_schur.h"
#include "hdsdp_linsolver.h"
#endif

typedef struct {
    
    /* Dual solver */
    hdsdp *HSolver;
    
    /* Basic conic information */
    int nRow;
    int nCones;
    
    /* Synchronized from dual solver, no memory needed */
    double *rowRHS;
    hdsdp_cone **HCones;
    hdsdp_kkt *HKKT;
    
    double *dRowDual;
    double *dRowDualStep;
    double *dPrimalAuxiVec1;
    double *dPrimalAuxiVec2;
    double dBarrierMu;
    
    /* Freshly allocated memory */
    double *dPrimalKKTRhs;
    double **dPrimalX;
    double **dPrimalScalX;
    double **dPrimalXStep;
    double **dPrimalMatBuffer;
    
    /* Each cone has its own primal ratio test */
    int iLanczos;
    hdsdp_lanczos **Lanczos;
    /* Each cone has its Primal Cholesky factor */
    hdsdp_linsys_fp **XFactors;
    
} hdsdp_psdp;


#endif /* def_hdsdp_psdp_h */
