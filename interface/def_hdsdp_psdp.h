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
    
    hdsdp *HSolver;
    
    int nRow;
    int nCol;
    int nCones;
    
    double *rowRHS;
    hdsdp_cone *HCone;
    hdsdp_kkt *HKKT;
    
    double *dRowDual;
    double *dRowDualStep;
    
    double *dPrimalAuxiVec1;
    double *dPrimalAuxiVec2;
    double *dPrimalKKTRhs;
    
    double dBarrierMu;
    
    double *dPrimalX;
    double *dPrimalScalX;
    double *dDualS;
    double *dPrimalXStep;
    double *dPrimalMatBuffer;
    
    hdsdp_lanczos *Lanczos;
    hdsdp_linsys_fp *XFactor;
    
} hdsdp_psdp;


#endif /* def_hdsdp_psdp_h */
