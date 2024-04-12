#ifndef def_hdsdp_h
#define def_hdsdp_h

#ifdef HEADERPATH
#include "interface/hdsdp.h"
#include "interface/hdsdp_conic.h"
#include "interface/hdsdp_utils.h"
#include "interface/hdsdp_user_data.h"
#include "interface/hdsdp_schur.h"
#include "interface/hdsdp_algo.h"
#else
#include "hdsdp.h"
#include "hdsdp_conic.h"
#include "hdsdp_utils.h"
#include "hdsdp_user_data.h"
#include "hdsdp_schur.h"
#include "hdsdp_algo.h"
#endif

#define NUM_INT_PARAM  20
#define NUM_DBL_PARAM  20

#define INT_FEATURE_I_NULLOBJ      0
#define INT_FEATURE_I_MANYCONES    1
#define INT_FEATURE_I_NOPINTERIOR  2
#define INT_FEATURE_I_NODINTERIOR  3
#define INT_FEATURE_I_VERYDENSE    4
#define INT_FEATURE_I_IMPTRACE     5
#define INT_FEATURE_I_IMPYBOUND    6
#define INT_FEATURE_N_SUMCONEDIMS  7
#define INT_FEATURE_N_MAXCONEDIM   8
#define INT_FEATURE_N_CONES        9
#define INT_FEATURE_N_ROWS         10
#define INT_FEATURE_N_SPSDPCONES   11
#define INT_FEATURE_N_DSSDPCONES   12
#define INT_FEATURE_N_LPCONES      13
#define INT_FEATURE_N_BNDCONES     14
#define INT_FEATURE_N_ZEORMATS     15
#define INT_FEATURE_N_SPMATS       16
#define INT_FEATURE_N_DSMATS       17
#define INT_FEATURE_N_SPR1MATS     18
#define INT_FEATURE_N_DSR1MATS     19

#define DBL_FEATURE_OBJFRONORM     0
#define DBL_FEATURE_OBJONENORM     1
#define DBL_FEATURE_RHSFRONORM     2
#define DBL_FEATURE_RHSONENORM     3
#define DBL_FEATURE_RHSINFNORM     4
#define DBL_FEATURE_OBJSCALING     5
#define DBL_FEATURE_RHSSCALING     6
#define DBL_FEATURE_DATAFRONORM    7
#define DBL_FEATURE_DATAONENORM    8
#define DBL_FEATURE_IMPYBOUNDUP    9
#define DBL_FEATURE_IMPYBOUNDLOW   10
#define DBL_FEATURE_IMPTRACEX      11

struct hdsdp_solver_internal {
    
    char coneModelName[100];
    char modelFeatures[200];
    
    /* Logging */
    int HPhaseA;
    
    /* User data */
    int nRows;
    double *rowRHS;
    
    /* Cones */
    int nCones;
    hdsdp_cone **HCones;
    /* Bound cone */
    hdsdp_cone *HBndCone;
    
    /* KKT solver */
    hdsdp_kkt *HKKT;
    
    /* Step and vector */
    double *dRowDual;
    double *dRowDualStep;
    double dBarHsdTau;
    double dBarHsdTauStep;
    double *dMinvRowRHS;
    double *dMinvASinv;
    double *dMinvASinvRdSinv;
    double *dMinvASinvCSinv;
    double *dHAuxiVec1;
    double *dHAuxiVec2;
    double dTraceSinv;
    
    /* Monitor */
    int whichMethod;
    int nIterCount;
    int nSmallStep;
    double dBarrierMu;
    double dProxNorm;
    double dPotentialVal;
    double dBarrierVal;
    double dDStep;
    double dResidual;
    double dPerturb;
    double *dPInfeasUpper;
    double *dPInfeasLower;
    
    /* Convergence criterion */
    double dAllConeDims;
    double dPotentialRho;
    double pObjVal;
    double dObjVal;
    double pObjInternal;
    double dObjInternal;
    double dObjImprove;
    double pInfeas;
    double dInfeas;
    double comp;
    
    /* Primal solution recovery */
    double *dAccRowDualMaker;
    double *dAccRowDualStepMaker;
    double dAccBarrierMaker;
    double *dInaccRowDualMaker;
    double *dInaccRowDualStepMaker;
    double dInaccBarrierMaker;
    
    /* Starting time */
    double dTimeBegin;
    
    hdsdp_status HStatus;
    
    /* Parameters */
    int HIntParams[20];
    double HDblParams[20];
    int HIntFeatures[20];
    double HDblFeatures[20];
    double dErrs[6];
};

#endif /* def_hdsdp_h */
