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

#define HDSDP_INTFEATURE_N_SUMCONEDIMS  0
#define HDSDP_INTFEATURE_N_CONES        1
#define HDSDP_INTFEATURE_N_ROWS         2
#define HDSDP_INTFEATURE_I_NULLOBJ      3
#define HDSDP_INTFEATURE_I_MANYBLOCKS   4
#define HDSDP_INTFEATURE_I_NOPINTERIOR  5
#define HDSDP_INTFEATURE_I_NODINTERIOR  6
#define HDSDP_INTFEATURE_I_VERYDENSE    7
#define HDSDP_INTFEATURE_I_IMPTRACE     8
#define HDSDP_INTFEATURE_I_IMPDBOUNDD   9

#define HDSDP_DBLFEATURE_OBJFRONORM     0
#define HDSDP_DBLFEATURE_OBJONENORM     1
#define HDSDP_DBLFEATURE_RHSFRONORM     2
#define HDSDP_DBLFEATURE_RHSONENORM     3

struct hdsdp_solver_internal {
    
    char *coneModelName[100];
    
    /* Logging */
    int HPhaseA;
    
    /* User data */
    int nRows;
    double *rowRHS;
    
    /* Cones */
    int nCones;
    hdsdp_cone **HCones;
    
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
    
    /* Monitor */
    int nIterCount;
    double dBarrierMu;
    double dProxNorm;
    double dPotentialVal;
    double dBarrierVal;
    double dPStep;
    double dDStep;
    double dResidual;
    
    /* Convergence criterion */
    double dPotentialRho;
    double pObjVal;
    double dObjVal;
    double pInfeas;
    double dInfeas;
    double comp;
    double pInfeasRel;
    double dInfeasRel;
    double compRel;
    
    double dTimeBegin;
    
    hdsdp_status HStatus;
    
    /* Parameters */
    int HIntParams[20];
    double HDblParams[20];
    int HIntFeatures[10];
    double HDblFeatures[10];
};

#endif /* def_hdsdp_h */
