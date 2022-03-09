#ifndef dsdpparam_h
#define dsdpparam_h

/* Implement the parameter interface for DSDP */
#include "dsdphsd.h"

#define DBL_PARAM_ASIGMA         0
#define DBL_PARAM_BSIGMA         1
#define DBL_PARAM_RHO            2
#define DBL_PARAM_INIT_POBJ      3
#define DBL_PARAM_INIT_BETA      4
#define DBL_PARAM_INIT_MU        5
#define DBL_PARAM_INIT_TAU       6
#define DBL_PARAM_INIT_KAPPA     7
#define DBL_PARAM_AALPHA         8
#define DBL_PARAM_NRM_THRESH     9
#define DBL_PARAM_INFEAS_THRESH 10
#define DBL_PARAM_ABS_OPTTOL    11
#define DBL_PARAM_REL_OPTTOL    12
#define DBL_PARAM_ABS_FEASTOL   13
#define DBL_PARAM_REL_FEASTOL   14


#define INT_PARAM_ACORRECTOR     0
#define INT_PARAM_BCORRECTOR     1
#define INT_PARAM_INITMETHOD     2
#define INIT_METHOD_FRO         20
#define INT_PARAM_AATTEMPT       3
#define AATEMPT_AGGRESSIVE      30
#define AATEMPT_MILD            31
#define AATEMT_CONSERVATIVE     32
#define INT_PARAM_CG_REUSE       4
#define INT_PARAM_PRESOLVE       5
#define PRESOLVE_AGGRESSIVE     50
#define PRESOLVE_CONSERVATIVE   51
#define INT_PARAM_AMAXITER       6
#define INT_PARAM_BMAXITER       7

#define NUM_DBL_PARAM  (DBL_PARAM_REL_FEASTOL + 1)
#define NUM_INT_PARAM  (INT_PARAM_BMAXITER   + 1)

typedef struct {
    
    DSDP_INT *intParams;
    double   *dblParams;
    
} hsdParam;

static DSDP_INT defaultIntParam[NUM_INT_PARAM] = {
    0,     // INT_PARAM_ACORRECTOR
    0,     // INT_PARAM_BCORRECTOR
    INIT_METHOD_FRO,
           // INT_PARAM_INITMETHOD
    AATEMPT_AGGRESSIVE,
           // INT_PARAM_AATTEMPT
    10,     // INT_PARAM_CG_REUSE
    PRESOLVE_AGGRESSIVE,
           // INT_PARAM_PRESOLVE
    100,   // INT_PARAM_AMAXITER
    500    // INT_PARAM_BMAXITER
};

static double defaultDblParam[NUM_DBL_PARAM] = {
    
    0.9,   // DBL_PARAM_ASIGMA
    0.1,   // DBL_PARAM_BSIGMA
    4.0,   // DBL_PARAM_RHO
    1e+05, // DBL_PARAM_INIT_POBJ
    1.0,   // DBL_PARAM_INIT_BETA
    1e+06, // DBL_PARAM_INIT_MU
    1.0,   // DBL_PARAM_INIT_TAU
    1.0,   // DBL_PARAM_INIT_KAPPA
    0.75,  // DBL_PARAM_AALPHA
    1e+08, // DBL_PARAM_NRM_THRESH
    1e+08, // DBL_PARAM_INFEAS_THRESH
    1e-06, // DBL_PARAM_ABS_OPTTOL
    1e-03, // DBL_PARAM_REL_OPTTOL
    1e-06, // DBL_PARAM_ABS_FEASTOL
    1e-03  // DBL_PARAM_REL_FEASTOL
};

static hsdParam defaultParam =
{
    defaultIntParam,
    defaultDblParam
};

extern DSDP_INT setDblParam ( hsdParam *param, DSDP_INT pName, double    dblVal );
extern DSDP_INT getDblParam ( hsdParam *param, DSDP_INT pName, double   *dblVal );
extern DSDP_INT setIntParam ( hsdParam *param, DSDP_INT pName, DSDP_INT  intVal );
extern DSDP_INT getIntParam ( hsdParam *param, DSDP_INT pName, DSDP_INT *intVal );
extern void     printParams ( hsdParam *param                                   );
extern void     DSDPParamPrint ( hsdParam *param );

#endif /* dsdpparam_h */
