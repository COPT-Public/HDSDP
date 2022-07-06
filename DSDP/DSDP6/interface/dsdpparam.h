#ifndef dsdpparam_h
#define dsdpparam_h

/* Implement the parameter interface for DSDP */
#include "dsdphsd.h"

#define DBL_PARAM_RHO            0
#define DBL_PARAM_RHON           1
#define DBL_PARAM_INIT_POBJ      2
#define DBL_PARAM_INIT_BETA      3
#define DBL_PARAM_INIT_MU        4
#define DBL_PARAM_INFEAS_THRESH  5
#define DBL_PARAM_ABS_OPTTOL     6
#define DBL_PARAM_REL_OPTTOL     7
#define DBL_PARAM_ABS_FEASTOL    8
#define DBL_PARAM_REL_FEASTOL    9
#define DBL_PARAM_PRLX_PENTALTY 10
#define DBL_PARAM_BOUND_X       11
#define DBL_PARAM_TIMELIMIT     12

#define INT_PARAM_ACORRECTOR     0
#define INT_PARAM_BCORRECTOR     1
#define INT_PARAM_CG_REUSE       2
#define INT_PARAM_AMAXITER       3
#define INT_PARAM_BMAXITER       4
#define INT_PARAM_GOLDSEARCH     5

#define NUM_DBL_PARAM  (DBL_PARAM_TIMELIMIT  + 1)
#define NUM_INT_PARAM  (INT_PARAM_GOLDSEARCH + 1)

typedef struct {
    DSDP_INT *intParams;
    double *dblParams;
} hsdParam;

static DSDP_INT defaultIntParam[NUM_INT_PARAM] = {
    4,     // INT_PARAM_ACORRECTOR
    0,     // INT_PARAM_BCORRECTOR
    100,   // INT_PARAM_CG_REUSE
    500,   // INT_PARAM_AMAXITER
    500,   // INT_PARAM_BMAXITER
    FALSE, // INT_PARAM_GOLDSEARCH
};

static double defaultDblParam[NUM_DBL_PARAM] = {
    4.0,   // DBL_PARAM_RHO
    4.0,   // DBL_PARAM_RHON
    1e+10, // DBL_PARAM_INIT_POBJ
    1e+05, // DBL_PARAM_INIT_BETA
    1e+15, // DBL_PARAM_INIT_MU
    1e+08, // DBL_PARAM_INFEAS_THRESH
    1e-06, // DBL_PARAM_ABS_OPTTOL
    1e-06, // DBL_PARAM_REL_OPTTOL
    1e-07, // DBL_PARAM_ABS_FEASTOL
    1e-07, // DBL_PARAM_REL_FEASTOL
    1e+07, // DBL_PARAM_PRLX_PENTALTY
    1e+08, // DBL_PARAM_BOUND_X
    15000, // DBL_PARAM_TIME_LIMIT
};

static hsdParam defaultParam =
{
    defaultIntParam,
    defaultDblParam
};

extern void setDblParam ( hsdParam *param, DSDP_INT pName, double    dblVal );
extern void getDblParam ( hsdParam *param, DSDP_INT pName, double   *dblVal );
extern void setIntParam ( hsdParam *param, DSDP_INT pName, DSDP_INT  intVal );
extern void getIntParam ( hsdParam *param, DSDP_INT pName, DSDP_INT *intVal );
extern void printParams ( hsdParam *param                                   );
extern void DSDPParamPrint ( hsdParam *param );

#endif /* dsdpparam_h */
