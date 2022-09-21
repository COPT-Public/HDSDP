#ifndef dsdpparam_h
#define dsdpparam_h

/* Implement the parameter interface for DSDP */
#include "dsdphsd.h"

#define NUM_DBL_PARAM  (DBL_PARAM_TIMELIMIT  + 1)
#define NUM_INT_PARAM  (INT_PARAM_GOLDSEARCH + 1)

typedef struct {
    DSDP_INT *intParams;
    double *dblParams;
} dsdpparam;

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

static dsdpparam defaultParam =
{
    defaultIntParam,
    defaultDblParam
};

extern void setDblParam    ( dsdpparam *param, DSDP_INT pName, double    dblVal );
extern void getDblParam    ( dsdpparam *param, DSDP_INT pName, double   *dblVal );
extern void setIntParam    ( dsdpparam *param, DSDP_INT pName, DSDP_INT  intVal );
extern void getIntParam    ( dsdpparam *param, DSDP_INT pName, DSDP_INT *intVal );
extern void printParams    ( dsdpparam *param                                   );
extern void DSDPParamPrint ( dsdpparam *param );

#endif /* dsdpparam_h */
