#ifndef pot_param_h
#define pot_param_h

#include "pot_def.h"

#define NUM_INT_PARAM (10)
#define NUM_DBL_PARAM (10)

static int defaultIntParam[NUM_INT_PARAM] = {
    10000,
    100
};

static double defaultDblParam[NUM_DBL_PARAM] = {
    1e-04,
    1e-04
};

#define POTLP_INFINITY 1e+30


#endif /* pot_param_h */
