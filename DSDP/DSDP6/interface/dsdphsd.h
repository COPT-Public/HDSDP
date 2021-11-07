#ifndef __DSDPHSD__
#define __DSDPHSD__

/*
    
    A primary implementation of the dual scaling algorithm for semi-definite programming
    using Homogeneous Self-dual Embedding
    
    Yinyu Ye, Stanford University
    Gwz,      Shanghai University of Finance and Economics
 
    Development log:
    
    Nov 5th, 2021.
    Ready to write up the interface for Single-block HSD
 
*/

#ifdef COMPMEX
#include "mex.h"
#endif

#ifdef mex_h
#define calloc mxCalloc
#define printf mexPrintf
#define free mxFree
// Define data type
typedef mwSignedIndex DSDP_INT;
#define KType double *
#else
#include <stdint.h>
typedef int32_t DSDP_INT;
#endif

typedef struct HSDSolver Solver;

// Memory handler
#define DSDPFREE(var) do {free((var)); (var) = NULL;} while (0)

// Boolean flag
#define TRUE                    1
#define FALSE                   0

// Solver return code
#define DSDP_RETCODE_OK         98
#define DSDP_RETCODE_FAILED     99

// Algorithm status
#define DSDP_OPTIMAL            100
#define DSDP_MAXITER            101
#define DSDP_INTERNAL_ERROR     102
#define DSDP_INACCURATE         103
#define DSDP_PD_INFEASIBLE      104

// Solver status
#define DSDP_STATUS_UNINIT      105
#define DSDP_STATUS_UNSOLVED    106
#define DSDP_STATUS_SOLVED      107
#define DSDP_STATUS_FAILED      108

// Version Number
#define MAJOR_VERSION            0
#define MINOR_VERSION            1


#endif /* dsdphsd_h */
