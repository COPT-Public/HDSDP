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

#define DSDP64

#ifdef COMPMEX
#include "mex.h"
#endif

#ifdef mex_h
#define calloc mxCalloc
#define printf mexPrintf
#define free mxFree
// Define data type
typedef mwSignedIndex DSDP_INT;
#else
#include <stdint.h>
typedef long int DSDP_INT;
#define ID "%ld"
#endif

typedef struct HSDSolver Solver;

// Memory handler
#define DSDP_FREE(var) do {free((var)); (var) = NULL;} while (0)

// Boolean flag
#define TRUE                    1
#define FALSE                   0

// Solver return code
#define DSDP_RETCODE_OK         2
#define DSDP_RETCODE_FAILED     3

// Algorithm status
#define DSDP_UNKNOWN            99
#define DSDP_OPTIMAL            100
#define DSDP_MAXITER            101
#define DSDP_INTERNAL_ERROR     102
#define DSDP_INACCURATE         103
#define DSDP_PD_INFEASIBLE      104

// Solver status
#define DSDP_STATUS_UNINIT      105
#define DSDP_STATUS_INIT_UNSET  106
#define DSDP_STATUS_SET         107
#define DSDP_STATUS_PRESOLVED   108
#define DSDP_STATUS_SOLVED      109
#define DSDP_STATUS_FAILED      110

// Memory strategy
#define DSDP_MEMORY_THRESHOLD   10000

// Version Number
#define MAJOR_VERSION            0
#define MINOR_VERSION            1

// Error log
#define error(etype, x)          printf("[%s]: %s", (etype), (x)); \
                                 retcode = DSDP_RETCODE_FAILED; return retcode;
#define checkCode                if (retcode != DSDP_RETCODE_OK) return retcode;
#define checkCodeFree            if (retcode != DSDP_RETCODE_OK) goto clean_up;

#endif /* dsdphsd_h */
