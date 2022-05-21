#ifndef __DSDPHSD__
#define __DSDPHSD__

/*
    
    A primary implementation of the dual scaling algorithm for semi-definite programming
    using Homogeneous Self-dual Embedding
    
    Yinyu Ye, Stanford University
    Gwz,      Shanghai University of Finance and Economics
 
*/
#ifdef superDebug
#undef superDebug
#endif

// #define superDebug

//#define compareMode

#ifdef RELEASE
#define assert(x) if (!(x)) { printf("Fatal error. \n"); exit(0); };
#else
#include <assert.h>
#endif

//#define superDebug
#ifndef MAXTIME
#define MAXTIME 7200
#endif

#ifdef SHOWALL
#undef SHOWALL
#endif

// #define DSDP64

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
#ifdef DSDP64
typedef long int DSDP_INT;
#define ID "%ld"
#else
typedef int DSDP_INT;
#define ID "%d"
#endif
#endif

//typedef HSDSolver Solver;

// Memory handler
#define DSDP_FREE(var) do {free((var)); (var) = NULL;} while (0)

// Boolean flag
#define TRUE                    1
#define FALSE                   0

// Constants
#define DSDP_INFINITY           1e+20

// Solver return code
#define DSDP_RETCODE_OK         2
#define DSDP_RETCODE_FAILED     3

// Algorithm status
#define DSDP_UNKNOWN            99
#define DSDP_OPTIMAL            100
#define DSDP_MAXITER            101
#define DSDP_INTERNAL_ERROR     102
#define DSDP_PD_FEASIBLE        103
#define DSDP_PINFEAS_DFEAS      104
#define DSDP_PUNKNOWN_DFEAS     105
#define DSDP_PFEAS_DINFEAS      106
#define DSDP_PUNKNOWN_DINFEAS   107

// Solver status
#define DSDP_STATUS_UNINIT      0
#define DSDP_STATUS_INIT_UNSET  1
#define DSDP_STATUS_SET         2
#define DSDP_STATUS_PRESOLVED   3
#define DSDP_STATUS_SOLVED      4
#define DSDP_STATUS_FAILED      5

// Memory strategy
#define DSDP_MEMORY_THRESHOLD   10000

// Error log
#define error(etype, x)          printf("[%s]: %s", (etype), (x)); assert( 0 ); \
                                 retcode = DSDP_RETCODE_FAILED; return retcode;
#define error_clean(etype, x)    printf("[%s]: %s", (etype), (x)); \
                                 retcode = DSDP_RETCODE_FAILED; \
                                 goto exit_cleanup;
#define checkCode                if (retcode != DSDP_RETCODE_OK) {error(etype, "\n") return retcode;}
#define checkCodeFree            if (retcode != DSDP_RETCODE_OK) {error(etype, "\n") goto clean_up;}

// Parameters

// Initialization
#define DSDP_INITMETHOD_FRO     105

// Phase A attempt
#define DSDP_ATTEMPT_NO         106
#define DSDP_ATTEMPT_CONSV      107
#define DSDP_ATTEMPT_AGG        108

#define VERSION_MAJOR           0
#define VERSION_MINOR           5
#define VERSION_TECHNICAL       0

#endif /* dsdphsd_h */
