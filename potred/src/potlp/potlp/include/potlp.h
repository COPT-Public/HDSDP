/** @file potlp.h
 *  @brief The header file for the user interface of potential reduction algorithm for linear programming
 *  @author Wenzhi Gao, Shanghai University of Finance and Economics
 *  @date Sep 22th, 2022
 *
 * @TODO: Add more detailed comments
 */
#ifndef potlp_h
#define potlp_h

#include <stdio.h>

#ifdef MATLAB_MEX_FILE
#include "mex.h"
#define calloc mxCalloc
#define printf mexPrintf
#define free mxFree
// Define data type
typedef mwSignedIndex DSDP_INT;
#else
#include <stdint.h>
#ifdef POTLP64
typedef long int potlp_int;
#define ID "%ld"
#else
typedef int potlp_int;
#define ID "%d"
#endif
#endif

// Memory handler
#define POTLP_FREE(var) do {free((var)); (var) = NULL;} while (0)

// Return code
#define POTLP_STATUS_OK     (0)
#define POTLP_STATUS_FAILED (1)

// Algorithm statuss
#define POTLP_UNKNOWN                (99)
#define POTLP_OPTIMAL               (100)
#define POTLP_MAXITER               (101)
#define POTLP_INFEAS_OR_UNBOUNDED   (102)
#define POTLP_TIMELIMIT             (103)


// Version information
#define VERSION_MAJOR           0
#define VERSION_MINOR           1
#define VERSION_TECHNICAL       0

#define BUILD_DATE_YEAR         2022
#define BUILD_DATE_MONTH        9
#define BUILD_DATE_DAY          22



#endif /* potlp_h */
