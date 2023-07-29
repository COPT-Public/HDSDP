/** @file hdsdp\_utils.h
 *  @brief Define utilities for HDSDP
 *  @date 11/24/2022
 */

#ifndef hdsdp_utils_h
#define hdsdp_utils_h

#ifdef HEADERPATH
#include "interface/hdsdp.h"
#else
#include "hdsdp.h"
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

/* Define macros */
#ifdef SILENT_SOLVER
#define hdsdp_printf
#else
#define hdsdp_printf printf
#endif

#define HDSDP_FREE(var) do {if (var) {free((var)); (var) = NULL;}} while (0)
#define HDSDP_INIT(var, type, size) (var) = (type *) calloc(size, sizeof(type))
#define HDSDP_REALLOC(var, type, size) (var) = (type *) realloc(var, sizeof(type) * (size))
#define HDSDP_MEMCPY(dst, src, type, size) memcpy(dst, src, sizeof(type) * (size))
#define HDSDP_ZERO(var, type, size) memset(var, 0, sizeof(type) * (size))
#define HDSDP_NULLCHECK(var) if(!(var)) {                     \
                                retcode = HDSDP_RETCODE_FAILED; \
                                goto exit_cleanup;              \
                             }
#define HDSDP_MEMCHECK(var) if (!(var)) {                       \
                                retcode = HDSDP_RETCODE_MEMORY; \
                                goto exit_cleanup;              \
                            }
#define HDSDP_ERROR_TRACE printf("File [%30s] Line [%d]\n", __FILE__, __LINE__)
#define HDSDP_CALL(func) retcode = (func);                      \
                         if (retcode != HDSDP_RETCODE_OK) {     \
                             goto exit_cleanup;                 \
                         }
#define HDSDP_MAX(x, y) ((x) > (y) ? (x) : (y))
#define HDSDP_MIN(x, y) ((x) < (y) ? (x) : (y))

#define PACK_NNZ(n) ((n) * ((n) + 1) / 2)
#define PACK_IDX(n, i, j) (int)((2 * (n) - (j) - 1) * (j) / 2) + (i)
#define FULL_IDX(n, i, j) ((j) * (n) + (i))
#define PACK_ENTRY(A, n, i, j) (A[(int)((2 * (n) - (j) - 1) * (j) / 2) + (i)])
#define FULL_ENTRY(A, n, i, j) (A[(j) * (n) + (i)])

#define HDSDP_PROFILER(func, ntest)                 \
double HLocalTProfiler = HUtilGetTimeStamp();       \
for ( int iTest = 0; iTest < (ntest); ++iTest ) {   \
    (func);                                         \
    printf("Run %d. Elapsed time: %f\n",                \
    iTest + 1, HUtilGetTimeStamp() - HLocalTProfiler);  \
}                                                   \
printf("Function Profiler: Line %d of %s by %d runs. "       \
"Average running time: %fs\n", __LINE__, __FILE__,  \
(ntest), (HUtilGetTimeStamp() - HLocalTProfiler) / ntest);

#define HDSDP_CODE_PROFILER_START double tHDSDPStart = HUtilGetTimeStamp()
#define HDSDP_CODE_PROFILER_END            \
    printf("Code Profiler Line %d of %s. "           \
    "Running time: %fs\n", __LINE__, __FILE__, \
    (HUtilGetTimeStamp() - tHDSDPStart))

#define set_func_pointer(A, B) (A = (typeof(A)) B)
#define set_int_param(hdsdp, param, val)  hdsdp->HIntParams[param] = val
#define set_dbl_param(hdsdp, param, val)  hdsdp->HDblParams[param] = val
#define get_int_param(hdsdp, param) hdsdp->HIntParams[param]
#define get_dbl_param(hdsdp, param) hdsdp->HDblParams[param]
#define set_int_feature(hdsdp, feature, val) hdsdp->HIntFeatures[feature] = val
#define set_dbl_feature(hdsdp, feature, val) hdsdp->HDblFeatures[feature] = val
#define get_int_feature(hdsdp, feature) hdsdp->HIntFeatures[feature]
#define get_dbl_feature(hdsdp, feature) hdsdp->HDblFeatures[feature]
#define print_int_param(hdsdp, param, name) printf("  %-20s: %d\n", name, hdsdp->HIntParams[param])
#define print_dbl_param(hdsdp, param, name) printf("  %-20s: %5.1e\n", name, hdsdp->HDblParams[param])
#define print_int_feature(hdsdp, feature, name) printf("  %10s: %d\n", name, hdsdp->HIntFeatures[feature])
#define print_dbl_feature(hdsdp, feature, name) printf("  %10s: %5.2e\n", name, hdsdp->HDblFeatures[feature])

#ifdef __cplusplus
extern "C" {
#endif

extern double HUtilGetTimeStamp( void );
extern void HUtilMatSymmetrize( int n, double *v );
extern int HUtilCheckIfAscending( int n, int *idx );
extern void HUtilDescendSortIntByInt( int *data, int *ref, int low, int up );
extern void HUtilSortIntbyDbl( int *data, double *ref, int low, int up );
extern void HUtilAscendSortDblByInt( double *data, int *ref, int low, int up );
extern void HUtilAscendSortDblByInt( double *data, int *ref, int low, int up );
extern void HUtilPrintDblContent( int n, double *d );
extern void HUtilPrintIntContent( int n, int *d );
extern double HUtilPrintDblSum( int n, double *d );
extern double HUtilPrintDblAbsSum( int n, double *d );
extern int HUtilCheckUserInterrupt( void );

extern void HUtilStartCtrlCCheck( void );
extern int HUtilCheckCtrlC( void );
extern void HUtilResetCtrl( void );

hdsdp_retcode HUtilKKTCheck( void *kkt );

#ifdef __cplusplus
}
#endif

#endif /* hdsdp_utils_h */
