/** @file hdsdp\_utils.c
 *  @brief Implement the utilities for HDSDP
 *  @date 11/24/2022
 */
#ifdef HEADERPATH
#include "interface/hdsdp_utils.h"
#else
#include "hdsdp_utils.h"
#endif

#include <math.h>
#include <sys/time.h>
#include <signal.h>

static int isCtrlC = 0;
static void monitorCtrlC( int sigNum ) {
    isCtrlC = 1;
    return;
}

static struct sigaction act;

/* TODO: Add compatibility for Windows platform */
static double my_clock( void ) {
    
    struct timeval t;
    gettimeofday(&t, NULL);
    return (1e-06 * t.tv_usec + t.tv_sec);
}

static int dpartitioni( int *ind, double *val, int l, int h ) {
    
    double tmp2 = 0.0, p = val[l];
    int tmp = l, tmp3 = 0;
    
    while ( l < h ) {
        
        while ( l < h && val[h] >= p ) { --h; }
        while ( l < h && val[l] <= p ) { ++l; }
        
        if ( l < h ) {
            tmp2 = val[l]; val[l] = val[h]; val[h] = tmp2;
            tmp3 = ind[l]; ind[l] = ind[h]; ind[h] = tmp3;
        }
    }
    
    tmp2 = val[l]; val[l] = val[tmp]; val[tmp] = tmp2;
    tmp3 = ind[l]; ind[l] = ind[tmp]; ind[tmp] = tmp3;
    
    return l;
}

static int ipartitiond( double *ind, int *val, int l, int h ) {
    
    int tmp2 = 0, p = val[l], tmp = l;
    double tmp3;
    
    while ( l < h ) {
        
        while ( l < h && val[h] >= p ) { --h; }
        while ( l < h && val[l] <= p ) { ++l; }
        
        if ( l < h ) {
            tmp2 = val[l]; val[l] = val[h]; val[h] = tmp2;
            tmp3 = ind[l]; ind[l] = ind[h]; ind[h] = tmp3;
        }
    }
    
    tmp2 = val[l]; val[l] = val[tmp]; val[tmp] = tmp2;
    tmp3 = ind[l]; ind[l] = ind[tmp]; ind[tmp] = tmp3;
    
    return l;
}

static int ipartitioni( int *ind, int *val, int l, int h ) {
    
    int tmp = l, tmp2 = 0, tmp3, p = val[l];
    
    while ( l < h ) {
        
        while ( l < h && val[h] <= p ) { --h; }
        while ( l < h && val[l] >= p ) { ++l; }
        
        if ( l < h ) {
            tmp2 = val[l]; val[l] = val[h]; val[h] = tmp2;
            tmp3 = ind[l]; ind[l] = ind[h]; ind[h] = tmp3;
        }
    }
    
    tmp2 = val[l]; val[l] = val[tmp]; val[tmp] = tmp2;
    tmp3 = ind[l]; ind[l] = ind[tmp]; ind[tmp] = tmp3;
    
    return l;
}

static int dpartitiond( double *ind, double *val, int l, int h ) {
    
    int tmp = l;
    double tmp2 = 0.0, tmp3, p = val[l];
    
    while ( l < h ) {
        
        while ( l < h && val[h] >= p ) { --h; }
        while ( l < h && val[l] <= p ) { ++l; }
        
        if ( l < h ) {
            tmp2 = val[l]; val[l] = val[h]; val[h] = tmp2;
            tmp3 = ind[l]; ind[l] = ind[h]; ind[h] = tmp3;
        }
    }
    
    tmp2 = val[l]; val[l] = val[tmp]; val[tmp] = tmp2;
    tmp3 = ind[l]; ind[l] = ind[tmp]; ind[tmp] = tmp3;
    
    return l;
}

extern double HUtilGetTimeStamp( void ) {
    
    return my_clock();
}

/** @brief Symmetrize an n by n matrix whose lower triangular is filled
 *
 */
extern void HUtilMatSymmetrize( int n, double *v ) {
    
    for ( int i = 0, j; i < n; ++i ) {
        for ( j = i + 1; j < n; ++j ) {
            FULL_ENTRY(v, n, i, j) = FULL_ENTRY(v, n, j, i);
        }
    }
    
    return;
}

/* Debugging */
extern void HUtilPrintDblContent( int n, double *d ) {
    
    for ( int i = 0; i < n; ++i ) {
        printf("%5.3e, ", d[i]);
    }
    printf("\n");
    return;
}

extern void HUtilPrintIntContent( int n, int *d ) {
    
    for ( int i = 0; i < n; ++i ) {
        printf("%5d, ", d[i]);
    }
    printf("\n");
    return;
}

extern void HUtilPrintDblSum( int n, double *d ) {
    
    double ds = 0.0;
    
    for ( int i = 0; i < n; ++i ) {
        ds += d[i];
    }
    
    printf("Sum = %20.10e \n", ds);
    return;
}

extern void HUtilPrintDblAbsSum( int n, double *d ) {
    
    double ds = 0.0;
    
    for ( int i = 0; i < n; ++i ) {
        ds += fabs(d[i]);
    }
    
    printf("Sum = %20.10e \n", ds);
    return;
}

/* Sorting */
extern int HUtilCheckIfAscending( int n, int *idx ) {
    /* Check is an integer array is ascending. */
    
    for ( int i = 0; i < n - 1; ++i ) {
        if ( idx[i] > idx[i + 1] ) {
            return 0;
        }
    }

    return 1;
}

extern void HUtilSortIntbyDbl( int *data, double *ref, int low, int up ) {
    
    if ( low < up ) {
        int p = dpartitioni(data, ref, low, up);
        HUtilSortIntbyDbl(data, ref, low, p - 1);
        HUtilSortIntbyDbl(data, ref, p + 1, up);
    }
    
    return;
}

extern void HUtilDescendSortIntByInt( int *data, int *ref, int low, int up ) {
    
    if ( low < up ) {
        int p = ipartitioni(data, ref, low, up);
        HUtilDescendSortIntByInt(data, ref, low, p - 1);
        HUtilDescendSortIntByInt(data, ref, p + 1, up);
    }
    
    return;
}

extern void HUtilAscendSortDblByInt( double *data, int *ref, int low, int up ) {
    
    if ( low < up ) {
        int p = ipartitiond(data, ref, low, up);
        HUtilAscendSortDblByInt(data, ref, low, p - 1);
        HUtilAscendSortDblByInt(data, ref, p + 1, up);
    }
    
    return;
}

extern void HUtilSortDblByDbl( double *data, double *ref, int low, int up ) {
    
    if ( low < up ) {
        int p = dpartitiond(data, ref, low, up);
        HUtilSortDblByDbl(data, ref, low, p - 1);
        HUtilSortDblByDbl(data, ref, p + 1, up);
    }
    
    return;
}

extern void HUtilStartCtrlCCheck( void ) {
    
    act.sa_handler = monitorCtrlC;
    sigaction(SIGINT, &act, NULL);
    
    return;
}

extern int HUtilCheckCtrlC( void ) {
    
    return isCtrlC;
}

extern void HUtilResetCtrl( void ) {
    
    isCtrlC = 0;
}
