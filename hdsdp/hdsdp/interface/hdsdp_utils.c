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

/* TODO: Add compatibility for Windows platform */
static double my_clock( void ) {
    
    struct timeval t;
    gettimeofday(&t, NULL);
    return (1e-06 * t.tv_usec + t.tv_sec);
}

static int partition( int *ind, double *val, int l, int h ) {
    
    double tmp2 = 0, tmp3, p = val[l];
    int tmp = l;
    
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
    
    printf("Sum = %10.6e \n", ds);
    return;
}

extern void HUtilSortbyDbl( int *ind, double *val, int low, int high ) {
    
    if ( low < high ) {
        int p = partition(ind, val, low, high);
        HUtilSortbyDbl(ind, val, low, p - 1);
        HUtilSortbyDbl(ind, val, p + 1, high);
    }
    
    return;
}
