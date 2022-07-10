#include "dsdpsort.h"

// Implement sorting in data extractor
static DSDP_INT partition( double *data, DSDP_INT *idxbase, DSDP_INT low, DSDP_INT high ) {
    
    DSDP_INT tmp = low, tmp2 = 0, p = idxbase[low]; double tmp3;
    while (low < high) {
        while (low < high && idxbase[high] >= p) { --high; }
        while (low < high && idxbase[low]  <= p) { ++low;  }
        if (low < high) {
            tmp2 = idxbase[low]; idxbase[low] = idxbase[high]; idxbase[high] = tmp2;
            tmp3 = data[low]; data[low] = data[high]; data[high] = tmp3;
        }
    }
    tmp2 = idxbase[low]; idxbase[low] = idxbase[tmp]; idxbase[tmp] = tmp2;
    tmp3 = data[low]; data[low] = data[tmp]; data[tmp] = tmp3;
    return low;
}

extern void dsdpSort( double *data, DSDP_INT *idxbase, DSDP_INT low, DSDP_INT high ) {
    if (low < high) {
        DSDP_INT pl = partition(data, idxbase, low, high);
        dsdpSort(data, idxbase, low, pl - 1); dsdpSort(data, idxbase, pl + 1, high);
    }
}

// Implement sorting in Schur block analysis
static DSDP_INT partition2( DSDP_INT *data, DSDP_INT *idxbase, DSDP_INT low, DSDP_INT high ) {
    
    DSDP_INT tmp = low, tmp2 = 0, p = idxbase[low];
    double tmp3;
    
    while (low < high) {
        while (low < high && idxbase[high] <= p) { --high; }
        while (low < high && idxbase[low]  >= p) { ++low;  }
        if (low < high) {
            tmp2 = idxbase[low]; idxbase[low] = idxbase[high]; idxbase[high] = tmp2;
            tmp3 = data[low]; data[low] = data[high]; data[high] = tmp3;
        }
    }
    
    tmp2 = idxbase[low]; idxbase[low] = idxbase[tmp]; idxbase[tmp] = tmp2;
    tmp3 = data[low]; data[low] = data[tmp]; data[tmp] = tmp3;
    return low;
}

extern void dsdpSort2( DSDP_INT *data, DSDP_INT *idxbase, DSDP_INT low, DSDP_INT high ) {
    if (low < high) {
        DSDP_INT pl = partition2(data, idxbase, low, high);
        dsdpSort2(data, idxbase, low, pl - 1); dsdpSort2(data, idxbase, pl + 1, high);
    }
}
