#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "dsdpfeast.h"
#include "dsdpeigfact.h"
#include "test.h"
#ifndef max
#define max(a, b) (a) < (b) ? (b): (a)
#endif

char etype[] = "Test feast";
DSDP_INT n = 10;

DSDP_INT feastAp[] = {0, 2, 3, 6, 6, 7, 8, 8, 9, 9, 9};
DSDP_INT feastAi[] = {4, 9, 6, 2, 5, 8, 6, 6, 9};
double feastAx[]  =  {-0.65681593,  1.26155072,  0.47542481, -0.33733642,  0.155489,
                      -1.48139907,  0.88415371,  0.12694707,  1.12873645};

double feastpackA[] = {145.0902, 107.4532, -2.4046, -0.6691, -11.3401, 9.7020, 137.9264,
                  188.1493, 22.3598, 81.2420, 495.2369, 32.9430, 25.1972, 4.0607,
                  2.8208, 190.7778, 122.7086, 80.0744, -2.2327, 29.9940, 11.3035,
                  14.2333, 12.8114, -0.2192, 1.9181, 23.7690, -2.8766, 87.7535,
                  71.1435, 5.6926, 4.3115, 5.4869, 3.2918, -0.2061, 254.9858, 44.6697,
                  1.8421, 43.1823, 0.0087, 7.7144, 70.8160, 3.0342, 67.2678, 1.4728,
                  16.5201, 172.6808, 188.3668, 24.0036, 41.3076, 406.5596, 7.6074,
                  214.3802, 61.4641, -22.6640, 349.4797};

DSDP_INT test_feast(void) {
    
    DSDP_INT retcode  = DSDP_RETCODE_OK;
    
    DSDP_INT packAnnz = feastAp[n];
    spsMat *data  = NULL;
    data    = (spsMat *) calloc(1, sizeof(spsMat));
    retcode = spsMatInit(data); checkCodeFree;
    retcode = spsMatAllocData(data, n, packAnnz);
    
    memcpy(data->p, feastAp, sizeof(DSDP_INT) * (n + 1));
    memcpy(data->i, feastAi, sizeof(DSDP_INT) * packAnnz);
    memcpy(data->x, feastAx, sizeof(double)   * packAnnz);
    
    double *eigvals = (double *) calloc(n, sizeof(double));
    double *eigvecs = (double *) calloc(n * n, sizeof(double));
    
    factorizeSparseData(data, -35.0, eigvals, eigvecs);
    
    for (DSDP_INT i = 0; i < n; ++i) {
        printf("lambda %d = %10.5e \n", i + 1, eigvals[i]);
    }
    
    dsMat *data2 = (dsMat *) calloc(1, sizeof(dsMat));
    retcode = denseMatAlloc(data2, n, FALSE); checkCodeFree;
    memcpy(data2->array, feastpackA, sizeof(double) * nsym(n));
    
    factorizeDenseData(data2, -1064.0, eigvals, eigvecs);
    for (DSDP_INT i = 0; i < n; ++i) {
        printf("lambda %d = %10.5e \n", i + 1, eigvals[i]);
    }

clean_up:
    
    DSDP_FREE(eigvals);
    DSDP_FREE(eigvecs);
    spsMatFree(data);
    DSDP_FREE(data);
    denseMatFree(data2);
    DSDP_FREE(data2);
    
    return 0;
}
