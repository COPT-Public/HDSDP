#include <stdio.h>
#include "dsdpreadsdpa.h"
#include "dsdppardiso.h"
#include "densemat.h"
#include "sparsemat.h"
// #define UNIT_TEST

#ifdef UNIT_TEST

#include "test.h"
int main (void)
{
     test_dsdp();
    
//    test_data();
//    test_sparse();
//    test_dense();
//    test_presolve();
//    test_vec();
//    test_pardiso();
//    test_feast();
    return 0;
}

#else

#define N 2000

void test(void) {
    DSDP_INT i, j, idx = 0;
    spsMat A; dsMat B; vec b;
    b.dim = N;
    
    double *x = (double *) calloc(N * N, sizeof(double));
    double *aux = (double *) calloc(N * N, sizeof(double));
    b.x = x;
    spsMatInit(&A); denseMatInit(&B);
    spsMatAlloc(&A, N); denseMatAlloc(&B, N, TRUE);
    
    for (i = 0; i < N; ++i) {
        A.p[i] = idx;
        for (j = i; j < N; ++j) {
            A.i[idx] = j;
            A.x[idx] += sqrt(i + j / N);
            if (i == j) {
                A.x[idx] = 1e+03;
            }
            idx += 1;
        }
    }
    
    A.p[N] = idx;
    
    for (i = 0; i < N; ++i) {
        for (j = 0; j < N; ++j) {
            B.array[i * N + j] += sqrt(i + j / N);
            if (i == j) {
                B.array[i * N + i] = 1e+03;
            }
        }
    }
    
    spsMatSymbolic(&A);
    spsMatIspd(&A, &idx);
    
    clock_t start = clock();
    for (i = 0; i < 2000; ++i) {
        spsMatVecSolve(&A, &b, aux);
    }
    printf("%f \n", (double) (clock() - start) / CLOCKS_PER_SEC);
    start = clock();
    char low = 'L';
    DSDP_INT dim = N, one = 1;
    double zero = 0.0, done = 1.0;
    for (i = 0; i < 2000; ++i) {
        dsymv(&low, &dim, &done, B.array, &dim, x, &one, &zero, aux, &one);
    }
    printf("%f \n", (double) (clock() - start) / CLOCKS_PER_SEC);
    spsMatFree(&A); denseMatFree(&B);
    DSDP_FREE(x); DSDP_FREE(aux);
    
}

int MKL_Get_Max_Threads(void);
int MKL_Set_Num_Threads(int nth);

int main (int argc, char **argv) {
//
//     test();
//    return 0;
    
    printf("---------------------------------------"
           "---------------------------------------"
            "----------------------\n");
    printf("| Checking running environment. \n");
    printf("| Number of threads available %d. ", MKL_Get_Max_Threads());
    int max_threads = MKL_Get_Max_Threads(), nthreads;
    if (max_threads >= 16) {
        nthreads = 16;
    } else if (max_threads >= 8) {
        nthreads = 8;
    } else {
        nthreads = max_threads;
    }
    
    MKL_Set_Num_Threads(nthreads);
    printf(" Optimizing over %d threads. \n", nthreads);
    
    
    // return DSDPSolveSDPA(argc, argv);
    int argc2 = 2;
    char *argv2[10];
    argv2[0] = argv[0];
    argv2[1] = "/Users/gaowenzhi/Desktop/gwz/DSDP/DSDP6/matlab/benchmark/sdplib/theta12.dat-s";
//    return DSDPSolveSDPA(argc2, argv2);
//
//    argv2[1] = "/Users/gaowenzhi/Desktop/gwz/DSDP/DSDP6/matlab/benchmark/sdplib/vibra3.dat-s";
    return DSDPSolveSDPA(argc2, argv2);
    
}

#endif
