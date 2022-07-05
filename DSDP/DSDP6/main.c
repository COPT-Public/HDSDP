#include <stdio.h>
#include "hijack/intel_cpu_feature_patch.c"
#include "hijack/intel_mkl_cpuid_patch.c"
#include "hijack/intel_mkl_feature_patch.c"
#include "dsdpreadsdpa.h"

int MKL_Get_Max_Threads(void);
int MKL_Set_Num_Threads(int nth);

void setUpMKL(void) {
    printf("---------------------------------------"
           "---------------------------------------"
            "----------------------\n");
    printf("| Checking running environment. \n");
    printf("| Number of threads available: %d. ", MKL_Get_Max_Threads());
    int max_threads = MKL_Get_Max_Threads(), nthreads;
    if (max_threads >= 16) {
        nthreads = 16;
    } else if (max_threads >= 8) {
        nthreads = 8;
    } else {
        nthreads = max_threads;
    }
    
    nthreads = 8;
    MKL_Set_Num_Threads(nthreads);
    printf(" Optimizing over %d threads. \n", nthreads);
}

int main (int argc, char **argv) {
    // Hijack Intel MKL
    intel_cpu_patch(); intel_mkl_patch();
    
    setUpMKL();
    // return DSDPSolveSDPA(argc, argv);
    int argc2 = 2; char *argv2[10];
    argv2[0] = argv[0];
    argv2[1] = "/Users/gaowenzhi/Desktop/gwz/DSDP/DSDP6/matlab/benchmark/sdplib/G48mc.dat-s";
    return DSDPSolveSDPA(argc2, argv2);
}
