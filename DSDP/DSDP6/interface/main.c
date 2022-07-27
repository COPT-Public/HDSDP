#include <stdio.h>
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
    if      (max_threads >= 16) { nthreads = 16; }
    else if (max_threads >= 12) { nthreads = 12; }
    else if (max_threads >= 8)  { nthreads = 8;  }
    else                        { nthreads = max_threads; }
    printf(" Optimizing over %d threads. \n", nthreads);
    MKL_Set_Num_Threads(nthreads);
}

int main (int argc, char **argv) {
    // Improve efficiency of Intel MKL
    // intel_cpu_patch(); intel_mkl_patch();
    setUpMKL();
//    int argc2 = 2; char *argv2[10];
//    argv2[0] = argv[0];
//    argv2[1] = "/Users/gaowenzhi/Desktop/gwz/DSDP/DSDP6/matlab/test/sdplib/biggs.dat-s";
#ifdef OPT_PRECOND
    printf("| Specially written for optimal diagonal pre-conditioning. \n");
    argv2[1] = "/Users/gaowenzhi/Desktop/opt_precond/datasets/suitesparse/sdp/494_bus.dat-s";
#endif
//    return DSDPSolveSDPA(argc2, argv2);
     return DSDPSolveSDPA(argc, argv);
}
