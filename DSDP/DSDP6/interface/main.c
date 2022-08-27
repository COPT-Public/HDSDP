#include <stdio.h>
#include "dsdpreadsdpa.h"
#include "linsysconfig.h"

int main (int argc, char **argv) {
    // Improve efficiency of Intel MKL
    // intel_cpu_patch(); intel_mkl_patch();
    setUpMKL();
//    int argc2 = 2; char *argv2[10];
//    argv2[0] = argv[0];
//    argv2[1] = "/Users/gaowenzhi/Desktop/EnergyInternet/code/sedumi_data/opt.dat-s";
#ifdef OPT_PRECOND
    printf("| Specially written for optimal diagonal pre-conditioning. \n");
    argv2[1] = "/Users/gaowenzhi/Desktop/opt_precond/datasets/suitesparse/sdp/494_bus.dat-s";
#endif
//    return DSDPSolveSDPA(argc2, argv2);
     return DSDPSolveSDPA(argc, argv);
}
