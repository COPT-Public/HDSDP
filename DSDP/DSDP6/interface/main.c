#include <stdio.h>
#include "dsdpreadsdpa.h"
#include "linsysconfig.h"

int main (int argc, char **argv) {
    // Improve efficiency of Intel MKL
    // intel_cpu_patch(); intel_mkl_patch();
    setUpMKL();
    int argc2 = 2; char *argv2[10];
    argv2[1] = "/Users/gaowenzhi/Desktop/gwz/DSDP/DSDP6/matlab/test/sdplib/trto3.dat-s";
    
    return DSDPSolveSDPA(argc2, argv2);
}
