#include <stdio.h>
#include "dsdpreadsdpa.h"
#include "linsysconfig.h"

int main (int argc, char **argv) {
    // Improve efficiency of Intel MKL
    // intel_cpu_patch(); intel_mkl_patch();
    setUpMKL();
    return DSDPSolveSDPA(argc, argv);
}
