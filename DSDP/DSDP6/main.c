#include <stdio.h>
#include "dsdpreadsdpa.h"

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

int main (int argc, char **argv) {
    return DSDPSolveSDPA(argc, argv);
//    int argc2 = 2;
//    char *argv2[10];
//    argv2[0] = argv[0];
//    argv2[1] = "/Users/gaowenzhi/Desktop/gwz/DSDP/DSDP6/matlab/benchmark/sdplib/hamming_8_3_4.dat-s";
//    return DSDPSolveSDPA(argc2, argv2);
}

#endif
