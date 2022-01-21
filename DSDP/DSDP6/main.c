#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#ifndef MEMWATCH
#define MEMWATCH
#endif
#include "memwatch.h"

#include "dsdppardiso.h"
#include "test.h"

int main (void)
{
    mwInit();
    test_dsdp();
//    test_data();
//    test_sparse();
//    test_dense();
//    test_presolve();
//    test_vec();
//    test_pardiso();
//    test_feast();
    mwAbort();
    
    return 0;
}
