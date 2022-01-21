#include <stdio.h>
#include <stdlib.h>
#include <math.h>


#include "dsdppardiso.h"
#include "test.h"

int main (void)
{
    
    for (int i = 0; i < 20; ++i ) {
        test_dsdp();
    }
    
//    test_data();
//    test_sparse();
//    test_dense();
//    test_presolve();
//    test_vec();
//    test_pardiso();
//    test_feast();
    
    return 0;
}
