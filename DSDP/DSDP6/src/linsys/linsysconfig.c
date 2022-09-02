#include "linsysconfig.h"
// Config the MKL package

int MKL_Get_Max_Threads(void);
int MKL_Set_Num_Threads(int nth);

void setUpMKL(void) {
    
    printf("---------------------------------------"
           "---------------------------------------"
            "--------------------\n");
    printf("| Checking running environment. \n");
    printf("| Number of threads available: %d. ", MKL_Get_Max_Threads());
    
    int max_threads = MKL_Get_Max_Threads(), nthreads;
    if      (max_threads >= 16) { nthreads = 16; }
    else if (max_threads >= 12) { nthreads = 12; }
    else if (max_threads >= 8)  { nthreads = 8;  }
    else                        { nthreads = max_threads; }
    nthreads = 1;
    printf(" Optimizing over %d threads. \n", nthreads);
    MKL_Set_Num_Threads(nthreads);
}
