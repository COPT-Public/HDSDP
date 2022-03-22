#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "sparsemat.h"
#include "dsdphsd.h"
#include "dsdppardiso.h"

// Define the format to printf MKL_INT values
#define IFORMAT ID
#define MKL_INT DSDP_INT

void set_eye(double *x) {
    memset(x, 0, sizeof(double) * 16);
    for (int i = 0; i < 4; ++i) {
        x[i * 4 + i] = 1.0;
    }
}



int test_pardiso(void) {
    
    int          n = 4;
    
    /*
     A =

        0.331829000000000                   0                   0  -0.000157649000000
                        0   0.033182900000000                   0  -0.022505800000000
                        0                   0   0.331829000000000   0.022470700000000
       -0.000157649000000  -0.022505800000000   0.022470700000000   0.580879000000000
    */
    int      Ap[5] = {0, 4, 7, 9, 10};
    int      Ai[10] = {0, 1, 2, 3, 1, 2, 3, 2, 3, 3};
    double   Ax[10] = {0.331829, 0.0, 0.0, -0.000157649, 0.03318289,
                      0.0, -0.022505758, 0.331829, 0.0224707, 0.580879};
    double  eye[16] = {1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0,
                       0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0};
    double  aux[16] = {1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0,
                       0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0};
    double Linv[16] = {0.0};
    int    perm[4] = {0, 1, 2, 3};
    
    int maxfct  = 1; // Maximum number of factors
    int mnum    = 1; // The matrix used for the solution phase
    int mtype   = 2; // Real and symmetric positive definite
    int msglvl  = 0; // Print information
    int idummy  = 0; // Dummy variable for taking up space
    int phs12   = 12;
    int phs331  = 331;
    int phs33   = 33;
    int phsfree = -1;
    int error   = 0;
    
    void *pdsWorker[64];
    memset(pdsWorker, 0, sizeof(void *) * 64);
    
    int params[64] = {
        
        1, /* Non-default value */ 3, /* P Nested dissection */ 0, /* Reserved          */
        0, /* No CG             */ 2, /* No user permutation */ 0, /* No overwriting    */
        0, /* Refinement report */ 0, /* Auto ItRef step     */ 0, /* Reserved          */
        8, /* Perturb           */ 0, /* Disable scaling     */ 0, /* No transpose      */
        0, /* Disable matching  */ 0, /* Report on pivots    */ 0, /* Output            */
        0, /* Output            */ 0, /* Output              */-1, /* No report         */
        0, /* No report         */ 0, /* Output              */ 0, /* Pivoting          */
        0, /* nPosEigVals       */ 0, /* nNegEigVals         */ 0, /* Classic factorize */
        0,                         0,                           1, /* Matrix checker    */
        0,                         0,                           0,
        0,                         0,                           0,
        0,                         1, /* 0-based solve       */ 0,
        0,                         0,                           0,
        0,                         0,                           0,
        0,                         0,                           0,
        0,                         0,                           0,
        0,                         0,                           0,
        0,                         0,                           0,
        0,                         0, /* No diagonal         */ 0,
        0,                         0,                           0,
        0,                         0,                           0,
        0
    };
    
    pardiso(pdsWorker, &maxfct, &mnum, &mtype, &phs12, &n,
            Ax, Ap, Ai, perm, &idummy, params, &msglvl,
            eye, aux, &error);
    
    pardiso(pdsWorker, &maxfct, &mnum, &mtype, &phs331, &n,
            Ax, Ap, Ai, &idummy, &n, params, &msglvl,
            eye, Linv, &error);
    
    
    
    for (int i = 0; i < n; ++i) {
        // printf("Row %d: ", i);
        for (int j = 0; j < n; ++j) {
            printf("%20.20e, ", Linv[i + n * j]);
        }
        printf("\n");
    }
    
    double tmp; printf("\n \n Linv^T * L^-T: \n");
    for (int i = 0, j, k; i < 4; ++i) {\
        // printf("Row %d: ", i);
        for (j = 0; j < 4; ++j) {
            tmp = 0.0;
            for (k = 0; k < 4; ++k) {
                tmp += Linv[i * 4 + k] * Linv[j * n + k];
            }
            printf("%20.20e, ", tmp);
        }
        printf("\n");
    }
    
    set_eye(eye);
    pardiso(pdsWorker, &maxfct, &mnum, &mtype, &phs33, &n,
            Ax, Ap, Ai, &idummy, &n, params, &msglvl,
            eye, Linv, &error);
    
    printf("\n \n Sinv: \n");
    for (int i = 0; i < n; ++i) {
        // printf("Row %d: ", i);
        for (int j = 0; j < n; ++j) {
            printf("%20.20e, ", Linv[i + n * j]);
        }
        printf("\n");
    }
    
    
    pardiso(pdsWorker, &maxfct, &mnum, &mtype, &phsfree, &n,
            NULL, NULL, NULL, &idummy, &idummy, params, &msglvl,
            NULL, NULL, &error);
    
    return error;
}
