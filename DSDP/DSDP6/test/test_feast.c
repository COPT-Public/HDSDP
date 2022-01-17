/*******************************************************************************
* Copyright 2005-2021 Intel Corporation.
*
* This software and the related documents are Intel copyrighted  materials,  and
* your use of  them is  governed by the  express license  under which  they were
* provided to you (License).  Unless the License provides otherwise, you may not
* use, modify, copy, publish, distribute,  disclose or transmit this software or
* the related documents without Intel's prior written permission.
*
* This software and the related documents  are provided as  is,  with no express
* or implied  warranties,  other  than those  that are  expressly stated  in the
* License.
*******************************************************************************/

/*
!   Content: Example for k Max/Min eigenvalue problem based on
!            Intel(R) Math Kernel Library (Intel(R) MKL) Extended
!            Eigensolver (CSR sparse format, double precision)
!
!*******************************************************************************
!
! The following routines are used in the example:
!          MKL_SPARSE_D_EV
!
! Consider the 4x4 matrix A
!
!                 |  6   2   0   0   |
!                 |  2   3   0   0   |
!     A   =       |  0   0   2  -1   |
!                 |  0   0  -1   2   |
!
! stored as sparse matrix.
!
!
!  The test calls mkl_sparse_d_ev routine to find several largest singular
!  values and corresponding right-singular vectors. Orthogonality of singular
!  vectors is tested using DGEMM routine
!
!*******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "dsdpfeast.h"
#include "test.h"
#ifndef max
#define max(a, b) (a) < (b) ? (b): (a)
#endif

DSDP_INT test_feast(void) {
    
    /* Matrix A of size N in CSR format */
    DSDP_INT     N = 4;                                   /* number of rows in matrix A */
    DSDP_INT ia[5] = {1,3,5,7,9};                         /* ia array from CSR format */
    DSDP_INT ja[8] = {1,2,1,2,3,4,3,4};                   /* ja array from CSR format */
    double   a[8] = {6.0,2.0,2.0,3.0,2.0,-1.0,-1.0,2.0}; /* val array from CSR format */

    double   Eig[4] = {1.0, 2.0, 3.0, 7.0}; /* Exact eigenvalues */

    /* mkl_sparse_d_ev input parameters */
    char         which = 'S'; /* Which eigenvalues to calculate. ('L' - largest (algebraic) eigenvalues, 'S' - smallest (algebraic) eigenvalues) */
    DSDP_INT      k0  = 3;     /* Desired number of max/min eigenvalues */

    /* mkl_sparse_d_ev output parameters */
    DSDP_INT      k;      /* Number of eigenvalues found (might be less than k0). */
    double       E[4];   /* Eigenvalues */
    double       X[4*4]; /* Eigenvectors */
    double       res[4]; /* Residual */

    /* Local variables */
    DSDP_INT      info;               /* Errors */
    DSDP_INT      i, j;
/* Leading dimension for destination array in GEMM */

    /* Sparse BLAS IE variables */
    sparse_matrix_t A = NULL; /* Handle containing sparse matrix in internal data structure */
    struct matrix_descr descr; /* Structure specifying sparse matrix properties */

    /* Create handle for matrix A stored in CSR format */
    descr.type = SPARSE_MATRIX_TYPE_GENERAL; /* Full matrix is stored */
    mkl_sparse_d_create_csr ( &A, SPARSE_INDEX_BASE_ONE, N, N, ia, ia+1, ja, a );

    /* Step 2. Call mkl_sparse_ee_init to define default input values */
    mkl_sparse_ee_init(pm);

    /* Step 3. Solve the standard Ax = ex eigenvalue problem. */
    info = mkl_sparse_d_ev(&which, pm, A, descr, k0, &k, E, X, res);

    printf("mkl_sparse_d_ev output info %d \n",info);
    if ( info != 0 )
    {
        printf("Routine mkl_sparse_d_ev returns code of ERROR: %i", (int)info);
        return 1;
    } else if (pm[9] != 0) {
        printf("Routine mkl_sparse_d_ev returns success but the reason for"
               " exiting iterations is %d != 0 (convergence)", (int)pm[9]);
        return 2;
    }

    printf("*************************************************\n");
    printf("************** REPORT ***************************\n");
    printf("*************************************************\n");
    printf("#mode found/subspace %d %d \n", k, k0);
    if (pm[2] == 1)
        printf("reason for iteration exit %d\n", pm[9]);
    printf("Index/Exact Eigenvalues/Estimated Eigenvalues/Residuals\n");
    for (i=0; i<k; i++)
    {
        printf("   %d  %.15e %.15e %.15e \n",i, Eig[i], E[i], res[i]);
    }

    mkl_sparse_destroy(A);
    return 0;
}
