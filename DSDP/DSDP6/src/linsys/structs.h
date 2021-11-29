#ifndef structs_h
#define structs_h

#include "cs.h"
#include "dsdphsd.h"
#include "dsdppardiso.h"

typedef struct {
    
    DSDP_INT dim;                       // Dimension of the sparse matrix
    DSDP_INT isFactorized;              // Whether the sparse matrix has been factorized
    cs_dl    *cscMat;                   // CSC matrix data using CXSparse
    void     *pdsWorker[PARDISOINDEX];  // Pardiso working array
    
} spsMat;


typedef struct {
    
    DSDP_INT dim;          // Dimension of the matrix
    DSDP_INT isFactorized; // Whether the dense matrix is factorized
    double   *array;       // The (dim + 1) * dim / 2 array
    double   *lfactor;     // The packed Cholesky factor returned by dppsv or other Lapack routines
    
} dsMat;

/* In DSDP rank-1 matrix is represented by 1 or -1 * a * a' and a is used to represent the matrix */
typedef struct {
    
    double   sign;   // The sign before the vector
    DSDP_INT dim;    // Dimension of the rank 1 matrix
    double   *x;     // Vector a
    DSDP_INT *nzIdx; // Index of nonzero elements
    DSDP_INT nnz;    // Number of nonzero elements in a
    
} r1Mat;

typedef struct {

    DSDP_INT dim;  // Dimension of the vectorß
    double   *x;   // Array storing the data
    
} vec;

#endif /* structs_h */
