#ifndef structs_h
#define structs_h

#include <string.h>
#include "cs.h"
#include "dsdphsd.h"
#include "dsdppardiso.h"

/* In DSDP rank-1 matrix is represented by d or -d * a * a' and a is used to represent the matrix */
typedef struct {
    
    double   sign;   // The sign before the vector
    DSDP_INT dim;    // Dimension of the rank 1 matrix
    double   *x;     // Vector a
    DSDP_INT *nzIdx; // Index of nonzero elements
    DSDP_INT nnz;    // Number of nonzero elements in a
    
} r1Mat;

/* In DSDP rank-k matrix is represented by the sum of multiple rank-1 matrices */
typedef struct {
    
    DSDP_INT dim;
    DSDP_INT isdata;
    DSDP_INT rank;
    r1Mat    **data;
    
} rkMat;

typedef struct {
    
    DSDP_INT dim;                       // Dimension of the sparse matrix
    DSDP_INT isFactorized;              // Whether the sparse matrix has been factorized
    DSDP_INT *p;                        // Column pointer
    DSDP_INT *i;                        // Row index
    double   *x;                        // Sparse data
    DSDP_INT nnz;                       // Nonzeros entries
    DSDP_INT *nzHash;                   // Position of non-zero entry in packed format
    
    rkMat    *factor;                   // Eigenvalue decomposition
    void     *pdsWorker[PARDISOINDEX];  // Pardiso working array
    
} spsMat;


typedef struct {
    
    DSDP_INT dim;          // Dimension of the matrix
    DSDP_INT isFactorized; // Whether the dense matrix is factorized
    DSDP_INT isillCond;    // Whether the matrix suffers from ill conditioning
    DSDP_INT lwork;        // Length of the working array
    
    rkMat    *factor;      // Eigenvalue decomposition
    double   *array;       // The schur matrix
    double   *lfactor;     // Cholesky factor
    DSDP_INT *ipiv;        // Array in case of indefinite Schur matrix
    double   *work;        // Working array for LDLT
    
} dsMat;


typedef struct {

    DSDP_INT dim;  // Dimension of the vectorß
    double   *x;   // Array storing the data
    
} vec;

typedef struct {
    
    DSDP_INT m;
    DSDP_INT stype;  // Type of the schur matrix
    dsMat  *denseM;  // Dense schur matrix
    spsMat *spsM;    // Sparse schur matrix
    
    DSDP_INT isillCond;    // Is the matrix ill-conditioned
    DSDP_INT isFactorized; // Is the matrix factorized
    
    double **diag;   // Store the address of the diagonal array
    
} schurMat;

#define packIdx(P, n, i, j) (P[(DSDP_INT)((2 * (n) - (j) - 1) * (j) / 2) + (i)])
#define fullIdx(P, n, i, j) (P[(j) * (n) + (i)])

#endif /* structs_h */
