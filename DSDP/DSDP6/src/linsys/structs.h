#ifndef structs_h
#define structs_h

#include "dsdphsd.h"
#include "pardiso.h"

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
    DSDP_INT *cidx;                     // Position of non-zero entry in packed format
    
    rkMat    *factor;                   // Eigenvalue decomposition
    void     *pdsWorker[PARDISOINDEX];  // Pardiso working array
    
    DSDP_INT nominalsps;                // Is matrix actually dense
    double   *Sinv;                     // Pointer towards Sinv
    
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

/** @brief The struct that implements the Lanczos routine for SDP step length computation
 * The struct contains the necessary workspace to compute the step length for a certain SDP block
 */
typedef struct {
    
    DSDP_INT n;  ///< Size of the matrix
    vec *v;      ///< Workspace array
    vec *w;      ///< Workspace array
    vec *z1;     ///< Workspace array
    vec *z2;     ///< Workspace array
    vec *va; ///< Workspace array
    
    double *V;   ///< Workspace array
    double *H;   ///< Workspace array
    double *Y;   ///< Workspace array
    double *d;   ///< Workspace array
    
    double *mataux; ///< Workspace array
    double *eigaux; ///< Workspace array
    
    DSDP_INT *eigintaux; ///< Workspace array
    DSDP_INT isuppz[4]; ///< Workspace array
    
    DSDP_INT iwork; ///< Workspace array
    DSDP_INT lwork; ///< Workspace array
    
} lczstep;

/** @brief The struct that implements pre-conditioned conjugate gradient solver
 * The struct contains the necessary workspace to solve a linear system using PCG method
 */
typedef struct {
    
    schurMat *M;        ///< LHS data
    vec      *r;        ///<  Residual
    vec      *rnew;     ///<  Workspace array
    vec      *d;        ///<  Workspace array
    vec      *Pinvr;    ///<  Workspace array
    vec      *Md;       ///<  Workspace array
    vec      *x;        ///<  CG solution vector
    vec      *aux;      ///<  CG auxiliary array
    
    DSDP_INT pType;     ///<  Pre-conditioner type
    schurMat *cholPre;  ///<  Cholesky pre-conditioner
    vec      *vecPre;   ///<  Diagonal pre-conditioner
        
    double   tol;       ///<  Relative tolerance of CG
    double   resinrm;   ///<  Residual norm
    DSDP_INT dim;       ///<  Dimension of linear system
    DSDP_INT niter;     ///<  Number of iterations
    DSDP_INT maxiter;   ///<  Maximum number of iterations
    DSDP_INT status;    ///<  Solution status
    DSDP_INT reuse;     ///<  Reuse Cholesky pre-conditioner
    DSDP_INT nused;     ///<  Number of  iterations current pre-conditioner is already used
    DSDP_INT nfailed;   ///<  Number of non-successfull solves
    
} cgsol;

/** @brief The struct that wraps up SPEIG package
 * The struct contains the necessary workspace to factorize a sequence of SDP coefficient matrices
 */
typedef struct {
    
    DSDP_INT  nmax;     ///<Maximum possible dimension of matrix to factorize
    DSDP_INT  lwork;    ///< Length of working space
    DSDP_INT  liwork;   ///< Length of working space
    
    double   *dwork;    ///< Double working space
    double   *dworkmat; ///< Double working space
    double   *dworkevc; ///< Double working space
    double   *dworkevl; ///< Double working space
    
    DSDP_INT *perm;     ///< Permutation vector
    DSDP_INT *pinv;     ///< Inverse of the permutation
    DSDP_INT *cnz;      ///< Column nnz counter
    DSDP_INT *iwork;    ///< Integer working space
    DSDP_INT *iworkup;  ///< Another integer working space
    
} speigfac;

/* Macros for matrix indexing */
#define packIdx(A, n, i, j) (A[(DSDP_INT)((2 * (n) - (j) - 1) * (j) / 2) + (i)])
#define fullIdx(A, n, i, j) (A[(j) * (n) + (i)])

#endif /* structs_h */
