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


/** @brief Working struct for the adaptive CG solver
 *
 * To make the CG solver more general, the operations on
 *  1. vector 2. matrix 3. matrix-vector
 * are presented in abstract function pointers and users can define their own implementations for various purposes
 *
 * More details: the following operations are expected
 * Vector:
 * 1. v_init(v) Initialize a vector struct
 * 2. v_alloc(v) Allocate internal memory for a vector struct
 * 3. v_free(v) Free the internal memory for a vector struct
 * 4. v_copy(s, t) Copy content from vector s to vector t
 * 5. v_reset(v) Set the content of vector 0
 * 6. v_norm(v, &nrm) Compute norm of v
 * 7. v_axpy(a, x, y)    Compute y = y + a * x
 * 8. v_zaxpby(z, a, x, b, y) Compute z = a * x + b * y
 * 9. v_dot(x, y, &dot) Compute x' * y
 *
 * Matrix:
 * 1. A_chol(A) Perform Cholesky decomposition
 * 2. is_illcond = A_cond(A) Get a boolean variable reflecting the conditioning of A
 *
 * Matrix-vector:
 * 1. diagpcd(diag, v) Perform diagonal pre-conditioning
 * 2. cholpcf(chol, v) Perform Cholesky pre-conditioning
 */
typedef struct {
    
    schurMat *A;         ///<  LHS data
    vec *r;         ///<  Residual
    vec *rnew;      ///<  Workspace array
    vec *d;         ///<  Workspace array
    vec *pinvr;     ///<  Workspace array
    vec *Ad;        ///<  Workspace array
    vec *x;         ///<  CG solution vector
    vec *aux;       ///<  CG auxiliary array
    vec *btmp;      ///<  Buffer for b
    
    DSDP_INT ptype;     ///<  Pre-conditioner type
    schurMat *chol;      ///<  Cholesky pre-conditioner
    vec *diag;      ///<  Diagonal pre-conditioner
        
    double tol;      ///< Relative tolerance of CG
    double rnrm;     ///< Residual norm
    double avgsvtime;  ///< Averate solution time of previous CG solves
    double avgfctime;  ///< Average factorization time of previous CG solves
    double currenttime; ///< Buffer of solution time in the current round
    double latesttime; ///< Time for the latest solve
    
    DSDP_INT  n;        ///<  Dimension of linear system
    DSDP_INT  niter;    ///<  Number of iterations in most recent solve
    DSDP_INT  maxiter;  ///<  Maximum number of iterations
    DSDP_INT  status;   ///<  Solution status
    DSDP_INT  reuse;    ///<  Reuse Cholesky pre-conditioner
    DSDP_INT  nused;    ///<  Number of rounds current Cholesky pre-conditioner is already used
    DSDP_INT  nmaxiter;  ///<  Number of non-successfull solves
    DSDP_INT  restart;  ///< Restart frequency
    DSDP_INT  nfactors; ///< Number of factorizes performed so far
    DSDP_INT  nrounds;  ///< Number of rounds of solves
    DSDP_INT  nsolverd; ///< Number of solves in a round
    DSDP_INT  nsolves;  ///< Number of linear systems solved
    
    /* Vector operations */
    void  (*v_init)   (vec *v); ///< Initialize vector
    DSDP_INT (*v_alloc)  (vec *v, DSDP_INT n); ///< Allocate memory for v
    void  (*v_free)   (vec *v);          ///< Free the internal memory of vector
    void  (*v_copy)   (vec *s, vec *t); ///< Copy s to t
    void  (*v_reset)  (vec *v);          ///< Reset vector to 0
    void  (*v_norm)   (vec *v, double *nrm); ///< nrm = norm(v)
    void  (*v_axpy)   (double a, vec *x, vec *y); ///< y = y + a * x
    void  (*v_axpby)  (double a, vec *x, double b, vec *y); ///< y = a * x + b * y
    void  (*v_zaxpby) (vec *z, double a, vec *x, double b, vec *y); ///< z = a * x + b * y
    void  (*v_dot)    (vec *x, vec *y, double *xTy); ///<  xTy = x' * y
    
    /* Matrix operations */
    DSDP_INT (*A_chol) (schurMat *A); ///< Compute Cholesky of A
    DSDP_INT (*A_cond) (schurMat *A); ///< Evaluate conditioning of A
    void  (*A_getdiag) (schurMat *A, vec *diag); ///< Compute diagonal pre-conditioner diag = diag(A)
    
    /* Matrix-vector operations */
    DSDP_INT (*diagpcd)(vec *diag, vec *v); ///< Diagonal Preconditioning operation v = P \ v
    DSDP_INT (*cholpcd)(schurMat *chol, vec *v, vec *aux); ///< Cholesky Preconditioning operation v = P \ v
    void  (*Av) (schurMat *A, vec *v, vec *Av); ///< Compute Av = A * v
    DSDP_INT (*Ainv) (schurMat *A, vec *v, vec *aux); ///< Solve Ainvv = A \ v
    
} adpcg;


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
