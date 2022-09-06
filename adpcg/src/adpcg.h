/** @file adpcg.h
 *  @brief Header for basic types and routine list
 *
 * Given a set of positive definite linear systems \f$ A^k x = b^k \f$, adpcg solves them adaptively with
 * pre-conditioning conjugate gradient method.
 *
 * The routine is employed in HDSDP.
 *
 * @author Wenzhi Gao, Shanghai University of Finance and Economics
 * @date Aug 29th, 2022
 *
 */

#ifndef adpcg_h
#define adpcg_h


#include <stddef.h>

#ifdef ADPCG_64
typedef int64_t cgint;
#define id "%lld"
#else
typedef int32_t cgint;
#define id "%d"
#endif

#define cgerr(x) printf(x);

/* Return code */
#define CG_OK                (0)
#define CG_ERR               (1)

/* Boolean */
#define CG_TRUE              (1)
#define CG_FALSE             (0)

/* Pre-conditioner */
#define CG_PRECOND_DIAG     (10)  ///< Use diagonal as the pre-conditioner
#define CG_PRECOND_CHOL     (11)  ///< Use Cholesky factor as the pre-conditioner
#define CG_NO_PRECOND       (12)  ///< Use direct solver

/* Solution status */
#define CG_STATUS_SOLVED    (100) ///< System is solved to desired accuracy within maximum iteration
#define CG_STATUS_MAXITER   (101) ///< System reaches maximum iteration
#define CG_STATUS_FAILED    (102) ///< CG fails to converge
#define CG_STATUS_UNKNOWN   (104) ///< CG status is not known. Solution not yet started
#define CG_STATUS_DIRECT    (105) ///< CG serves as a wrapper for direct solver

/* Auxiliary macros */
#ifndef MIN
#define MIN(a, b) (((a)<(b))?(a):(b))
#endif /* MIN */
#ifndef MAX
#define MAX(a, b) (((a)>(b))?(a):(b))
#endif  /* MAX */

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
    
    void *A;         ///<  LHS data
    void *r;         ///<  Residual
    void *rnew;      ///<  Workspace array
    void *d;         ///<  Workspace array
    void *pinvr;     ///<  Workspace array
    void *Ad;        ///<  Workspace array
    void *x;         ///<  CG solution vector
    void *aux;       ///<  CG auxiliary array
    void *btmp;      ///<  CG temporary RHS array
    
    cgint ptype;     ///<  Pre-conditioner type
    void *chol;      ///<  Cholesky pre-conditioner
    void *diag;      ///<  Diagonal pre-conditioner
        
    double tol;      ///< Relative tolerance of CG
    double rnrm;     ///< Residual norm
    double avgsvtime;  ///< Averate solution time of previous CG solves
    double avgfctime;  ///< Average factorization time of previous CG solves
    double currenttime; ///< Buffer of solution time in the current round
    double latesttime; ///< Time for the latest solve
    
    cgint  n;        ///<  Dimension of linear system
    cgint  niter;    ///<  Number of iterations in most recent solve
    cgint  maxiter;  ///<  Maximum number of iterations
    cgint  status;   ///<  Solution status
    cgint  reuse;    ///<  Reuse Cholesky pre-conditioner
    cgint  nused;    ///<  Number of rounds current Cholesky pre-conditioner is already used
    cgint  nmaxiter;  ///<  Number of non-successfull solves
    cgint  restart;  ///< Restart frequency
    cgint  nfactors; ///< Number of factorizes performed so far
    cgint  nrounds;  ///< Number of rounds of solves
    cgint  nsolverd; ///< Number of solves in a round
    cgint  nsolves;  ///< Number of linear systems solved
    
    /* Vector operations */
    void  (*v_init)   (void *v); ///< Initialize vector
    cgint (*v_alloc)  (void *v, cgint n); ///< Allocate memory for v
    void  (*v_free)   (void *v);          ///< Free the internal memory of vector
    void  (*v_copy)   (void *s, void *t); ///< Copy s to t
    void  (*v_reset)  (void *v);          ///< Reset vector to 0
    void  (*v_norm)   (void *v, double *nrm); ///< nrm = norm(v)
    void  (*v_axpy)   (double a, void *x, void *y); ///< y = y + a * x
    void  (*v_axpby)  (double a, void *x, double b, void *y); ///< y = a * x + b * y
    void  (*v_zaxpby) (void *z, double a, void *x, double b, void *y); ///< z = a * x + b * y
    void  (*v_dot)    (void *x, void *y, double *xTy); ///<  xTy = x' * y
    
    /* Matrix operations */
    cgint (*A_chol) (void *A); ///< Compute Cholesky of A
    cgint (*A_cond) (void *A); ///< Evaluate conditioning of A
    void  (*A_getdiag) (void *A, void *diag); ///< Compute diagonal pre-conditioner diag = diag(A)
    
    /* Matrix-vector operations */
    cgint (*diagpcd)(void *diag, void *v); ///< Diagonal Preconditioning operation v = P \ v
    cgint (*cholpcd)(void *chol, void *v, void *aux); ///< Cholesky Preconditioning operation v = P \ v
    void  (*Av) (void *A, void *v, void *Av); ///< Compute Av = A * v
    cgint (*Ainv) (void *A, void *v, void *aux); ///< Solve Ainvv = A \ v
    
} adpcg;


#endif /* adpcg_h */
