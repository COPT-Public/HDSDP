/** @file def\_hdsdp\_linsolver.h
 *  @brief HDSDP linear system solver
 */
#ifndef def_hdsdp_linsolver_h
#define def_hdsdp_linsolver_h

#include "interface/hdsdp.h"

typedef enum {
    /* Direct solver supports both Schur complement and matrix variable */
    HDSDP_LINSYS_DENSE_DIRECT,
    HDSDP_LINSYS_SMALL_DIRECT, /* Tailored for extremely small cones of size <= 3 */
    HDSDP_LINSYS_SPARSE_DIRECT,
    
    /* Iterative solver is only used for Schur complement */
    HDSDP_LINSYS_DENSE_ITERATIVE,
    HDSDP_LINSYS_SPARSE_ITERATIVE
    
} linsys_type;

/* Define linear system solver */
typedef struct {
    
    int nCol;
    void *chol;
    linsys_type LinType;
    
    hdsdp_retcode (*cholCreate) ( void **, int );
    void (*cholSetParam) ( void *, void * );
    hdsdp_retcode (*cholSymbolic) ( void *, int *, int * );
    hdsdp_retcode (*cholNumeric) ( void *, int *, int *, double * );
    hdsdp_retcode (*cholPsdCheck) ( void *, int *, int *, double *, int * );
    
    void (*cholFSolve) ( void *, int, double *, double * );
    void (*cholBSolve) ( void *, int, double *, double * );
    hdsdp_retcode (*cholSolve) ( void *, int, double * );
    hdsdp_retcode (*cholGetDiag) ( void *, double * );
    void (*cholInvert) ( void *, double *, double * );
    
    void (*cholDestroy) ( void ** );
    
} hdsdp_linsys;

/* Sparse direct */
typedef struct {
    
    int nCol;
    int *colMatBeg;
    int *colMatIdx;
    double *colMatElem;
    
    double *dWork;
    
    void *pt[64];
    int iparm[64];
    
} pardiso_linsys;

/* Dense direct */
typedef struct {
    
    int nCol;
    double *dFullMatElem;
    
} lapack_linsys;

typedef struct {
    
    int nCol;
    double dSmallMatElem[9];
    
} small_linsys;

typedef struct {
    
    int maxIter;
    int nRestartFreq;
    double absTol;
    double relTol;
    
} iterative_params;

typedef enum {
    
    ITERATIVE_STATUS_OK,
    ITERATIVE_STATUS_NUMERICAL,
    ITERATIVE_STATUS_MAXITER,
    ITERATIVE_STATUS_FAILED
    
} iter_status;

/* Dense iterative */
typedef struct {
    
    int nCol;
    
    double *pFullMatElem;
    double *iterResi;
    double *iterResiNew;
    double *iterDirection;
    double *preInvResi;
    double *MTimesDirection;
    double *iterVec;
    double *iterVecAuxi;
    double *rhsBuffer;
    
    /* Pre-conditioner */
    int useJacobi;
    double *JacobiPrecond;
    lapack_linsys *lap;
    
    /* Statistics */
    double iterResiNorm;
    double avgSolveTime;
    double avgFactorTime;
    
    int nIters;
    iter_status solStatus;
    int nFactors;
    int nSolves;
    
    iterative_params params;
    
} iterative_linsys;

/* Sparse direct */
#define PARDISO_RET_OK          ( 0)
#define PARDISO_RET_INDEFINITE  (-4)
#define PARDISO_SYM_POSDEFINITE ( 2)
#define PARDISO_SYM_INDEFINITE  (-2)
#define PARDISO_PHASE_SYM       (11)      // Pardiso symbolic analysis
#define PARDISO_PHASE_SYM_FAC   (12)      // Pardiso symbolic analysis
#define PARDISO_PHASE_FAC       (22)      // Pardiso numerical factorization
#define PARDISO_PHASE_FORWARD  (331)
#define PARDISO_PHASE_BACKWARD (333)
#define PARDISO_PHASE_SOLVE     (33)      // Solve linear system
#define PARDISO_PHASE_FREE      (-1)      // Free internal data structure

#define PARDISO_PARAM_NONDEFAULT    (0)
#define PARDISO_PARAM_SYMBOLIC      (1)
#define PARDISO_PARAM_SYMBOLIC_MMD  (0)
#define PARDISO_PARAM_SYMBOLIC_ND   (2)
#define PARDISO_PARAM_REFINEMENT    (7)
#define PARDISO_PARAM_INPLACE       (5)
#define PARDISO_PARAM_PERTURBATION  (9)
#define PARDISO_PARAM_SCALING      (10)
#define PARDISO_PARAM_MATCHING     (12)
#define PARDISO_PARAM_FACNNZ       (17)
#define PARDISO_PARAM_FACFLOP      (18)
#define PARDISO_PARAM_THREADS      (33)
#define PARDISO_PARAM_INDEX        (34)
#define PARDISO_PARAM_INDEX_C       (1)
#define PARDISO_PARAM_DIAGONAL     (55)
#define PARDISO_PARAM_DIAGONAL_ON   (1)

#define set_pardiso_param(iparm, param, val) iparm[param] = val
#define get_pardiso_output(iparm, param) iparm[param]

extern void pardisoinit ( void *, int *, int * );
extern void pardiso     ( void     *, int    *, int *, int *, int *, int *,
                          double   *, int    *, int *, int *, int *, int *,
                          int *, double      *, double   *, int * );
extern void pardiso_getdiag ( const void *, void *, void *, const int *, int * );


/* Dense direct */
#define LAPACK_RET_OK    ( 0 )
#define LAPACK_UPLOW_LOW ('L')
#define LAPACK_NOTRANS   ('N')
#define LAPACK_TRANS     ('T')
#define LAPACK_SIDE_LEFT ('L')
#define LAPACK_DIAG_NONUNIT ('N')

extern void dtrsv( const char *uplo, const char *trans, const char *diag,
                   const int *n, const double *a, const int *lda, double *x,
                   const int *incx );
void dtrsm( const char *side, const char *uplo, const char *transa,
            const char *diag, const int *m, const int *n, const double *alpha,
            const double *a, const int *lda, double *b, const int *ldb );
void dpotrf ( const char *uplo, const int *n, double *a, const int *lda, int *info );
void dpotri( const char *uplo, const int *n, double *a, const int *lda, int *info );
void dpotrs( const char *uplo, const int *n, const int *nrhs, const double *a,
             const int *lda, double *b, const int *ldb, int *info );




#endif /* def_hdsdp_linsolver_h */
