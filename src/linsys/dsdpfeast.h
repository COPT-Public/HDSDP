#ifndef dsdpfeast_h
#define dsdpfeast_h

/* A wrapper of MKL feast algorithm to compute the max/min eigenvalue of a sparse matrix */
#include "dsdphsd.h"

struct  sparse_matrix;
typedef struct sparse_matrix *sparse_matrix_t;

typedef enum {
    SPARSE_STATUS_SUCCESS           = 0,    /* the operation was successful */
    SPARSE_STATUS_NOT_INITIALIZED   = 1,    /* empty handle or matrix arrays */
    SPARSE_STATUS_ALLOC_FAILED      = 2,    /* internal error: memory allocation failed */
    SPARSE_STATUS_INVALID_VALUE     = 3,    /* invalid input value */
    SPARSE_STATUS_EXECUTION_FAILED  = 4,    /* e.g. 0-diagonal element for triangular solver, etc. */
    SPARSE_STATUS_INTERNAL_ERROR    = 5,    /* internal error */
    SPARSE_STATUS_NOT_SUPPORTED     = 6     /* e.g. operation for double precision doesn't support other types */
} sparse_status_t;

typedef enum {
    SPARSE_INDEX_BASE_ZERO  = 0,           /* C-style */
    SPARSE_INDEX_BASE_ONE   = 1            /* Fortran-style */
} sparse_index_base_t;

typedef enum {
    SPARSE_MATRIX_TYPE_GENERAL            = 20,   /*    General case                    */
    SPARSE_MATRIX_TYPE_SYMMETRIC          = 21,   /*    Triangular part of              */
    SPARSE_MATRIX_TYPE_HERMITIAN          = 22,   /*    the matrix is to be processed   */
    SPARSE_MATRIX_TYPE_TRIANGULAR         = 23,
    SPARSE_MATRIX_TYPE_DIAGONAL           = 24,   /* diagonal matrix; only diagonal elements will be processed */
    SPARSE_MATRIX_TYPE_BLOCK_TRIANGULAR   = 25,
    SPARSE_MATRIX_TYPE_BLOCK_DIAGONAL     = 26    /* block-diagonal matrix; only diagonal blocks will be processed */
} sparse_matrix_type_t;

typedef enum {
    SPARSE_FILL_MODE_LOWER  = 40,           /* lower triangular part of the matrix is stored */
    SPARSE_FILL_MODE_UPPER  = 41,            /* upper triangular part of the matrix is stored */
    SPARSE_FILL_MODE_FULL   = 42            /* upper triangular part of the matrix is stored */
} sparse_fill_mode_t;

typedef enum {
    SPARSE_DIAG_NON_UNIT    = 50,           /* triangular matrix with non-unit diagonal */
    SPARSE_DIAG_UNIT        = 51            /* triangular matrix with unit diagonal */
} sparse_diag_type_t;

struct matrix_descr {
        sparse_matrix_type_t  type;       /* matrix type: general, diagonal or triangular / symmetric / hermitian */
        sparse_fill_mode_t    mode;       /* upper or lower triangular part of the matrix ( for triangular / symmetric / hermitian case) */
        sparse_diag_type_t    diag;       /* unit or non-unit diagonal ( for triangular / symmetric / hermitian case) */
};

sparse_status_t mkl_sparse_d_create_csc( sparse_matrix_t           *A,
                                         const sparse_index_base_t indexing, /* indexing: C-style or Fortran-style */
                                         const DSDP_INT            rows,
                                         const DSDP_INT            cols,
                                         DSDP_INT                  *cols_start,
                                         DSDP_INT                  *cols_end,
                                         DSDP_INT                  *row_indx,
                                         double                    *values );

sparse_status_t mkl_sparse_d_create_csr( sparse_matrix_t           *A,
                                         const sparse_index_base_t indexing, /* indexing: C-style or Fortran-style */
                                         const DSDP_INT            rows,
                                         const DSDP_INT            cols,
                                         DSDP_INT                  *rows_start,
                                         DSDP_INT                  *rows_end,
                                         DSDP_INT                  *col_indx,
                                         double                    *values );

sparse_status_t mkl_sparse_d_ev( char                *which,
                                 DSDP_INT            *pm,
                                 sparse_matrix_t     A,
                                 struct matrix_descr descrA,
                                 DSDP_INT            k0,
                                 DSDP_INT            *k,
                                 double              *E,
                                 double              *X,
                                 double              *res );

sparse_status_t mkl_sparse_ee_init( DSDP_INT *pm);

sparse_status_t mkl_sparse_destroy( sparse_matrix_t  A );

static struct matrix_descr dsdp_descr = {
    SPARSE_MATRIX_TYPE_SYMMETRIC,
    SPARSE_FILL_MODE_UPPER,
    SPARSE_DIAG_NON_UNIT
};

static char MAX_EIG = 'L';
static char MIN_EIG = 'S';
static DSDP_INT k0  = 1; // We only need one extremal eigenvalue
static DSDP_INT k   = 0; // Used to catch MKL output
static double resi  = 0.0; // Used to catch residual
static DSDP_INT pm[128] = {
    0, 6 /* Stopping tolerance 1e-05 */,
          0 /* Decide algorithm one run time */,
             0, 0, 0, 0 /* No eigen vector (memory still needs allocation */,
                         0,
    0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0,
};

static DSDP_INT fpm[128] = {
    0,  /*  0: Runtime status */
    8,  /*  1: Countour points */
    12, /*  2: Error trace precision */
    20, /*  3: Max refinement */
    0,  /*  4: No user initial subspace */
    0,  /*  5: Stopping test */
    5,  /*  6: Error trace single precision */
    0, 0, 0, 0, 0, 0,
    0,  /* 13: Standard use */
    0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0,
    0,  /* 26: No input check */
    0,  /* 27: Check definite B */
    0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0,
    0,  /* 63: Default Pardiso iparam */
    1, /* Non-default value */ 3, /* P Nested dissection */ 0, /* Reserved          */
    0, /* No CG             */ 0, /* No user permutation */ 0, /* No overwriting    */
    0, /* Refinement report */ 0, /* Auto ItRef step     */ 0, /* Reserved          */
    6, /* Perturb           */ 0, /* Disable scaling     */ 0, /* No transpose      */
    0, /* Disable matching  */ 0, /* Report on pivots    */ 0, /* Output            */
    0, /* Output            */ 0, /* Output              */-1, /* No report         */
    0, /* No report         */ 0, /* Output              */ 1, /* Pivoting          */
    0, /* nPosEigVals       */ 0, /* nNegEigVals         */ 0, /* Classic factorize */
    0,                         0,                           0, /* Matrix checker    */
    0,                         0,                           0,
    0,                         0,                           0,
    0,                         0, /* 1-based solve       */ 0,
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

void dfeast_scsrev ( const char     *uplo,
                     const DSDP_INT *n,
                     const double   *sa,
                     const DSDP_INT *isa,
                     const DSDP_INT *jsa,
                     DSDP_INT       *fpm,
                     double         *epsout,
                     DSDP_INT       *loop,
                     const double   *emin,
                     const double   *emax,
                     DSDP_INT       *m0,
                     double         *e,
                     double         *x,
                     DSDP_INT       *m,
                     double         *res,
                     DSDP_INT       *info );

void dfeast_syev ( const char       *uplo,
                   const DSDP_INT   *n,
                   const double     *a,
                   const DSDP_INT   *lda,
                   DSDP_INT         *fpm,
                   double           *epsout,
                   DSDP_INT         *loop,
                   const double     *emin,
                   const double     *emax,
                   DSDP_INT         *m0,
                   double           *e,
                   double           *x,
                   DSDP_INT         *m,
                   double           *res,
                   DSDP_INT         *info );

#endif /* dsdpfeast_h */
