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

#endif /* dsdpfeast_h */
