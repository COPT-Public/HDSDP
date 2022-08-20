#ifndef spinfo_h
#define spinfo_h

/*
 Define other related macros for speigs
*/

/* Return code */
#define SP_EIGS_OK           (0)
#define SP_EIGS_ERR          (1)

/* Boolean */
#define TRUE                 (1)
#define FALSE                (0)

/* Matrix type */
#define MATRIX_TYPE_ZERO     (1)
#define MATRIX_TYPE_SPARSE   (2)
#define MATRIX_TYPE_GENERAL  (3)
#define MATRIX_TYPE_RANKONE  (4)
#define MATRIX_TYPE_DIAG     (5)
#define MATRIX_TYPE_TWOTWO   (6)

/* Some constants */
#define ROOT                 (7.0710678118654757273731092936941422522068e-01)
#define LAPACK_IWORK         (12)
#define LAPACK_LWORK         (30)


#endif /* spinfo_h */
