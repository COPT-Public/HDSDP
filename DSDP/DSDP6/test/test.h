#ifndef test_h
#define test_h

#include "memwatch.h"
#include "vec.h"
#include "sparsemat.h"
#include "densemat.h"
#include "rankonemat.h"

#ifndef MEMWATCH
#define MEMWATCH
#endif

#include "memwatch.h"
// #define printf mwTrace

#define passed(test) printf("Test %s passed. \n", (test))

#ifdef __cplusplus
extern "C" {
#endif

/* DSDP interface */
DSDP_INT test_data(void);

/* Vector utility */
DSDP_INT test_vec(void);
DSDP_INT test_sparse(void);
DSDP_INT test_dense(void);
DSDP_INT test_pardiso(void);
DSDP_INT test_presolve(void);
DSDP_INT test_feast(void);

#ifdef __cplusplus
}
#endif

#endif /* test_h */
