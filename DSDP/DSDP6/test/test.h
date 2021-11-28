#ifndef test_h
#define test_h

#include "vec.h"
#include "sparsemat.h"
#include "densemat.h"
#include "rankonemat.h"

#define passed(test) printf("Test %s passed. \n", (test))

#ifdef __cplusplus
extern "C" {
#endif

/* Vector utility */
DSDP_INT test_vec(void);
DSDP_INT test_dense(void);
DSDP_INT test_pardiso(void);
DSDP_INT test_presolve(void);

#ifdef __cplusplus
}
#endif

#endif /* test_h */
