#include <stdio.h>
#include "vec.h"

extern DSDP_INT test_vec(void) {
    
    printf("Starting test \n");
    
    DSDP_INT retcode = DSDP_RETCODE_OK;
    
    // Create
    vec x;
    vec y;
    vec z;
    
    vec *px = &x;
    vec *py = &y;
    vec *pz = &z;
    
    DSDP_INT n = 10;
    
    vec_init(px);
    vec_init(py);
    vec_init(pz);
    
    // px = py = pz = 0.0
    vec_alloc(px, n);
    vec_alloc(py, n);
    vec_alloc(pz, n);
    
    vec_print(pz);
    
    // px = py = 2.0
    vec_set(px, 2.0);
    vec_copy(px, py);
    vec_print(px);
    
    // pz = 0.5
    vec_inv(pz, py);
    vec_print(pz);
    
    // px = 0.0
    vec_reset(px);
    vec_print(px);
    
    // px = 4.0
    vec_invsqr(px, pz);
    vec_print(px);
    
    // py = 3.2
    vec_axpy(0.3, px, py);
    vec_print(py);
    
    // pz = 3.6
    vec_zaxpby(pz, 0.5, px, 0.5, py);
    vec_print(pz);
    
    // px = 9.2
    vec_axpby(0.5, px, 2.0, pz);
    vec_print(pz);
    
exit_cleanup:
    
    vec_free(px);
    vec_free(py);
    vec_free(pz);
    
    printf("Test complete \n");
    
    return retcode;
    
}
