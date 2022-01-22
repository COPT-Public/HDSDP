#include "dsdppfeascheck.h"
#include "dsdputils.h"

/* Implement DSDP primal feasibility check */
extern DSDP_INT dsdpCheckPhaseAPfeas( HSDSolver *dsdpSolver, double dtaudelta,
                                     vec *dydelta, DSDP_INT *ispfeas ) {
    // Phase A feasibility check
    // - Ry + C * (tau - dtaudelta) - dsdpgetATy(A, y - dydelta)
    DSDP_INT retcode = DSDP_RETCODE_OK;
    
    vec_axpby(1.0, dsdpSolver->y, -1.0, dydelta);
    retcode = getPhaseACheckerS(dsdpSolver, dydelta->x, dsdpSolver->tau - dtaudelta);
    retcode = dsdpCheckerInCone(dsdpSolver, ispfeas);
    
    return retcode;
}

extern DSDP_INT dsdpCheckBackwardNewton( HSDSolver *dsdpSolver, DSDP_INT *ispfeas ) {
    
    // Check backward newton step C - dsdpgetATy(A, y - dymuprimal)
    DSDP_INT retcode = DSDP_RETCODE_OK;
    
    vec_zaxpby(dsdpSolver->b2, 1.0, dsdpSolver->y,
               -1.0, dsdpSolver->b1);
    
    retcode = getPhaseBCheckerS(dsdpSolver, dsdpSolver->b2->x);
    retcode = dsdpCheckerInCone(dsdpSolver, ispfeas);
    
    return retcode;
}
