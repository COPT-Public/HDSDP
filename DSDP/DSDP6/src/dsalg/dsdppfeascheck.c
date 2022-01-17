#include "dsdppfeascheck.h"
#include "dsdputils.h"
/* Implement DSDP primal feasibility check */

extern DSDP_INT dsdpCheckPhaseAPfeas( HSDSolver *dsdpSolver, double dtaudelta, vec *dydelta, DSDP_INT *ispfeas ) {
    // Phase A feasibility check
    // - Ry + C * (tau - dtaudelta) - dsdpgetATy(A, y - dydelta
    DSDP_INT retcode = DSDP_RETCODE_OK;
    
    vec_axpby(1.0, dsdpSolver->y, -1.0, dydelta);
    getPhaseAS(dsdpSolver, dydelta->x, dsdpSolver->tau - dtaudelta);
    dsdpInCone(dsdpSolver, ispfeas);
    
    return retcode;
}
