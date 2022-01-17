#include "dsdpdinfeas.h"
#include "dsdputils.h"
#include "residualsetup.h"

// Implement the phase A of the dual scaling algorithm: Dual infeasibility elimination
static char etype[] = "DSDP Dual infeasibility elimination";
/* Driver routine for dual infeasibility elimination
 This routine
 
 
 */
extern DSDP_INT DSDPDInfeasEliminator( HSDSolver *dsdpSolver ) {
    
    DSDP_INT retcode = DSDP_RETCODE_OK;
    
    // dObj = b' * y;
    retcode = getDualObj(dsdpSolver);
    
    // [M, u, asinv, ~, ~, ~, csinv, csinvcsinv, asinvrysinv, csinvrysinv] = ...
    // dsdpgetSchur(A, S, C, Rd, initstrategy);
    
    
    
    
    
    
    
    
    return retcode;
}


