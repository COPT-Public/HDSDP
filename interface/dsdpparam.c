#include <stdio.h>
#include "dsdpparam.h"
#include "dsdplog.h"

/* Implement DSDP parameter interface */
static void printDblParam( const double *dblParams, DSDP_INT pName ) {
    
    printf("| ");
    switch (pName) {
        case DBL_PARAM_RHO          : printf("Rho (0.0, inf]: "); break;
        case DBL_PARAM_RHON         : printf("Rhon [1.0, 10.0]: "); break;
        case DBL_PARAM_INIT_POBJ    : printf("Initial pObj (-inf, inf): "); break;
        case DBL_PARAM_INIT_BETA    : printf("Initial beta [1.0, inf): "); break;
        case DBL_PARAM_INIT_MU      : printf("Initial mu (0.0, inf): "); break;
        case DBL_PARAM_INFEAS_THRESH: printf("Infeasibility threshold (0.0, inf): "); break;
        case DBL_PARAM_ABS_OPTTOL   : printf("Absolute optimality tolerance (0.0, inf): "); break;
        case DBL_PARAM_REL_OPTTOL   : printf("Relative optimality tolerance (0.0, inf): "); break;
        case DBL_PARAM_ABS_FEASTOL  : printf("Absolute feasibility tolerance (0.0, inf): "); break;
        case DBL_PARAM_REL_FEASTOL  : printf("Relative feasibixlity tolerance (0.0, inf): "); break;
        case DBL_PARAM_PRLX_PENTALTY: printf("Primal relaxation penalty (0.0, inf): "); break;
        case DBL_PARAM_BOUND_X      : printf("Primal variable bound (0.0, inf): "); break;
        case DBL_PARAM_TIMELIMIT    : printf("Time limit (0.0, inf): "); break;
        default: printf("Invalid parameter code %d \n", pName); return; break;
    }
    
    printf("%g \n", dblParams[pName]);
}

static void printIntParam( const DSDP_INT *intParams, DSDP_INT pName ) {
    
    printf("| ");
    switch (pName) {
        case INT_PARAM_ACORRECTOR: printf("A corrector step (0, 20]: "); break;
        case INT_PARAM_BCORRECTOR: printf("B corrector step [0, 20]: "); break;
        case INT_PARAM_CG_REUSE  : printf("CG pre-conditioner reuse [0, 100]: "); break;
        case INT_PARAM_AMAXITER  : printf("A maximum iteration (0, inf): "); break;
        case INT_PARAM_BMAXITER  : printf("B maximum iteration (0, inf): "); break;
        case INT_PARAM_GOLDSEARCH: printf("Golden linesearch {0, 1}: "); break;
        default: printf("Invalid parameter code %d \n", pName); return;
    }
    
    printf("%d \n", intParams[pName]);
}

static DSDP_INT checkDblParam( DSDP_INT pName, double dblVal ) {

    DSDP_INT retcode = DSDP_RETCODE_OK;
    
    switch (pName) {
        case DBL_PARAM_RHO          : if (dblVal < 0.0) { retcode = DSDP_RETCODE_FAILED; } break;
        case DBL_PARAM_RHON         : if (dblVal > 10.0 || dblVal < 1.0) { retcode = DSDP_RETCODE_FAILED; } break;
        case DBL_PARAM_INIT_POBJ    : if (dblVal >= DSDP_INFINITY || dblVal <= -DSDP_INFINITY) { retcode = DSDP_RETCODE_FAILED; } break;
        case DBL_PARAM_INIT_BETA    : if (dblVal < 0.0 || dblVal >= DSDP_INFINITY) { retcode = DSDP_RETCODE_FAILED; } break;
        case DBL_PARAM_INIT_MU      : if (dblVal < 0.0 || dblVal >= DSDP_INFINITY) { retcode = DSDP_RETCODE_FAILED; } break;
        case DBL_PARAM_INFEAS_THRESH: if (dblVal < 0.0 || dblVal >= DSDP_INFINITY) { retcode = DSDP_RETCODE_FAILED; } break;
        case DBL_PARAM_ABS_OPTTOL   : if (dblVal < 0.0 || dblVal >= DSDP_INFINITY) { retcode = DSDP_RETCODE_FAILED; } break;
        case DBL_PARAM_REL_OPTTOL   : if (dblVal < 0.0 || dblVal >= DSDP_INFINITY) { retcode = DSDP_RETCODE_FAILED; } break;
        case DBL_PARAM_ABS_FEASTOL  : if (dblVal < 0.0 || dblVal >= DSDP_INFINITY) { retcode = DSDP_RETCODE_FAILED; } break;
        case DBL_PARAM_REL_FEASTOL  : if (dblVal < 0.0 || dblVal >= DSDP_INFINITY) { retcode = DSDP_RETCODE_FAILED; } break;
        case DBL_PARAM_PRLX_PENTALTY: if (dblVal < 0.0 || dblVal >= DSDP_INFINITY) { retcode = DSDP_RETCODE_FAILED; } break;
        case DBL_PARAM_BOUND_X      : if (dblVal < 0.0 || dblVal >= DSDP_INFINITY) { retcode = DSDP_RETCODE_FAILED; } break;
        case DBL_PARAM_TIMELIMIT    : if (dblVal < 0.0 || dblVal >= DSDP_INFINITY) { retcode = DSDP_RETCODE_FAILED; } break;
        default                     : printf("| Invalid parameter code %d \n", pName); retcode = DSDP_RETCODE_FAILED; break;
    }

    return retcode;
}

static DSDP_INT checkIntParam( DSDP_INT pName, DSDP_INT intVal ) {

    DSDP_INT retcode = DSDP_RETCODE_OK;
    
    switch (pName) {
        case INT_PARAM_ACORRECTOR: if (intVal <= 0 || intVal > 20) { retcode = DSDP_RETCODE_FAILED; } break;
        case INT_PARAM_BCORRECTOR: if (intVal < 0 || intVal > 20) { retcode = DSDP_RETCODE_FAILED; } break;
        case INT_PARAM_CG_REUSE  : if (intVal < 0 || intVal > 5) { retcode = DSDP_RETCODE_FAILED; } break;
        case INT_PARAM_AMAXITER  : if (intVal <= 0 || intVal > DSDP_INFINITY) { retcode = DSDP_RETCODE_FAILED; } break;
        case INT_PARAM_BMAXITER  : if (intVal <= 0 || intVal > DSDP_INFINITY) { retcode = DSDP_RETCODE_FAILED; } break;
        case INT_PARAM_GOLDSEARCH: if (intVal != 0 && intVal != 1) { retcode = DSDP_RETCODE_FAILED; } break;
        default: printf("Invalid parameter code %d \n", pName); retcode = DSDP_RETCODE_FAILED; break;
    }
    
    return retcode;
}

extern void setDblParam( dsdpparam *param, DSDP_INT pName, double dblVal ) {
    
    if (pName > NUM_DBL_PARAM || pName < 0) {
        printf("| Invalid Parameter Code %d \n", pName); return;
    }
    if (checkDblParam(pName, dblVal) == DSDP_RETCODE_OK) {
        param->dblParams[pName] = dblVal;
    } else {
        printf("| Invalid parameter range. Parameter unchanged. \n"); return;
    }
    param->dblParams[pName] = dblVal;
}

extern void getDblParam( dsdpparam *param, DSDP_INT pName, double *dblVal ) {
    
    if (pName > NUM_DBL_PARAM || pName < 0) {
        printf("| Invalid Parameter Code %d \n", pName); return;
    }
    *dblVal = param->dblParams[pName];
}

extern void setIntParam( dsdpparam *param, DSDP_INT pName, DSDP_INT intVal ) {
    
    if (pName > NUM_INT_PARAM || pName < 0) {
        printf("| Invalid Parameter Code %d \n", pName); return;
    }
    if (checkIntParam(pName, intVal) == DSDP_RETCODE_OK) {
        param->intParams[pName] = intVal;
    } else {
        printf("| Invalid parameter range. Parameter unchanged. \n");
    }
    return;
}

extern void getIntParam( dsdpparam *param, DSDP_INT pName, DSDP_INT *intVal ) {
    
    if (pName > NUM_INT_PARAM || pName < 0) {
        printf("Invalid Parameter Code %d \n", pName); return;
    }
    *intVal = param->intParams[pName];
}

extern void printParams( dsdpparam *param ) {
    
    printf("| Parameter Summary \n");
    for (DSDP_INT pName = 0; pName < NUM_DBL_PARAM; ++pName) {
        printDblParam(param->dblParams, pName);
    }
    for (DSDP_INT pName = 0; pName < NUM_INT_PARAM; ++pName ) {
        printIntParam(param->intParams, pName);
    }
}

/* DSDP Summary parameters printer */
extern void DSDPParamPrint( dsdpparam *param ) {
    // Print DSDP Parameters at the beginning. Invoked before solution
    printf("| Parameter Summary: \n");
    showBeautifulDashlines();
    printDblParam(param->dblParams, DBL_PARAM_RHON);
    printIntParam(param->intParams, INT_PARAM_GOLDSEARCH);
    printDblParam(param->dblParams, DBL_PARAM_PRLX_PENTALTY);
    printDblParam(param->dblParams, DBL_PARAM_TIMELIMIT);
}
