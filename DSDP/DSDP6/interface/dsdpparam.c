#include <stdio.h>
#include "dsdphsd.h"
#include "dsdpparam.h"
#include "dsdplog.h"

/* Implement DSDP parameter interface */
static void printDblParam( const double *dblParams, DSDP_INT pName ) {
    
    printf("| ");
    
    if (pName > NUM_DBL_PARAM || pName < 0) {
        printf("Invalid Double Parameter Code %d \n", pName);
        return;
    }
    
    switch (pName) {
            
        case DBL_PARAM_ASIGMA:
            printf("Double Parameter A sigma (0.0, 1.0): ");
            break;
        case DBL_PARAM_BSIGMA:
            printf("Double Parameter B sigma (0.0, 1.0): ");
            break;
        case DBL_PARAM_RHO:
            printf("Double Parameter rho [1.0, 10.0]: ");
            break;
        case DBL_PARAM_INIT_POBJ:
            printf("Double Parameter initial pObj (-inf, inf): ");
            break;
        case DBL_PARAM_INIT_BETA:
            printf("Double Parameter initial beta [1.0, inf): ");
            break;
        case DBL_PARAM_INIT_MU:
            printf("Double Parameter initial mu (0.0, inf): ");
            break;
        case DBL_PARAM_INIT_TAU:
            printf("Double Parameter initial tau (0.0, inf): ");
            break;
        case DBL_PARAM_INIT_KAPPA:
            printf("Double Parameter initial kappa (0.0, inf): ");
            break;
        case DBL_PARAM_AALPHA:
            printf("Double Parameter A alpha (0.0, 1.0): ");
            break;
        case DBL_PARAM_NRM_THRESH:
            printf("Double Parameter norm threshold (0.0, inf): ");
            break;
        case DBL_PARAM_INFEAS_THRESH:
            printf("Double Parameter infeasibility threshold (0.0, inf): ");
            break;
        case DBL_PARAM_ABS_OPTTOL:
            printf("Double Parameter absolute optimality tolerance (0.0, inf): ");
            break;
        case DBL_PARAM_REL_OPTTOL:
            printf("Double Parameter relative optimality tolerance (0.0, inf): ");
            break;
        case DBL_PARAM_ABS_FEASTOL:
            printf("Double Parameter absolute feasibility tolerance (0.0, inf): ");
            break;
        case DBL_PARAM_REL_FEASTOL:
            printf("Double Parameter relative feasibility tolerance (0.0, inf): ");
            break;
        case DBL_PARAM_PRLX_PENTALTY:
            printf("Double Parameter primal relaxation penalty (0.0, inf): ");
            break;
        default:
            printf("Invalid Parameter Code %d \n", pName);
            return;
            break;
            
    }
    
    printf("%g \n", dblParams[pName]);
    return;
}

static void printIntParam( const DSDP_INT *intParams, DSDP_INT pName ) {
    
    
    printf("| ");
    if (pName > NUM_INT_PARAM || pName < 0) {
        printf("Invalid Integer Parameter Code %d \n", pName);
        return;
    }
    
    switch (pName) {
        case INT_PARAM_ACORRECTOR:
            printf("Integer Parameter A corrector step (0, 20): ");
            break;
        case INT_PARAM_BCORRECTOR:
            printf("Integer Parameter B corrector step (0, 20): ");
            break;
        case INT_PARAM_INITMETHOD:
            if (intParams[pName] == INIT_METHOD_FRO) {
                printf("Integer Parameter initialization method: Frobenius Norm \n");
            } else {
                printf("Invalid initialization method \n");
            }
            return; break;
        case INT_PARAM_AATTEMPT:
            if (intParams[pName] == AATEMPT_AGGRESSIVE) {
                printf("Integer Parameter A attempt mode: Aggressive \n");
            } else if (intParams[pName] == AATEMPT_MILD) {
                printf("Integer Parameter A attempt mode: Mild \n");
            } else if (intParams[pName] == AATEMT_CONSERVATIVE) {
                printf("Integer Parameter A attempt mode: Conservative \n");
            } else {
                printf("Invalid A attempt parameter \n");
            }
            return; break;
        case INT_PARAM_CG_REUSE:
            printf("Integer Parameter CG pre-conditioner reuse [0, 5]: ");
            break;
        case INT_PARAM_PRESOLVE:
            if (intParams[pName] == PRESOLVE_AGGRESSIVE) {
                printf("Integer Parameter presolve mode: Aggressive \n");
            } else if (intParams[pName] == PRESOLVE_CONSERVATIVE) {
                printf("Integer Parameter presolve mode: Conservative \n");
            } else {
                printf("Invalid presolve mode \n");
            }
            return; break;
        case INT_PARAM_AMAXITER:
            printf("Integer Parameter A maximum iteration: ");
            break;
        case INT_PARAM_BMAXITER:
            printf("Integer Parameter B maximum iteration: ");
            break;
        case INT_PARAM_PRELAX:
            printf("Integer Parameter Primal Relaxation: ");
            break;
        default:
            printf("Invalid Parameter Code %d \n", pName);
            return; break;
    }
    
    printf("%d \n", intParams[pName]);
    return;
}

static DSDP_INT checkDblParam( DSDP_INT pName, double dblVal ) {
    // Check double parameter
    DSDP_INT retcode = DSDP_RETCODE_OK;
    // Currently do nothing
    return retcode;
}

static DSDP_INT checkIntParam( DSDP_INT pName, DSDP_INT intVal ) {
    // Check integer parameter
    DSDP_INT retcode = DSDP_RETCODE_OK;
    // Currently do nothing
    return retcode;
}

extern DSDP_INT setDblParam( hsdParam *param, DSDP_INT pName, double dblVal ) {
    // Set double parameter
    DSDP_INT retcode = DSDP_RETCODE_OK;
    
    if (pName > NUM_DBL_PARAM || pName < 0) {
        printf("Invalid Parameter Code %d \n", pName);
        return DSDP_RETCODE_FAILED;
    }
    
    // printf("| Before: \n");
    // printDblParam(param->dblParams, pName);
    retcode = checkDblParam(pName, dblVal);
    
    if (retcode == DSDP_RETCODE_OK) {
        param->intParams[pName] = dblVal;
    } else {
        printf("| Invalid parameter range. Parameter unchanged. \n");
    }
    
    param->dblParams[pName] = dblVal;
    // printf("| After: \n");
    // printDblParam(param->dblParams, pName);
    
    return retcode;
}

extern DSDP_INT getDblParam( hsdParam *param, DSDP_INT pName, double *dblVal ) {
    
    // Get double parameter
    DSDP_INT retcode = DSDP_RETCODE_OK;
    
    if (pName > NUM_DBL_PARAM || pName < 0) {
        printf("Invalid Parameter Code %d \n", pName);
        return DSDP_RETCODE_FAILED;
    }
    
    *dblVal = param->dblParams[pName];
    return retcode;
}

extern DSDP_INT setIntParam( hsdParam *param, DSDP_INT pName, DSDP_INT intVal ) {
    // Set int parameter
    DSDP_INT retcode = DSDP_RETCODE_OK;
    
    if (pName > NUM_INT_PARAM || pName < 0) {
        printf("Invalid Parameter Code %d \n", pName);
        return DSDP_RETCODE_FAILED;
    }
    
    printf("| Before: \n");
    printIntParam(param->intParams, pName);
    retcode = checkIntParam(pName, intVal);
    
    if (retcode == DSDP_RETCODE_OK) {
        param->intParams[pName] = intVal;
    } else {
        printf("| Invalid parameter range. Parameter unchanged. \n");
    }
    
    printf("| After: \n");
    printIntParam(param->intParams, pName);
    return retcode;
}

extern DSDP_INT getIntParam( hsdParam *param, DSDP_INT pName, DSDP_INT *intVal ) {
    // Get int parameter
    DSDP_INT retcode = DSDP_RETCODE_OK;
    
    if (pName > NUM_INT_PARAM || pName < 0) {
        printf("Invalid Parameter Code %d \n", pName);
        return DSDP_RETCODE_FAILED;
    }
    
    *intVal = param->intParams[pName];
    return retcode;
}

extern void printParams( hsdParam *param ) {
    // Print all the parameters
    printf("| Parameter Summary \n");
    for (DSDP_INT pName = 0; pName < NUM_DBL_PARAM; ++pName) {
        printDblParam(param->dblParams, pName);
    }
    
    for (DSDP_INT pName = 0; pName < NUM_INT_PARAM; ++pName ) {
        printIntParam(param->intParams, pName);
    }
    
    return;
}

/* DSDP Summary parameters printer */
extern void DSDPParamPrint( hsdParam *param ) {
    // Print DSDP Parameters at the beginning. Invoked before solution
    
    dsdpshowdash();
    printf("| Parameter Summary: \n");
    dsdpshowdash();
    printDblParam(param->dblParams, DBL_PARAM_RHO);
    printIntParam(param->intParams, INT_PARAM_PRESOLVE);
    
    DSDP_INT prelax;
    getIntParam(param, INT_PARAM_PRELAX, &prelax);
    printIntParam(param->intParams, INT_PARAM_PRELAX);
    
    if (prelax) {
        printDblParam(param->dblParams, DBL_PARAM_PRLX_PENTALTY);
    }
    
    printIntParam(param->intParams, INT_PARAM_AATTEMPT);
    
    return;
}
