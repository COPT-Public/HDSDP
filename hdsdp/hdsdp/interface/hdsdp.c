/** @file hdsdp.c
 *  @brief A homogeneous dual-scaling interior point solver
 */

#ifdef HEADERPATH
#include "interface/def_hdsdp.h"
#include "interface/hdsdp.h"
#include "interface/hdsdp_conic.h"
#include "interface/hdsdp_utils.h"
#include "interface/hdsdp_user_data.h"
#include "interface/hdsdp_schur.h"
#include "interface/hdsdp_algo.h"
#include "linalg/dense_opts.h"
#else
#include "def_hdsdp.h"
#include "hdsdp.h"
#include "hdsdp_conic.h"
#include "hdsdp_utils.h"
#include "def_hdsdp_user_data.h"
#include "hdsdp_user_data.h"
#include "hdsdp_schur.h"
#include "hdsdp_algo.h"
#include "dense_opts.h"
#include "vec_opts.h"
#endif

#include <math.h>

static void HDSDPIGetStatistics( hdsdp *HSolver ) {
    /* Collect statistics */
    hdsdp_printf("  Collecting statistcs \n");
    
    /* Get sum of conic dimensions. Used to convert the identity matrix to Frobenius norm */
    int sumConeDims = 0;
    int maxConeDim = 0;
    for ( int iCone = 0; iCone < HSolver->nCones; ++iCone ) {
        if ( HSolver->HCones[iCone]->cone == HDSDP_CONETYPE_DENSE_SDP ||
             HSolver->HCones[iCone]->cone == HDSDP_CONETYPE_SPARSE_SDP ) {
            maxConeDim = HDSDP_MAX(maxConeDim, HConeGetDim(HSolver->HCones[iCone]));
        }
        sumConeDims += HConeGetDim(HSolver->HCones[iCone]);
    }
    
    set_int_feature(HSolver, INT_FEATURE_N_MAXCONEDIM, maxConeDim);
    set_int_feature(HSolver, INT_FEATURE_N_SUMCONEDIMS, sumConeDims);
    set_int_feature(HSolver, INT_FEATURE_N_ROWS, HSolver->nRows);
    set_int_feature(HSolver, INT_FEATURE_N_CONES, HSolver->nCones);
    
    HSolver->dAllConeDims = (double) sumConeDims + 2 * HSolver->nRows;
    
    /* Get norms */
    double objOneNorm = 0.0;
    double objFroNorm = 0.0;
    double dataOneNorm = 0.0;
    double dataFroNorm = 0.0;
    double dFroNormTmp = 0.0;
    double rhsOneNorm = 0.0;
    double rhsFroNorm = 0.0;
    double rhsInfNorm = 0.0;
    
    for ( int iCone = 0; iCone < HSolver->nCones; ++iCone ) {
        objOneNorm += HConeGetObjNorm(HSolver->HCones[iCone], ABS_NORM);
        dataOneNorm += HConeGetCoeffNorm(HSolver->HCones[iCone], ABS_NORM);
        dFroNormTmp = HConeGetObjNorm(HSolver->HCones[iCone], FRO_NORM);
        objFroNorm += dFroNormTmp * dFroNormTmp;
        dFroNormTmp = HConeGetCoeffNorm(HSolver->HCones[iCone], FRO_NORM);
        dataFroNorm += dFroNormTmp * dFroNormTmp;
    }
    
    objFroNorm = sqrt(objFroNorm);
    dataFroNorm = sqrt(dataFroNorm);
    
    set_dbl_feature(HSolver, DBL_FEATURE_OBJFRONORM, objFroNorm);
    set_dbl_feature(HSolver, DBL_FEATURE_DATAFRONORM, dataFroNorm);
    set_dbl_feature(HSolver, DBL_FEATURE_OBJONENORM, objOneNorm);
    set_dbl_feature(HSolver, DBL_FEATURE_DATAONENORM, dataOneNorm);
    
    for ( int iRow = 0; iRow < HSolver->nRows; ++iRow ) {
        rhsOneNorm += fabs(HSolver->rowRHS[iRow]);
        rhsFroNorm += HSolver->rowRHS[iRow] * HSolver->rowRHS[iRow];
        rhsInfNorm = HDSDP_MAX(rhsInfNorm, fabs(HSolver->rowRHS[iRow]));
    }
    
    rhsFroNorm = sqrt(rhsFroNorm);
    set_dbl_feature(HSolver, DBL_FEATURE_RHSONENORM, rhsOneNorm);
    set_dbl_feature(HSolver, DBL_FEATURE_RHSFRONORM, rhsFroNorm);
    set_dbl_feature(HSolver, DBL_FEATURE_RHSINFNORM, rhsInfNorm);
    
    /* Is objective NULL ? */
    if ( objFroNorm == 0.0 ) {
        set_int_feature(HSolver, INT_FEATURE_I_NULLOBJ, 1);
    }
    
    /* Are there many cones ? */
    if ( HSolver->nCones >= 100 ) {
        set_int_feature(HSolver, INT_FEATURE_I_MANYCONES, 1);
    }
    
    /* How many sparse cones ? */
    int coneStats[7] = {0};
    for ( int iCone = 0; iCone < HSolver->nCones; ++iCone ) {
        coneStats[(int) HSolver->HCones[iCone]->cone] += 1;
    }
    
    set_int_feature(HSolver, INT_FEATURE_N_DSSDPCONES, coneStats[HDSDP_CONETYPE_DENSE_SDP]);
    set_int_feature(HSolver, INT_FEATURE_N_SPSDPCONES, coneStats[HDSDP_CONETYPE_SPARSE_SDP]);
    set_int_feature(HSolver, INT_FEATURE_N_LPCONES, coneStats[HDSDP_CONETYPE_LP]);
    
    return;
}

static void HDSDPIPrintStatistics( hdsdp *HSolver ) {
    
    hdsdp_printf("\nStatistics \n");
    
    print_int_feature(HSolver, INT_FEATURE_N_ROWS, "Number of rows");
    print_int_feature(HSolver, INT_FEATURE_N_CONES, "Number of cones");
    print_int_feature(HSolver, INT_FEATURE_N_SPSDPCONES, "Number of sparse SDP cones");
    print_int_feature(HSolver, INT_FEATURE_N_DSSDPCONES, "Number of dense SDP cones");
    print_int_feature(HSolver, INT_FEATURE_N_LPCONES, "Number of LP cones");
    print_int_feature(HSolver, INT_FEATURE_N_SUMCONEDIMS, "Cone dimensions");
    print_int_feature(HSolver, INT_FEATURE_N_MAXCONEDIM, "Max matrix cone dimension");
    print_dbl_feature(HSolver, DBL_FEATURE_OBJONENORM, "Norm of objective");
    print_dbl_feature(HSolver, DBL_FEATURE_DATAONENORM, "Norm of SDP data");
    print_dbl_feature(HSolver, DBL_FEATURE_RHSONENORM, "Norm of RHS");
    
    return;
}

static void HDSDPIAdjustConeParams( hdsdp *HSolver ) {
    
    /* Detect features when there is one cone */
    int isHit = 0;
    int isOneCone = 1;
    
    if ( get_int_feature(HSolver, INT_FEATURE_N_DSSDPCONES) + \
         get_int_feature(HSolver, INT_FEATURE_N_SPSDPCONES) > 1 ) {
        isOneCone = 0;
    }
    
    if ( isOneCone ) {
        HConeDetectFeature(HSolver->HCones[0], HSolver->rowRHS, HSolver->HIntFeatures, HSolver->HDblFeatures);
    }
    
    if ( get_int_feature(HSolver, INT_FEATURE_N_LPCONES) &&
         get_int_feature(HSolver, INT_FEATURE_N_CONES) < 10 ) {
        HConeDetectFeature(HSolver->HCones[HSolver->nCones - 1], HSolver->rowRHS, HSolver->HIntFeatures, HSolver->HDblFeatures);
    }
    
    int isImpliedTrace = get_int_feature(HSolver, INT_FEATURE_I_IMPTRACE);
    int isNoPrimalInterior = get_int_feature(HSolver, INT_FEATURE_I_NOPINTERIOR);
    int isNoDualInterior = get_int_feature(HSolver, INT_FEATURE_I_NODINTERIOR);
    int isExtremelyDense = get_int_feature(HSolver, INT_FEATURE_I_VERYDENSE);
    int isImpliedDual = get_int_feature(HSolver, INT_FEATURE_I_IMPYBOUND);
    int isNoObj = get_int_feature(HSolver, INT_FEATURE_I_NULLOBJ);
    
    HDSDP_ZERO(HSolver->modelFeatures, char, 200);
    
    if ( (isImpliedTrace + isNoPrimalInterior + isNoDualInterior + isExtremelyDense + isNoObj) || isImpliedDual ) {
        isHit = 1;
        strcat(HSolver->modelFeatures, "This is a ");
    }
    
    if ( isExtremelyDense ) {
        set_int_param(HSolver, INT_PARAM_CORRECTORA, 4);
        set_dbl_param(HSolver, DBL_PARAM_DUALSTART, 1.0);
        set_dbl_param(HSolver, DBL_PARAM_DUALBOX_UP, 1e+04);
        set_dbl_param(HSolver, DBL_PARAM_DUALBOX_LOW, -1e+04);
        strcat(HSolver->modelFeatures, "dense ");
    }
    
    if ( isImpliedTrace ) {
        set_dbl_param(HSolver, DBL_PARAM_DUALSTART, 1e+03);
        set_dbl_param(HSolver, DBL_PARAM_TRXESTIMATE, get_dbl_feature(HSolver, DBL_FEATURE_IMPTRACEX));
        set_dbl_param(HSolver, DBL_PARAM_POBJSTART, 1e+08);
        set_dbl_param(HSolver, DBL_PARAM_POTRHOVAL, 5.0);
        set_dbl_param(HSolver, DBL_PARAM_DUALBOX_UP, 1e+06);
        set_dbl_param(HSolver, DBL_PARAM_DUALBOX_LOW, -1e+06);
        strcat(HSolver->modelFeatures, "trace-implied ");
    }
    
    if ( isNoPrimalInterior ) {
        set_dbl_param(HSolver, DBL_PARAM_DUALBOX_UP, 1e+04);
        set_dbl_param(HSolver, DBL_PARAM_DUALBOX_LOW, -1e+04);
        set_dbl_param(HSolver, DBL_PARAM_DUALSTART, 1e+03);
        set_dbl_param(HSolver, DBL_PARAM_PRECORDACC, 1e-07);
        strcat(HSolver->modelFeatures, "no-primal interior ");
    }
    
    if ( isImpliedDual ) {
        
        double dBoxUpper = get_dbl_param(HSolver, DBL_PARAM_DUALBOX_UP);
        double dBoxLower = get_dbl_param(HSolver, DBL_PARAM_DUALBOX_LOW);
        int isUpper = 0;
        int isLower = 0;
        
        if ( get_dbl_feature(HSolver, DBL_FEATURE_IMPYBOUNDUP) ) {
            isUpper = 1;
            dBoxUpper = HDSDP_MIN(dBoxUpper, get_dbl_feature(HSolver, DBL_FEATURE_IMPYBOUNDUP));
            set_dbl_param(HSolver, DBL_PARAM_DUALBOX_UP, dBoxUpper);
        }
        if ( get_dbl_feature(HSolver, DBL_FEATURE_IMPYBOUNDLOW) ) {
            isLower = 1;
            dBoxLower = HDSDP_MAX(dBoxLower, get_dbl_feature(HSolver, DBL_FEATURE_IMPYBOUNDLOW));
            set_dbl_param(HSolver, DBL_PARAM_DUALBOX_LOW, dBoxLower);
        }
        
        if ( isUpper && isLower ) {
            set_dbl_param(HSolver, DBL_PARAM_DUALSTART, 1e+02);
            set_dbl_param(HSolver, DBL_PARAM_POBJSTART, 1e+05);
        } else {
            set_dbl_param(HSolver, DBL_PARAM_DUALSTART, 1e+05);
            set_dbl_param(HSolver, DBL_PARAM_POBJSTART, 1e+10);
            set_int_param(HSolver, INT_PARAM_CORRECTORA, 12);
            set_int_param(HSolver, INT_PARAM_CORRECTORB, 12);
        }
        
        set_dbl_param(HSolver, DBL_PARAM_ABSOPTTOL, 1e-01);
        set_dbl_param(HSolver, DBL_PARAM_RELOPTTOL, 1e-04);
        set_dbl_param(HSolver, DBL_PARAM_PRECORDACC, 1e-05);
        
        strcat(HSolver->modelFeatures, "dual-bounded ");
    }
    
    if ( isNoDualInterior ) {
        set_dbl_param(HSolver, DBL_PARAM_DUALBOX_UP, 1.0);
        set_dbl_param(HSolver, DBL_PARAM_DUALBOX_LOW, -1.0);
        
        if ( HSolver->dAllConeDims > 100000 ) {
            set_dbl_param(HSolver, DBL_PARAM_DUALSTART, 1e+00);
            set_dbl_param(HSolver, DBL_PARAM_ABSFEASTOL, 1e-04);
            set_dbl_param(HSolver, DBL_PARAM_RELFEASTOL, 1e-05);
        } else {
            set_dbl_param(HSolver, DBL_PARAM_DUALBOX_UP, 1e+01);
            set_dbl_param(HSolver, DBL_PARAM_DUALBOX_LOW, -1e+01);
            set_dbl_param(HSolver, DBL_PARAM_ABSFEASTOL, 1e-05);
            set_dbl_param(HSolver, DBL_PARAM_RELFEASTOL, 1e-07);
        }
        
        set_dbl_param(HSolver, DBL_PARAM_PRECORDACC, 1e-05);
        
        strcat(HSolver->modelFeatures, "no-dual interior ");
    }
    
    if ( isNoObj ) {
        set_dbl_param(HSolver, DBL_PARAM_DUALSTART, 1.0);
        set_dbl_param(HSolver, DBL_PARAM_DUALBOX_UP, 1.0);
        set_dbl_param(HSolver, DBL_PARAM_DUALBOX_LOW, -1.0);
        strcat(HSolver->modelFeatures, "no objective ");
    }
    
    if ( isHit ) {
        strcat(HSolver->modelFeatures, "SDP problem\n");
    }
    
    return;
}

static void HDSDPIAdjustParams( hdsdp *HSolver ) {
    
    hdsdp_printf("  Making adjustments\n");
    
    /* Scale cone objective and RHS */
    double objOneNorm = get_dbl_feature(HSolver, DBL_FEATURE_OBJONENORM);
    double rhsInfNorm = get_dbl_feature(HSolver, DBL_FEATURE_RHSINFNORM);
    
    double objScal = 0.0;
    double rhsScal = 0.0;
    
    objScal = 1.0;
    
    if ( objOneNorm > 1e+10 ) {
        objScal = 1e-08;
    } else if ( objOneNorm > 1e+08 ) {
        objScal = 1e-06;
    } else if ( objOneNorm > 1e+05 ) {
        objScal = 1e-05;
    }
    
    if ( rhsInfNorm > 1e+10 ) {
        rhsScal = 1e-08;
    } else if ( rhsInfNorm > 1e+08 ) {
        rhsScal = 1e-06;
    } else {
        rhsScal = 1.0;
    }
    
    set_dbl_feature(HSolver, DBL_FEATURE_OBJSCALING, objScal);
    set_dbl_feature(HSolver, DBL_FEATURE_RHSSCALING, rhsScal);
    
    for ( int iCone = 0; iCone < HSolver->nCones; ++iCone ) {
        HConeScalByConstant(HSolver->HCones[iCone], objScal);
    }
    
    for ( int iRow = 0; iRow < HSolver->nRows; ++iRow ) {
        HSolver->rowRHS[iRow] = rhsScal * HSolver->rowRHS[iRow];
    }
    
    hdsdp_printf("    Scale cone objective by %5.1e \n", objScal);
    hdsdp_printf("    Scale rhs by %5.1e \n", rhsScal);
    
    int nMaxThreads = HUtilGetGlobalMKLThreads();
    int nTargetThreads = get_int_param(HSolver, INT_PARAM_THREADS);

    if ( nMaxThreads > nTargetThreads ) {
        HUtilSetGlobalMKLThreads(nTargetThreads);
    } else {
        HUtilSetGlobalMKLThreads(nMaxThreads);
        hdsdp_printf("    Hardware has %d thread(s)\n", nMaxThreads);
        set_int_param(HSolver, INT_PARAM_THREADS, nMaxThreads);
    }
    
    /* Determine correctors */
    int nCorrA = 0;
    int nCorrB = 0;
    
    nCorrA = (int) (HSolver->nRows - 2) / get_int_feature(HSolver, INT_FEATURE_N_MAXCONEDIM);
    
    if ( get_int_feature(HSolver, INT_FEATURE_N_SUMCONEDIMS) < 100 && nCorrA == 0 ) {
        nCorrA = 1;
    }
    
    if ( nCorrA >= 1 ) {
        nCorrA += 1;
    }
    
    nCorrA = nCorrA * nCorrA;
    
    if ( HSolver->nRows < 2000 && nCorrA > 10 ) {
        nCorrA = 10;
    }
    
    nCorrB = nCorrA;
    
    if ( get_int_feature(HSolver, INT_FEATURE_N_MAXCONEDIM) >= 5 * HSolver->nRows ) {
        nCorrB = 0;
        nCorrA = 2;
    } else if ( get_int_feature(HSolver, INT_FEATURE_N_MAXCONEDIM) >= HSolver->nRows ) {
        nCorrB = HDSDP_MIN(nCorrB, 2);
        nCorrA = 4;
    } else {
        nCorrA = 6;
    }
    
    if ( HSolver->nRows > 20 * get_int_feature(HSolver, INT_FEATURE_N_MAXCONEDIM) ) {
        nCorrB = HDSDP_MAX(nCorrB, 12);
        nCorrA = 12;
    } else if ( HSolver->nRows > 5 * get_int_feature(HSolver, INT_FEATURE_N_MAXCONEDIM) ) {
        nCorrB = HDSDP_MAX(nCorrB, 10);
        nCorrA = 10;
    } else if ( HSolver->nRows > 2 * get_int_feature(HSolver, INT_FEATURE_N_MAXCONEDIM) ) {
        nCorrB = HDSDP_MAX(nCorrB, 8);
        nCorrA = 8;
    }
    
    nCorrB = HDSDP_MIN(nCorrB, 12);
    nCorrA = HDSDP_MAX(nCorrA, 2);

    set_int_param(HSolver, INT_PARAM_CORRECTORA, nCorrA);
    set_int_param(HSolver, INT_PARAM_CORRECTORB, nCorrB);
    
    if ( get_int_feature(HSolver, INT_FEATURE_I_MANYCONES) ) {
        set_int_param(HSolver, INT_PARAM_CORRECTORA, 6);
        set_int_param(HSolver, INT_PARAM_CORRECTORB, 0);
        set_dbl_param(HSolver, DBL_PARAM_DUALSTART, 1.0);
        set_dbl_param(HSolver, DBL_PARAM_POBJSTART, 1e+10);
    }
    
    HDSDPIAdjustConeParams(HSolver);
    
    return;
}

static void HDSDPIGetDefaultParams( hdsdp *HSolver ) {
    
    set_int_param(HSolver, INT_PARAM_MAXITER, 500);
    set_int_param(HSolver, INT_PARAM_CORRECTORA, 12);
    set_int_param(HSolver, INT_PARAM_CORRECTORB, 12);
    set_int_param(HSolver, INT_PARAM_THREADS, 12);
    
    set_dbl_param(HSolver, DBL_PARAM_ABSOPTTOL, 1e-08);
    set_dbl_param(HSolver, DBL_PARAM_ABSFEASTOL, 1e-08);
    set_dbl_param(HSolver, DBL_PARAM_RELOPTTOL, 1e-10);
    set_dbl_param(HSolver, DBL_PARAM_RELFEASTOL, 1e-10);
    set_dbl_param(HSolver, DBL_PARAM_TIMELIMIT, 3600.0);
    set_dbl_param(HSolver, DBL_PARAM_POTRHOVAL, 4.0);
    set_dbl_param(HSolver, DBL_PARAM_HSDGAMMA, 0.5);
    set_dbl_param(HSolver, DBL_PARAM_DUALBOX_UP, 1e+07);
    set_dbl_param(HSolver, DBL_PARAM_DUALBOX_LOW, -1e+07);
    set_dbl_param(HSolver, DBL_PARAM_BARMUSTART, 1e+05);
    set_dbl_param(HSolver, DBL_PARAM_POBJSTART, 1e+10);
    set_dbl_param(HSolver, DBL_PARAM_DUALSTART, 1e+05);
    set_dbl_param(HSolver, DBL_PARAM_TRXESTIMATE, 1e+08);
    set_dbl_param(HSolver, DBL_PARAM_PRECORDACC, 1e-08);
    
    return;
}

static void HDSDPIPrintParams( hdsdp *HSolver ) {
    
    hdsdp_printf("\nParameters\n");
    print_int_param(HSolver, INT_PARAM_MAXITER, "Maximum iteration");
    print_int_param(HSolver, INT_PARAM_CORRECTORA, "Infeasible corrector");
    print_int_param(HSolver, INT_PARAM_CORRECTORB, "Feasible corrector");
    print_int_param(HSolver, INT_PARAM_THREADS, "Threads");
    
    print_dbl_param(HSolver, DBL_PARAM_ABSOPTTOL, "Abs optimality");
    print_dbl_param(HSolver, DBL_PARAM_RELOPTTOL, "Rel optimality");
    print_dbl_param(HSolver, DBL_PARAM_ABSFEASTOL, "Abs feasibility");
    print_dbl_param(HSolver, DBL_PARAM_RELFEASTOL, "Rel feasibility");
    print_dbl_param(HSolver, DBL_PARAM_TIMELIMIT, "Time limit");
    print_dbl_param(HSolver, DBL_PARAM_POTRHOVAL, "Potential param");
    print_dbl_param(HSolver, DBL_PARAM_HSDGAMMA, "Mu reduction rate");
    print_dbl_param(HSolver, DBL_PARAM_DUALBOX_LOW, "Dual box low");
    print_dbl_param(HSolver, DBL_PARAM_DUALBOX_UP, "Dual box up");
    print_dbl_param(HSolver, DBL_PARAM_BARMUSTART, "Starting barrier mu");
    print_dbl_param(HSolver, DBL_PARAM_DUALSTART, "Starting residual");
    print_dbl_param(HSolver, DBL_PARAM_POBJSTART, "Starting primal");
    
    hdsdp_printf("\n");
    
    return;
}

static void HDSDPIPrintSolutionStats( hdsdp *HSolver ) {
    
    if ( HSolver->HStatus == HDSDP_UNKNOWN ) {
        hdsdp_printf("\nSDP Status: %s \n", "Unknown");
    } else if ( HSolver->HStatus == HDSDP_DUAL_FEASIBLE ) {
        hdsdp_printf("\nSDP Status: %s \n", "Dual feasible");
    } else if ( HSolver->HStatus == HDSDP_DUAL_OPTIMAL ) {
        hdsdp_printf("\nSDP Status: %s \n", "Dual optimal ");
    } else if ( HSolver->HStatus == HDSDP_PRIMAL_DUAL_OPTIMAL ) {
        hdsdp_printf("\nSDP Status: %s \n", "Primal dual optimal");
    } else if ( HSolver->HStatus == HDSDP_MAXITER ) {
        hdsdp_printf("\nSDP Status: %s \n", "Maximum iteration");
    } else if ( HSolver->HStatus == HDSDP_SUSPECT_INFEAS_OR_UNBOUNDED ) {
        hdsdp_printf("\nSDP Status: %s \n", "Suspected infeasible or unbounded");
    } else if ( HSolver->HStatus == HDSDP_INFEAS_OR_UNBOUNDED ) {
        hdsdp_printf("\nSDP Status: %s \n", "Infeasible or unbounded");
    } else if ( HSolver->HStatus == HDSDP_TIMELIMIT ) {
        hdsdp_printf("\nSDP Status: %s \n", "Time limit");
    } else if ( HSolver->HStatus == HDSDP_USER_INTERRUPT ) {
        hdsdp_printf("\nSDP Status: %s \n", "User interrupt");
    } else if ( HSolver->HStatus == HDSDP_INTERNAL_ERROR ) {
        hdsdp_printf("\nSDP Status: %s \n", "Internal error");
    } else if ( HSolver->HStatus == HDSDP_NUMERICAL ) {
        hdsdp_printf("\nSDP Status: %s \n", "Numerical error");
    } else {
        assert( 0 );
    }
    
    if ( HSolver->HStatus != HDSDP_SUSPECT_INFEAS_OR_UNBOUNDED &&
         HSolver->HStatus != HDSDP_INFEAS_OR_UNBOUNDED ) {
        hdsdp_printf("  pObj %+15.10e\n", HSolver->pObjVal);
        hdsdp_printf("  dObj %+15.10e\n", HSolver->dObjVal);
        hdsdp_printf("PD Gap %+15.10e\n", HSolver->pObjVal - HSolver->dObjVal);
        hdsdp_printf("  Time %3.1f seconds\n", HUtilGetTimeStamp() - HSolver->dTimeBegin);
    }
    
    hdsdp_printf("\n");

    return;
}

/* Define HDSDP Solver interface */
extern hdsdp_retcode HDSDPCreate( hdsdp **pHSolver ) {
    
    hdsdp_retcode retcode = HDSDP_RETCODE_OK;
    
    if ( !pHSolver ) {
        retcode = HDSDP_RETCODE_FAILED;
        return retcode;
    }
    
    hdsdp *HSolver = NULL;
    HDSDP_INIT(HSolver, hdsdp, 1);
    HDSDP_MEMCHECK(HSolver);
    
    *pHSolver = HSolver;
    
    HUtilStartCtrlCCheck();
    
exit_cleanup:
    return retcode;
}

extern hdsdp_retcode HDSDPInit( hdsdp *HSolver, int nRows, int nCones ) {
    /* Initialize the dual potential reduction solver */
    hdsdp_retcode retcode = HDSDP_RETCODE_OK;
    
    HSolver->nRows = nRows;
    HSolver->nCones = nCones;
    
    HDSDPIGetDefaultParams(HSolver);
    
    HDSDP_INIT(HSolver->rowRHS, double, nRows);
    HDSDP_MEMCHECK(HSolver->rowRHS);
        
    HDSDP_INIT(HSolver->HCones, hdsdp_cone *, nCones);
    HDSDP_MEMCHECK(HSolver->HCones);
    
    for ( int iCone = 0; iCone < nCones; ++iCone ) {
        HDSDP_CALL(HConeCreate(&HSolver->HCones[iCone]));
    }
    
    HDSDP_CALL(HConeCreate(&HSolver->HBndCone));
    HDSDP_CALL(HKKTCreate(&HSolver->HKKT));
    
    /* Allocate vectors */
    HDSDP_INIT(HSolver->dRowDual, double, nRows);
    HDSDP_MEMCHECK(HSolver->dRowDual);
    
    HDSDP_INIT(HSolver->dRowDualStep, double, nRows);
    HDSDP_MEMCHECK(HSolver->dRowDualStep);
    
    HDSDP_INIT(HSolver->dMinvRowRHS, double, nRows);
    HDSDP_MEMCHECK(HSolver->dMinvRowRHS);
    
    HDSDP_INIT(HSolver->dMinvASinv, double, nRows);
    HDSDP_MEMCHECK(HSolver->dMinvASinv);
    
    HDSDP_INIT(HSolver->dMinvASinvRdSinv, double, nRows);
    HDSDP_MEMCHECK(HSolver->dMinvASinvRdSinv);
    
    HDSDP_INIT(HSolver->dMinvASinvCSinv, double, nRows);
    HDSDP_MEMCHECK(HSolver->dMinvASinvCSinv);
    
    HDSDP_INIT(HSolver->dHAuxiVec1, double, nRows);
    HDSDP_MEMCHECK(HSolver->dHAuxiVec1);
    
    HDSDP_INIT(HSolver->dHAuxiVec2, double, nRows);
    HDSDP_MEMCHECK(HSolver->dHAuxiVec2);
    
    HDSDP_INIT(HSolver->dPInfeasUpper, double, nRows);
    HDSDP_MEMCHECK(HSolver->dPInfeasUpper);
    
    HDSDP_INIT(HSolver->dPInfeasLower, double, nRows);
    HDSDP_MEMCHECK(HSolver->dPInfeasLower);
    
    HDSDP_INIT(HSolver->dAccRowDualMaker, double, nRows);
    HDSDP_MEMCHECK(HSolver->dAccRowDualMaker);
    
    HDSDP_INIT(HSolver->dAccRowDualStepMaker, double, nRows);
    HDSDP_MEMCHECK(HSolver->dAccRowDualStepMaker);
    
    HDSDP_INIT(HSolver->dInaccRowDualMaker, double, nRows);
    HDSDP_MEMCHECK(HSolver->dInaccRowDualMaker);
    
    HDSDP_INIT(HSolver->dInaccRowDualStepMaker, double, nRows);
    HDSDP_MEMCHECK(HSolver->dInaccRowDualStepMaker);
    
    HSolver->dBarrierMu = 1e+10;
    HSolver->comp = HDSDP_INFINITY;
    
    HSolver->dAccBarrierMaker = -1.0;
    HSolver->dInaccBarrierMaker = -1.0;
    
    HSolver->HStatus = HDSDP_UNKNOWN;
    
exit_cleanup:
    return retcode;
}

extern hdsdp_retcode HDSDPSetCone( hdsdp *HSolver, int iCone, void *userCone ) {
    
    hdsdp_retcode retcode = HDSDP_RETCODE_OK;
    
    HDSDP_CALL(HConeSetData(HSolver->HCones[iCone], userCone));
    
exit_cleanup:
    return retcode;
}

extern void HDSDPSetDualObjective( hdsdp *HSolver, double *dObj ) {
    
    HDSDP_MEMCPY(HSolver->rowRHS, dObj, double, HSolver->nRows);
    
    return;
}

extern void HDSDPSetDualStart( hdsdp *HSolver, double *dStart ) {
    
    if ( dStart ) {
        HDSDP_MEMCPY(HSolver->dRowDual, dStart, double, HSolver->nRows);
    }
    
    return;
}

extern void HDSDPSetIntParam( hdsdp *HSolver, int intParam, int intParamVal ) {
    
    if ( intParam >= 3 || intParam < 0 ) {
        return;
    }
    
    HSolver->HIntParams[intParam] = intParamVal;
    return;
}


extern void HDSDPSetDblParam( hdsdp *HSolver, int dblParam, double dblParamVal ) {
    
    if ( dblParam >= 13 || dblParam < 0 ) {
        return;
    }
    HSolver->HDblParams[dblParam] = dblParamVal;
    return;
    
}

extern hdsdp_retcode HDSDPOptimize( hdsdp *HSolver, int dOptOnly ) {
    
    hdsdp_retcode retcode = HDSDP_RETCODE_OK;
    
    /* Start optimization */
    hdsdp_printf("\nHDSDP: software for semi-definite programming \n\n");
    hdsdp_printf("Wenzhi Gao, Dongdong Ge, Yinyu Ye, 2023\n");
    hdsdp_printf("---------------------------------------------\n");
    HSolver->dTimeBegin = HUtilGetTimeStamp();
    
    /* Process conic data */
    hdsdp_printf("Pre-solver starts \n");
    hdsdp_printf("  Processing the cones \n");
    for ( int iCone = 0; iCone < HSolver->nCones; ++iCone ) {
        HDSDP_CALL(HConeProcData(HSolver->HCones[iCone]));
        HDSDP_CALL(HConePresolveData(HSolver->HCones[iCone]));
    }
    
    hdsdp_printf("  Starting KKT analysis \n");
    /* Prepare KKT solver */
    HDSDP_CALL(HKKTInit(HSolver->HKKT, HSolver->nRows, HSolver->nCones, HSolver->HCones));
    
    /* Get statistics */
    HDSDPIGetStatistics(HSolver);
    
    /* Adjust parameters */
    HDSDPIAdjustParams(HSolver);
    
    hdsdp_printf("  Creating range constraint \n");
    user_data udata;
    udata.coneMatBeg = NULL;
    udata.coneMatIdx = NULL;
    udata.cone = HDSDP_CONETYPE_SCALAR_BOUND;
    udata.nConicCol = 0;
    udata.nConicRow = HSolver->nRows;
    double boxRange[2] = {get_dbl_param(HSolver, DBL_PARAM_DUALBOX_LOW),
                          get_dbl_param(HSolver, DBL_PARAM_DUALBOX_UP)};
    udata.coneMatElem = boxRange;
    HConeSetData(HSolver->HBndCone, &udata);
    HDSDP_CALL(HConeProcData(HSolver->HBndCone));
    HDSDP_CALL(HConePresolveData(HSolver->HBndCone));
    
    hdsdp_printf("    Dual box low: %+5.1e\n", boxRange[0]);
    hdsdp_printf("    Dual box up : %+5.1e\n", boxRange[1]);
    
    /* End pre-solver */
    hdsdp_printf("Pre-solver ends. Elapsed time: %3.1f seconds \n", HUtilGetTimeStamp() - HSolver->dTimeBegin);
    
    HDSDPIPrintStatistics(HSolver);
    HDSDPIPrintParams(HSolver);
    
    hdsdp_printf("%s", HSolver->modelFeatures);
    hdsdp_printf("Optimizing over %d thread(s) \n", get_int_param(HSolver, INT_PARAM_THREADS));
    
    /* Invoke solver */
    retcode = HDSDPSolve(HSolver, dOptOnly);
    
    hdsdp_printf("\nOptimization time: %3.1f seconds\n", HUtilGetTimeStamp() - HSolver->dTimeBegin);
    
    if ( HSolver->HStatus != HDSDP_INFEAS_OR_UNBOUNDED &&
         HSolver->HStatus != HDSDP_SUSPECT_INFEAS_OR_UNBOUNDED ) {
        retcode = HDSDPCheckSolution(HSolver, HSolver->dErrs);
    }
    
    hdsdp_printf("DIMACS error metric:\n    %5.2e %5.2e %5.2e %5.2e %5.2e %5.2e \n",
                 HSolver->dErrs[0], HSolver->dErrs[1], HSolver->dErrs[2],
                 HSolver->dErrs[3], HSolver->dErrs[4], HSolver->dErrs[5]);
    
    HDSDPIPrintSolutionStats(HSolver);
    
exit_cleanup:
    return retcode;
}

extern hdsdp_retcode HDSDPGetRowDual( hdsdp *HSolver, double *pObjVal, double *dObjVal, double *dualVal ) {
    
    hdsdp_retcode retcode = HDSDP_RETCODE_OK;
    
    if ( pObjVal ) {
        *pObjVal = HSolver->pObjVal;
    }
    
    if ( dObjVal ) {
        *dObjVal = HSolver->dObjVal;
    }
    
    if ( dualVal ) {
        HDSDP_MEMCPY(dualVal, HSolver->dRowDual, double, HSolver->nRows);
    }
    
exit_cleanup:
    return retcode;
}

extern void HDSDPGetConeValues( hdsdp *HSolver, int iCone, double *conePrimal, double *coneDual, double *coneAuxi ) {
    
    double dBarrierMaker = HSolver->dAccBarrierMaker;
    double *dRowDualMaker = HSolver->dAccRowDualMaker;
    double *dRowDualStepMaker = HSolver->dAccRowDualStepMaker;
    
#if 0
#include "debug_data.h"
    dBarrierMaker = globalmu;
    dRowDualMaker = globaly;
    dRowDualStepMaker = globaldy;
    
    for ( int iCone = 0; iCone < HSolver->nCones; ++iCone ) {
        HConeReduceResi(HSolver->HCones[iCone], 0.0);
        HConeSetPerturb(HSolver->HCones[iCone], 7.3018532792337075E-8);
    }
#endif
    
    if ( dBarrierMaker <= 0.0 ) {
        dBarrierMaker = HSolver->dInaccBarrierMaker;
        dRowDualMaker = HSolver->dInaccRowDualMaker;
        dRowDualStepMaker = HSolver->dInaccRowDualStepMaker;
    }
        
    if ( conePrimal ) {
        HConeGetPrimal(HSolver->HCones[iCone], dBarrierMaker,
                       dRowDualMaker, dRowDualStepMaker, conePrimal, coneAuxi);
    }
    
    if ( coneDual ) {
        HConeGetDual(HSolver->HCones[iCone], coneDual, coneAuxi);
    }
    
    return;
}

#define DIMACS_ERROR_1 0
#define DIMACS_ERROR_2 1
#define DIMACS_ERROR_3 2
#define DIMACS_ERROR_4 3
#define DIMACS_ERROR_5 4
#define DIMACS_ERROR_6 5
extern hdsdp_retcode HDSDPCheckSolution( hdsdp *HSolver, double dErrs[6] ) {
    /* Compute and report 6 DIMACS error */
    
    hdsdp_retcode retcode = HDSDP_RETCODE_OK;
        
    dErrs[DIMACS_ERROR_1] = 1.0;
    dErrs[DIMACS_ERROR_2] = 1.0;
    dErrs[DIMACS_ERROR_3] = 1.0;
    dErrs[DIMACS_ERROR_4] = 1.0;
    dErrs[DIMACS_ERROR_5] = 1.0;
    dErrs[DIMACS_ERROR_6] = 1.0;
    
    if ( HSolver->dInaccBarrierMaker < 0.0 && HSolver->dAccBarrierMaker < 0.0 ) {
        HSolver->HStatus = HDSDP_NUMERICAL;
        return HDSDP_RETCODE_OK;
    }
    
    double dPrimalInfeas = HDSDP_INFINITY;
    /* We guarantee dual infeasibility up to a diagonal perturbation */
    double dDualInfeas = HSolver->dPerturb * sqrt(get_int_feature(HSolver, INT_FEATURE_N_SUMCONEDIMS));
    double dMinPrimalEVal = HDSDP_INFINITY;
    /* We guarantee dual interior point solution */
    double dCompl = 0.0;
    double dPrimalObj = 0.0;
    double dDualObj = 0.0;
    
    double dPrimalScal = get_dbl_feature(HSolver, DBL_FEATURE_RHSSCALING);
    double dDualScal = get_dbl_feature(HSolver, DBL_FEATURE_OBJSCALING);
    
    int nMaxDim = get_int_feature(HSolver, INT_FEATURE_N_MAXCONEDIM);
    int nMaxBufferDim = 0;
    int nConeDimSqr = 0;
    
    for ( int iCone = 0; iCone < HSolver->nCones; ++iCone ) {
        nMaxBufferDim = HDSDP_MAX(nMaxBufferDim, HConeGetVarBufferDim(HSolver->HCones[iCone]));
    }
    
    double *dPrimalMatBuffer = NULL;
    double *dDualMatBuffer = NULL;
    double *dAuxiMat = NULL;
    double dEValAuxi[2] = {0.0, 0.0};
    
    int lWork = nMaxDim * 30;
    int liWork = nMaxDim * 12;
    double *dWork = NULL;
    int *iWork = NULL;
    
    /* Prepare space */
    HDSDP_INIT(dPrimalMatBuffer, double, nMaxBufferDim);
    HDSDP_MEMCHECK(dPrimalMatBuffer);
    HDSDP_INIT(dDualMatBuffer, double, nMaxBufferDim);
    HDSDP_MEMCHECK(dDualMatBuffer);
    HDSDP_INIT(dAuxiMat, double, nMaxBufferDim);
    HDSDP_MEMCHECK(dAuxiMat);
    HDSDP_INIT(iWork, int, liWork);
    HDSDP_MEMCHECK(iWork);
    HDSDP_INIT(dWork, double, lWork);
    HDSDP_MEMCHECK(dWork);
    
    /* Compute dual objective */
    for ( int iRow = 0; iRow < HSolver->nRows; ++iRow ) {
        dDualObj += HSolver->rowRHS[iRow] * HSolver->dRowDual[iRow];
    }
        
    /* Check DIMACS errors of the SDP solution */
    HDSDP_ZERO(HSolver->dHAuxiVec1, double, HSolver->nRows);
    for ( int iCone = 0; iCone < HSolver->nCones; ++iCone ) {
        HDSDPGetConeValues(HSolver, iCone, dPrimalMatBuffer, dDualMatBuffer, dAuxiMat);
        /* Get A * X */
        HConeComputeATimesX(HSolver->HCones[iCone], dPrimalMatBuffer, HSolver->dHAuxiVec1);
        /* Compute complementarity */
        nConeDimSqr = HConeGetVarBufferDim(HSolver->HCones[iCone]);
        dCompl += dot(&nConeDimSqr, dPrimalMatBuffer, &HIntConstantOne, dDualMatBuffer, &HIntConstantOne);
        dPrimalObj += HConeComputeTraceCX(HSolver->HCones[iCone], dPrimalMatBuffer);
        
        /* Eigen-decompose X to get the minimum eigenvalue */
        if ( HSolver->HCones[iCone]->cone == HDSDP_CONETYPE_LP ) {
            dEValAuxi[0] = HUtilGetDblMinimum(nConeDimSqr, dPrimalMatBuffer);
        } else {
            int nConeDim = HConeGetDim(HSolver->HCones[iCone]);
            fds_syev(nConeDim, dPrimalMatBuffer, dEValAuxi, dAuxiMat, 1, dWork, iWork, lWork, liWork);
        }
        
        dMinPrimalEVal = HDSDP_MIN(dMinPrimalEVal, dEValAuxi[0]);
    }
    
    dDualObj = dDualObj / (dPrimalScal * dDualScal);
    dPrimalObj = dPrimalObj / (dPrimalScal * dDualScal);
    
    for ( int iRow = 0; iRow < HSolver->nRows; ++iRow ) {
        HSolver->dHAuxiVec1[iRow] -= HSolver->rowRHS[iRow];
    }
    
    dPrimalInfeas = nrm2(&HSolver->nRows, HSolver->dHAuxiVec1, &HIntConstantOne);
    dPrimalInfeas = dPrimalInfeas / get_dbl_feature(HSolver, DBL_FEATURE_RHSSCALING);
    dDualInfeas = dDualInfeas / get_dbl_feature(HSolver, DBL_FEATURE_OBJSCALING);
    
    /* Collect DIMACS errors */
    /* DIMACS error 1: pInfeas / (1 + ||b||_1)*/
    dErrs[DIMACS_ERROR_1] = dPrimalInfeas / (1.0 + get_dbl_feature(HSolver, DBL_FEATURE_RHSONENORM));
    /* DIMACS error 2: -lambda_min(X) / (1 + ||b||_1) */
    if ( dMinPrimalEVal < 0.0 ) {
        dErrs[DIMACS_ERROR_2] = - dMinPrimalEVal / (1.0 + get_dbl_feature(HSolver, DBL_FEATURE_RHSONENORM));
    } else {
        dErrs[DIMACS_ERROR_2] = 0.0;
    }
    /* DIMACS error 3: dInfeas / (1 + ||C||) */
    dErrs[DIMACS_ERROR_3] = dDualInfeas / (1.0 + get_dbl_feature(HSolver, DBL_FEATURE_OBJONENORM));
    /* DIMACS error 4: -lambda_min(S) / (1 + ||C||_1) */
    dErrs[DIMACS_ERROR_4] = 0.0;
    /* DIMACS error 5 */
    dErrs[DIMACS_ERROR_5] = (dPrimalObj - dDualObj) / (fabs(dPrimalObj) + fabs(dDualObj) + 1.0);
    /* DIMACS error 6 */
    dErrs[DIMACS_ERROR_6] = dCompl / (fabs(dPrimalObj) + fabs(dDualObj) + 1.0);
    
    double dMaxDimacsErr = 0.0;
    
    for ( int iElem = 0; iElem < 6; ++iElem ) {
        dMaxDimacsErr = HDSDP_MAX(dMaxDimacsErr, fabs(dErrs[iElem]));
    }
    
    HSolver->pObjVal = dPrimalObj;
    HSolver->dObjVal = dDualObj;
    
    // hdsdp_printf("Maximum error: %5.2e \n", dMaxDimacsErr);
    if ( dMaxDimacsErr > 1e-02 ) {
        if ( HSolver->dAccBarrierMaker < 0.0 ) {
            HSolver->HStatus = HDSDP_NUMERICAL;
        } else {
            /* The primal solution is not good. Switch to the other */
            hdsdp_printf("\nDealing with primal solution\n");
            HSolver->dAccBarrierMaker = -1.0;
            HDSDP_FREE(dPrimalMatBuffer);
            HDSDP_FREE(dDualMatBuffer);
            HDSDP_FREE(dAuxiMat);
            HDSDP_FREE(dWork);
            HDSDP_FREE(iWork);
            return HDSDPCheckSolution(HSolver, dErrs);
        }
        
    } else {
        HSolver->HStatus = HDSDP_PRIMAL_DUAL_OPTIMAL;
    }
    
exit_cleanup:
    
    HDSDP_FREE(dPrimalMatBuffer);
    HDSDP_FREE(dDualMatBuffer);
    HDSDP_FREE(dAuxiMat);
    HDSDP_FREE(dWork);
    HDSDP_FREE(iWork);

    return retcode;
}

extern void HDSDPClear( hdsdp *HSolver ) {
    
    if ( !HSolver ) {
        return;
    }
    
    HDSDP_FREE(HSolver->rowRHS);
    
    for ( int iCone = 0; iCone < HSolver->nCones; ++iCone ) {
        HConeDestroy(&HSolver->HCones[iCone]);
    }
    
    HDSDP_FREE(HSolver->HCones);
    HConeDestroy(&HSolver->HBndCone);
    
    HKKTDestroy(&HSolver->HKKT);
    
    HDSDP_FREE(HSolver->dRowDual);
    HDSDP_FREE(HSolver->dRowDualStep);
    HDSDP_FREE(HSolver->dMinvRowRHS);
    HDSDP_FREE(HSolver->dMinvASinv);
    HDSDP_FREE(HSolver->dMinvASinvRdSinv);
    HDSDP_FREE(HSolver->dMinvASinvCSinv);
    HDSDP_FREE(HSolver->dHAuxiVec1);
    HDSDP_FREE(HSolver->dHAuxiVec2);
    HDSDP_FREE(HSolver->dPInfeasUpper);
    HDSDP_FREE(HSolver->dPInfeasLower);
    HDSDP_FREE(HSolver->dInaccRowDualMaker);
    HDSDP_FREE(HSolver->dAccRowDualMaker);
    HDSDP_FREE(HSolver->dInaccRowDualStepMaker);
    HDSDP_FREE(HSolver->dAccRowDualStepMaker);
    
    HDSDP_ZERO(HSolver, hdsdp, 1);
    
    return;
}

extern void HDSDPDestroy( hdsdp **pHSolver ) {
    
    if ( !pHSolver ) {
        return;
    }
    
    HDSDPClear(*pHSolver);
    HDSDP_FREE(*pHSolver);
    
    hdsdp_printf("HDSDP ends. Exiting \n");
    
    return;
}
