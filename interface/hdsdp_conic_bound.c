#ifdef HEADERPATH
#include "interface/hdsdp_conic_sdp.h"
#include "interface/def_hdsdp_user_data.h"
#include "interface/hdsdp_utils.h"
#include "interface/def_hdsdp_schur.h"
#include "linalg/hdsdp_sdpdata.h"
#include "linalg/sparse_opts.h"
#include "linalg/dense_opts.h"
#include "linalg/vec_opts.h"
#include "linalg/hdsdp_linsolver.h"
#include "external/hdsdp_cs.h"
#else
#include "hdsdp_conic_sdp.h"
#include "def_hdsdp_user_data.h"
#include "hdsdp_utils.h"
#include "def_hdsdp_schur.h"
#include "hdsdp_sdpdata.h"
#include "sparse_opts.h"
#include "dense_opts.h"
#include "vec_opts.h"
#include "hdsdp_linsolver.h"
#include "hdsdp_cs.h"
#endif

#include <math.h>

/* Implement the scalar two-sided bound
 
    l * e <= I * y <= u * e
 
   and we rewrite
     - I * y <= -l * e
       I * y <=  u * e
 */

static inline void sBoundConeIUpdateBuffer( hdsdp_cone_bound_scalar *cone, double dCCoef, double dACoefScal, double *dACoef, int whichBuffer ) {
    
    double *ltarget = NULL;
    double *utarget = NULL;
    
    switch (whichBuffer) {
            
        case BUFFER_DUALVAR:
            ltarget = cone->dualLower;
            utarget = cone->dualUpper;
            break;
        case BUFFER_DUALCHECK:
            ltarget = cone->dualLowerChecker;
            utarget = cone->dualUpperChecker;
            break;
        case BUFFER_DUALSTEP:
            ltarget = cone->dualLowerStep;
            utarget = cone->dualUpperStep;
            break;
        default:
            break;
    }
    
    for ( int iRow = 0; iRow < cone->nRow; ++iRow ) {
        utarget[iRow] = dCCoef * cone->dBoundUp + dACoefScal * dACoef[iRow];
        ltarget[iRow] = -dCCoef * cone->dBoundLow - dACoefScal * dACoef[iRow];
    }
    
    return;
}

extern hdsdp_retcode sBoundConeCreateImpl( hdsdp_cone_bound_scalar **pCone ) {
    
    hdsdp_retcode retcode = HDSDP_RETCODE_OK;
    HDSDP_NULLCHECK(pCone);
    
    hdsdp_cone_bound_scalar *cone = NULL;
    
    HDSDP_INIT(cone, hdsdp_cone_bound_scalar, 1);
    HDSDP_MEMCHECK(cone);
    HDSDP_ZERO(cone, hdsdp_cone_bound_scalar, 1);
    *pCone = cone;
    
exit_cleanup:
    return retcode;
}

extern hdsdp_retcode sBoundConeProcDataImpl( hdsdp_cone_bound_scalar *cone, int nRow, int nCol, int *coneMatBeg, int *coneMatIdx, double *coneMatElem ) {
    
    hdsdp_retcode retcode = HDSDP_RETCODE_OK;
    
    cone->nRow = nRow;
    
    HDSDP_INIT(cone->dualLower, double, nRow);
    HDSDP_MEMCHECK(cone->dualLower);
    
    HDSDP_INIT(cone->dualUpper, double, nRow);
    HDSDP_MEMCHECK(cone->dualUpper);
    
    HDSDP_INIT(cone->dualLowerInverse, double, nRow);
    HDSDP_MEMCHECK(cone->dualLowerInverse);
    
    HDSDP_INIT(cone->dualUpperInverse, double, nRow);
    HDSDP_MEMCHECK(cone->dualUpperInverse);
    
    HDSDP_INIT(cone->dualLowerStep, double, nRow);
    HDSDP_MEMCHECK(cone->dualLowerStep);
    
    HDSDP_INIT(cone->dualUpperStep, double, nRow);
    HDSDP_MEMCHECK(cone->dualUpperStep);
    
    HDSDP_INIT(cone->dualUpperChecker, double, nRow);
    HDSDP_MEMCHECK(cone->dualUpperChecker);
    
    HDSDP_INIT(cone->dualLowerChecker, double, nRow);
    HDSDP_MEMCHECK(cone->dualLowerChecker);
    
    cone->dBoundLow = coneMatElem[0];
    cone->dBoundUp = coneMatElem[1];
    
    assert( cone->dBoundLow < cone->dBoundUp );
    
exit_cleanup:
    return retcode;
}

extern hdsdp_retcode sBoundConePresolveDummyImpl( hdsdp_cone_bound_scalar *cone ) {
    
    (void) cone;
    return HDSDP_RETCODE_OK;
}

extern void sBoundConeSetStartDummyImpl( hdsdp_cone_bound_scalar *cone, double dummy ) {
    
    (void) dummy;
    return;
}

extern double sBoundConeGetObjNormDummyImpl( hdsdp_cone_bound_scalar *cone, int dummy ) {
    /* Scalar bound cone is not counted in the objective */
    return 0.0;
}

extern double sBoundConeGetCoeffNormDummyImpl( hdsdp_cone_bound_scalar *cone, int dummy ) {
    /* Scalar bound cone is not counted in the coefficient */
    return 0.0;
}

extern void sBoundConeScalDummyImpl( hdsdp_cone_bound_scalar *cone ) {
    
    return;
}

extern void sBoundConeUpdateImpl( hdsdp_cone_bound_scalar *cone, double barHsdTau, double *rowDual ) {
    
    sBoundConeIUpdateBuffer(cone, barHsdTau, -1.0, rowDual, BUFFER_DUALVAR);
    return;
}

extern hdsdp_retcode sBoundConeRatioTestImpl( hdsdp_cone_bound_scalar *cone, double barHsdTauStep, double *rowDualStep, double dAdaRatio, int whichBuffer, double *maxStep ) {
    
    double *ltarget = NULL;
    double *utarget = NULL;
    
    sBoundConeIUpdateBuffer(cone, barHsdTauStep, -1.0, rowDualStep, BUFFER_DUALSTEP);
    
    if ( whichBuffer == BUFFER_DUALVAR ) {
        ltarget = cone->dualLower;
        utarget = cone->dualUpper;
    } else {
        ltarget = cone->dualLowerChecker;
        utarget = cone->dualUpperChecker;
    }
    
    /* Ensure that sl + dsl >= 0 and su + dsu >= 0 */
    double dStepTmp = HDSDP_INFINITY;
    double dMaxStep = 0.0;
    
    for ( int iRow = 0; iRow < cone->nRow; ++iRow ) {
        dStepTmp = cone->dualLowerStep[iRow] / ltarget[iRow];
        dMaxStep = HDSDP_MIN(dMaxStep, dStepTmp);
        dStepTmp = cone->dualUpperStep[iRow] / utarget[iRow];
        dMaxStep = HDSDP_MIN(dMaxStep, dStepTmp);
    }
    
    dMaxStep = 1.0 / dMaxStep;
    
    if ( dMaxStep > 0.0 ) {
        dMaxStep = 100.0;
    } else {
        dMaxStep = -dMaxStep;
    }
    
    *maxStep = dMaxStep;
    
    return HDSDP_RETCODE_OK;
}

extern int sBoundConeGetDimImpl( hdsdp_cone_bound_scalar *cone ) {
    
    return cone->nRow;
}

extern hdsdp_retcode sBoundConeGetKKT( hdsdp_cone_bound_scalar *cone, void *kkt, int typeKKT ) {
    
    hdsdp_kkt *Hkkt = (hdsdp_kkt *) kkt;
    
    for ( int iRow = 0; iRow < cone->nRow; ++iRow ) {
        cone->dualLowerInverse[iRow] = 1.0 / cone->dualLower[iRow];
        cone->dualUpperInverse[iRow] = 1.0 / cone->dualUpper[iRow];
    }
    
    for ( int iRow = 0; iRow < Hkkt->nRow; ++iRow ) {
        Hkkt->dASinvVec[iRow] -= cone->dualLowerInverse[iRow];
        Hkkt->dASinvVec[iRow] += cone->dualUpperInverse[iRow];
    }
    
    if ( typeKKT == KKT_TYPE_CORRECTOR ) {
        return HDSDP_RETCODE_OK;
    }
    
    /* Set up KKT system */
    for ( int iRow = 0; iRow < cone->nRow; ++iRow ) {
        *Hkkt->kktDiag[iRow] += cone->dualLowerInverse[iRow] * cone->dualLowerInverse[iRow] + \
                                    cone->dualUpperInverse[iRow] * cone->dualUpperInverse[iRow];
    }
    
    if ( typeKKT == KKT_TYPE_HOMOGENEOUS ) {
        /* Get kkt->dCSinv and kkt->dCSinvCSinv */
        double dSinvsqrVal = 0.0;
        for ( int iRow = 0; iRow < cone->nRow; ++iRow ) {
            Hkkt->dCSinv += cone->dBoundUp * cone->dualUpperInverse[iRow];
            dSinvsqrVal = cone->dualUpperInverse[iRow] * cone->dualUpperInverse[iRow];
            Hkkt->dASinvCSinvVec[iRow] += cone->dBoundUp * dSinvsqrVal;
            Hkkt->dCSinvCSinv += cone->dBoundUp * cone->dBoundUp * dSinvsqrVal;
            /* -I * y + s_l = - l
             Hkkt->dCSinv += -(-cone->dBoundLow) * cone->dualUpperInverse[iRow]; */
            Hkkt->dCSinv -= cone->dBoundLow * cone->dualLowerInverse[iRow];
            dSinvsqrVal = cone->dualLowerInverse[iRow] * cone->dualLowerInverse[iRow];
            Hkkt->dASinvCSinvVec[iRow] += cone->dBoundLow * dSinvsqrVal;
            Hkkt->dCSinvCSinv += cone->dBoundLow * cone->dBoundLow * dSinvsqrVal;
        }
    }
    
    return HDSDP_RETCODE_OK;
}

extern hdsdp_retcode sBoundConeGetKKTByFixedStrategy( hdsdp_cone_bound_scalar *cone, void *kkt, int typeKKT, int ikktStrategy ) {
    
    (void) ikktStrategy;
    return sBoundConeGetKKT(cone, kkt, typeKKT);
}

extern int64_t sBoundConeGetSymNnzImpl( hdsdp_cone_bound_scalar *cone ) {
    
    return cone->nRow;
}

extern void sBoundConeAddSymNnzImpl( hdsdp_cone_bound_scalar *cone, int iCol, int *schurMatCol ) {
    
    schurMatCol[iCol] = 1;
    return;
}

extern void sBoundConeGetSymMappingImpl( hdsdp_cone_bound_scalar *cone, int dummy1, int *dummy2 ) {
    
    (void) dummy1;
    (void) dummy2;
    
    return;
}

extern hdsdp_retcode sBoundConeInteriorCheck( hdsdp_cone_bound_scalar *cone, double barHsdTau, double *rowDual, int *isInterior ) {

    *isInterior = 0;
    sBoundConeUpdateImpl(cone, barHsdTau, rowDual);
    
    for ( int iRow = 0; iRow < cone->nRow; ++iRow ) {
        if ( cone->dualLower[iRow] <= 0.0 ) {
            return HDSDP_RETCODE_OK;
        }
    }
    
    for ( int iRow = 0; iRow < cone->nRow; ++iRow ) {
        if ( cone->dualUpper[iRow] <= 0.0 ) {
            return HDSDP_RETCODE_OK;
        }
    }
    
    *isInterior = 1;
    
exit_cleanup:
    return HDSDP_RETCODE_OK;
}

extern hdsdp_retcode sBoundConeInteriorCheckExpert( hdsdp_cone_bound_scalar *cone, double dCCoef, double dACoefScal, double *dACoef, double dEyeCoef, int whichBuffer, int *isInterior ) {
    
    *isInterior = 0;
    sBoundConeIUpdateBuffer(cone, dCCoef, dACoefScal, dACoef, whichBuffer);
    
    double *ltarget = NULL;
    double *utarget = NULL;
    
    if ( whichBuffer == BUFFER_DUALVAR ) {
        ltarget = cone->dualLower;
        utarget = cone->dualUpper;
    } else {
        ltarget = cone->dualLowerChecker;
        utarget = cone->dualUpperChecker;
    }
    
    for ( int iRow = 0; iRow < cone->nRow; ++iRow ) {
        if ( ltarget[iRow] <= 0.0 ) {
            return HDSDP_RETCODE_OK;
        }
    }
    
    for ( int iRow = 0; iRow < cone->nRow; ++iRow ) {
        if ( utarget[iRow] <= 0.0 ) {
            return HDSDP_RETCODE_OK;
        }
    }
    
    *isInterior = 1;
    
exit_cleanup:
    return HDSDP_RETCODE_OK;
}

extern void sBoundConeReduceResidual( hdsdp_cone_bound_scalar *cone, double dummy ) {
    
    (void) dummy;
    return;
}

extern void sBoundConeSetPerturb( hdsdp_cone_bound_scalar *cone, double dummy ) {
    
    (void) dummy;
    return;
}

extern hdsdp_retcode sBoundConeGetBarrier( hdsdp_cone_bound_scalar *cone, double barHsdTau, double *rowDual, int whichBuffer, double *logdet ) {
    
    double dLogDeterminant = 0.0;
    
    double *ltarget = NULL;
    double *utarget = NULL;
    
    if ( whichBuffer == BUFFER_DUALVAR ) {
        ltarget = cone->dualLower;
        utarget = cone->dualUpper;
    } else {
        ltarget = cone->dualLowerChecker;
        utarget = cone->dualUpperChecker;
    }
    
    if ( rowDual ) {
        sBoundConeUpdateImpl(cone, barHsdTau, rowDual);
    }
    
    for ( int iCol = 0; iCol < cone->nRow; ++iCol ) {
        dLogDeterminant += log(ltarget[iCol]) + log(utarget[iCol]);
    }
    
    if ( dLogDeterminant != dLogDeterminant ) {
        hdsdp_printf("Bound constraint is violated\n");
        return HDSDP_RETCODE_FAILED;
    }
    
    *logdet = dLogDeterminant;
 
exit_cleanup:
    return HDSDP_RETCODE_OK;
}

extern hdsdp_retcode sBoundConeAddStepToBufferAndCheck( hdsdp_cone_bound_scalar *cone, double dStep, int whichBuffer, int *isInterior ) {
    
    *isInterior = 0;
    
    double *ltarget = NULL;
    double *utarget = NULL;
    
    if ( whichBuffer == BUFFER_DUALVAR ) {
        ltarget = cone->dualLower;
        utarget = cone->dualUpper;
    } else {
        ltarget = cone->dualLowerChecker;
        utarget = cone->dualUpperChecker;
        HDSDP_MEMCPY(cone->dualLowerChecker, cone->dualLower, double, cone->nRow);
        HDSDP_MEMCPY(cone->dualUpperChecker, cone->dualUpper, double, cone->nRow);
    }
    
    for ( int iRow = 0; iRow < cone->nRow; ++iRow ) {
        ltarget[iRow] += dStep * cone->dualLowerStep[iRow];
        utarget[iRow] += dStep * cone->dualUpperStep[iRow];
        
        if ( ltarget[iRow] <= 0.0 || utarget[iRow] <= 0.0 ) {
            *isInterior = 0;
            return HDSDP_RETCODE_OK;
        }
    }
    
    *isInterior = 1;
    
    return HDSDP_RETCODE_OK;
}

extern void sBoundConeGetPrimal( hdsdp_cone_bound_scalar *cone, double dBarrierMu, double *dRowDual, double *dRowDualStep, double *dBoundLowerPrimal, double *dBoundUpperPrimal ) {
    
    
    double dRowDualStepElem = 0.0;
    double dualLowerElem = 0.0;
    double dualUpperElem = 0.0;
    
    for ( int iRow = 0; iRow < cone->nRow; ++iRow ) {
        dRowDualStepElem = - dRowDualStep[iRow];
        dualUpperElem = cone->dBoundUp - dRowDual[iRow];
        dualLowerElem = dRowDual[iRow] - cone->dBoundLow;
        dBoundLowerPrimal[iRow] = 1.0 / dualLowerElem - dRowDualStepElem / ( dualLowerElem * dualLowerElem );
        dBoundUpperPrimal[iRow] = 1.0 / dualUpperElem + dRowDualStepElem / ( dualUpperElem * dualUpperElem );
        dBoundLowerPrimal[iRow] *= dBarrierMu;
        dBoundUpperPrimal[iRow] *= dBarrierMu;
    }
    
    return;
}

extern void sBoundConeClearImpl( hdsdp_cone_bound_scalar *cone ) {
    
    if ( !cone ) {
        return;
    }
    
    HDSDP_FREE(cone->dualLower);
    HDSDP_FREE(cone->dualUpper);
    HDSDP_FREE(cone->dualLowerInverse);
    HDSDP_FREE(cone->dualUpperInverse);
    HDSDP_FREE(cone->dualLowerStep);
    HDSDP_FREE(cone->dualUpperStep);
    HDSDP_FREE(cone->dualLowerChecker);
    HDSDP_FREE(cone->dualUpperChecker);
    
    HDSDP_ZERO(cone, hdsdp_cone_bound_scalar, 1);
    
    return;
}

extern void sBoundConeDestroyImpl( hdsdp_cone_bound_scalar **pCone ) {
    
    if ( !pCone ) {
        return;
    }
    
    sBoundConeClearImpl(*pCone);
    HDSDP_FREE(*pCone);
    
    return;
}

extern void sBoundConeViewImpl( hdsdp_cone_bound_scalar *cone ) {
    
    printf("- Scaler bound cone of %d rows. \n", cone->nRow);
    printf("- LB: %10.6e UB: %10.6e \n", cone->dBoundLow, cone->dBoundUp);
    
    return;
}
