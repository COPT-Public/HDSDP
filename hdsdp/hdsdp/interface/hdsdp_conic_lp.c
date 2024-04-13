#ifdef HEADERPATH
#include "interface/def_hdsdp.h"
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
#include "def_hdsdp.h"
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

/* Implement the LP cone
   A' * y + s = c
*/

static inline int LPConeICheckInterior( int nCol, double *colDual ) {
    
    for ( int iCol = 0; iCol < nCol; ++iCol ) {
        if ( colDual[iCol] <= 0.0 ) {
            return 0;
        }
    }
    
    return 1;
}

static inline void LPConeIUpdateBuffer( hdsdp_cone_lp *cone, double dCCoef, double dACoefScal, double *dACoef, double dEyeCoef, int whichBuffer ) {
    
    double *target = NULL;
    
    switch (whichBuffer) {
        case BUFFER_DUALVAR:
            target = cone->colDual;
            break;
        case BUFFER_DUALCHECK:
            target = cone->colDualChecker;
            break;
        case BUFFER_DUALSTEP:
            target = cone->colDualStep;
            break;
        default:
            break;
    }
    
    HDSDP_ZERO(target, double, cone->nCol);
    
    csp_Axpby(cone->nRow, cone->rowMatBeg, cone->rowMatIdx, cone->rowMatElem, dACoefScal, dACoef, target);
    
    for ( int iCol = 0; iCol < cone->nCol; ++iCol ) {
        target[iCol] += dCCoef * cone->colObj[iCol];
    }
    
    if ( whichBuffer != BUFFER_DUALSTEP ) {
        dEyeCoef += cone->dualPerturb;
    }
    
    if ( dEyeCoef != 0.0 ) {
        for ( int iCol = 0; iCol < cone->nCol; ++iCol ) {
            target[iCol] += dEyeCoef;
        }
    }
    
    return;
}

extern hdsdp_retcode LPConeCreateImpl( hdsdp_cone_lp **pCone ) {
    
    hdsdp_retcode retcode = HDSDP_RETCODE_OK;
    
    HDSDP_NULLCHECK(pCone);
    hdsdp_cone_lp *cone = NULL;
    
    HDSDP_INIT(cone, hdsdp_cone_lp, 1);
    HDSDP_MEMCHECK(cone);
    HDSDP_ZERO(cone, hdsdp_cone_lp, 1);
    *pCone = cone;
    
exit_cleanup:
    return retcode;
}

extern hdsdp_retcode LPConeProcDataImpl( hdsdp_cone_lp *cone, int nRow, int nCol, int *coneMatBeg, int *coneMatIdx, double *coneMatElem ) {
    
    hdsdp_retcode retcode = HDSDP_RETCODE_OK;
    
    cone->nRow = nRow;
    cone->nCol = nCol;
    
    HDSDP_INIT(cone->colObj, double, nCol);
    HDSDP_MEMCHECK(cone->colObj);
    HDSDP_INIT(cone->colDual, double, nCol);
    HDSDP_MEMCHECK(cone->colDual);
    HDSDP_INIT(cone->colDualChecker, double, nCol);
    HDSDP_MEMCHECK(cone->colDualChecker);
    HDSDP_INIT(cone->colDualInverse, double, nCol);
    HDSDP_MEMCHECK(cone->colDualInverse);
    HDSDP_INIT(cone->colDualStep, double, nCol);
    HDSDP_MEMCHECK(cone->colDualStep);
    HDSDP_INIT(cone->colBuffer, double, nCol);
    HDSDP_MEMCHECK(cone->colBuffer);
    
    int nObjNnz = coneMatBeg[1];
    
    HDSDP_INIT(cone->rowMatBeg, int, nRow + 1);
    HDSDP_MEMCHECK(cone->rowMatBeg);
    HDSDP_INIT(cone->rowMatIdx, int, coneMatBeg[nRow + 1] - nObjNnz);
    HDSDP_MEMCHECK(cone->rowMatIdx);
    HDSDP_INIT(cone->rowMatElem, double, coneMatBeg[nRow + 1] - nObjNnz);
    HDSDP_MEMCHECK(cone->rowMatElem);
    
    /* Copy LP coefficient */
    for ( int iElem = 0; iElem < coneMatBeg[1]; ++iElem ) {
        cone->colObj[coneMatIdx[iElem]] = coneMatElem[iElem];
    }
    
    /* Copy LP constraint matrix */
    HDSDP_MEMCPY(cone->rowMatBeg, coneMatBeg + 1, int, nRow + 1);
    HDSDP_MEMCPY(cone->rowMatIdx, coneMatIdx + nObjNnz, int, coneMatBeg[nRow + 1] - nObjNnz);
    HDSDP_MEMCPY(cone->rowMatElem, coneMatElem + nObjNnz, double, coneMatBeg[nRow + 1] - nObjNnz);
    
    /* Remove the shift of objective */
    for ( int iCol = 0; iCol < nRow + 1; ++iCol ) {
        cone->rowMatBeg[iCol] -= nObjNnz;
    }
    
    int nMaxKKTNz = 0;
    nMaxKKTNz = nRow * nRow;
    cone->coneKKTNnz = nMaxKKTNz;
    
#if 0
    for ( int iCol = 0; iCol < cone->nCol; ++iCol ) {
        nMaxKKTNz += (cone->rowMatBeg[iCol + 1] - cone->rowMatBeg[iCol]) * \
                        (cone->rowMatBeg[iCol + 1] - cone->rowMatBeg[iCol]);
    }
    HDSDP_INIT(cone->kktMapping, int, nMaxKKTNz);
    HDSDP_MEMCHECK(cone->kktMapping);
#endif
    
exit_cleanup:
    return retcode;
}

extern hdsdp_retcode LPConePresolveImpl( hdsdp_cone_lp *cone ) {
    
    return HDSDP_RETCODE_OK;
}

extern void LPConeSetStartImpl( hdsdp_cone_lp *cone, double rResi ) {
    
    cone->dualResidual = rResi;
    return;
}

extern double LPConeGetObjNorm( hdsdp_cone_lp *cone, int whichNorm ) {
    
    double dNorm = 0.0;
    if ( whichNorm == ABS_NORM ) {
        for ( int iCol = 0; iCol < cone->nCol; ++iCol ) {
            dNorm += fabs(cone->colObj[iCol]);
        }
    } else {
        for ( int iCol = 0; iCol < cone->nCol; ++iCol ) {
            dNorm += cone->colObj[iCol] * cone->colObj[iCol];
        }
        
        dNorm = sqrt(dNorm);
    }
    
    return dNorm;
}

extern double LPConeGetCoeffNorm( hdsdp_cone_lp *cone, int whichNorm ) {
    
    if ( whichNorm == ABS_NORM ) {
        return csp_sum_abs(cone->nRow, cone->rowMatBeg, cone->rowMatIdx, cone->rowMatElem);
    } else {
        return csp_fro_norm(cone->nRow, cone->rowMatBeg, cone->rowMatIdx, cone->rowMatElem);
    }
    
    assert( 0 );
}

extern void LPConeScal( hdsdp_cone_lp *cone, double dScal ) {
    
    rscl(&cone->nCol, &dScal, cone->colObj, &HIntConstantOne);
    return;
}

extern void LPConeUpdateImpl( hdsdp_cone_lp *cone, double barHsdTau, double *rowDual ) {
    
    LPConeIUpdateBuffer(cone, barHsdTau, -1.0, rowDual, -cone->dualResidual, BUFFER_DUALVAR);
    
    return;
}

extern hdsdp_retcode LPConeRatioTestImpl( hdsdp_cone_lp *cone, double barHsdTauStep, double *rowDualStep, double dAdaRatio,
                                         int whichBuffer, double *maxStep ) {
    
    
    double *target = NULL;
    LPConeIUpdateBuffer(cone, barHsdTauStep, -1.0, rowDualStep, dAdaRatio * cone->dualResidual, BUFFER_DUALSTEP);
    
    if ( whichBuffer == BUFFER_DUALVAR ) {
        target = cone->colDual;
    } else {
        target = cone->colDualChecker;
    }
    
    double dStepTmp = HDSDP_INFINITY;
    double dMaxStep = 0.0;
    
    for ( int iCol = 0; iCol < cone->nCol; ++iCol ) {
        dStepTmp = cone->colDualStep[iCol] / target[iCol];
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

extern int LPConeGetDim( hdsdp_cone_lp *cone ) {
    
    return cone->nCol;
}

extern hdsdp_retcode LPConeGetKKT( hdsdp_cone_lp *cone, int iCone, void *kkt, int typeKKT ) {
    
    hdsdp_kkt *Hkkt = (hdsdp_kkt *) kkt;
    
    /* We currently implement the dense version */
    assert( !Hkkt->isKKTSparse );
    
    if ( typeKKT == KKT_TYPE_PRIMAL ) {
        if ( !Hkkt->dPrimalX || !Hkkt->dPrimalX[iCone] ) {
            return HDSDP_RETCODE_FAILED;
        }
        for ( int iCol = 0; iCol < cone->nCol; ++iCol ) {
            cone->colDualInverse[iCol] = Hkkt->dPrimalX[iCone][iCol];
        }
    } else {
        for ( int iCol = 0; iCol < cone->nCol; ++iCol ) {
            cone->colDualInverse[iCol] = 1.0 / cone->colDual[iCol];
        }
    }
    
    csp_ATxpby(cone->nRow, cone->rowMatBeg, cone->rowMatIdx,
              cone->rowMatElem, 1.0, cone->colDualInverse, Hkkt->dASinvVec);
    
    if ( cone->dualResidual ) {
        
        for ( int iCol = 0; iCol < cone->nCol; ++iCol ) {
            Hkkt->dTraceSinv += cone->colDualInverse[iCol];
        }
        
        for ( int iCol = 0; iCol < cone->nCol; ++iCol ) {
            cone->colBuffer[iCol] = cone->dualResidual * cone->colDualInverse[iCol] * cone->colDualInverse[iCol];
        }
        csp_ATxpby(cone->nRow, cone->rowMatBeg, cone->rowMatIdx, cone->rowMatElem,
                   1.0, cone->colBuffer, Hkkt->dASinvRdSinvVec);
    }
    
    if ( typeKKT == KKT_TYPE_CORRECTOR ) {
        return HDSDP_RETCODE_OK;
    }
    
    /* Set up KKT system M <- M + A * S^-2 * A^T */
    if ( Hkkt->isKKTSparse ) {
        assert( 0 );
    } else {
        for ( int iRow = 0; iRow < cone->nRow; ++iRow ) {
            HDSDP_ZERO(cone->colBuffer, double, cone->nCol);
            for ( int iElem = cone->rowMatBeg[iRow]; iElem < cone->rowMatBeg[iRow + 1]; ++iElem ) {
                int iCol = cone->rowMatIdx[iElem];
                cone->colBuffer[iCol] = \
                    cone->rowMatElem[iElem] * (cone->colDualInverse[iCol] * cone->colDualInverse[iCol]);
            }
            for ( int iCol = 0; iCol <= iRow; ++iCol ) {
                double dMElem = 0.0;
                for ( int iElem = cone->rowMatBeg[iCol]; iElem < cone->rowMatBeg[iCol + 1]; ++iElem ) {
                    dMElem += cone->rowMatElem[iElem] * cone->colBuffer[cone->rowMatIdx[iElem]];
                }
                FULL_ENTRY(Hkkt->kktMatElem, cone->nRow, iRow, iCol) += dMElem;
            }
        }
    }
    
    if ( typeKKT == KKT_TYPE_HOMOGENEOUS ) {
        for ( int iCol = 0; iCol < cone->nCol; ++iCol ) {
            double dcsinv = cone->colObj[iCol] * cone->colDualInverse[iCol];
            Hkkt->dCSinv += dcsinv;
            Hkkt->dCSinvCSinv += dcsinv * dcsinv;
        }
        
        for ( int iCol = 0; iCol < cone->nCol; ++iCol ) {
            cone->colBuffer[iCol] = cone->colObj[iCol] * cone->colDualInverse[iCol] * cone->colDualInverse[iCol];
        }
        
        csp_ATxpby(cone->nRow, cone->rowMatBeg, cone->rowMatIdx, cone->rowMatElem, 1.0, cone->colBuffer, Hkkt->dASinvCSinvVec);
    }
    
    return HDSDP_RETCODE_OK;
}

extern void LPConeBuildPrimalXSXDirection( hdsdp_cone_lp *cone, void *kkt, double *dPrimalScalMatrix, double *dPrimalXSXBuffer, int iDualMat ) {
    
    if ( iDualMat ) {
        for ( int iCol = 0; iCol < cone->nCol; ++iCol ) {
            dPrimalXSXBuffer[iCol] = \
            dPrimalScalMatrix[iCol] * dPrimalScalMatrix[iCol] * cone->colDual[iCol];
        }
    } else {
        for ( int iCol = 0; iCol < cone->nCol; ++iCol ) {
            dPrimalXSXBuffer[iCol] = \
            dPrimalScalMatrix[iCol] * dPrimalScalMatrix[iCol] * cone->colDualStep[iCol];
        }
    }
    
    return;
}

extern double LPConeXDotS( hdsdp_cone_lp *cone, double *dConePrimal ) {
    
    return dot(&cone->nCol, cone->colDual, &HIntConstantOne, dConePrimal, &HIntConstantOne);
}

extern hdsdp_retcode LPConeGetKKTFixedStrategy( hdsdp_cone_lp *cone, int iCone, void *kkt, int typeKKT, int ikktStrategy ) {
    
    (void) ikktStrategy;
    return LPConeGetKKT(cone, iCone, kkt, typeKKT);
}

extern int64_t LPConeGetSymNnzImpl( hdsdp_cone_lp *cone ) {
    
    return cone->coneKKTNnz;
}

extern void LPConeAddSymNnzImpl( hdsdp_cone_sdp_dense *cone, int iCol, int *schurMatCol ) {
    
    assert( 0 );
    return;
}

extern void LPConeGetSymMapping( hdsdp_cone_lp *cone, int iCol, int *schurMatCol ) {
    
    assert( 0 );
    return;
}

extern hdsdp_retcode LPConeInteriorCheck( hdsdp_cone_lp *cone, double barHsdTau, double *rowDual, int *isInterior ) {
    
    LPConeUpdateImpl(cone, barHsdTau, rowDual);
    *isInterior = LPConeICheckInterior(cone->nCol, cone->colDual);
    
    return HDSDP_RETCODE_OK;
}

extern hdsdp_retcode LPConeInteriorCheckExpert( hdsdp_cone_lp *cone, double dCCoef, double dACoefScal, double *dACoef,
                                                double dEyeCoef, int whichBuffer, int *isInterior ) {
    
    LPConeIUpdateBuffer(cone, dCCoef, dACoefScal, dACoef, dEyeCoef, whichBuffer);
    
    if ( whichBuffer == BUFFER_DUALVAR ) {
        *isInterior = LPConeICheckInterior(cone->nCol, cone->colDual);
    } else {
        *isInterior = LPConeICheckInterior(cone->nCol, cone->colDualChecker);
    }
    
    return HDSDP_RETCODE_OK;
}

extern void LPConeReduceResidual( hdsdp_cone_lp *cone, double resiReduction ) {
    
    cone->dualResidual = resiReduction;
    return;
}

extern void LPConeSetPerturb( hdsdp_cone_lp *cone, double dDualPerturb ) {
    
    assert( dDualPerturb >= 0.0 );
    cone->dualPerturb = dDualPerturb;
    return;
}

extern hdsdp_retcode LPConeGetBarrier( hdsdp_cone_lp *cone, double barHsdTau, double *rowDual, int whichBuffer, double *logdet ) {
    
    double dLogDeterminant = 0.0;
    double *target = NULL;
    
    if ( whichBuffer == BUFFER_DUALVAR ) {
        target = cone->colDual;
    } else {
        target = cone->colDualChecker;
    }
    
    if ( rowDual ) {
        LPConeUpdateImpl(cone, barHsdTau, rowDual);
    }
    
    for ( int iCol = 0; iCol < cone->nCol; ++iCol ) {
        dLogDeterminant += log(target[iCol]);
    }
    
    *logdet = dLogDeterminant;
    return HDSDP_RETCODE_OK;
}

extern hdsdp_retcode LPConeAddStepToBufferAndCheck( hdsdp_cone_lp *cone, double dStep, int whichBuffer, int *isInterior ) {
    
    double *target = NULL;
    
    if ( whichBuffer == BUFFER_DUALVAR ) {
        target = cone->colDual;
    } else {
        target = cone->colDualChecker;
        HDSDP_MEMCPY(cone->colDualChecker, cone->colDual, double, cone->nCol);
    }
    
    axpy(&cone->nCol, &dStep, cone->colDualStep, &HIntConstantOne, target, &HIntConstantOne);
    *isInterior = LPConeICheckInterior(cone->nCol, target);
    
    return HDSDP_RETCODE_OK;
}

extern void LPConeGetPrimal( hdsdp_cone_lp *cone, double dBarrierMu, double *dRowDual, double *dRowDualStep, double *dConePrimal, double *dAuxiMat ) {
    
    int isInterior = 0;
    LPConeInteriorCheckExpert(cone, 1.0, -1.0, dRowDual, 0.0, BUFFER_DUALCHECK, &isInterior);
    
    if ( !isInterior ) {
        hdsdp_printf("Recovery step is infeasible\n");
        return;
    }
    
    LPConeIUpdateBuffer(cone, 0.0, 1.0, dRowDualStep, 0.0, BUFFER_DUALSTEP);
    
    HDSDP_ZERO(dConePrimal, double, cone->nCol);
    
    for ( int iCol = 0; iCol < cone->nCol; ++iCol ) {
        dConePrimal[iCol] = dBarrierMu * ( cone->colDualChecker[iCol] + cone->colDualStep[iCol] ) / \
                            (cone->colDualChecker[iCol] * cone->colDualChecker[iCol]);
    }
    
    return;
}

extern void LPConeGetDual( hdsdp_cone_lp *cone, double *dConeDual, double *ddummy ) {
    
    HDSDP_MEMCPY(dConeDual, cone->colDual, double, cone->nCol);
    return;
}

extern double LPConeTraceCX( hdsdp_cone_lp *cone, double *dConePrimal ) {
    
    double dTraceCX = 0.0;
    
    for ( int iCol = 0; iCol < cone->nCol; ++iCol ) {
        dTraceCX += cone->colObj[iCol] * dConePrimal[iCol];
    }
    
    return dTraceCX;
}

extern void LPConeATimesX( hdsdp_cone_lp *cone, double *dPrimalX, double *dATimesX ) {
    
    csp_ATxpby(cone->nRow, cone->rowMatBeg, cone->rowMatIdx,
               cone->rowMatElem, 1.0, dPrimalX, dATimesX);
    
    return;
}

extern void LPConeClearImpl( hdsdp_cone_lp *cone ) {
    
    if ( !cone ) {
        return;
    }
    
    HDSDP_FREE(cone->colObj);
    HDSDP_FREE(cone->colDual);
    HDSDP_FREE(cone->colDualInverse);
    HDSDP_FREE(cone->colDualChecker);
    HDSDP_FREE(cone->colDualStep);
    HDSDP_FREE(cone->colBuffer);
    HDSDP_FREE(cone->rowMatBeg);
    HDSDP_FREE(cone->rowMatIdx);
    HDSDP_FREE(cone->rowMatElem);
    HDSDP_FREE(cone->kktMapping);
    
    HDSDP_ZERO(cone, hdsdp_cone_lp, 1);
    
    return;
}

extern void LPConeDestroyImpl( hdsdp_cone_lp **pCone ) {
    
    if ( !pCone ) {
        return;
    }
    
    LPConeClearImpl(*pCone);
    HDSDP_FREE(*pCone);
    
    return;
}

extern void LPConeViewImpl( hdsdp_cone_lp *cone ) {
    
    hdsdp_printf("LP Cone of %d variables and %d constraints \n", cone->nRow, cone->nCol);
    
    return;
}

extern void LPConeGetStatsImpl( hdsdp_cone_lp *cone, double *rowRHS, int coneIntFeatures[20], double coneDblFeatures[20] ) {
    
    /* Detect two special structures for LPs:
     
     1. Implied dual equality (free primal variable)
     2. Implied dual box constraint
    
     */
    
    if ( cone->nCol % 2 != 0 || cone->nCol < 100 ) {
        return;
    }
    
    /* Use two auxiliary arrays from the LP datas structure */
    double *dDualBoundUpperTmp = NULL;
    double *dDualBoundLowerTmp = NULL;
    
    HDSDP_INIT(dDualBoundLowerTmp, double, cone->nRow);
    HDSDP_INIT(dDualBoundUpperTmp, double, cone->nRow);
    HDSDP_ZERO(dDualBoundLowerTmp, double, cone->nRow);
    HDSDP_ZERO(dDualBoundUpperTmp, double, cone->nRow);
    
    if ( !dDualBoundLowerTmp || !dDualBoundUpperTmp ) {
        HDSDP_FREE(dDualBoundLowerTmp);
        HDSDP_FREE(dDualBoundUpperTmp);
        return;
    }
    
    int nHalfCols = (int) cone->nCol / 2;
    
    int isImpliedDual = 1;
    int isImpliedUpper = 0;
    int isImpliedLower = 0;
    
    double dMaxDualBoundUpper = 1.0;
    double dMinDualBoundLower = -1.0;
    
    
    /* First detect if l <= y <= u */
    for ( int iCol = 0; iCol < cone->nRow; ++iCol ) {
        for ( int iElem = cone->rowMatBeg[iCol]; iElem < cone->rowMatBeg[iCol + 1]; ++iElem ) {
            
            if ( cone->rowMatBeg[iCol + 1] - cone->rowMatBeg[iCol] > 2 ) {
                isImpliedDual = 0;
                break;
            }
            
            if ( cone->rowMatElem[iElem] > 0.0 ) {
                if ( dDualBoundUpperTmp[iCol] ) {
                    isImpliedDual = 0;
                    break;
                }
                isImpliedUpper = 1;
                double dUpperBound = cone->colObj[cone->rowMatIdx[iElem]] / cone->rowMatElem[iElem];
                dDualBoundUpperTmp[iCol] = HDSDP_MAX(dDualBoundUpperTmp[iCol], dUpperBound);
            } else {
                if ( dDualBoundLowerTmp[iCol] ) {
                    isImpliedDual = 0;
                    break;
                }
                isImpliedLower = 1;
                double dLowerBound = cone->colObj[cone->rowMatIdx[iElem]] / cone->rowMatElem[iElem];
                dDualBoundLowerTmp[iCol] = HDSDP_MIN(dDualBoundLowerTmp[iCol], dLowerBound);
            }
        }
        
        if ( !isImpliedDual ) {
            break;
        }
    }
    
    if ( isImpliedDual ) {
        
        coneIntFeatures[INT_FEATURE_I_IMPYBOUND] = 1;
        
        if ( isImpliedUpper ) {
            for ( int iRow = 0; iRow < cone->nRow; ++iRow ) {
                dMaxDualBoundUpper = HDSDP_MAX(dMaxDualBoundUpper, dDualBoundUpperTmp[iRow]);
            }
            
            if ( dMaxDualBoundUpper <= 0.0 ) {
                dMaxDualBoundUpper = 1.0;
            }
            
            coneDblFeatures[DBL_FEATURE_IMPYBOUNDUP] = dMaxDualBoundUpper;
        }
        
        if ( isImpliedLower ) {
            for ( int iRow = 0; iRow < cone->nRow; ++iRow ) {
                dMinDualBoundLower = HDSDP_MIN(dMinDualBoundLower, dDualBoundLowerTmp[iRow]);
            }
            
            if ( dMinDualBoundLower >= 0.0 ) {
                dMinDualBoundLower = -1.0;
            }
            
            coneDblFeatures[DBL_FEATURE_IMPYBOUNDLOW] = dMinDualBoundLower;
        }
    }
    
    HDSDP_FREE(dDualBoundLowerTmp);
    HDSDP_FREE(dDualBoundUpperTmp);
    
    for ( int iCol = 0; iCol < nHalfCols; ++iCol ) {
        if ( cone->colObj[iCol] + cone->colObj[iCol + nHalfCols] != 0.0 ) {
            return;
        }
    }
    
    for ( int iCol = 0; iCol < cone->nRow; ++iCol ) {
        int colNnz = cone->rowMatBeg[iCol + 1] - cone->rowMatBeg[iCol];
        int nHalfNnz = (int) colNnz / 2;
        if ( colNnz % 2 != 0 ) {
            return;
        }
        
        for ( int iRow = 0; iRow < nHalfNnz; ++iRow ) {
            if ( cone->rowMatElem[cone->rowMatBeg[iCol] + iRow] + \
                cone->rowMatElem[cone->rowMatBeg[iCol] + iRow + nHalfNnz] != 0.0 ) {
                return;
            }
        }
    }
    
    coneIntFeatures[INT_FEATURE_I_NODINTERIOR] = 1;
    
    return;
}
