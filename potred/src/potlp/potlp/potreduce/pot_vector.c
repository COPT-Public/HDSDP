/** @file pot\_vector.c
 *  @brief Implement the vector methods
 *
 * @TODO: Add more detailed comments
 */

#include "pot_vector.h"
#include "vec_mat.h"

#include <math.h>

extern pot_int potVecCreate( pot_vec **ppVec ) {
    
    pot_vec *pVec = NULL;
    POTLP_INIT(pVec, pot_vec, 1);
    
    if ( !pVec ) {
        return RETCODE_FAILED;
    }
    
    *ppVec = pVec;
    return RETCODE_OK;
}

extern pot_int potVecInit( pot_vec *pVec, pot_int vDim, pot_int vConeDim ) {
    
    pot_int retcode = RETCODE_OK;
    
    if ( !pVec ) {
        retcode = RETCODE_FAILED;
        goto exit_cleanup;
    }
    
    if ( pVec->n > 0 || pVec->x ) {
        retcode = RETCODE_FAILED;
        goto exit_cleanup;
    }
    
    if ( vConeDim > vDim || vConeDim < 0 ) {
        retcode = RETCODE_FAILED;
        goto exit_cleanup;
    }
    
    POTLP_INIT(pVec->x, double, vDim);
    
    if ( !pVec->x ) {
        retcode = RETCODE_FAILED;
        goto exit_cleanup;
    }
    
    pVec->n = vDim; pVec->ncone = vConeDim;
    pVec->nrm = -1.0;
    
exit_cleanup:
    return retcode;
}

extern void potVecDiff( pot_vec *pVecOut, pot_vec *pVecInPrev, pot_vec *pVecInPres ) {
    
    assert( pVecOut->n == pVecInPres->n );
    assert( pVecInPres->n == pVecInPrev->n );
    
    for ( int i = 0; i < pVecOut->n; ++i ) {
        pVecOut->x[i] = pVecInPres->x[i] - pVecInPrev->x[i];
    }
    
    pVecOut->nrm = -1.0;
}

extern void potVeczAxpby( pot_vec *pVecZ, double alpha, pot_vec *pVecX, double beta, pot_vec *pVecY ) {
    
    assert( pVecZ->n == pVecX->n );
    assert( pVecX->n == pVecY->n );
    
    for ( int i = 0; i < pVecX->n; ++i ) {
        pVecZ->x[i] = alpha * pVecX->x[i] + beta * pVecY->x[i];
    }
    
    pVecZ->nrm = -1.0;
    return;
}

extern double potVecNormalize( pot_vec *pVec ) {
    
    assert( pVec->nrm == -1.0 );
    pVec->nrm = nrm2(&pVec->n, pVec->x, &potIntConstantOne);
    
    if ( pVec->nrm ) {
        rscl(&pVec->n, &pVec->nrm, pVec->x, &potIntConstantOne);
    }
    
    return pVec->nrm;
}

extern void potVecScal( pot_vec *pVexX, double sVal ) {
    
    scal(&pVexX->n, &sVal, pVexX->x, &potIntConstantOne);
    
    if ( pVexX ->nrm != -1.0 ) {
        pVexX->nrm = pVexX->nrm * fabs(sVal);
    }
    
    return;
}

extern void potVecArrScal( pot_vec *pVecX, double *dArray ) {
    
    for ( int i = 0; i < pVecX->n; ++i ) {
        pVecX->x[i] /= dArray[i];
    }
    
    return;
}

extern double potVecScaledDot( pot_vec *pVecX, pot_vec *pVecY, pot_vec *pVecZ ) {
    
    double snrm = 0.0, tmp1 = 0.0, tmp2;
    for ( int i = pVecX->n - pVecX->ncone; i < pVecX->n; ++i ) {
        tmp1 = pVecY->x[i] / pVecX->x[i];
        tmp2 = pVecZ->x[i] / pVecX->x[i];
        snrm += tmp1 * tmp2;
    }
    
    return snrm;
}

extern void potVecConeNormalize( pot_vec *pVec ) {
    
    pot_int nCone = pVec->ncone;
    double *coneStart = pVec->x + pVec->n - nCone;
    double coneNorm = nrm2(&nCone, coneStart, &potIntConstantOne);
    rscl(&nCone, &coneNorm, coneStart, &potIntConstantOne);
    
    return;
}

extern double potVecDot( pot_vec *pVecX, pot_vec *pVecY ) {
    
    assert( pVecX->n == pVecY->n );
    return dot(&pVecX->n, pVecX->x, &potIntConstantOne, pVecY->x, &potIntConstantOne);
}

/** @brief Utility for setting up the gradient of the potential function
 *
 */
extern void potVecAxinvpBy( double alpha1, double alpha2, pot_vec *pVecX, double beta, pot_vec *pVecY ) {
    
    assert( pVecX->ncone == pVecY->ncone );
    
    double *pX = pVecX->x, *pY = pVecY->x;
    int iCol = 0;
    for ( iCol = 0; iCol < pVecX->n - pVecX->ncone; ++iCol ) {
        pY[iCol] = beta * pY[iCol];
    }
    
    for ( ; iCol < pVecX->n; ++iCol ) {
        pY[iCol] = alpha1 / (pX[iCol] * alpha2) + beta * pY[iCol];
    }
    
    return;
}

extern void potVecAxpy( double alpha, pot_vec *pVecX, pot_vec *pVecY ) {
    
    assert( pVecX->n == pVecY->n );
    axpy(&pVecX->n, &alpha, pVecX->x, &potIntConstantOne, pVecY->x, &potIntConstantOne);
    pVecY->nrm = -1.0;
    
    return;
}

extern void potVecConeAxpy( double alpha, pot_vec *pVecX, pot_vec *pVecY ) {
    
    assert( pVecX->n == pVecY->n );
    pot_int cShift = pVecX->n - pVecX->ncone;
    axpy(&pVecX->ncone, &alpha, pVecX->x + cShift, &potIntConstantOne, pVecY->x + cShift, &potIntConstantOne);
    pVecY->nrm = -1.0;
    
    return;
}

/** @brief Used for setting up the the Hessian vector product
 * alpha is expected to be a small coefficient that cancels the dimensions in x close to 0.0
 */
extern void potVecConeAxinvsqrVpy( double alpha, pot_vec *pVecX, pot_vec *pVecV, pot_vec *pVecY ) {
    
    double alphaInv = 1.0 / alpha;
    double xinvalpha = 0.0;
    
    if ( alpha > 0.0 ) {
        for ( int i = pVecX->n - pVecX->ncone; i < pVecX->n; ++i ) {
            xinvalpha = alphaInv * pVecX->x[i];
            pVecY->x[i] += pVecV->x[i] / (xinvalpha * xinvalpha);
        }
    } else {
        for ( int i = pVecX->n - pVecX->ncone; i < pVecX->n; ++i ) {
            xinvalpha = alphaInv * pVecX->x[i];
            pVecY->x[i] -= pVecV->x[i] / (xinvalpha * xinvalpha);
        }
    }
    
    return;
}

extern void potVecConeScal( pot_vec *pVecX, pot_vec *pVecY ) {
    
    assert( pVecX->ncone == pVecY->ncone );
    pot_int cShift = pVecX->n - pVecX->ncone;
    vvscl(&pVecX->ncone, pVecX->x + cShift, pVecY->x + cShift);
    
    return;
}

extern void potVecSimplexProj( pot_vec *pVecX ) {
    
    double coneSumVal =potVecSumCone(pVecX);
    rscl(&pVecX->n, &coneSumVal, pVecX->x, &potIntConstantOne);
    
    return;
}

extern double potVecConeMin( pot_vec *pVecX ) {
    
    int ishift = pVecX->n - pVecX->ncone;
    int idmin = idamin(&pVecX->ncone, pVecX->x + ishift, &potIntConstantOne);
    
    return pVecX->x[ishift + idmin];
}

extern double potVecConeMax( pot_vec *pVecX ) {
    
    int ishift = pVecX->n - pVecX->ncone;
    int idmin = idamax(&pVecX->ncone, pVecX->x + ishift, &potIntConstantOne);
    
    return pVecX->x[ishift + idmin];
}

extern double potVecLogDet( pot_vec *pVecX ) {
    
    return sumlogdet(&pVecX->ncone, pVecX->x + pVecX->n - pVecX->ncone);
}

extern double potVecSumCone( pot_vec *pVecX ) {
    
    double eTx = 0.0;
    for ( int i = pVecX->n - pVecX->ncone; i < pVecX->n; ++i ) {
        eTx += pVecX->x[i];
    }
    
    return eTx;
}

extern double potVecSumScalCone( pot_vec *pVecX, pot_vec *pVecY ) {
    
    double xTy = 0.0;
    for ( int i = pVecX->n - pVecX->ncone; i < pVecX->n; ++i ) {
        xTy += pVecX->x[i] * pVecY->x[i];
    }
    
    return xTy;
}

extern double potVecRatioTest( pot_vec *pVecX, pot_vec *pVecdX, double dCone ) {
    
    double ratio = 0.0;
    
    if ( dCone == 0.0 ) {
        for ( int i = pVecX->n - pVecX->ncone; i < pVecX->n; ++i ) {
            ratio = POTLP_MIN(ratio, pVecdX->x[i] / pVecX->x[i]);
            
        }
    } else {
        double d = 0.0;
        for ( int i = pVecX->n - pVecX->ncone; i < pVecX->n; ++i ) {
            d = pVecX->x[i] - dCone;
            if ( d < 0.0 ) {
                return 0.0;
            }
            ratio = POTLP_MIN(ratio, pVecdX->x[i] / d);
        }
    }
    
    return (- 1.0 / ratio);
}

extern void potVecConeAddConstant( pot_vec *pVecX, double cVal ) {
    
    for ( int i = pVecX->n - pVecX->ncone; i < pVecX->n; ++i ) {
        pVecX->x[i] += cVal;
    }
    
    pVecX->nrm = -1.0;
    return;
}

extern void potVecCopy( pot_vec *srcVec, pot_vec *dstVec ) {
    
    assert( srcVec->n == dstVec->n );
    
    POTLP_MEMCPY(dstVec->x, srcVec->x, double, srcVec->n);
    
    dstVec->nrm = srcVec->nrm;
    
    return;
}

extern void potVecReset( pot_vec *pVec ) {
    
    POTLP_ZERO(pVec->x, double, pVec->n);
    pVec->nrm = -1.0;
    
    return;
}

extern void potVecImport( pot_vec *pVec, double *dVec ) {
    
    POTLP_MEMCPY(pVec->x, dVec, double, pVec->n);
    pVec->nrm = -1.0;
    
    return;
}

extern void potVecExport( pot_vec *pVec, double *dVec ) {
    
    POTLP_MEMCPY(dVec, pVec->x, double, pVec->n);
    return;
}

extern void potVecClear( pot_vec *pVec ) {
    
    if (!pVec) {
        return;
    }
    
    pVec->n = 0; pVec->ncone = 0; pVec->nrm = -1.0;
    POTLP_FREE(pVec->x);
    
    return;
}

extern void potVecDestroy( pot_vec **ppVec ) {
    
    if (!ppVec) {
        return;
    }
    
    potVecClear(*ppVec);
    POTLP_FREE(*ppVec);
    
    return;
}
