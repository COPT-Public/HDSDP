#include "pot_lanczos.h"
#include "pot_vector.h"
#include "pot_solver.h"
#include "vec_mat.h"

#include <math.h>

#define SYEV_WORK 30
#define SYEV_IWORK 12

static void dArrSymmetrize( pot_int n, double *dArray ) {
    
    double Aij, Aji;
    for ( int i = 0, j; i < n; ++i ) {
        for ( j = i + 1; j < n; ++j ) {
            Aij = dArray[j * n + i];
            Aji = dArray[i * n + j];
            dArray[j * n + i] = dArray[i * n + j] = (Aij + Aji) * 0.5;
        }
    }
    
    return;
}

static void potLanczosIInitV( pot_vec *vVec ) {
    
    srand(vVec->ncone);
    for ( int i = 0; i < vVec->x[i]; ++i ) {
        srand(rand());
        vVec->x[i] = sqrt(sqrt((rand() % 1627))) * (rand() % 2 - 0.5);
    }
    
    return;
}

extern pot_int potLanczosCreate( pot_lanczos **ppotLanczos ) {
    
    pot_int retcode = RETCODE_OK;
    
    if ( !ppotLanczos ) {
        retcode = RETCODE_FAILED;
        goto exit_cleanup;
    }
    
    pot_lanczos *potLanczos = NULL;
    POTLP_INIT(potLanczos, pot_lanczos, 1);
    
    if ( !potLanczos ) {
        retcode = RETCODE_FAILED;
        goto exit_cleanup;
    }
    
    POTLP_ZERO(potLanczos, pot_lanczos, 1);
    *ppotLanczos = potLanczos;
    
exit_cleanup:
    return retcode;
}

extern pot_int potLanczosInit( pot_lanczos *potLanczos, pot_int nCols, pot_int nCones ) {
    
    pot_int retcode = RETCODE_OK;
    potLanczos->n = nCols;
    
    if ( potLanczos->vVec   || potLanczos->wVec     || potLanczos->z1Vec ||
         potLanczos->z2Vec  || potLanczos->vaVec    || potLanczos->VMat  ||
         potLanczos->HMat   || potLanczos->YMat     || potLanczos->UMat  ||
         potLanczos->dArray || potLanczos->eiDblMat || potLanczos->eiIntMat ) {
        
        retcode = RETCODE_FAILED;
        goto exit_cleanup;
    }
    
    potLanczos->maxIter = POTLP_MIN(nCols, 1000);
    
    POT_CALL(potVecCreate(&potLanczos->vVec));
    POT_CALL(potVecInit(potLanczos->vVec, nCols, nCones));
    
    POT_CALL(potVecCreate(&potLanczos->wVec));
    POT_CALL(potVecInit(potLanczos->wVec, nCols, nCones));
    
    POT_CALL(potVecCreate(&potLanczos->z1Vec));
    POT_CALL(potVecInit(potLanczos->z1Vec, nCols, nCones));
    
    POT_CALL(potVecCreate(&potLanczos->z2Vec));
    POT_CALL(potVecInit(potLanczos->z2Vec, nCols, nCones));
    
    POT_CALL(potVecCreate(&potLanczos->vaVec));
    POT_CALL(potVecInit(potLanczos->vaVec, nCols, nCones));
    
    POTLP_INIT(potLanczos->VMat, double, nCols * (potLanczos->maxIter + 1));
    POTLP_INIT(potLanczos->HMat, double, (potLanczos->maxIter + 1) * (potLanczos->maxIter + 1));
    POTLP_INIT(potLanczos->YMat, double, potLanczos->maxIter * 2);
    POTLP_INIT(potLanczos->dArray, double, potLanczos->maxIter);
    POTLP_INIT(potLanczos->UMat, double, potLanczos->maxIter * potLanczos->maxIter);
    POTLP_INIT(potLanczos->eiDblMat, double, potLanczos->maxIter * SYEV_WORK);
    POTLP_INIT(potLanczos->eiIntMat, pot_int, potLanczos->maxIter * SYEV_IWORK);
    
    if ( !potLanczos->YMat || !potLanczos->HMat || !potLanczos->YMat ||
         !potLanczos->dArray || !potLanczos->UMat || !potLanczos->eiDblMat ||
         !potLanczos->eiIntMat ) {
        retcode = RETCODE_FAILED;
        goto exit_cleanup;
    }
    
exit_cleanup:
    
    return retcode;
}

extern void potLanczosInitData( pot_lanczos *potLanczos, void *MMat, void (*lczMatVec) (void *, pot_vec *, pot_vec *) ) {
    
    if ( potLanczos->MMat || potLanczos->lczMatVec ) {
        return;
    }
    
    potLanczos->MMat = MMat;
    potLanczos->lczMatVec = lczMatVec;
    
    return;
}

#define LANCZOS_DEBUG(format, info)

/** @brief
 * Implement an abstract Lanczos method to find the negative curvature of the Hessian
 */
extern pot_int potLanczosSolve( pot_lanczos *potLanczos, pot_vec *lczStart, pot_vec *nCurv ) {
    
    pot_int retcode = RETCODE_OK;
    
    pot_lanczos *potl = potLanczos;
    
    pot_vec *vVec = potl->vVec;
    pot_vec *wVec = potl->wVec;
    pot_vec *z1Vec = potl->z1Vec;
    pot_vec *z2Vec = potl->z2Vec;
    
    double  *wVecData = wVec->x;
    double  *z1VecData = z1Vec->x;
    double  *z2VecData = z2Vec->x;
    
    double *VMat = potl->VMat;
    double *HMat = potl->HMat;
    double *YMat = potl->YMat;
    double *UMat = potl->UMat;
    double *dArray = potl->dArray;
    
    double  *eiDblMat = potl->eiDblMat;
    pot_int *eiIntMat = potl->eiIntMat;
    
    int maxIter = potl->maxIter;
    
    pot_int nVRow = vVec->n;
    pot_int nHRow = maxIter + 1;
    
#define H(i, j) HMat[nHRow * (j) + (i)]
#define V(i, j) VMat[nVRow * (j) + (i)]
    
    pot_int ldWork = maxIter * SYEV_WORK;
    pot_int liWork = maxIter * SYEV_IWORK;
    
    /* Initialize */
    if ( lczStart ) {
        potVecCopy(lczStart, vVec);
    } else {
        potLanczosIInitV(vVec);
    }
    
    potVecNormalize(vVec);
    potVecExport(vVec, VMat);
    
    int k;
    int checkFreq = (int) potl->maxIter / 10;
    double finalEigs = 0.0;
    
    /* Start Lanczos iterations */
    for ( k = 0; k < maxIter; ++k ) {
        
        potl->lczMatVec(potl->MMat, vVec, wVec);
        
        double normPrev = nrm2(&nVRow, wVecData, &potIntConstantOne);
        
        if ( k > 0 ) {
            double negHElem = - H(k, k - 1);
            axpy(&nVRow, &negHElem, &V(0, k - 1), &potIntConstantOne, wVecData, &potIntConstantOne);
        }
        
        double vAlp = - dot(&nVRow, wVecData, &potIntConstantOne, &V(0, k), &potIntConstantOne);
        axpy(&nVRow, &vAlp, &V(0, k), &potIntConstantOne, wVecData, &potIntConstantOne);
        double normPres = nrm2(&nVRow, wVecData, &potIntConstantOne);
        
        H(k, k) = - vAlp;
        
        LANCZOS_DEBUG("Lanczos Alp value: %f \n", -vAlp);
        
        /* Refinement */
        if ( normPres < 0.99 * normPrev || 1 ) {
            wVec->nrm = -1.0; // TODO: Not good. Try rewriting this
            // TODO: Use dgemv
            for ( int i = 0; i < k; ++i ) {
                double alpha = - dot(&nVRow, &V(0, i), &potIntConstantOne, wVecData, &potIntConstantOne);
                axpy(&nVRow, &alpha, &V(0, i), &potIntConstantOne, wVecData, &potIntConstantOne);
                H(i, k) -= alpha;
            }
        }
        
        potVecCopy(wVec, vVec);
        normPres = potVecNormalize(vVec);
        potVecExport(vVec, &V(0, k + 1));
        H(k + 1, k) = H(k, k + 1) = normPres;
        
        if ( (k + 1) % checkFreq == 0 || k >= maxIter - 1 ) {
            
            int kPlus1 = k + 1;
            
            for ( int i = 0; i < kPlus1; ++i ) {
                POTLP_MEMCPY(UMat + kPlus1 * i, &H(0, i), double, kPlus1);
            }
            
            dArrSymmetrize(kPlus1, UMat);
            
            POT_CALL(psyev(kPlus1, UMat, dArray, YMat, eiDblMat, eiIntMat, ldWork, liWork));
            
            double resiVal = fabs( H(kPlus1, k) * YMat[kPlus1 + k] );
            LANCZOS_DEBUG("Lanczos outer resi value: %f \n", resiVal);
            
            if ( resiVal < 1e-04 || k >= maxIter - 1 ) {
                
                LANCZOS_DEBUG("Lanczos iteration %d \n", k);
                
                double eigMin1 = dArray[1];
                double eigMin2 = dArray[0];
                
                pgemv(nVRow, kPlus1, VMat, YMat + kPlus1, z1VecData);
                potl->lczMatVec(potl->MMat, z1Vec, z2Vec);
                potVecCopy(z2Vec, nCurv);
                
                double negEig = -eigMin1;
                axpy(&nVRow, &negEig, z1VecData,
                     &potIntConstantOne, z2VecData, &potIntConstantOne);
                
                double resiVal1 = nrm2(&nVRow, z2VecData, &potIntConstantOne);
                
                pgemv(nVRow, kPlus1, VMat, YMat, z2VecData);
                potl->lczMatVec(potl->MMat, z2Vec, z1Vec);
                
                axpy(&nVRow, &negEig, z2VecData,
                     &potIntConstantOne, z1VecData, &potIntConstantOne);
                
                double resiVal2 = nrm2(&nVRow, z1VecData, &potIntConstantOne);
                double resiDiff = eigMin1 - eigMin2 - resiVal2;
                double valGamma = ( resiDiff > 0 ) ? resiDiff : 1e-16;
                
                double resiVal1sqr = resiVal1 * resiVal1 / valGamma;
                double deltaVal = POTLP_MIN(resiVal1, resiVal1sqr);
                
                finalEigs = -eigMin1;
                
                LANCZOS_DEBUG("Lanczos eigenvalue %f \n", finalEigs);
                
                if ( eigMin1 - fabs(deltaVal) > 0 ) {
                    break;
                }
            }
        }
    }
    
    LANCZOS_DEBUG("Lanczos finished in %d iterations \n", k);
    LANCZOS_DEBUG("Lanczos provides eigenvalue of %f \n", finalEigs);
    
    if ( finalEigs > 0.0 || k == maxIter - 1) {
        retcode = RETCODE_FAILED;
    }
    
exit_cleanup:
    return retcode;
}

extern void potLanczosClear( pot_lanczos *potLanczos ) {
    
    if ( !potLanczos ) {
        return;
    }
    
    potLanczos->n = 0;
    potLanczos->maxIter = 0;
    
    potLanczos->MMat = NULL;
    
    potVecDestroy(&potLanczos->vVec);
    potVecDestroy(&potLanczos->wVec);
    potVecDestroy(&potLanczos->z1Vec);
    potVecDestroy(&potLanczos->z2Vec);
    potVecDestroy(&potLanczos->vaVec);
    
    POTLP_FREE(potLanczos->VMat);
    POTLP_FREE(potLanczos->HMat);
    POTLP_FREE(potLanczos->YMat);
    POTLP_FREE(potLanczos->dArray);
    POTLP_FREE(potLanczos->UMat);
    POTLP_FREE(potLanczos->eiDblMat);
    POTLP_FREE(potLanczos->eiIntMat);
    
    potLanczos->lczMatVec = NULL;
    
    return;
}

extern void potLanczosDestroy( pot_lanczos **ppotLanczos ) {
    
    if ( !ppotLanczos ) {
        return;
    }
    
    potLanczosClear(*ppotLanczos);
    POTLP_FREE(*ppotLanczos);
    
    return;
}
