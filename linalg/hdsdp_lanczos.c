#ifdef HEADERPATH
#include "linalg/def_hdsdp_lanczos.h"
#include "linalg/hdsdp_lanczos.h"
#include "linalg/vec_opts.h"
#include "linalg/dense_opts.h"
#include "interface/hdsdp_utils.h"
#else
#include "def_hdsdp_lanczos.h"
#include "hdsdp_lanczos.h"
#include "vec_opts.h"
#include "dense_opts.h"
#include "hdsdp_utils.h"
#endif

#include <math.h>

#define SYEV_WORK  (30)
#define SYEV_IWORK (12)

static void dArrSymmetrize( int n, double *dArray ) {
    
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

static void HLanczosIPrepare( int n, double *vVec ) {
    
    srand(n);
    for ( int i = 0; i < n; ++i ) {
        srand(rand());
        vVec[i] = sqrt(sqrt((rand() % 1627))) * (rand() % 2 - 0.5);
    }
    
    return;
}

static void HLanczosICleanUp( hdsdp_lanczos *HLanczos ) {
    
    HDSDP_ZERO(HLanczos->wVec, double, HLanczos->nCol);
    HDSDP_ZERO(HLanczos->z1Vec, double, HLanczos->nCol);
    HDSDP_ZERO(HLanczos->z2Vec, double, HLanczos->nCol);
    HDSDP_ZERO(HLanczos->vaVec, double, HLanczos->nCol);
    HDSDP_ZERO(HLanczos->VMat, double, HLanczos->nCol * (HLanczos->nMaxSpaceDim + 1));
    HDSDP_ZERO(HLanczos->HMat, double, (HLanczos->nMaxSpaceDim + 1) * (HLanczos->nMaxSpaceDim + 1));
    HDSDP_ZERO(HLanczos->YMat, double, HLanczos->nMaxSpaceDim * 2);
    HDSDP_ZERO(HLanczos->dArray, double, HLanczos->nMaxSpaceDim);
    HDSDP_ZERO(HLanczos->UMat, double, HLanczos->nMaxSpaceDim * HLanczos->nMaxSpaceDim);
    HDSDP_ZERO(HLanczos->eigDblMat, double, HLanczos->nMaxSpaceDim * SYEV_WORK);
    HDSDP_ZERO(HLanczos->eigIntMat, int, HLanczos->nMaxSpaceDim * SYEV_IWORK);
    
    return;
}

extern hdsdp_retcode HLanczosCreate( hdsdp_lanczos **pHLanczos ) {
    
    hdsdp_retcode retcode = HDSDP_RETCODE_OK;
    
    if ( !pHLanczos ) {
        retcode = HDSDP_RETCODE_FAILED;
        goto exit_cleanup;
    }
    
    hdsdp_lanczos *HLanczos = NULL;
    HDSDP_INIT(HLanczos, hdsdp_lanczos, 1);
    HDSDP_MEMCHECK(HLanczos);
    HDSDP_ZERO(HLanczos, hdsdp_lanczos, 1);
    *pHLanczos = HLanczos;
    
exit_cleanup:
    return retcode;
}

extern hdsdp_retcode HLanczosInit( hdsdp_lanczos *HLanczos, int nCol, int nSpaceDim ) {
    
    hdsdp_retcode retcode = HDSDP_RETCODE_OK;
    
    HLanczos->nCol = nCol;
    HLanczos->nMaxSpaceDim = nSpaceDim;
    
    HDSDP_INIT(HLanczos->vVec, double, nCol);
    HDSDP_MEMCHECK(HLanczos->vVec);
    
    HDSDP_INIT(HLanczos->wVec, double, nCol);
    HDSDP_MEMCHECK(HLanczos->wVec);
    
    HDSDP_INIT(HLanczos->z1Vec, double, nCol);
    HDSDP_MEMCHECK(HLanczos->z1Vec);
    
    HDSDP_INIT(HLanczos->z2Vec, double, nCol);
    HDSDP_MEMCHECK(HLanczos->z2Vec);
    
    HDSDP_INIT(HLanczos->vaVec, double, nCol);
    HDSDP_MEMCHECK(HLanczos->vaVec);
    
    HDSDP_INIT(HLanczos->VMat, double, nCol * (HLanczos->nMaxSpaceDim + 1));
    HDSDP_MEMCHECK(HLanczos->VMat);
    
    HDSDP_INIT(HLanczos->HMat, double, (HLanczos->nMaxSpaceDim + 1) * (HLanczos->nMaxSpaceDim + 1));
    HDSDP_MEMCHECK(HLanczos->HMat);
    
    HDSDP_INIT(HLanczos->YMat, double, HLanczos->nMaxSpaceDim * 2);
    HDSDP_MEMCHECK(HLanczos->YMat);
    
    HDSDP_INIT(HLanczos->dArray, double, HLanczos->nMaxSpaceDim);
    HDSDP_MEMCHECK(HLanczos->dArray);
    
    HDSDP_INIT(HLanczos->dLanczosWarmStart, double, HLanczos->nCol);
    HDSDP_MEMCHECK(HLanczos->dLanczosWarmStart);
    
    HDSDP_INIT(HLanczos->UMat, double, HLanczos->nMaxSpaceDim * HLanczos->nMaxSpaceDim);
    HDSDP_MEMCHECK(HLanczos->UMat);
    
    HDSDP_INIT(HLanczos->eigDblMat, double, HLanczos->nMaxSpaceDim * SYEV_WORK);
    HDSDP_MEMCHECK(HLanczos->eigDblMat);
    
    HDSDP_INIT(HLanczos->eigIntMat, int, HLanczos->nMaxSpaceDim * SYEV_IWORK);
    HDSDP_MEMCHECK(HLanczos->eigIntMat);
    
exit_cleanup:
    return retcode;
}

extern void HLanczosSetData( hdsdp_lanczos *HLanczos, void *MMat, void (*Mvec) (void *, double *, double *) ) {
    
    if ( HLanczos->MMat || HLanczos->Mvec ) {
        return;
    }
    
    HLanczos->MMat = MMat;
    HLanczos->Mvec = Mvec;
    
    return;
}

#ifdef HDSDP_LANCZOS_DEBUG
#undef HDSDP_LANCZOS_DEBUG
#define HDSDP_LANCZOS_DEBUG(format, info) printf(format, info)
#else
#define HDSDP_LANCZOS_DEBUG(format, info)
#endif
#define H(i, j) HLanczos->HMat[nHRow * (j) + (i)]
#define V(i, j) HLanczos->VMat[nVRow * (j) + (i)]
extern hdsdp_retcode HLanczosSolve( hdsdp_lanczos *HLanczos, double *LanczosStart, double *dMaxStep ) {
    
    hdsdp_retcode retcode = HDSDP_RETCODE_OK;
    
    /* Prepare the initial point */
    if ( LanczosStart ) {
        HDSDP_LANCZOS_DEBUG("Loaded user Lanczos warm-start. %s\n", "");
        HDSDP_MEMCPY(HLanczos->vVec, LanczosStart, double, HLanczos->nCol);
    } else {
        if ( HLanczos->nComputed == 0 ) {
            HDSDP_LANCZOS_DEBUG("Starting Lanczos from scratch. %s\n", "");
            HLanczosIPrepare(HLanczos->nCol, HLanczos->vVec);
        } else {
            HLanczosICleanUp(HLanczos);
            HDSDP_LANCZOS_DEBUG("Loaded Lanczos warm start from previous iteration. %s", "");
            HDSDP_MEMCPY(HLanczos->vVec, HLanczos->dLanczosWarmStart, double, HLanczos->nCol);
        }
    }
    
    /* Normalize the starting vector use it as the initial basis */
    normalize(&HLanczos->nCol, HLanczos->vVec);
    HDSDP_MEMCPY(HLanczos->VMat, HLanczos->vVec, double, HLanczos->nCol);
    
    /* Configure and start Lanczos iteration */
    int k = 0;
    int LCheckFrequency = (int) HLanczos->nMaxSpaceDim / 5;
    LCheckFrequency = HDSDP_MIN(LCheckFrequency, 5);
    
    int nVRow = HLanczos->nCol;
    int nHRow = HLanczos->nMaxSpaceDim + 1;
    
    int ldWork = HLanczos->nMaxSpaceDim * SYEV_WORK;
    int liWork = HLanczos->nMaxSpaceDim * SYEV_IWORK;
    
    for ( k = 0; k < HLanczos->nMaxSpaceDim; ++k ) {
        
        HLanczos->Mvec(HLanczos->MMat, HLanczos->vVec, HLanczos->wVec);
        // double normPrev = nrm2(&HLanczos->nCol, HLanczos->wVec, &HIntConstantOne);
        
        if ( k > 0 ) {
            double negHElem = - H(k, k - 1);
            axpy(&nVRow, &negHElem, &V(0, k - 1), &HIntConstantOne, HLanczos->wVec, &HIntConstantOne);
        }
        
        double vAlp = -dot(&nVRow, HLanczos->wVec, &HIntConstantOne, &V(0, k), &HIntConstantOne);
        axpy(&nVRow, &vAlp, &V(0, k), &HIntConstantOne, HLanczos->wVec, &HIntConstantOne);
        double normPres = nrm2(&nVRow, HLanczos->wVec, &HIntConstantOne);
        
        H(k, k) = - vAlp;
        HDSDP_LANCZOS_DEBUG("Lanczos Alp value: %f \n", -vAlp);
        
        HDSDP_MEMCPY(HLanczos->vVec, HLanczos->wVec, double, HLanczos->nCol);
        normPres = normalize(&HLanczos->nCol, HLanczos->vVec);
        HDSDP_MEMCPY(&V(0, k + 1), HLanczos->vVec, double, HLanczos->nCol);
        H(k + 1, k) = H(k, k + 1) = normPres;
        
        /* Frequently check subspace */
        if ( ( k + 1 ) % LCheckFrequency == 0 || k > HLanczos->nMaxSpaceDim - 1 ) {
            
            HDSDP_LANCZOS_DEBUG("Entering Lanczos internal check at iteration %d.\n", k);
            
            int kPlus1 = k + 1;
            for ( int i = 0; i < kPlus1; ++i ) {
                HDSDP_MEMCPY(HLanczos->UMat + kPlus1 * i, &H(0, i), double, kPlus1);
            }
            
            dArrSymmetrize(kPlus1, HLanczos->UMat);
            HDSDP_CALL(fds_syev(kPlus1, HLanczos->UMat, HLanczos->dArray, HLanczos->YMat,
                                HLanczos->eigDblMat, HLanczos->eigIntMat, ldWork, liWork));
            
            double resiVal = fabs( H(kPlus1, k) * HLanczos->YMat[kPlus1 + k] );
            HDSDP_LANCZOS_DEBUG("Lanczos outer resi value: %f \n", resiVal);
            
            if ( resiVal < 1e-04 || k >= HLanczos->nMaxSpaceDim - 1 ) {
                HDSDP_LANCZOS_DEBUG("Lanczos inner iteration %d \n", k);
                
                double eigMin1 = HLanczos->dArray[1];
                double eigMin2 = HLanczos->dArray[0];
                
                fds_gemv(nVRow, kPlus1, HLanczos->VMat,
                         HLanczos->YMat + kPlus1, HLanczos->z1Vec);
                HLanczos->Mvec(HLanczos->MMat, HLanczos->z1Vec, HLanczos->z2Vec);
                
                /* Record warm-start */
                HDSDP_MEMCPY(HLanczos->dLanczosWarmStart, HLanczos->z2Vec, double, HLanczos->nCol);
                
                double negEig = -eigMin1;
                axpy(&nVRow, &negEig, HLanczos->z1Vec,
                     &HIntConstantOne, HLanczos->z2Vec, &HIntConstantOne);
                
                double resiVal1 = nrm2(&nVRow, HLanczos->z2Vec, &HIntConstantOne);
                fds_gemv(nVRow, kPlus1, HLanczos->VMat,
                         HLanczos->YMat, HLanczos->z2Vec);
                HLanczos->Mvec(HLanczos->MMat, HLanczos->z2Vec, HLanczos->z1Vec);
                axpy(&nVRow, &negEig, HLanczos->z2Vec,
                     &HIntConstantOne, HLanczos->z1Vec, &HIntConstantOne);
                
                /* Compute bound on the stepsize */
                double resiVal2 = nrm2(&nVRow, HLanczos->z1Vec, &HIntConstantOne);
                double resiDiff = eigMin1 - eigMin2 - resiVal2;
                double valGamma = ( resiDiff > 0 ) ? resiDiff : 1e-16;
                double resiVal1sqr = resiVal1 * resiVal1 / valGamma;
                valGamma = HDSDP_MIN(resiVal1, resiVal1sqr);
                
                if ( valGamma < 1e-03 || valGamma + eigMin1 <= 0.5 ) {
                    if ( valGamma + eigMin1 <= 0.0 ) {
                        *dMaxStep = HDSDP_INFINITY;
                    } else {
                        *dMaxStep = 1.0 / ( valGamma + eigMin1 );
                    }
                    break;
                    
                } else {
                    if ( normPres == 0.0 ) {
                        retcode = HDSDP_RETCODE_FAILED;
                        goto exit_cleanup;
                    }
                    *dMaxStep = 1.0 / ( valGamma + eigMin1 );
                }
            }
        }
    }
    
    HLanczos->nComputed += 1;
    
exit_cleanup:
    return retcode;
}

extern void HLanczosClear( hdsdp_lanczos *HLanczos ) {
    
    if ( !HLanczos ) {
        return;
    }
    
    HDSDP_FREE(HLanczos->vVec);
    HDSDP_FREE(HLanczos->wVec);
    HDSDP_FREE(HLanczos->z1Vec);
    HDSDP_FREE(HLanczos->z2Vec);
    HDSDP_FREE(HLanczos->vaVec);
    
    HDSDP_FREE(HLanczos->VMat);
    HDSDP_FREE(HLanczos->HMat);
    HDSDP_FREE(HLanczos->YMat);
    HDSDP_FREE(HLanczos->UMat);
    
    HDSDP_FREE(HLanczos->dArray);
    HDSDP_FREE(HLanczos->dLanczosWarmStart);
    HDSDP_FREE(HLanczos->eigDblMat);
    HDSDP_FREE(HLanczos->eigIntMat);
    
    HDSDP_ZERO(HLanczos, hdsdp_lanczos, 1);
    
    return;
}

extern void HLanczosDestroy( hdsdp_lanczos **pHLanczos ) {
    
    if ( !pHLanczos ) {
        return;
    }
    
    HLanczosClear(*pHLanczos);
    HDSDP_FREE(*pHLanczos);
    
    return;
}
