#ifdef HEADERPATH
#include "linalg/def_hdsdp_mineig.h"
#include "linalg/hdsdp_mineig.h"
#include "linalg/hdsdp_lanczos.h"
#include "linalg/sparse_opts.h"
#include "linalg/vec_opts.h"
#include "linalg/dense_opts.h"
#include "interface/hdsdp_utils.h"
#else
#include "def_hdsdp_mineig.h"
#include "hdsdp_mineig.h"
#include "hdsdp_lanczos.h"
#include "sparse_opts.h"
#include "vec_opts.h"
#include "dense_opts.h"
#include "hdsdp_utils.h"
#endif

#include <math.h>

#ifndef SYEV_WORK
#define SYEV_WORK  (30)
#endif

#ifndef SYEV_IWORK
#define SYEV_IWORK (12)
#endif

/** @brief Symmetrize an array
 */
static void dArrSymmetrize( int nCol, double *dArray ) {
    
    double Aij, Aji;
    for ( int i = 0, j; i < nCol; ++i ) {
        for ( j = i + 1; j < nCol; ++j ) {
            Aij = dArray[j * nCol + i];
            Aji = dArray[i * nCol + j];
            dArray[j * nCol + i] = dArray[i * nCol + j] = (Aij + Aji) * 0.5;
        }
    }
    return;
}

static void HMinEigIGetStart( int nCol, double *vVec ) {
    
    srand((unsigned int) nCol);
    for ( int i = 0; i < nCol; ++i ) {
        srand((unsigned int) rand());
        vVec[i] = sqrt(sqrt((rand() % 1627))) * (rand() % 2 - 0.5);
        vVec[i] = 1.0;
    }
    
    return;
}

static void HMinEigIMatVec( hdsdp_mineig *eig, double *dXVec, double *dMatXVec ) {
    
    HDSDP_ZERO(dMatXVec, double, eig->nCol);
    
    const int *Ap = eig->colMatBeg;
    const int *Ai = eig->colMatIdx;
    const double *Ax = eig->colMatElem;
    double *x = dXVec;
    double *y = dMatXVec;
    
    for ( int i = 0, j; i < eig->nCol; ++i ) {
        y[Ai[Ap[i]]] -= x[i] * Ax[Ap[i]];
        for ( j = Ap[i] + 1; j < Ap[i + 1]; ++j ) {
            y[Ai[j]] -= x[i] * Ax[j];
            y[i] -= x[Ai[j]] * Ax[j];
        }
    }
    
    return;
}

#ifdef HDSDP_LANCZOS_DEBUG
#undef HDSDP_LANCZOS_DEBUG
#define HDSDP_LANCZOS_DEBUG(format, info) printf(format, info)
#else
#define HDSDP_LANCZOS_DEBUG(format, info)
#endif
#define H(i, j) eig->HMat[nHRow * (j) + (i)]
#define V(i, j) eig->VMat[nVRow * (j) + (i)]
/** @brief Apply Lanczos iteration to approximately find the minimum eigenvalue of a symmertic matrix A
 * The implementation is adopted from SDPT3 https://github.com/sqlp/sdpt3/blob/master/Solver/steplength.m
 */
static hdsdp_retcode HMinEigILanczos( hdsdp_mineig *eig, double *dPerturb ) {
    
    hdsdp_retcode retcode = HDSDP_RETCODE_OK;
    
    /* Get initial vector */
    HMinEigIGetStart(eig->nCol, eig->vVec);
    
    /* Normalize the starting vector use it as the initial basis */
    normalize(&eig->nCol, eig->vVec);
    
    /* Copy the basis to VMat */
    HDSDP_MEMCPY(eig->VMat, eig->vVec, double, eig->nCol);
    
    /* Configure and start Lanczos iteration */
    int k = 0;
    int LCheckFrequency = (int) eig->nKrylovDim / 5;
    
    int nVRow = eig->nCol;
    int nHRow = eig->nKrylovDim + 1;
    int ldWork = eig->nKrylovDim * SYEV_WORK;
    int liWork = eig->nKrylovDim * SYEV_IWORK;
    
    for ( k = 0; k < eig->nKrylovDim; ++k ) {
            
        HMinEigIMatVec(eig, eig->vVec, eig->wVec);
            
        if ( k > 0 ) {
            double negHElem = - H(k, k - 1);
            axpy(&nVRow, &negHElem, &V(0, k - 1), &HIntConstantOne, eig->wVec, &HIntConstantOne);
        }
            
        double vAlp = -dot(&nVRow, eig->wVec, &HIntConstantOne, &V(0, k), &HIntConstantOne);
        axpy(&nVRow, &vAlp, &V(0, k), &HIntConstantOne, eig->wVec, &HIntConstantOne);
        double normPres = nrm2(&nVRow, eig->wVec, &HIntConstantOne);
            
        H(k, k) = - vAlp;
        HDSDP_LANCZOS_DEBUG("Lanczos Alp value: %f \n", -vAlp);
            
        HDSDP_MEMCPY(eig->vVec, eig->wVec, double, eig->nCol);
        normPres = normalize(&eig->nCol, eig->vVec);
            
        if ( normPres > 0.0 ) {
            HDSDP_MEMCPY(&V(0, k + 1), eig->vVec, double, eig->nCol);
            H(k + 1, k) = H(k, k + 1) = normPres;
        }
        
        /* Frequently check subspace */
        if ( ( k + 1 ) % LCheckFrequency == 0 || k > eig->nKrylovDim - 1 || normPres == 0.0 ) {
            
            HDSDP_LANCZOS_DEBUG("Entering Lanczos internal check at iteration %d.\n", k);
            
            int kPlus1 = k + 1;
            for ( int i = 0; i < kPlus1; ++i ) {
                HDSDP_MEMCPY(eig->UMat + kPlus1 * i, &H(0, i), double, kPlus1);
            }
            
            dArrSymmetrize(kPlus1, eig->UMat);
            HDSDP_CALL(fds_syev(kPlus1, eig->UMat, eig->dArray, eig->YMat, 2,
                                eig->eigDblMat, eig->eigIntMat, ldWork, liWork));
            
            double resiVal = fabs( H(kPlus1, k) * eig->YMat[kPlus1 + k] );
            HDSDP_LANCZOS_DEBUG("Lanczos outer resi value: %f \n", resiVal);
            
            if ( resiVal < 1e-04 || k >= eig->nKrylovDim - 1 ) {
                
                HDSDP_LANCZOS_DEBUG("Lanczos inner iteration %d \n", k);
                
                double eigMin1 = eig->dArray[1];
                double eigMin2 = eig->dArray[0];
                
                fds_gemv(nVRow, kPlus1, eig->VMat,
                            eig->YMat + kPlus1, eig->z1Vec);
                HMinEigIMatVec(eig, eig->z1Vec, eig->z2Vec);
                
                double negEig = -eigMin1;
                axpy(&nVRow, &negEig, eig->z1Vec, &HIntConstantOne, eig->z2Vec, &HIntConstantOne);
                
                double resiVal1 = nrm2(&nVRow, eig->z2Vec, &HIntConstantOne);
                fds_gemv(nVRow, kPlus1, eig->VMat, eig->YMat, eig->z2Vec);
                HMinEigIMatVec(eig, eig->z2Vec, eig->z1Vec);
                axpy(&nVRow, &negEig, eig->z2Vec, &HIntConstantOne, eig->z1Vec, &HIntConstantOne);
                    
                /* Compute bound on the stepsize */
                double resiVal2 = nrm2(&nVRow, eig->z1Vec, &HIntConstantOne);
                double resiDiff = eigMin1 - eigMin2 - resiVal2;
                double valGamma = ( resiDiff > 0 ) ? resiDiff : 1e-16;
                double resiVal1sqr = resiVal1 * resiVal1 / valGamma;
                valGamma = HDSDP_MIN(resiVal1, resiVal1sqr);
                    
                if ( valGamma < 1e-03 ) {
                    *dPerturb = - (valGamma + eigMin1);
                    break;
                } else {
                    
                    if ( normPres == 0.0 ) {
                        retcode = HDSDP_RETCODE_FAILED;
                        goto exit_cleanup;
                    }
                    
                    *dPerturb = valGamma + eigMin1;
                }
            }
        }
    }
    
exit_cleanup:
    return retcode;
}


/** @brief Overwrite the diagonal by the original matrix plus multiple of identity and try factorizing the perturbed matrix
 *
 */
static hdsdp_retcode HMinEigPsdCheck( hdsdp_mineig *eig, double dPerturb, int *isPsd ) {
    
    hdsdp_retcode retcode = HDSDP_RETCODE_OK;
    
    /* Overwrite the diagonal by the original matrix plus multiple of identity */
    for ( int iCol = 0; iCol < eig->nCol; ++iCol ) {
        eig->colMatElemPlusEye[eig->colMatBeg[iCol]] = eig->colMatElem[eig->colMatBeg[iCol]] + dPerturb;
    }
    
    /* Check PSD */
    HDSDP_CALL(eig->ispsd(eig->chol, eig->colMatBeg, eig->colMatIdx, eig->colMatElemPlusEye, isPsd));
    
exit_cleanup:
    return retcode;
}


#ifndef NO_BISECTION
//#define NO_BISECTION
#endif
/** @brief Apply bisection to find a perturbation so that the perturbed matrix becomes positive semidefinite
 *
 */
static hdsdp_retcode HMinEigIBisection( hdsdp_mineig *eig, double dPerturbGuess, double *dPerturb ) {
    
    hdsdp_retcode retcode = HDSDP_RETCODE_OK;
    *dPerturb = dPerturbGuess;
#ifdef NO_BISECTION
    return retcode;
#endif
    
    int isPsd = 0;
    
    double dPerturbUpper = -dPerturbGuess * 1.01;
    double dPerturbLower = 0.0;
    double dBisecTol = 1e-12;
    
    /* Step 1. Finding an upper bound of perturbation starting from the guess */
    while ( !isPsd ) {
        
        HDSDP_CALL(HMinEigPsdCheck(eig, dPerturbUpper, &isPsd));
        if ( !isPsd ) {
            dPerturbLower = dPerturbUpper;
        } else {
            break;
        }
        
        if ( dPerturbUpper < 1e-05 ) {
            dPerturbUpper *= 2.0;
        } else {
            dPerturbUpper *= 10.0;
        }
    }
    
    /* Step 2. Apply bisection to find the minimum perturbation */
    while ( 1 ) {
        
        double dPerturbMid = dPerturbLower + 0.5 * (dPerturbUpper - dPerturbLower);
        HDSDP_CALL(HMinEigPsdCheck(eig, dPerturbMid, &isPsd));
        
        if ( isPsd ) {
            dPerturbUpper = dPerturbMid;
        } else {
            dPerturbLower = dPerturbMid;
        }
        
        if ( fabs(dPerturbUpper - dPerturbLower) <= dBisecTol * (fabs(dPerturbUpper) + fabs(dPerturbUpper) + 1.0) ) {
            *dPerturb = dPerturbMid;
            break;
        }
    }
    
exit_cleanup:
    return retcode;
}


extern hdsdp_retcode HMinEigCreate( hdsdp_mineig **peig ) {
    
    hdsdp_retcode retcode = HDSDP_RETCODE_OK;
    
    if ( !peig ) {
        retcode = HDSDP_RETCODE_FAILED;
        goto exit_cleanup;
    }
    
    hdsdp_mineig *eig = NULL;
    HDSDP_INIT(eig, hdsdp_mineig, 1);
    HDSDP_MEMCHECK(eig);
    HDSDP_ZERO(eig, hdsdp_mineig, 1);
    
    *peig = eig;
    
exit_cleanup:
    return retcode;
}

extern hdsdp_retcode HMinEigInit( hdsdp_mineig *eig, int nCol, int nKrylovDim ) {
    
    hdsdp_retcode retcode = HDSDP_RETCODE_OK;
    
    if ( nCol <= 0 || nKrylovDim <= 0 ) {
        retcode = HDSDP_RETCODE_FAILED;
        goto exit_cleanup;
    }
    
    /* Special case of 1D matrix */
    if ( nCol == 1 ) {
        goto exit_cleanup;
    }
    
    eig->nCol = nCol;
    eig->nKrylovDim = nKrylovDim;
    
    /* Allocate Lanczos working arrays */
    HDSDP_INIT(eig->vVec, double, nCol);
    HDSDP_MEMCHECK(eig->vVec);

    HDSDP_INIT(eig->wVec, double, nCol);
    HDSDP_MEMCHECK(eig->wVec);

    HDSDP_INIT(eig->z1Vec, double, nCol);
    HDSDP_MEMCHECK(eig->z1Vec);

    HDSDP_INIT(eig->z2Vec, double, nCol);
    HDSDP_MEMCHECK(eig->z2Vec);

    HDSDP_INIT(eig->vaVec, double, nCol);
    HDSDP_MEMCHECK(eig->vaVec);

    HDSDP_INIT(eig->VMat, double, nCol * (eig->nKrylovDim + 1));
    HDSDP_MEMCHECK(eig->VMat);

    HDSDP_INIT(eig->HMat, double, (eig->nKrylovDim + 1) * (eig->nKrylovDim + 1));
    HDSDP_MEMCHECK(eig->HMat);

    HDSDP_INIT(eig->YMat, double, eig->nKrylovDim * 2);
    HDSDP_MEMCHECK(eig->YMat);

    HDSDP_INIT(eig->dArray, double, eig->nKrylovDim);
    HDSDP_MEMCHECK(eig->dArray);

    HDSDP_INIT(eig->UMat, double, eig->nKrylovDim * eig->nKrylovDim);
    HDSDP_MEMCHECK(eig->UMat);

    HDSDP_INIT(eig->eigDblMat, double, eig->nKrylovDim * SYEV_WORK);
    HDSDP_MEMCHECK(eig->eigDblMat);

    HDSDP_INIT(eig->eigIntMat, int, eig->nKrylovDim * SYEV_IWORK);
    HDSDP_MEMCHECK(eig->eigIntMat);
    
exit_cleanup:
    return retcode;
}

extern hdsdp_retcode HMinEigComputeMinEig( hdsdp_mineig *eig, int *colMatBeg, int *colMatIdx, double *colMatElem, void *chol,
                                          int (*ispsd) ( void *, int *, int *, double *, int * ), mineig_algo algo, double *dPerturb ) {
    
    hdsdp_retcode retcode = HDSDP_RETCODE_OK;
    
    if ( !colMatBeg || !colMatIdx || !colMatElem ) {
        retcode = HDSDP_RETCODE_FAILED;
        goto exit_cleanup;
    }
    
    /* Copy pointers */
    eig->colMatBeg = colMatBeg;
    eig->colMatIdx = colMatIdx;
    eig->colMatElem = colMatElem;
    eig->nElem = colMatBeg[eig->nCol];
    
    /* Allocate memory for trial decomposition */
    HDSDP_INIT(eig->colMatElemPlusEye, double, eig->nElem);
    HDSDP_MEMCHECK(eig->nElem);
    HDSDP_MEMCPY(eig->colMatElemPlusEye, eig->colMatElem, double, eig->nElem);
    
    eig->chol = chol;
    eig->ispsd = ispsd;
    
    double dPerturbGuess = 0.0;
    
    if ( algo == MIN_EIG_LANCZOS ) {
        HDSDP_CALL(HMinEigILanczos(eig, &dPerturbGuess));
    }
    
    HDSDP_CALL(HMinEigIBisection(eig, dPerturbGuess, dPerturb));
    printf("Minimum eigenvalue: %5.3e \n", *dPerturb);
    
exit_cleanup:
    return retcode;
}

extern void HMinEigClear( hdsdp_mineig *eig ) {
    
    if ( !eig ) {
       return;
    }
    
    HDSDP_FREE(eig->colMatElemPlusEye);

    HDSDP_FREE(eig->vVec);
    HDSDP_FREE(eig->wVec);
    HDSDP_FREE(eig->z1Vec);
    HDSDP_FREE(eig->z2Vec);
    HDSDP_FREE(eig->vaVec);

    HDSDP_FREE(eig->VMat);
    HDSDP_FREE(eig->HMat);
    HDSDP_FREE(eig->YMat);
    HDSDP_FREE(eig->UMat);

    HDSDP_FREE(eig->dArray);
    HDSDP_FREE(eig->eigDblMat);
    HDSDP_FREE(eig->eigIntMat);

    HDSDP_ZERO(eig, hdsdp_lanczos, 1);
    
    return;
}

extern void HMinEigDestroy( hdsdp_mineig **peig ) {
    
    if ( !peig ) {
        return;
    }
    
    HMinEigClear(*peig);
    HDSDP_FREE(*peig);
    
    return;
}
