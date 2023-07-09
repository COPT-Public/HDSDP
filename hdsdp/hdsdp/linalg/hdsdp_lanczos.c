#ifdef HEADERPATH
#include "linalg/def_hdsdp_lanczos.h"
#include "linalg/hdsdp_lanczos.h"
#include "interface/hdsdp_utils.h"
#else
#include "def_hdsdp_lanczos.h"
#include "hdsdp_lanczos.h"
#include "hdsdp_utils.h"
#endif

#include <math.h>

#define SYEV_WORK 30
#define SYEV_IWORK 12

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

extern hdsdp_retcode HLanczosSolve( hdsdp_lanczos *HLanczos, double *dMinEVal ) {
    
    hdsdp_retcode retcode = HDSDP_RETCODE_OK;
    
    /* TODO: Move Lanczos implementation from LP here */
    
    
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
