#ifndef def_hdsdp_lanczos_h
#define def_hdsdp_lanczos_h

#ifdef HEADERPATH
#include "interface/hdsdp.h"
#else
#include "hdsdp.h"
#endif

/* Define the Lanczos data structure for HDSDP */
typedef struct {
    
    /* Matrix dimension */
    int nCol;
    
    /* Maximum subspace dimension*/
    int nMaxSpaceDim;
    
    /* Matrix data */
    void *MMat;
    
    /* Auxiliary data */
    double *vVec;
    double *wVec;
    double *z1Vec;
    double *z2Vec;
    double *vaVec;
    
    double *VMat;
    double *HMat;
    double *YMat;
    double *UMat;
    
    double *dLanczosWarmStart;
    double *dArray;
    double *eigDblMat;
    int    *eigIntMat;
    
    /* Matrix vector multiplication */
    void (*Mvec) (void *, double *, double *);
    
    /* Statistics */
    int nComputed;
    
} hdsdp_lanczos;

#endif /* def_hdsdp_lanczos_h */
