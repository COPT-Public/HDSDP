#ifndef def_hdsdp_lanczos_h
#define def_hdsdp_lanczos_h

#include "interface/hdsdp.h"

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
    
    double *dArray;
    double *eigDblMat;
    int    *eigIntMat;
    
    /* Matrix vector multiplication */
    void (*Mvec) (void *, double *, double *);
    
} hdsdp_lanczos;

#endif /* def_hdsdp_lanczos_h */
