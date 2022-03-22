#ifndef dsdplanczos_h
#define dsdplanczos_h

#include "dsdphsd.h"
#include "structs.h"

// Implementation of Lanczos iteration for SDP stepsize
#define MAXITER 10

typedef struct {
    
    DSDP_INT n;
    vec *v;
    vec *w;
    vec *z1;
    vec *z2;
    vec *vecaux;
    
    double *V;
    double *H;
    double *Y;
    double *d;
    
    double *mataux;
    double *eigaux;
    
    DSDP_INT *eigintaux;
    DSDP_INT isuppz[4];
    
    DSDP_INT iwork;
    DSDP_INT lwork;
    
} DSDPLanczos;

#ifdef __cplusplus
extern "C" {
#endif

extern DSDP_INT dsdpLanczosInit ( DSDPLanczos *lczSolver );
extern DSDP_INT dsdpLanczosAlloc( DSDPLanczos *lczSolver, DSDP_INT n );
extern DSDP_INT dsdpLanczosStep ( DSDPLanczos *lczSolver, spsMat *S, spsMat *dS, double *lbd, double *delta );
extern void     dsdpLanczosFree ( DSDPLanczos *lczSolver );

#ifdef __cplusplus
}
#endif

#endif /* dsdplanczos_h */
