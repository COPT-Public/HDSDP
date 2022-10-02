#ifndef pot_objfunc_h
#define pot_objfunc_h

#include "pot_structs.h"


#ifdef __cplusplus
extern "C" {
#endif

extern pot_int potObjFInit( pot_fx *potObjF, pot_int nCols );
extern void potObjFSetX( pot_fx *potObjF, pot_vec *xVec );
extern pot_int potObjFInitData( pot_fx *potObjF, void *inputData );
extern double potObjFVal( pot_fx *potObjF, pot_vec *xVec );
extern void potObjFGrad( pot_fx *potObjF, pot_vec *xVec, pot_vec *fGrad );
extern void potObjFHess( pot_fx *potObjF, pot_vec *xVec, double *fHess );
extern void potObjFHVec( pot_fx *potObjF, pot_vec *xVec, pot_vec *fHVec );
extern void potObjFMonitor( pot_fx *potObjF, pot_int *info );
extern void potObjFDestroy( pot_fx *potObjF );

#ifdef __cplusplus
}
#endif

#endif /* pot_objfunc_h */
