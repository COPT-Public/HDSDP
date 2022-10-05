#ifndef pot_constr_mat_h
#define pot_constr_mat_h

#include "pot_structs.h"

#ifdef __cplusplus
extern "C" {
#endif

extern void potConstrMatPrepareX( pot_constr_mat *potConstrMat, pot_vec *xVec );
extern void potConstrMatProj( pot_constr_mat *potConstrMat, pot_vec *xVec, pot_vec *xVecP );
extern void potConstrMatScalProj( pot_constr_mat *potConstrMat, pot_vec *xVec, pot_vec *yVec, pot_vec *yVecP );
extern void potConstrMatMonitor( pot_constr_mat *potConstrMat, void *info );

#ifdef __cplusplus
}
#endif

#endif /* pot_constr_mat_h */
