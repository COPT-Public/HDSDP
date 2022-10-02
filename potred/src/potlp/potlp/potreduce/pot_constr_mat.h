#ifndef pot_constr_mat_h
#define pot_constr_mat_h

#include "pot_structs.h"

#ifdef __cplusplus
extern "C" {
#endif
extern pot_int potConstrMatInit( pot_constr_mat *potConstrMat, pot_int nRows, pot_int nCols );
extern pot_int potConstrMatInitData( pot_constr_mat *potConstrMat, void *inputData );
extern void potConstrMatPrepareX( pot_constr_mat *potConstrMat, pot_vec *xVec );
extern void potConstrMatProj( pot_constr_mat *potConstrMat, pot_vec *xVec, pot_vec *yVec );
extern void potConstrMatMonitor( pot_constr_mat *potConstrMat );
extern void potConstrMatDestroy( pot_constr_mat *potConstrMat );

#ifdef __cplusplus
}
#endif

#endif /* pot_constr_mat_h */
