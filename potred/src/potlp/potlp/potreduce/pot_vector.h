#ifndef pot_vector_h
#define pot_vector_h

#include "pot_structs.h"

#ifdef __cplusplus
extern "C" {
#endif

extern pot_int potVecInit( pot_vec *pVec, pot_int vDim, pot_int vConeDim );
extern void potVecClear( pot_vec *pVec );
extern void potVecDestroy( pot_vec **pVec );
extern void potVecCopy( pot_vec *srcVec, pot_vec *dstVec );

#ifdef __cplusplus
}
#endif

#endif /* pot_vector_h */
