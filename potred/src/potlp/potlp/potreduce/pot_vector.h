#ifndef pot_vector_h
#define pot_vector_h

#include "pot_structs.h"

extern pot_int potVecInit( pot_vec *pVec, pot_int vDim, pot_int vConeDim );
extern void potVecDestroy( pot_vec *pVec );
extern void potVecCopy( pot_vec *srcVec, pot_vec *dstVec );

#endif /* pot_vector_h */
