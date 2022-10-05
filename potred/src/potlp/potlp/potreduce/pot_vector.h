#ifndef pot_vector_h
#define pot_vector_h

#include "pot_structs.h"

#ifdef __cplusplus
extern "C" {
#endif

extern pot_int potVecCreate( pot_vec **ppVec );
extern pot_int potVecInit( pot_vec *pVec, pot_int vDim, pot_int vConeDim );
extern void potVecDiff( pot_vec *pVecOut, pot_vec *pVecInPrev, pot_vec *pVecInPres );
extern void potVeczAxpby( pot_vec *pVecZ, double alpha, pot_vec *pVecX, double beta, pot_vec *pVecY );
extern double potVecNormalize( pot_vec *pVec );
extern double potVecScaledDot( pot_vec *pVecX, pot_vec *pVecY, pot_vec *pVecZ );
extern double potVecDot( pot_vec *pVecX, pot_vec *pVecY );
extern void potVecAxinvpBy( double alpha1, double alpha2, pot_vec *pVecX, double beta, pot_vec *pVecY );
extern void potVecAxpy( double alpha, pot_vec *pVecX, pot_vec *pVecY );
extern void potVecConeAxpy( double alpha, pot_vec *pVecX, pot_vec *pVecY );
extern double potVecLogDet( pot_vec *pVecX );
extern double potVecSumCone( pot_vec *pVecX );
extern double potVecSumScalCone( pot_vec *pVecX, pot_vec *pVecY );
extern void potVecConeAddConstant( pot_vec *pVecX, double cVal );
extern void potVecClear( pot_vec *pVec );
extern void potVecDestroy( pot_vec **pVec );
extern void potVecCopy( pot_vec *srcVec, pot_vec *dstVec );
extern void potVecReset( pot_vec *pVec );

#ifdef __cplusplus
}
#endif

#endif /* pot_vector_h */
