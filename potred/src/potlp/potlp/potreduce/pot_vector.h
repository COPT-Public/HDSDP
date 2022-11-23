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
extern void potVecScal( pot_vec *pVexX, double sVal );
extern void potVecArrScal( pot_vec *pVecX, double *dArray );
extern double potVecScaledDot( pot_vec *pVecX, pot_vec *pVecY, pot_vec *pVecZ );
extern void potVecConeNormalize( pot_vec *pVec );
extern double potVecDot( pot_vec *pVecX, pot_vec *pVecY );
extern void potVecAxinvpBy( double alpha1, double alpha2, pot_vec *pVecX, double beta, pot_vec *pVecY );
extern void potVecAxpy( double alpha, pot_vec *pVecX, pot_vec *pVecY );
extern void potVecConeAxpy( double alpha, pot_vec *pVecX, pot_vec *pVecY );
extern void potVecConeAxinvsqrVpy( double alpha, pot_vec *pVecX, pot_vec *pVecV, pot_vec *pVecY );
extern void potVecConeScal( pot_vec *pVecX, pot_vec *pVecY );
extern void potVecSimplexProj( pot_vec *pVecX );
extern double potVecConeMin( pot_vec *pVecX );
extern double potVecConeMax( pot_vec *pVecX );
extern double potVecLogDet( pot_vec *pVecX );
extern double potVecSumCone( pot_vec *pVecX );
extern double potVecSumScalCone( pot_vec *pVecX, pot_vec *pVecY );
extern double potVecRatioTest( pot_vec *pVecX, pot_vec *pVecdX, double dCone );
extern void potVecConeAddConstant( pot_vec *pVecX, double cVal );
extern void potVecClear( pot_vec *pVec );
extern void potVecDestroy( pot_vec **pVec );
extern void potVecCopy( pot_vec *srcVec, pot_vec *dstVec );
extern void potVecReset( pot_vec *pVec );
extern void potVecExport( pot_vec *pVec, double *dVec );

#ifdef __cplusplus
}
#endif

#endif /* pot_vector_h */
