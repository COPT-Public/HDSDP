/** @file speigs.h
 *  @brief Header for constants and routine list
 *
 * Implement the eigen-decomposition algorithm from DSDP5.8 by Steve Benson.
 *
 * Given a real symmetric matrix A, the routine explores special structures
 * within and computes the full eigen-decomposition of the matrix.
 * In the backend the routine calls Lapack dsyev to decompose the pre-processed system.
 *
 * This implementation is less regirous since it is called consecutively many times in HDSDP.
 * We therefore refer to https://github.com/leavesgrp/SPEIGS for a more careful implementation of the SPEIG routine library
 *
 *  @author Wenzhi Gao, Shanghai University of Finance and Economics
 *  @date Aug, 27th, 2022
 *
 */

#ifndef speigs_h
#define speigs_h

/* Include headers */
#include "structs.h"

/* Some constants */
#define LWORK   (30)
#define IWORK   (12)

#ifdef __cplusplus
extern "C" {
#endif

extern void     speigInit     ( speigfac *eigfac );
extern DSDP_INT speigAlloc    ( speigfac *eigfac, DSDP_INT nmax );
extern DSDP_INT speigSpFactor ( speigfac *eigfac, spsMat *A, double *eigvals, double *eigvecs );
extern void     speigFree     ( speigfac *eigfac );

#ifdef __cplusplus
}
#endif

#endif /* speigs_h */
