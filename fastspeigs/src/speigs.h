#ifndef speigs_h
#define speigs_h

#include <stddef.h>

/*
 Implement the eigen-decomposition algorithm from DSDP5.8 by Steve Benson.
 
 Given a real symmetric matrix A, the routine explores special structures
 within and computes the full eigen-decomposition of the matrix.
 In the backend the routine calls Lapack dsyev or Netlib Eispack to
 decompose the pre-processed system.

 This routine is also employed in HDSDP solver for SDP.
*/

#ifndef BACKEND_LAPACK
#define BACKEND_LAPACK
#endif

#ifdef SPEIG_64
typedef int32_t spint;
#else
typedef int64_t spint;
#endif

#ifdef BACKEND_LAPACK
#define eig_dsyevr
#else
#endif

#define sperr printf


#define SPEIG_VER  (1) // Version number

#endif /* speigs_h */
