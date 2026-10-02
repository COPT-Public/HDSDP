#ifndef lapack_names_h
#define lapack_names_h

#ifdef HEADERPATH
#include "interface/hdsdp.h"
#else
#include "hdsdp.h"
#endif

#ifdef UNDERBLAS
    #define ddot    ddot_
    #define daxpy   daxpy_
    #define dsymm   dsymm_
    #define dtrsm   dtrsm_
    #define dscal   dscal_
    #define dsyr    dsyr_
    #define dgemv   dgemv_
    #define dnrm2   dnrm2_
    #define dspmv   dspmv_
    #define dger    dger_
    #define drscl   drscl_
    #define dsymv   dsymv_
    #define dsytrs  dsytrs_
    #define dsyevr  dsyevr_
    #define dsytrf  dsytrf_
    #define dpotri  dpotri_
    #define dpotrs  dpotrs_
    #define dpotrf  dpotrf_
    #define dgemm   dgemm_
#endif

#ifdef CAPBLAS
    #define ddot    DDOT
    #define daxpy   DAXPY
    #define dsymm   DSYMM
    #define dtrsm   DTRSM
    #define dscal   DSCAL
    #define dsyr    DSYR
    #define dgemv   DGEMV
    #define dnrm2   DNRM2
    #define dspmv   DSPMV
    #define dger    DGER
    #define drscl   DRSCL
    #define dsymv   DSYMV
    #define dsytrs  DSYTRS
    #define dsyevr  DSYEVR
    #define dsytrf  DSYTRF
    #define dpotri  DPOTRI
    #define dpotrs  DPOTRS
    #define dpotrf  DPOTRF
    #define dgemm   DGEMM
#endif

#ifdef UNDERCAPBLAS
    #define ddot    DDOT_
    #define daxpy   DAXPY_
    #define dsymm   DSYMM_
    #define dtrsm   DTRSM_
    #define dscal   DSCAL_
    #define dsyr    DSYR_
    #define dgemv   DGEMV_
    #define dnrm2   DNRM2_
    #define dspmv   DSPMV_
    #define dger    DGER_
    #define drscl   DRSCL_
    #define dsymv   DSYMV_
    #define dsytrs  DSYTRS_
    #define dsyevr  DSYEVR_
    #define dsytrf  DSYTRF_
    #define dpotri  DPOTRI_
    #define dpotrs  DPOTRS_
    #define dpotrf  DPOTRF_
    #define dgemm   DGEMM_
#endif

#endif
