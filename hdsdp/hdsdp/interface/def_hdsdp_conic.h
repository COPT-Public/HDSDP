/** @file def\_hdsdp\_conic
 *  @brief Define the basic conic structure
 *
 */
#ifndef def_hdsdp_conic_h
#define def_hdsdp_conic_h

#ifdef HEADERPATH
#include "interface/hdsdp.h"
#include "linalg/def_hdsdp_sdpdata.h"
#include "linalg/def_hdsdp_linsolver.h"
#else
#include "hdsdp.h"
#include "def_hdsdp_sdpdata.h"
#include "def_hdsdp_linsolver.h"
#endif
#include <stdint.h>

/* Define conic type */
typedef enum {
    
    HDSDP_CONETYPE_UNKNOWN,
    HDSDP_CONETYPE_LP,         /* A' * y <= c */
    HDSDP_CONETYPE_BOUND,      /*      y <= u */
    HDSDP_CONETYPE_DENSE_SDP,
    HDSDP_CONETYPE_SPARSE_SDP,
    HDSDP_CONETYPE_SOCP
    
} cone_type;

/** @struct hdsdp\_cone
 *  @brief Define the HDSDP general conic interface
 *
 * In HDSDP, the general conic interface supports the following functionalities
 *
 *  Connection with interface:
 *  1. Set data
 *  2. Process data
 *  3. Destroy data
 *
 *  Algorithm:
 *  1. Conic initialize
 *  2. Conic iterate maintenance
 *  3. Conic numeric assembly of important quantities
 *  4. Conic symbolic assembly of important quantities
 *  5. Conic ratio test
 *  6. Conic barrier computation
 *  7. Conic primal projection
 *  8. Conic primal variable recovery
 *  9. Conic scal
 *
 */
typedef struct {
    
    cone_type cone; ///< What cone is it?
    
    void  *usrData;
    void  *coneData;
    
    /* Conic data interface */
    hdsdp_retcode (*coneCreate)    ( void ** );
    hdsdp_retcode (*coneProcData)    ( void *, int, int, int *, int *, double * );
    hdsdp_retcode (*conePresolveData) ( void * );
    void          (*coneDestroyData) ( void ** );
    
    /* Conic algorithm interface */
    void   (*coneSetStart)    ( void *, double );
    void   (*coneUpdate)      ( void *, double, double * );
    double (*coneRatioTest)   ( void *, double, double * );
    
    /* Schur complement and algorithm iterates */
    int64_t (*coneGetSymNnz)   ( void * );
    void    (*coneAddSymNz)    ( void *, int * );
    void    (*coneBuildSchur)  ( void *, int, void * );
    
    /* Barrier, projection and recovery */
    double  (*coneGetBarrier)  ( void *, double, double * );
    int     (*conePFeasCheck)  ( void *, double, double * );
    void    (*conePRecover)    ( void *, double * );
    
    double  (*coneGetCoeffNorm) ( void *, int );
    double  (*coneGetObjNorm)  ( void *, int );
    void    (*coneScal)        ( void *, double );
    
    /* Debugging */
    void    (*coneView)        ( void * );
    
} hdsdp_cone;

/* A dense SDP block */
typedef struct {
    
    int   nRow;
    int   nCol;
    
    int  isDualSparse;
    
    double dualResidual;
    
    /* Dual symbolic structure */
    int  *dualMatBeg;
    int  *dualMatIdx;
    int  *dualPosToElemMap;
    
    /* Dual matrix */
    double *dualMatElem;
    /* Dual buffer for checking positive definiteness */
    double *dualCheckerElem;
    /* Dual step buffer */
    double *dualStep;
    /* Dual factorization */
    hdsdp_linsys_fp *dualFactor;
    /* Dual Lanczos*/
    
    sdp_coeff **sdpRow;
    sdp_coeff  *sdpObj;
    
    int *sdpConePerm;
    
    /* SDP block statistics */
    int sdpConeStats[5]; ///< Number of coefficients of each type
    
} hdsdp_cone_sdp_dense;

/* A sparse SDP block */
typedef struct {
    
    int   nRow;
    int   nCol;
    
    int isDualSparse;
    
    double dualResidual;
    
    /* Dual symbolic structure */
    int  *dualMatBeg;
    int  *dualMatIdx;
    int  *dualPosToElemMap;
    
    /* Dual matrix */
    double *dualMatElem;
    /* Dual checker */
    double *dualCheckerElem;
    /* Dual step buffer */
    double *dualStep;
    /* Dual factorization */
    hdsdp_linsys_fp *dualFactor;
    
    int nRowElem;
    int *rowIdx;
    sdp_coeff **sdpRow;
    sdp_coeff  *sdpObj;
    
    int sdpConeStats[5];
    
} hdsdp_cone_sdp_sparse;

/* TODO: The structures for the following cones are not yet formally set up */
/* An LP cone */
typedef struct {
    
    int     nRow;
    int     nCol;
    
    double *colDual;
    double *colDualChecker;
    double *colDualStep;
    
    int *rowMatBeg;
    int *rowMatIdx;
    double *rowMatElem;
    
} hdsdp_cone_lp;

/* A bound cone */
typedef struct {
    
    int nRow;
    
    double *colDualUpSlack;
    double *colDualLowSlack;
    
    int    *dBoundLowIdx;
    double *dBoundLow;
    int    *dBoundUpIdx;
    double *dBoundUp;
    
} hdsdp_cone_bound;

#endif /* def_hdsdp_conic_h */
