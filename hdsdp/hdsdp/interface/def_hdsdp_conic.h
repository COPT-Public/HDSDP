/** @file def\_hdsdp\_conic
 *  @brief Define the basic conic structure
 *
 */
#ifndef def_hdsdp_conic_h
#define def_hdsdp_conic_h

#include "hdsdp.h"

/* Define conic type */
typedef enum {
    
    HDSDP_CONIC_LP,
    HDSDP_CONIC_SDP,
    HDSDP_CONIC_SOCP
    
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
    void   *coneData;
    
    /* Conic data interface */
    hdsdp_retcode (*coneInit)        ( void *, cone_type );
    hdsdp_retcode (*coneSetData)     ( void *, void * );
    hdsdp_retcode (*coneProcData)    ( void * );
    void          (*coneDestroyData) ( void * );
    
    /* Conic algorithm interface */
    void   (*coneSetStart)    ( void *, double );
    void   (*coneUpdate)      ( void *, double, double * );
    double (*coneRatioTest)   ( void *, double, double * );
    
    /* Schur complement and algorithm iterates */
    void   (*coneGetSymNnz)   ( void * );
    void   (*coneAddSymNz)    ( void *, int * );
    void   (*coneBuildM)      ( void *, void * );
    void   (*coneBuildAs)     ( void *, double * );
    void   (*coneBuildAscs)   ( void *, double * );
    void   (*coneBuildCs2c)   ( void *, double * );
    
    /* Barrier, projection and recovery */
    double (*coneGetBarrier)  ( void *, double, double * );
    int    (*conePFeasCheck)  ( void *, double, double * );
    void   (*conePRecover )   ( void *, double * );
    
    void   (*coneScal)        ( void *, double );
    
} hdsdp_cone;


#endif /* def_hdsdp_conic_h */
