/** @file potstructs.c
 *  @brief Define the data structure for potential reduction
 *
 * @TODO: Add more detailed comments
 */

#ifndef pot_structs_h
#define pot_structs_h

#include "potlp.h"
#include "pot_param.h"


/** @brief Struture that contains vector implementations
 *
 */
typedef struct {
    pot_int  n;
    pot_int  ncone;
    double   nrm;
    double    *x;
} pot_vec;

/** @brief Struture that contains sparse matrix implementation
 *
 */
typedef struct {
    pot_int  m;
    pot_int  n;
    
    void *AMatData;
    
    /* Abstract implementation */
    pot_int (*AMatInit) ( void **, const void * );
    void (*AMatPrepareX) ( void *, pot_vec * );
    void (*AMatProject) ( void *, pot_vec * );
    void (*AMatMonitor) ( void *, pot_int * );
    void (*AMatDestroy) ( void ** );
    
} pot_constr_mat;

/** @brief Struture that hosts the objective methods
 *
 */
typedef struct {
    
    pot_int n;    ///< Dimension of x
    
    pot_vec *xVec; ///< Pointer to x
    void *objFData; ///< Data that formulates \f$ f(x) \f$
    
    /* Abstract implementation */
    pot_int (*objFInit) ( void **, const void * ); ///< Method of initialization
    double (*objFVal) ( void *, pot_vec * ); ///< Method of computing objective
    void (*objFGrad)  ( void *, pot_vec *, pot_vec * ); ///< Method of computing gradient
    void (*objFHess)  ( void *, pot_vec *, double * ); ///< Method of computing Hessian
    void (*objFHVec)  ( void *, pot_vec *, double * ); ///< Method of Hessian-vector product
    void (*objFMonitor) ( void *, pot_int * ); ///< Method of internal progress monitor
    void (*objFDestroy) ( void ** ); ///< Method of destroy
    
} pot_fx;


/** @brief Struture that hosts the potential reduction framwork
 *
 */
typedef struct {
    
    pot_int n;     ///< Dimension of the problem

    pot_fx *objFunc; ///< Objective function pointer
    
    double   fVal; ///< Objective
    double   zVal; ///< Lower bound
    pot_vec *xVec; ///< Variable
    pot_vec *xVecOld; ///< Old x
    pot_vec *gVec; ///< Gradient
    pot_vec *mVec; ///< Momentum
    pot_vec *xStepVec; ///< Step direction
    double  *HessMat; ///< Hessian
    
    pot_constr_mat *AMat;
    void (*AMatProj) ( pot_constr_mat *, pot_vec * );
    
    double  rhoVal;
    double  potVal;
    double  betaRadius;
    
    double  *projHessMat[4];
    double  *projGMat[4];
    double  *projgVec[2];
    
    pot_vec *auxVec1;
    pot_vec *auxVec2;
    
    pot_int *intParams[NUM_INT_PARAM];
    pot_int *dblParams[NUM_DBL_PARAM];
    
    
} pot_solver;


#endif /* pot_structs_h */
