/** @file potstructs.c
 *  @brief Define the data structure for potential reduction
 *
 * @TODO: Add more detailed comments
 */

#ifndef pot_structs_h
#define pot_structs_h

#include "pot_solver.h"
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
    void (*AMatPrepareX) ( void *, pot_vec * );
    void (*AMatProject) ( void *, pot_vec * );
    void (*AMatScalProject) ( void *, pot_vec *, pot_vec * );
    void (*AMatMonitor) ( void *, void * );
    
} pot_constr_mat;

/** @brief Struture that hosts the objective methods
 *
 */
typedef struct {
    
    pot_int n;    ///< Dimension of x
    
    void *objFData; ///< Data that formulates \f$ f(x) \f$
    
    /* Abstract implementation */
    double (*objFVal) ( void *, pot_vec * ); ///< Method of computing objective
    void (*objFGrad)  ( void *, pot_vec *, pot_vec * ); ///< Method of computing gradient
    void (*objFHess)  ( void *, pot_vec *, double * ); ///< Method of computing Hessian
    void (*objFHVec)  ( void *, pot_vec *, pot_vec * ); ///< Method of Hessian-vector product
    void (*objFMonitor) ( void *, void * ); ///< Method of internal progress monitor
    
} pot_fx;


/** @brief Structure that hosts Lanczos computation
 *
 *
 */
typedef struct {
    
    pot_int n; ///< Dimension of x
    pot_int maxIter; /// Maximum subspace dimension
     
    void *MMat;  ///< Abstract Matrix
    
    pot_vec *vVec; ///< Auxiliary array v
    pot_vec *wVec; ///< Auxiliary array w
    pot_vec *z1Vec; ///< Auxiliary array z1
    pot_vec *z2Vec; ///< Auxiliary array z2
    pot_vec *vaVec; ///< Auxiliary array va
    
    double *VMat; ///< Auxiliary matrix V
    double *HMat; ///< Auxiliary matrix H
    double *YMat; ///< Auxiliary matrix Y
    double *UMat; ///< Auxiliary matrix U
    double *dArray; ///< Auxiliary array d
    
    double  *eiDblMat; ///< Eigen double workspace
    pot_int *eiIntMat; ///< Auxiliary int workspace
    
    void (*lczMatVec) (void *, pot_vec *, pot_vec *);
    
} pot_lanczos;


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
    pot_vec *gVec; ///< Gradient of f
    pot_vec *gkVec; ///< Gradient of potential function
    pot_vec *mkVec; ///< Momentum
    pot_vec *xStepVec; ///< Step direction
    double  *HessMat; ///< Hessian
    
    pot_constr_mat *AMat;
    
    pot_lanczos *lczTool;
    
    double  rhoVal;
    double  potVal;
    double  betaRadius;
    
    double  projHessMat[4];
    double  projGMat[4];
    double  projgVec[2];
    
    pot_vec *auxVec1;
    pot_vec *auxVec2;
    
    int intParams[NUM_INT_PARAM];
    int dblParams[NUM_DBL_PARAM];
    
    pot_int useCurvature;
    
} pot_solver;


#endif /* pot_structs_h */
