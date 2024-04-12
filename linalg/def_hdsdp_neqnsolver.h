#ifndef def_hdsdp_neqnsolver_h
#define def_hdsdp_neqnsolver_h

/* Implement normal equation solver for LPs with dense column
Assume that we solve an LP in standard form
 
 min           c' * x
 subject to     A * x == b
                    x >= 0
 
 and the normal equation is given by
 
                    A * X^2 * A' l = r
 
 When A contains dense column, we partition columns of A into A = [A_S, A_D]
 and solve
 
        ((A * X_S^2 A' + (A_D * X_D) * (A_D * X_D)') + R * I) l = r,
 
 where R is the regularization term for the normal equation.
*/

#ifdef HEADERPATH
#include "linalg/hdsdp_linsolver.h"
#else
#include "hdsdp_linsolver.h"
#endif

typedef struct {
    
    int nCol;
    int nRow;
    
    double dReg;
    int *colSparsity;
    int *colSparsitySortedId;
    
    /* A */
    int *colMatBeg;
    int *colMatIdx;
    double *colMatElem;
    
    /* A_S */
    int nColSparse;
    int *spColIdx;
    int *colMatSparseBeg;
    int *colMatSparseIdx;
    double *colMatSparseElem;
    
    /* A_D */
    int nColDense;
    int *dsColIdx;
    double *colMatDenseElem;
    
    /* A_S * X_S^2 * A_S' */
    int *schurMatBeg;
    int *schurMatIdx;
    double *schurMatElem;
    
    double *dRowBuffer;
    double *dColBuffer;
    
    /* Internal linear system solver */
    hdsdp_linsys_fp *spchol;
    
} hdsdp_normal_linsys;


/* Iterative method */
typedef struct {
    
    int nCol;
    
    double *iterResi;
    double *iterResiNew;
    double *iterDirection;
    double *preInvResi;
    double *MTimesDirection;
    double *iterVec;
    double *rhsBuffer;
    
    /* Matrix data */
    void *MMat;
    void (*Mvec) (void *, double *, double *);
    
    /* Preconditioner */
    void *Mcond;
    void (*Mprecond) (void *, double *, double *);
    
    /* Statistics */
    double iterResiNorm;
    double solveTime;
    
    int nIters;
    iter_status solStatus;
    int nSolves;
    
    iterative_params params;
    
} hdsdp_abstract_iterative;


#endif /* def_hdsdp_neqnsolver_h */
