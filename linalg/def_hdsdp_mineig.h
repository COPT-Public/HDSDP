#ifndef def_hdsdp_mineig_h
#define def_hdsdp_mineig_h

/** @struct hdsdp\_mineig */
typedef struct {
    
    /* Basic information and data */
    
    int nCol;           ///< Dimension of the matrix
    int nElem;          ///< Nnz of the matrix
    void *chol;         ///< Cholesky data structure
    int (*ispsd) ( void *, int *, int *, double *, int * );         ///< Cholesky certificate
    int *colMatBeg;     ///< CSC beg
    int *colMatIdx;     ///< CSC idx
    double *colMatElem; ///< CSC elem
    double *colMatElemPlusEye; ///< Elements of the perturbed matrix 
    double dPerturb;    ///< Recent perturbation value
    
    /* Lanczos Krylov solver */
    int nKrylovDim; ///< Dimension of Krylov subspace. Set as a parameter
    double *vVec; ///< Krylov solver auxiliary vector of size nCol
    double *wVec; ///< Krylov solver auxiliary vector of size nCol
    double *z1Vec; ///< Krylov solver auxiliary vector of size nCol
    double *z2Vec; ///< Krylov solver auxiliary vector of size nCol
    double *vaVec; ///< Krylov solver auxiliary vector of size nCol
    
    double *VMat; ///< Krylov tridiagonal matrix
    double *HMat;
    double *YMat;
    double *UMat;
    
    double *dArray;
    double *eigDblMat;
    int    *eigIntMat;
    
} hdsdp_mineig;

typedef enum {
    
    MIN_EIG_LANCZOS,    ///< Lanczos iteration
    MIN_EIG_BISECTION   ///< Bisection
    
} mineig_algo;


#endif /* def_hdsdp_mineig_h */
