#ifndef linsys_h
#define linsys_h

#include "pot_def.h"

typedef struct {
    
    int   nCol;
    void *solver;
    
    int  (*LCreate)  ( void **, int );
    void (*LDestroy) ( void ** );
    
    int  (*LSFac)  ( void *, int *, int * );
    int  (*LNFac)  ( void *, int *, int *, double * );
    int  (*LSolve) ( void *, double * );

} pot_linsys;

extern pot_int potLinsysCreate( pot_linsys **ppotLinsys );
extern pot_int potLinsysInit( pot_linsys *potLinsys, pot_int nCol );
extern pot_int potLinsysSymFactorize( pot_linsys *potLinsys, pot_int *colMatBeg, pot_int *colMatIdx );
extern pot_int potLinsysNumFactorize( pot_linsys *potLinsys, int *colMatBeg, int *colMatIdx, double *colMatElem );
extern pot_int potLinsysSolve( pot_linsys *potLinsys, double *rhsVec, double *solVec );
extern void potLinsysClear( pot_linsys *potLinsys );
extern void potLinsysDestroy( pot_linsys **ppotLinsys );

#endif /* linsys_h */
