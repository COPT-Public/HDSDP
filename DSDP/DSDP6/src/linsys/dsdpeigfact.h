#ifndef dsdpeigfact_h
#define dsdpeigfact_h

#include "dsdphsd.h"
#include "structs.h"

#ifdef __cplusplus
extern "C" {
#endif

extern DSDP_INT factorizeSparseData( spsMat *A, double elow, double *eigvals, double *eigvecs );
extern DSDP_INT factorizeDenseData ( dsMat *A, double elow, double *eigvals, double *eigvecs  );

#ifdef __cplusplus
}
#endif
#endif /* dsdpeigfact_h */
