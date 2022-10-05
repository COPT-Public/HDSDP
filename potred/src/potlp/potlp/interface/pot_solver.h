#ifndef pot_solver_h
#define pot_solver_h

#include "pot_structs.h"

extern pot_int potLPCreate( pot_solver **ppot );
extern pot_int potLPInit( pot_solver *pot, pot_int vDim, pot_int vConeDim );
extern pot_int potLPSetObj( pot_solver *pot, pot_fx *objFunc );
extern pot_int potLPSetLinearConstrs( pot_solver *pot, pot_constr_mat *AMat );
extern void potLPClear( pot_solver *pot );
extern void potLPDestroy( pot_solver **ppot );

#endif /* pot_solver_h */
