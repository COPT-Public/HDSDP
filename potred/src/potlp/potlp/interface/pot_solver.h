#ifndef pot_solver_h
#define pot_solver_h

#include "pot_def.h"
#include "pot_structs.h"

#define POT_CALL(func) if ((func) != RETCODE_OK) {     \
                            retcode = RETCODE_FAILED;  \
                            goto exit_cleanup;         \
                        }

extern pot_int potLPCreate( pot_solver **ppot );
extern pot_int potLPInit( pot_solver *pot, pot_int vDim, pot_int vConeDim );
extern pot_int potLPSetObj( pot_solver *pot, pot_fx *objFunc );
extern pot_int potLPSetLinearConstrs( pot_solver *pot, pot_constr_mat *AMat );
extern pot_int potReductionSolve( pot_solver *pot );
extern void potLPClear( pot_solver *pot );
extern void potLPDestroy( pot_solver **ppot );

#endif /* pot_solver_h */
