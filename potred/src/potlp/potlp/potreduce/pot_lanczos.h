#ifndef pot_lanczos_h
#define pot_lanczos_h

#include "pot_structs.h"

extern pot_int potLanczosCreate( pot_lanczos **ppotLanczos );
extern pot_int potLanczosInit( pot_lanczos *potLanczos, pot_int nCols );
extern void potLanczosInitData( pot_lanczos *potLanczos, void *MMat, void (*lczMatVec) (void *, pot_vec *, pot_vec *) );

extern pot_int potLanczosSolve( pot_lanczos *potLanczos, pot_vec *nCurv );
extern void potLanczosClear( pot_lanczos *potLanczos );
extern void potLanczosDestroy( pot_lanczos **ppotLanczos );

#endif /* pot_lanczos_h */
