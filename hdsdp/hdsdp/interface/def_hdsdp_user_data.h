#ifndef def_hdsdp_user_data_h
#define def_hdsdp_user_data_h

#ifdef HEADERPATH
#include "interface/def_hdsdp_conic.h"
#else
#include "def_hdsdp_conic.h"
#endif
/* Interface of user data */

/** @struct hdsdp\_user\_data
 *  @brief HDSDP user conic data for SDP and LP
 *
 * Conic data is defined differently for different cones
 *
 * SDP cone: CSC representation of an [(n + 1) \* n / 2] by [m + 1] matrix, containing SDP matrix
 * coefficients and objective coefficients in each column, only lower triangular is stored
 *
 * LP cone and bound cone: CSC representation of an [n] by [m] matrix, containg LP data
 * 
 */
struct hdsdp_user_data {
    
    cone_type cone;
    
    int     nConicRow; /* Number of constraints */
    int     nConicCol; /* For SDP cone: conic dimension, for LP cone: number of LP columns */
    int    *coneMatBeg;
    int    *coneMatIdx;
    double *coneMatElem;
    
};

#endif /* def_hdsdp_user_data_h */
