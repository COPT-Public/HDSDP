#ifndef def_hdsdp_user_data_h
#define def_hdsdp_user_data_h

#include "def_hdsdp_conic.h"
/* Interface of user data */

/** @struct hdsdp\_user\_data
 *  @brief HDSDP user conic data for SDP and LP
 *
 */
struct hdsdp_user_data {
    
    cone_type cone;
    
    int     nConicRow;
    int     nConicCol;
    int    *coneMatBeg;
    int    *coneMatIdx;
    double *coneMatElem;
    
};

#endif /* def_hdsdp_user_data_h */
