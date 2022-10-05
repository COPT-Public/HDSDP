/** @file pot\_utils.c
 *  @brief Implement the potential reduction utilities
 *
 * @TODO: Add more detailed comments
 */

#ifndef pot_utils_h
#define pot_utils_h

#include "pot_def.h"

#define POT_CALL(func) if ((func) != RETCODE_OK) {     \
                            retcode = RETCODE_FAILED;  \
                            goto exit_cleanup;         \
                        }
#endif /* pot_utils_h */
