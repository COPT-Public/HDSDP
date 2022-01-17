#include <stdio.h>
#include "dsdplog.h"
#define dsdplog printf

extern void dsdpshowdash(void) {
    // --------------------------------------------------------------------------------
    for (DSDP_INT i = 0; i < 50; ++i) {
        printf("-");
    }
    printf("\n");
}
