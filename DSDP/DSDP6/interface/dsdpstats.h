#ifndef dsdpstats_h
#define dsdpstats_h

/* Implement the statistics interface of DSDP */
#include "dsdphsd.h"

#define NUM_STATISTICS       66

typedef struct {
    double stats[NUM_STATISTICS];
} DSDPStats;


#ifdef __cplusplus
extern "C" {
#endif

extern void DSDPStatInit    ( DSDPStats *stats );
extern void DSDPStatUpdate  ( DSDPStats *stat, DSDP_INT sName, double  val );
extern void DSDPGetStats    ( DSDPStats *stat, DSDP_INT sName, double *val );
extern void DSDPDataStatPrint   ( DSDPStats *stat );
extern void DSDPDIMACErrorPrint ( DSDPStats *stat );
extern void DSDPBProfilerPrint  ( DSDPStats *stat );

#ifdef __cplusplus
}
#endif

#endif /* dsdpstats_h */
