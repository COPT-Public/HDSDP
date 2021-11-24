#ifndef dsdpparam_h
#define dsdpparam_h

/* Implement the parameter interface for DSDP */
#include <stdio.h>
#include <dsdphsd.h>

typedef struct {
    
    DSDP_INT maxIter;    // Algorithm maximum iteration
    DSDP_INT doPresolve; // Whether to do presolve
    
    // Tolerance
    double absOptTol;    // Algorithm absolute optimality tolerance
    double relOptTol;    // Algorithm relative optimality tolerance
    double absFeasTol;   // Algorithm absolute feasibility tolerance
    double relFeasTol;   // Algorithm relative feasibility tolerance

} hsdParam;

static hsdParam defaultParam =
{
    TRUE,
    1000,
    1e-08,
    1e-06,
    1e-08,
    1e-06
};

#endif /* dsdpparam_h */
