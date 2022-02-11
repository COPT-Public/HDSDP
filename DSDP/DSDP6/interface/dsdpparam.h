#ifndef dsdpparam_h
#define dsdpparam_h

/* Implement the parameter interface for DSDP */
#include <stdio.h>
#include "dsdphsd.h"

typedef struct {
    
    DSDP_INT AmaxIter;    // Algorithm maximum iteration
    DSDP_INT BmaxIter;
    DSDP_INT doPresolve; // Whether to do presolve
    
    double Asigma;
    double Bsigma;
    
    
    double rho;
    
    DSDP_INT initMethod;
    double initpObj;
    double initBeta;
    double initMu;
    double initTau;
    double initKappa;
    
    double Aalpha;
    double Balpha;
    
    DSDP_INT Ancorr;
    DSDP_INT Bncorr;
    double corrDelta;
    double corrAlpha;
    
    DSDP_INT Aattempt;
    
    double AnrmThresh;
    double BpInfeasSusp;
    
    // Tolerance
    double absOptTol;    // Algorithm absolute optimality tolerance
    double relOptTol;    // Algorithm relative optimality tolerance
    double absFeasTol;   // Algorithm absolute feasibility tolerance
    double relFeasTol;   // Algorithm relative feasibility tolerance
    
    // Others
    DSDP_INT CGreuse;
    
} hsdParam;

static hsdParam defaultParam =
{
    100,
    200,
    TRUE,
    
    0.1,
    0.1,
    
    5.0,
    
    DSDP_INITMETHOD_FRO,
    1e+05,
    1.0,
    1e+06,
    1.0,
    1.0,
    
    0.75,
    0.97,
    
    0,
    4,
    
    0.0,
    10.0,
    
    4,
    
    1e+20,
    1e+08,
    
    1e-04,
    1e-04,
    1e-04,
    1e-04,
    
    2,
};

#endif /* dsdpparam_h */
