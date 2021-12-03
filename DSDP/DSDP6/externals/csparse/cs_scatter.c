#include "cs.h"
/* x = x + beta * A(:,j), where x is a dense vector and A(:,j) is sparse */
CS_INT cs_scatter (const cs *A, CS_INT j, CS_ENTRY beta, CS_INT *w, CS_ENTRY *x, CS_INT mark,
    cs *C, CS_INT nz)
{
    CS_INT i, p, *Ap, *Ai, *Ci ;
    CS_ENTRY *Ax ;
    if (!CS_CSC (A) || !w || !CS_CSC (C)) return (-1) ;     /* check inputs */
    Ap = A->p ; Ai = A->i ; Ax = A->x ; Ci = C->i ;
    for (p = Ap [j] ; p < Ap [j+1] ; p++)
    {
        i = Ai [p] ;                            /* A(i,j) is nonzero */
        if (w [i] < mark)
        {
            w [i] = mark ;                      /* i is new entry in column j */
            Ci [nz++] = i ;                     /* add i to pattern of C(:,j) */
            if (x) x [i] = beta * Ax [p] ;      /* x(i) = beta*A(i,j) */
        }
        else if (x) x [i] += beta * Ax [p] ;    /* i exists in C(:,j) already */
    }
    return (nz) ;
}

CS_INT cs_scatter2 (const cs *A, DSDP_INT i, double *Arow) {
    // Scatter some column of the CSC matrix
    DSDP_INT retcode = DSDP_RETCODE_OK;
    
    double *Ax   = A->x;
    DSDP_INT *Ai = A->i;
    
    for (DSDP_INT j = A->p[i]; j < A->p[i + 1]; ++j ) {
        Arow[Ai[j]] = Ax[j];
    }
    
    return retcode;
}
