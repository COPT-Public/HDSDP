#include <math.h>
#include "vec_mat.h"

#ifdef MYBLAS
#undef MYBLAS
#endif

#define MYBLAS

/* Blas functions */
extern double dnrm2( int *n, double *x, int *incx );
extern void daxpy( int *n, double *a, double *x, int *incx, double *y, int *incy );
extern double ddot( int *n, double *x, int *incx, double *y, int *incy );
extern void dscal( int *n, double *sa, double *sx, int *incx );
extern void drscl( int *n, double *sa, double *sx, int *incx );
extern int idamax( int *n, double *x, int *incx );
extern void dgemv( char *trans, pot_int *m, pot_int *n, double *alpha,
                   double *a, pot_int *lda, double *x, pot_int *incx,
                   double *beta, double *y, pot_int *incy );

extern double nrm2( pot_int *n, double *x, pot_int *incx ) {
#ifdef MYBLAS
    assert( *incx == 1 );
    
    double nrm = 0.0;
    
    for ( int i = 0; i < *n; ++i ) {
        nrm += x[i] * x[i];
    }
    
    return sqrt(nrm);
#else
    return dnrm2(n, x, incx);
#endif
}

extern void axpy( pot_int *n, double *a, double *x, pot_int *incx, double *y, pot_int *incy ) {
#ifdef MYBLAS
    assert( *incx == 1 && *incy == 1 );
    
    for ( int i = 0; i < *n; ++i ) {
        y[i] += (*a) * x[i];
    }
#else
    daxpy(n, a, x, incx, y, incy);
#endif
    return;
}

extern double dot( pot_int *n, double *x, pot_int *incx, double *y, pot_int *incy ) {
#ifdef MYBLAS
    assert( *incx == 1 && *incy == 1 );
    
    double dres = 0.0;
    
    for ( int i = 0; i < *n; ++i ) {
        dres += x[i] * y[i];
    }
    
    return dres;
#else
    return ddot(n, x, incx, y, incy);
#endif
}

extern void scal( pot_int *n, double *sa, double *sx, pot_int *incx ) {
#ifdef MYBLAS
    assert( *incx == 1 );
    double a = *sa;
    
    if ( a == 1.0 ) {
        return;
    }
    
    for ( int i = 0; i < *n; ++i ) {
        sx[i] = sx[i] * a;
    }
#else
    dscal(n, sa, sx, incx);
#endif
    return;
}

/* Use standard Blas for this sensitive operation */
extern void rscl( pot_int *n, double *sa, double *sx, pot_int *incx ) {
#if 0
    assert( *incx == 1 );
    double a = *sa;
    
    assert( a != 0.0 );
    assert( a > 0.0 );
    
    if ( a == 1.0 ) {
        return;
    }
    
    if ( fabs(a) < 1e-16 ) {
        a = (a > 0) ? 1e-16 : -1e-16;
    }
    
    for ( int i = 0; i < *n; ++i ) {
        sx[i] = sx[i] / a;
    }
#else
    drscl(n, sa, sx, incx);
#endif
    
    return;
}

extern pot_int idamax( pot_int *n, double *x, pot_int *incx ) {
    
    pot_int idmax = 0;
    double damax = 0.0;
    
    for ( int i = 0; i < *n; ++i ) {
        double ax = fabs(x[i]);
        if ( ax > damax ) {
            damax = ax; idmax = i;
        }
    }
    
    return idmax;
}

#ifdef UNDERBLAS
#define eig dsyevr_
#else
#define eig dsyevr
#endif

extern void eig( const char     *jobz,
                 const char     *range,
                 const char     *uplo,
                 const pot_int  *n,
                 double         *a,
                 const pot_int  *lda,
                 const double   *vl,
                 const double   *vu,
                 const pot_int  *il,
                 const pot_int  *iu,
                 const double   *abstol,
                 pot_int        *m,
                 double         *w,
                 double         *z,
                 const pot_int  *ldz,
                 pot_int        *isuppz,
                 double         *work,
                 const pot_int  *lwork,
                 pot_int        *iwork,
                 const pot_int  *liwork,
                 pot_int        *info );

extern void gemv( char *trans, pot_int *m, pot_int *n, double *alpha,
                  double *a, pot_int *lda, double *x, pot_int *incx,
                  double *beta, double *y, pot_int *incy ) {
    
    dgemv(trans, m, n, alpha, a, lda, x, incx, beta, y, incy);
    
    return;
}

extern double sumlogdet( pot_int *n, double *x ) {
    
    double logdet = 0.0;
    for ( int i = 0; i < *n; ++i ) {
        logdet += log(x[i]);
    }
    
    return logdet;
}

extern void vvscl( pot_int *n, double *s, double *x ) {
    
    for ( int i = 0; i < *n; ++i ) {
        x[i] = x[i] * s[i];
    }
    
    return;
}

extern void vvrscl( pot_int *n, double *s, double *x ) {
    
    for ( int i = 0; i < *n; ++i ) {
        x[i] = x[i] / s[i];
    }
    
    return;
}

extern pot_int psyev( pot_int n, double *U, double *d, double *Y,
                      double *work, pot_int *iwork, pot_int lwork, pot_int liwork ) {
    
    pot_int retcode = RETCODE_OK;
    
    char jobz = 'V', range = 'I', uplo = 'U';
    int isuppz[4] = {0};
    int il = n - 1, iu = n;
    int m = 2;
    int info = 0;
    
    dsyevr(&jobz, &range, &uplo, &n, U, &n,
           &potDblConstantZero, &potDblConstantZero,
           &il, &iu, &potDblConstantZero, &m, d, Y,
           &n, isuppz, work, &lwork, iwork, &liwork, &info);
    
    if ( info != 0 ) {
        retcode = RETCODE_FAILED;
    }
    
    return retcode;
}

extern void pgemv( pot_int m, pot_int n, double *M, double *v, double *y ) {
    
    char trans = 'N';
    dgemv(&trans, &m, &n, &potDblConstantOne, M, &m, v,
         &potIntConstantOne, &potDblConstantZero, y, &potIntConstantOne);
    
    return;
}

extern void spMatAxpy( int n, int *Ap, int *Ai, double *Ax, double a, double *x, double *y ) {
    
    for ( int i = 0, j; i < n; ++i ) {
        for ( j = Ap[i]; j < Ap[i + 1]; ++j ) {
            y[Ai[j]] += a * x[i] * Ax[j];
        }
    }
    
    return;
}

extern void spMatATxpy( int n, int *Ap, int *Ai, double *Ax, double a, double *x, double *y ) {
    
    for ( int i = 0, j; i < n; ++i ) {
        double aTy = 0.0;
        for ( j = Ap[i]; j < Ap[i + 1]; ++j ) {
            aTy += x[Ai[j]] * Ax[j];
        }
        y[i] += a * aTy;
    }
    
    return;
}

extern void spMatMaxRowAbs( int n, int *Ap, int *Ai, double *Ax, double *row ) {
    
    double x = 0.0;
    for ( int i = 0, j; i < n; ++i ) {
        for ( j = Ap[i]; j < Ap[i + 1]; ++j ) {
            x = fabs(Ax[j]);
            row[Ai[j]] = ( x > row[Ai[j]] ) ? x : row[Ai[j]];
        }
    }
    
    return;
}

extern void spMatMaxColAbs( int n, int *Ap, int *Ai, double *Ax, double *col ) {
    
    double cmax = 0.0, x = 0.0;
    for ( int i = 0, j; i < n; ++i ) {
        cmax = 0.0;
        for ( j = Ap[i]; j < Ap[i + 1]; ++j ) {
            x = fabs(Ax[j]);
            cmax = ( x > cmax ) ? x : cmax;
        }
        col[i] = ( cmax > col[i] ) ? cmax : col[i];
    }
    
    return;
}

extern void spMatRowScal( int n, int *Ap, int *Ai, double *Ax, double *row ) {
    
    for ( int i = 0, j; i < n; ++i ) {
        for ( j = Ap[i]; j < Ap[i + 1]; ++j ) {
            Ax[j] = Ax[j] / row[Ai[j]];
        }
    }
    
    return;
}

extern void spMatColScal( int n, int *Ap, int *Ai, double *Ax, double *col ) {
    
    for ( int i = 0, j; i < n; ++i ) {
        for ( j = Ap[i]; j < Ap[i + 1]; ++j ) {
            Ax[j] = Ax[j] / col[i];
        }
    }

    return;
}

extern int spMatBuildQMat( int qm, int qn, int *Qp, int *Qi, double *Qx,
                           int am, int an, int *Ap, int *Ai, double *Ax,
                           double *b, double *c ) {
    
    pot_int retcode = RETCODE_OK;
    
    /* Auxiliary array */
    int *amaux = NULL;
    int *anaux = NULL;
    
    int nzA = Ap[an];
    
    POTLP_INIT(amaux, int, am + 1);
    POTLP_INIT(anaux, int, an + 1);
    
    if ( !amaux || !anaux ) {
        retcode = RETCODE_FAILED;
        goto exit_cleanup;
    }
    
    for ( int i = 1; i < am + 1; ++i ) {
        amaux[i] = 1;
    }
    
    /* Build the first column block
       [  0 ]
       [ -A']
       [  b']
     */
    
    
    for ( int i = 0; i < nzA; ++i ) {
        amaux[Ai[i] + 1] += 1;
    }
    
    for ( int i = 0; i < am; ++i ) {
        amaux[i + 1] += amaux[i];
    }
    
    POTLP_MEMCPY(Qp, amaux, int, am + 1);
    for ( int i = 0, j, k; i < an; ++i ) {
        for ( j = Ap[i]; j < Ap[i + 1]; ++j ) {
            k = amaux[Ai[j]];
            Qi[k] = i + am; Qx[k] = -Ax[j];
            amaux[Ai[j]] += 1;
        }
    }
    
    for ( int i = 0; i < am; ++i ) {
        Qi[amaux[i]] = am + an;
        Qx[amaux[i]] = b[i];
    }
    
    /* Build the second column block
       [  A ]
       [  0 ]
       [ -c']
     */
    
    anaux[0] = Qp[am];
    for ( int i = 1; i < an + 1; ++i ) {
        anaux[i] = Ap[i] + anaux[0] + i; /* i reserved for c' */
    }
    
    int *pQp = Qp + am;
    int *pQi = Qi + nzA + am;
    int *pAi = Ai;
    double *pQx = Qx + nzA + am;
    double *pAx = Ax;

    POTLP_MEMCPY(pQp, anaux, int, an + 1);
    pQp += an;
    for ( int i = 0, k; i < an; ++i ) {
        /* Copy A */
        k = Ap[i + 1] - Ap[i];
        POTLP_MEMCPY(pQi, pAi, int, k);
        POTLP_MEMCPY(pQx, pAx, double, k);
        /* Copy -c */
        pQi += k; pAi += k;
        *pQi = am + an; pQi += 1;
        pQx += k; pAx += k;
        *pQx = -c[i]; pQx += 1;
    }
    
    /* Build the third column block and alternatively the fourth block
       [  0 (  0 ) ]
       [ -I (  0 ) ]
       [  0 ( -1 ) ]
     */
    
    /* No kappa if blkn = an */
    int blkn = an + 1;
    for ( int i = 0; i < blkn; ++i ) {
        pQp[i + 1] = pQp[i] + 1;
        pQi[i] = am + i;
        pQx[i] = -1.0;
    }
    
    pQp += blkn; pQi += blkn; pQx += blkn;
    
    /* Build the last column block
       [ -b ]
       [  c ]
       [  0 ]
     */
    pQp[1] = pQp[0] + am + an;
    for ( int i = 0; i < am; ++i ) {
        pQi[i] = i; pQx[i] = -b[i];
    }
    
    pQi += am; pQx += am;
    for ( int i = 0; i < an; ++i ) {
        pQi[i] = i + am; pQx[i] = c[i];
    }
    
    /* Done */
    assert( pQp[1] == 2 * nzA + 2 * am + 3 * an + 1 );
    
exit_cleanup:
    
    POTLP_FREE(amaux);
    POTLP_FREE(anaux);
    return retcode;
}

#define RUIZ_DEBUG(format, info) // printf(format, info);
extern int spMatRuizScal( int m, int n, int *Ap, int *Ai, double *Ax, double *D, double *E, int maxIter ) {
    
    pot_int retcode = RETCODE_OK;
    
    pot_int nRow = m;
    pot_int nCol = n;
    
    double *ruizScalDiagRow = D;
    double *ruizScalDiagCol = E;
    double *ruizWorkDiagRow = NULL;
    double *ruizWorkDiagCol = NULL;
    
    /* Allocate workspace */
    POTLP_INIT(ruizWorkDiagRow, double, nRow);
    POTLP_INIT(ruizWorkDiagCol, double, nCol);
    
    if ( !ruizWorkDiagRow || !ruizWorkDiagCol ) {
        retcode = RETCODE_FAILED;
        goto exit_cleanup;
    }
    
    /* Initialize scalers */
    for ( int i = 0; i < nRow; ++i ) {
        ruizScalDiagRow[i] = 1.0;
    }
    
    for ( int i = 0; i < nCol; ++i ) {
        ruizScalDiagCol[i] = 1.0;
    }
    
    RUIZ_DEBUG("Start Ruiz-scaling %s\n", "");
    
    for ( int i = 0; i < maxIter; ++i ) {
        
        POTLP_ZERO(ruizWorkDiagRow, double, nRow);
        spMatMaxRowAbs(nCol, Ap, Ai, Ax, ruizWorkDiagRow);
        POTLP_ZERO(ruizWorkDiagCol, double, nCol);
        spMatMaxColAbs(nCol, Ap, Ai, Ax, ruizWorkDiagCol);
        
        /* sqrt operation */
        double maxRuizDiagDeviate = 0.0;
        double ruizDiagDeviate = 0.0;
        
        for ( int j = 0; j < nRow; ++j ) {
            ruizWorkDiagRow[j] = sqrtl(ruizWorkDiagRow[j]);
            ruizScalDiagRow[j] = ruizScalDiagRow[j] * ruizWorkDiagRow[j];
            ruizDiagDeviate = fabs(ruizWorkDiagRow[j] - 1.0);
            maxRuizDiagDeviate = POTLP_MAX(maxRuizDiagDeviate, ruizDiagDeviate);
        }
        
        for ( int j = 0; j < nCol; ++j ) {
            ruizWorkDiagCol[j] = sqrtl(ruizWorkDiagCol[j]);
            ruizScalDiagCol[j] = ruizScalDiagCol[j] * ruizWorkDiagCol[j];
            ruizDiagDeviate = fabs(ruizWorkDiagCol[j] - 1.0);
            maxRuizDiagDeviate = POTLP_MAX(maxRuizDiagDeviate, ruizDiagDeviate);
        }
        
        RUIZ_DEBUG("Ruiz Deviation %e \n", maxRuizDiagDeviate);
        
        if ( maxRuizDiagDeviate < 1e-08 ) {
            RUIZ_DEBUG("Ruiz Successfully Ends in %d iterations \n", i);
            break;
        }
        
        /* Scaling */
        spMatRowScal(nCol, Ap, Ai, Ax, ruizWorkDiagRow);
        spMatColScal(nCol, Ap, Ai, Ax, ruizWorkDiagCol);
    }
    
    RUIZ_DEBUG("Ruiz-scaling Ends %s\n", "");
    
    for ( int i = 0; i < nRow; ++i ) {
        ruizScalDiagRow[i] = 1.0 / ruizScalDiagRow[i];
    }
    
    for ( int i = 0; i < nCol; ++i ) {
        ruizScalDiagCol[i] = 1.0 / ruizScalDiagCol[i];
    }
    
    
exit_cleanup:
    
    POTLP_FREE(ruizWorkDiagRow);
    POTLP_FREE(ruizWorkDiagCol);
    
    return retcode;
}

extern int spMatL2Scal( int m, int n, int *Ap, int *Ai, double *Ax, double *D, double *E ) {
    
    pot_int retcode = RETCODE_OK;
    
    pot_int nRow = m;
    pot_int nCol = n;
    
    double *L2WorkDiagRow = NULL;
    double *L2WorkDiagCol = NULL;
    
    /* Allocate workspace */
    POTLP_INIT(L2WorkDiagRow, double, nRow);
    POTLP_INIT(L2WorkDiagCol, double, nCol);
    
    for ( int i = 0, j; i < nCol; ++i ) {
        for ( j = Ap[i]; j < Ap[i + 1]; ++j ) {
            L2WorkDiagRow[Ai[j]] += Ax[j] * Ax[j];
        }
    }
    
    for ( int i = 0; i < nRow; ++i ) {
        L2WorkDiagRow[i] = sqrt(L2WorkDiagRow[i]);
        if ( L2WorkDiagRow[i] == 0.0 ) {
            L2WorkDiagRow[i] = 1.0;
        }
    }
    
    spMatRowScal(nCol, Ap, Ai, Ax, L2WorkDiagRow);
    
    for ( int i = 0, j; i < nCol; ++i ) {
        for ( j = Ap[i]; j < Ap[i + 1]; ++j ) {
            L2WorkDiagCol[i] += Ax[j] * Ax[j];
        }
    }
    
    for ( int i = 0; i < nCol; ++i ) {
        L2WorkDiagCol[i] = sqrt(L2WorkDiagCol[i]);
        if ( L2WorkDiagCol[i] == 0.0 ) {
            L2WorkDiagCol[i] = 1.0;
        }
    }
    
    spMatColScal(nCol, Ap, Ai, Ax, L2WorkDiagCol);;
    
    /* Update scaling */
    vvrscl(&nRow, L2WorkDiagRow, D);
    vvrscl(&nCol, L2WorkDiagCol, E);
    
exit_cleanup:
    
    POTLP_FREE(L2WorkDiagRow);
    POTLP_FREE(L2WorkDiagCol);
    
    return retcode;
}
