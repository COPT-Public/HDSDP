#include <math.h>
#include "vec_mat.h"

extern double nrm2( pot_int *n, double *x, pot_int *incx ) {
    
    assert( *incx == 1 );
    
    double nrm = 0.0;
    
    for ( int i = 0; i < *n; ++i ) {
        nrm += x[i] * x[i];
    }
    
    return sqrt(nrm);
}

extern void axpy( pot_int *n, double *a, double *x, pot_int *incx, double *y, pot_int *incy ) {
    
    assert( *incx == 1 && *incy == 1 );
    
    for ( int i = 0; i < *n; ++i ) {
        y[i] += (*a) * x[i];
    }
    
    return;
}

extern double dot( pot_int *n, double *x, pot_int *incx, double *y, pot_int *incy ) {
    
    assert( *incx == 1 && *incy == 1 );
    
    double dres = 0.0;
    
    for ( int i = 0; i < *n; ++i ) {
        dres += x[i] * y[i];
    }
    
    return dres;
}

extern void scl( pot_int *n, double *sa, double *sx, pot_int *incx ) {
    
    assert( *incx == 1 );
    double a = *sa;
    
    for ( int i = 0; i < *n; ++i ) {
        sx[i] = sx[i] * a;
    }
    
    return;
}


extern void rscl( pot_int *n, double *sa, double *sx, pot_int *incx ) {
    
    assert( *incx == 1 );
    double a = *sa;
    
    assert( a != 0.0 );
    assert( a > 0.0 );
    
    if ( fabs(a) < 1e-16 ) {
        a = (a > 0) ? 1e-16 : -1e-16;
    }
    
    for ( int i = 0; i < *n; ++i ) {
        sx[i] = sx[i] / a;
    }
    
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

extern void dgemv( char *trans, pot_int *m, pot_int *n, double *alpha,
                   double *a, pot_int *lda, double *x, pot_int *incx,
                   double *beta, double *y, pot_int *incy );

extern void gemv( char *trans, pot_int *m, pot_int *n, double *alpha,
                  double *a, pot_int *lda, double *x, pot_int *incx,
                  double *beta, double *y, pot_int *incy ) {
    
    dgemv(trans, m, n, alpha, a, lda, x, incx, beta, y, incy);
    
    return;
}

extern double sumlogdet( pot_int *n, double *x ) {
    
    double logdet = 0.0;
    for ( int i = 0; i < *n; ++i ) {
        logdet += logl(x[i]);
    }
    
    return logdet;
}

extern void vvscl( pot_int *n, double *s, double *x ) {
    
    for ( int i = 0; i < *n; ++i ) {
        assert( s[i] > 0.0 );
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
