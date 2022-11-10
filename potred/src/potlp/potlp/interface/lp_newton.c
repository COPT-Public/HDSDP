/** @file lp\_newton.c
 *  @brief Add in the third direction of Newton to the potential reduction algorithm
 * This direction resembles the second-order interior point method and solves the augmented system for the directions
 * Only one iteration is implemented.
 *
 */

#include "vec_mat.h"
#include "linsys.h"
#include "lp_newton.h"

#include <math.h>

static double LpNewtonIBarrier( int nCol, double *x, double *s, double kappa, double tau ) {
    
    double mu = kappa * tau;
    for ( int i = 0; i < nCol; ++i ) {
        mu += x[i] * s[i];
    }
    
    return mu / (nCol + 1);
}

static double LpNewtonIRatioTest( int nCol, double *x, double *dx, double *s, double *ds,
                                  double kappa, double dkappa, double tau, double dtau ) {
    
    /* 1 / abs(min([dx./x; ds./s; dkappa / kappa; dtau./tau])) */
    
    double alphaTmp = POTLP_INFINITY;
    double ratio;
    
    for ( int i = 0; i < nCol; ++i ) {
        ratio = dx[i] / x[i];
        alphaTmp = POTLP_MIN(ratio, alphaTmp);
    }
    
    for ( int i = 0; i < nCol; ++i ) {
        ratio = ds[i] / s[i];
        alphaTmp = POTLP_MIN(ratio, alphaTmp);
    }
    
    ratio = dkappa / kappa;
    alphaTmp = POTLP_MIN(ratio, alphaTmp);
    
    ratio = dtau / tau;
    alphaTmp = POTLP_MIN(ratio, alphaTmp);
    
    return 1.0 / alphaTmp;
}

extern pot_int LpNewtonCreate( lp_newton **pnewton ) {
    
    pot_int retcode = RETCODE_OK;
    
    if ( !pnewton ) {
        retcode = RETCODE_FAILED;
        goto exit_cleanup;
    }
    
    lp_newton *nt = NULL;
    POTLP_INIT(nt, lp_newton, 1);
    
    if ( !nt ) {
        retcode = RETCODE_FAILED;
        goto exit_cleanup;
    }
    
    POTLP_ZERO(nt, lp_newton, 1);
    
    nt->beta = 0.995;
    nt->gamma = 0.7;
    
    *pnewton = nt;
    
exit_cleanup:
    return retcode;
}

extern pot_int LpNewtonInit( lp_newton *newton, pot_int nCol, pot_int nRow, pot_int *colMatBeg, pot_int *colMatIdx, double *colMatElem ) {
    
    pot_int retcode = RETCODE_OK;
    
    newton->nCol = nCol;
    newton->nRow = nRow;
    
    int nzA = colMatBeg[nCol];
    int ntCol = nRow + nCol;
    int ntNnz = nzA + nCol;
        
    POTLP_INIT(newton->AugBeg, int, ntCol + 1);
    POTLP_INIT(newton->AugIdx, int, ntNnz);
    POTLP_INIT(newton->AugElem, double, ntNnz);
    POTLP_INIT(newton->colBackup, double, nzA);
    
    /*
       [ I  AD']
       [ AD  0 ]
     */
    
    /* rowAux is auxiliary for matrix transposition */
    int *rowAux = NULL;
    POTLP_INIT(rowAux, int, nRow + 1);
    
    if ( !newton->AugBeg || !newton->AugIdx || !newton->AugElem ||
         !newton->colBackup || !rowAux ) {
        retcode = RETCODE_FAILED;
        goto exit_cleanup;
    }
    
    int *Ap = newton->AugBeg;
    int *Ai = newton->AugIdx;
    double *Ax = newton->AugElem;
    
    /* Assemble the identity */
    for ( int i = 0; i < nCol; ++i ) {
        Ap[i] = i; Ai[i] = i; Ax[i] = 1.0;
    }
    
    Ap += nCol; Ai += nCol; Ax += nCol;
    
    /* Transposition */
    for ( int i = 0; i < nzA; ++i ) {
        rowAux[colMatIdx[i] + 1] += 1;
    }
    
    for ( int i = 0; i < nRow; ++i ) {
        rowAux[i + 1] += rowAux[i];
    }
    
    POTLP_MEMCPY(Ap, rowAux, int, nRow + 1);
    for ( int i = 0; i < nRow + 1; ++i ) {
        Ap[i] += nCol;
    }
    
    for ( int i = 0, j, k; i < nCol; ++i ) {
        for ( j = colMatBeg[i]; j < colMatBeg[i + 1]; ++j ) {
            k = rowAux[colMatIdx[j]];
            Ai[k] = i; newton->colBackup[k] = colMatElem[j];
            rowAux[colMatIdx[j]] += 1;
        }
    }
    
    POT_CALL(potLinsysCreate(&newton->lpLinsys));
    POT_CALL(potLinsysInit(newton->lpLinsys, ntCol));
    POT_CALL(potLinsysSymFactorize(newton->lpLinsys, newton->AugBeg, newton->AugIdx));
    
    POTLP_INIT(newton->dd, double, nCol);
    POTLP_INIT(newton->xse, double, nCol);
    POTLP_INIT(newton->d1, double, nCol + nRow);
    POTLP_INIT(newton->d2, double, nCol + nRow);
    POTLP_INIT(newton->daux, double, nCol + nRow);
    
    POTLP_INIT(newton->dx, double, 2 * nCol + nRow);
    newton->dy = newton->dx + nCol;
    newton->ds = newton->dy + nRow;
    
    if ( !newton->dd || !newton->d1 || !newton->d2 || !newton->daux || !newton->dx ) {
        retcode = RETCODE_FAILED;
        goto exit_cleanup;
    }
    
exit_cleanup:
    POTLP_FREE(rowAux);
    return retcode;
}

/** @brief Implement one Newton's step
 * Currently no centering step is taken. On exit, colVal, rowDual, colDual, kappa, tau are modified by a Newton's step
 */
extern pot_int LpNewtonOneStep( lp_newton *newton, double *lpObj, double *lpRHS, int *colMatBeg, int *colMatIdx, double *colMatElem,
                                double *colVal, double *rowDual, double *colDual, double *kappa, double *tau,
                                double *pRes, double *dRes, double pObjVal, double dObjVal, double spxSize ) {
    
    pot_int retcode = RETCODE_OK;
    
    int nCol = newton->nCol;
    int nRow = newton->nRow;
    int nNt = nCol + nRow;
    
    /* Prepare iteration */
    double *b = lpRHS, *c = lpObj, kovert = *kappa / *tau;
    double *x = colVal, *y = rowDual, *s = colDual;
    double alpha = 0.0, beta = newton->beta, gamma = newton->gamma;
    double *D = newton->dd, *XSe = newton->xse, *d1 = newton->d1, *d2 = newton->d2, *daux = newton->daux;
    double *dx = newton->dx, *dy = newton->dy, *ds = newton->ds;
    double dkappa = 0.0, dtau = 0.0;
    double mu = LpNewtonIBarrier(nCol, x, s, kovert, 1.0);
    
    newton->mu = mu;
    
    /* Prepare Newton's system */
    int *ADBeg = newton->AugBeg + nCol;
    double *ADElem = newton->AugElem + nCol;
    
    /*
     XSe = sqrt(x .* s);
     D = sqrt(s) ./ sqrt(x);
     */
    
    for ( int i = 0; i < nCol; ++i ) {
        XSe[i] = x[i] * s[i];
        D[i] = sqrtl(s[i]) / sqrtl(x[i]);
    }
    
    /*
     ADinv = A * diag(Dinv);
     M = [speye(n), ADinv';
          ADinv, sparse(m, m)];
     */
    POTLP_MEMCPY(ADElem, newton->colBackup, double, colMatBeg[nCol]);
    for ( int i = 0, j; i < nCol; ++i ) {
        for ( j = ADBeg[i]; j < ADBeg[i + 1]; ++j ) {
            ADElem[j] /= D[i];
        }
    }
    
    POT_CALL(potLinsysNumFactorize(newton->lpLinsys, newton->AugBeg, newton->AugIdx, newton->AugElem));
    /*
     Dinvc = c ./ D;
     rhs1 = [-Dinvc; b];
     aux = [Dinvc; b];
     
     DinvXinvrmu1 = XSe - (mu * sigma) ./ XSe;
     rhs2 = [rd ./ D + DinvXinvrmu1; rp];
     */
    double coverd;
    double mugamma = mu * gamma;
    for ( int i = 0; i < nCol; ++i ) {
        coverd = c[i] / D[i];
        d1[i] = -coverd; daux[i] = coverd;
    }
    
    POTLP_MEMCPY(d1 + nCol, b, double, nRow);
    POTLP_MEMCPY(daux + nCol, b, double, nRow);
    
    for ( int i = 0; i < nCol; ++i ) {
        d2[i] = dRes[i] / D[i] + XSe[i] - mugamma / XSe[i];
    }
    
    POTLP_MEMCPY(d2, pRes, double, nRow);
    
    /* Solve the linear system
       d1 = M \ rhs1;
       d2 = M \ rhs2;
     */
    potLinsysSolve(newton->lpLinsys, d1, NULL);
    potLinsysSolve(newton->lpLinsys, d2, NULL);
    
    /*
     dtau = aux' * d2 + dobj - pobj - sigma * mu / tau;
     dtau = - dtau / (kappa / tau - aux' * d1);
     */
    double auxTd1 = dot(&nNt, daux, &potIntConstantOne, d1, &potIntConstantOne);
    double auxTd2 = dot(&nNt, daux, &potIntConstantOne, d2, &potIntConstantOne);
    dtau = auxTd2 + dObjVal - pObjVal - mugamma;
    dtau = - dtau / (kovert - auxTd1);
    
    /*
     dxdy = d1 * dtau - d2;
     dx = dxdy(1:n) ./ D;
     dy = - dxdy(n+1:end);
     dkappa = - kappa * dtau / tau - kappa + sigma * mu / tau;
     ds = -D.* dxdy(1:n) - s + (mu * sigma) ./ x;
     */
    
    for ( int i = 0; i < nNt; ++i ) {
        daux[i] = d1[i] * dtau - d2[i];
    }
    
    double auxoverd;
    for ( int i = 0; i < nCol; ++i ) {
        auxoverd = daux[i] / D[i];
        daux[i] = auxoverd;
        ds[i] = - auxoverd - s[i] + mugamma / x[i];
    }
    
    for ( int i = 0; i < nRow; ++i ) {
        dy[i] = -daux[i + nCol];
    }
    
    dkappa = -kovert * dtau - kovert + mugamma;
    alpha = LpNewtonIRatioTest(nCol, x, dx, s, ds, kovert, dkappa, 1.0, dtau);
    alpha = alpha * beta;
    
    /*
     alpha = 0.995 / abs(min([dx./x; ds./s; dkappa / kappa; dtau./tau]));
     x = x + alpha * dx;
     y = y + alpha * dy;
     s = s + alpha * ds;
     kappa = kappa + alpha * dkappa;
     tau = tau + alpha * dtau;
     */
    axpy(&nCol, &alpha, dx, &potIntConstantOne, x, &potIntConstantOne);
    axpy(&nRow, &alpha, dy, &potIntConstantOne, y, &potIntConstantOne);
    axpy(&nCol, &alpha, ds, &potIntConstantOne, s, &potIntConstantOne);
    kovert += alpha * dkappa;
    double newtau = 1.0 + alpha * dtau;
    
    /* Back to the simplex */
    double simp = kovert + newtau;
    for ( int i = 0; i < nCol; ++i ) {
        simp += x[i];
        simp += s[i];
    }
    
    simp = spxSize / simp;
    scal(&nCol, &simp, x, &potIntConstantOne);
    scal(&nRow, &simp, y, &potIntConstantOne);
    scal(&nCol, &simp, s, &potIntConstantOne);
    kovert *= simp; newtau *= simp;
    *kappa = kovert; *tau = newtau;
    
exit_cleanup:
    return retcode;
}

extern void LpNewtonClear( lp_newton *newton ) {
    
    if ( !newton ) {
        return;
    }
    
    potLinsysDestroy(&newton->lpLinsys);
    
    POTLP_FREE(newton->AugBeg);
    POTLP_FREE(newton->AugIdx);
    POTLP_FREE(newton->AugElem);
    POTLP_FREE(newton->colBackup);
    
    POTLP_FREE(newton->dd);
    POTLP_FREE(newton->xse);
    POTLP_FREE(newton->d1);
    POTLP_FREE(newton->d2);
    POTLP_FREE(newton->daux);
    POTLP_FREE(newton->dx);
    POTLP_FREE(newton->dy);
    POTLP_FREE(newton->ds);
    
    POTLP_ZERO(newton, lp_newton, 1);
    return;
}

extern void LpNewtonDestroy( lp_newton **pnewton ) {
    
    if ( !pnewton ) {
        return;
    }
    
    LpNewtonClear(*pnewton);
    POTLP_FREE(*pnewton);
    
    return;
}
