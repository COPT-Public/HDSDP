/** @file hdsdp\_sdpdata.h
 *  @brief Implenent the SDP Coeffcient data operations
 *  @date 11/24/2022
 */

#include "hdsdp_sdpdata.h"
#include "vec_opts.h"
#include "hdsdp_utils.h"

#include <math.h>
#include <assert.h>

/** @brief A = alpha \* A for sparse matrix A
 *  @param[in] A sdpSparseData pointer
 *  @param[in] alpha scale parameter
 */
extern void dataMatScalSparseImpl( void *A, double alpha ) {
    
    sdpSparseData *spA = (sdpSparseData *) A;
    int incx = 1;
    
    scal(&spA->nTriMatElem, &alpha, spA->triMatElem, &incx);
    
    return;
}

/** @brief A = alpha \* A for dense matrix A
 *  @param[in] A sdpDenseData pointer
 *  @param[in] alpha scale parameter
 */
extern void dataMatScalDenseImpl( void *A, double alpha ) {
    
    sdpDenseData *dsA = (sdpDenseData *) A;
    int nElem = PACK_NNZ(dsA->nSDPCol);
    int incx = 1;
    
    scal(&nElem, &alpha, dsA->dsMatElem, &incx);
    
    return;
}

/** @brief A = alpha \* A for rank one sparse matrix A
 *  @param[in] A sdpRankOneSparseData pointer
 *  @param[in] alpha scale parameter
 */
extern void dataMatScalRankOneSparseImpl( void *A, double alpha ) {
    
    /* avoid alpha = 0.0 case */
    assert( fabs(alpha) >= 1e-8 );
    
    sdpRankOneSparseData *spR1A = (sdpRankOneSparseData *) A;
    int incx = 1;
    
    /* check whether sign changes */
    if ( alpha < 0.0 ) {
        spR1A->spR1FactorSign = -spR1A->spR1FactorSign;
        alpha = -alpha;
    }
    alpha = sqrt(alpha);
    
    scal(&spR1A->nSpR1FactorElem, &alpha, spR1A->spR1MatElem, &incx);
    
    return;
}

/** @brief A = alpha \* A for rank one dense matrix A
 *  @param[in] A sdpRankOneDenseData pointer
 *  @param[in] alpha scale parameter
 */
extern void dataMatScalRankOneDenseImpl( void *A, double alpha ) {
    
    /* avoid alpha = 0.0 case */
    assert( fabs(alpha) >= 1e-8 );
    
    sdpRankOneDenseData *dsR1A = (sdpRankOneDenseData *) A;
    int incx = 1;
    
    /* check whether sign changes */
    if ( alpha < 0.0 ) {
        dsR1A->r1FactorSign = -dsR1A->r1FactorSign;
        alpha = -alpha;
    }
    alpha = sqrt(alpha);
    
    scal(&dsR1A->nSDPCol, &alpha, dsR1A->r1MatFactor, &incx);
    
    return;
}

/** @brief Calculate norm of sparse matrix A
 *  @param[in] A sdpSparseData pointer
 *  @param[in] type type of norm
 *  @return norm of A
 */
extern double dataMatNormSparseImpl( void *A, int type ) {
    
    sdpSparseData *spA = (sdpSparseData *) A;
    double nrmA = 0.0;
    
    if ( type == F_NORM ) {
        for ( int i = 0; i < spA->nTriMatElem; ++i ) {
            nrmA += ( spA->triMatCol[i] == spA->triMatRow[i] ) ?
                      spA->triMatElem[i] * spA->triMatElem[i] :
                      spA->triMatElem[i] * spA->triMatElem[i] * 2;
        }
        nrmA = sqrt(nrmA);
    } else if ( type == ABS_NORM ) {
        for ( int i = 0; i < spA->nTriMatElem; ++i ) {
            nrmA += ( spA->triMatCol[i] == spA->triMatRow[i] ) ?
                      fabs(spA->triMatElem[i]) :
                      fabs(spA->triMatElem[i]) * 2;
        }
    }
    
    return nrmA;
}

/** @brief Calculate norm of dense matrix A
 *  @param[in] A sdpDenseData pointer
 *  @param[in] type type of norm
 *  @return norm of A
 */
extern double dataMatNormDenseImpl( void *A, int type ) {
    
    sdpDenseData *dsA = (sdpDenseData *) A;
    double nrmA = 0.0;
    int nCol = dsA->nSDPCol;
    int incx = 1;
    int colLen; ///< exclude diagonal
    double colNrm; ///< exclude diagonal
    double *p = dsA->dsMatElem;
    
    if ( type == F_NORM ) {
        for ( int i = 0; i < nCol; ++i ) {
            nrmA += (*p) * (*p);
            colLen = nCol - i - 1;
            colNrm = nrm2(&colLen, p + 1, &incx);
            nrmA += colNrm * colNrm * 2;
            p = p + colLen + 1;
        }
        nrmA = sqrt(nrmA);
    } else if ( type == ABS_NORM ) {
        for ( int i = 0; i < nCol; ++i ) {
            nrmA += fabs(*p);
            colLen = nCol - i - 1;
            nrmA += nrm1(&colLen, p + 1, &incx) * 2;
            p = p + colLen + 1;
        }
    }
    
    return nrmA;
}

/** @brief Calculate norm of rank one sparse matrix A
 *  @param[in] A sdpRankOneSparseData pointer
 *  @param[in] type type of norm
 *  @return norm of A
 */
extern double dataMatNormRankOneSparseImpl( void *A, int type ) {
    
    sdpRankOneSparseData *spR1A = (sdpRankOneSparseData *) A;
    double nrmA = 0.0;
    int incx = 1;
    
    if ( type == F_NORM ) {
        nrmA = nrm2(&spR1A->nSpR1FactorElem, spR1A->spR1MatElem, &incx);
        nrmA = nrmA * nrmA;
    } else if ( type == ABS_NORM ) {
        nrmA = nrm1(&spR1A->nSpR1FactorElem, spR1A->spR1MatElem, &incx);
        nrmA = nrmA * nrmA;
    }
    
    return nrmA;
}

/** @brief Calculate norm of rank one densematrix A
 *  @param[in] A sdpRankOneDenseData pointer
 *  @param[in] type type of norm
 *  @return norm of A
 */
extern double dataMatNormRankOneDenseImpl( void *A, int type ) {
    
    sdpRankOneDenseData *dsR1A = (sdpRankOneDenseData *) A;
    double nrmA = 0.0;
    int incx = 1;
    
    if ( type == F_NORM ) {
        nrmA = nrm2(&dsR1A->nSDPCol, dsR1A->r1MatFactor, &incx);
        nrmA = nrmA * nrmA;
    } else if ( type == ABS_NORM ) {
        nrmA = nrm1(&dsR1A->nSDPCol, dsR1A->r1MatFactor, &incx);
        nrmA = nrmA * nrmA;
    }
    
    return nrmA;
}

/** @brief Calculate number of nonzero elements in sparse matrix A
 *  @param[in] A sdpSparseData pointer
 *  @return number of nonzero elements in A
 */
extern int dataMatGetNnzSparseImpl( void *A ) {
    
    sdpSparseData *spA = (sdpSparseData *) A;
    
    return spA->nTriMatElem;
}

/** @brief Calculate number of nonzero elements in dense matrix A
 *  @param[in] A sdpDenseData pointer
 *  @return number of nonzero elements in A
 */
extern int dataMatGetNnzDenseImpl( void *A ) {
    
    sdpDenseData *dsA = (sdpDenseData *) A;
    
    return PACK_NNZ(dsA->nSDPCol);
}

/** @brief Calculate number of nonzero elements in rank one sparse matrix A
 *  @param[in] A sdpRankOneSparseData pointer
 *  @return number of nonzero elements in A
 */
extern int dataMatGetNnzRankOneSparseImpl( void *A ) {
    
    sdpRankOneSparseData *spR1A = (sdpRankOneSparseData *) A;
    
    return PACK_NNZ(spR1A->nSpR1FactorElem);
}

/** @brief Calculate number of nonzero elements in rank one dense matrix A
 *  @param[in] A sdpRankOneDenseData pointer
 *  @return number of nonzero elements in A
 */
extern int dataMatGetNnzRankOneDenseImpl( void *A ) {
    
    sdpRankOneDenseData *dsR1A = (sdpRankOneDenseData *) A;
    
    return PACK_NNZ(dsR1A->nSDPCol);
}

/** @brief Construct full matrix from sparse matrix A
 *  @param[in] A sdpSparseData pointer
 *  @param[out] v store full matrix
 */
extern void dataMatDumpSparseImpl( void *A, double *v ) {
    
    sdpSparseData *spA = (sdpSparseData *) A;
    HDSDP_ZERO(v, double, spA->nSDPCol * spA->nSDPCol);
    int nCol = spA->nSDPCol;
    
    for ( int i = 0; i < spA->nTriMatElem; ++i ) {
        v[spA->triMatCol[i] * nCol + spA->triMatRow[i]] = spA->triMatElem[i];
        v[spA->triMatRow[i] * nCol + spA->triMatCol[i]] = spA->triMatElem[i];
    }
    
    return;
}

/** @brief Construct full matrix from dense matrix A
 *  @param[in] A sdpDenseData pointer
 *  @param[out] v store full matrix
 */
extern void dataMatDumpDenseImpl( void *A, double *v ) {
    
    sdpDenseData *dsA = (sdpDenseData *) A;
    int nCol = dsA->nSDPCol;
    int dsIdx; // index for (col, row) = (i, j) in dsMatElem
    
    for ( int i = 0, j; i < nCol; ++i) { // colIdx
        for ( j = i; j < nCol; ++j ) { // rowIdx
            dsIdx = i * (nCol * 2 - i - 1) / 2 + j;
            v[i * nCol + j] = dsA->dsMatElem[dsIdx];
            v[j * nCol + i] = v[i * nCol + j];
        }
    }
    
    return;
}

/** @brief Construct full matrix from rank one sparse matrix A
 *  @param[in] A sdpRankOneSparseData pointer
 *  @param[out] v store full matrix
 */
extern void dataMatDumpRankOneSparseImpl( void *A, double *v ) {
    
    sdpRankOneSparseData *spR1A = (sdpRankOneSparseData *) A;
    HDSDP_ZERO(v, double, spR1A->nSDPCol * spR1A->nSDPCol);
    int nCol = spR1A->nSDPCol;
    int *spR1MatIdx = spR1A->spR1MatIdx;
    
    for ( int i = 0, j; i < spR1A->nSpR1FactorElem; ++i ) { // colIdx
        for ( j = i; j < spR1A->nSpR1FactorElem; ++j ) { // rowIdx
            v[spR1MatIdx[i] * nCol + spR1MatIdx[j]] = \
            spR1A->spR1FactorSign * spR1A->spR1MatElem[i] * spR1A->spR1MatElem[j];
            if ( j > i ) {
                v[spR1MatIdx[j] * nCol + spR1MatIdx[i]] = \
                v[spR1MatIdx[i] * nCol + spR1MatIdx[j]];
            }
        }
    }
    
    return;
}

/** @brief Construct full matrix from rank one dense matrix A
 *  @param[in] A sdpRankOneDenseData pointer
 *  @param[out] v store full matrix
 */
extern void dataMatDumpRankOneDenseImpl( void *A, double *v ) {
    
    sdpRankOneDenseData *dsR1A = (sdpRankOneDenseData *) A;
    
    for ( int i = 0, j; i < dsR1A->nSDPCol; ++i ) { // colIdx
        for ( j = i; j < dsR1A->nSDPCol; ++j ) { // rowIdx
            v[i * dsR1A->nSDPCol + j] = \
            dsR1A->r1FactorSign * dsR1A->r1MatFactor[i] * dsR1A->r1MatFactor[j];
            if ( j > i ) {
                v[j * dsR1A->nSDPCol + i] = v[i * dsR1A->nSDPCol + j];
            }
        }
    }
    
    return;
}
