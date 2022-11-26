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
extern void dataMatScalSparse( void *A, double alpha ) {
    
    sdpSparseData *spA = (sdpSparseData *)A;
    int incx = 1;
    
    scal(&spA->nTriMatElem, &alpha, spA->triMatElem, &incx);
    
    return;
}

/** @brief A = alpha \* A for dense matrix A
 *  @param[in] A sdpDenseData pointer
 *  @param[in] alpha scale parameter
 */
extern void dataMatScalDense( void *A, double alpha ) {
    
    sdpDenseData *dnsA = (sdpDenseData *)A;
    int nElem = dnsA->nSDPCol * (dnsA->nSDPCol + 1) / 2, incx = 1;
    
    scal(&nElem, &alpha, dnsA->dsMatElem, &incx);
    
    return;
}

/** @brief A = alpha \* A for rank one sparse matrix A
 *  @param[in] A sdpRankOneSparseData pointer
 *  @param[in] alpha scale parameter
 */
extern void dataMatScalRankOneSparse( void *A, double alpha ) {
    
    assert( fabs(alpha) >= 1e-8 );
    
    sdpRankOneSparseData *spR1A = (sdpRankOneSparseData *)A;
    int incx = 1;
    
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
extern void dataMatScalRankOneDense( void *A, double alpha ) {
    
    assert( fabs(alpha) >= 1e-8 );
    
    sdpRankOneDenseData *dnsR1A = (sdpRankOneDenseData *)A;
    int incx = 1;
    
    if ( alpha < 0.0 ) {
        dnsR1A->r1FactorSign = -dnsR1A->r1FactorSign;
        alpha = -alpha;
    }
    alpha = sqrt(alpha);
    
    scal(&dnsR1A->nSDPCol, &alpha, dnsR1A->r1MatFactor, &incx);
    
    return;
}

/** @brief Calculate norm of sparse matrix A
 *  @param[in] A sdpSparseData pointer
 *  @param[in] type type of norm
 *  @return norm of A
 */
extern double dataMatNormSparse( void *A, int type ) {
    
    sdpSparseData *spA = (sdpSparseData *)A;
    double nrmA = 0.0;
    int i;
    
    if (type == 2) {
        for ( i = 0; i < spA->nTriMatElem; ++i ) {
            if ( spA->triMatCol[i] == spA->triMatRow[i] ) {
                nrmA += spA->triMatElem[i] * spA->triMatElem[i];
            } else {
                nrmA += spA->triMatElem[i] * spA->triMatElem[i] * 2;
            }
        }
        nrmA = sqrt(nrmA);
    } else if (type == 1) {
        for ( i = 0; i < spA->nTriMatElem; ++i ) {
            if ( spA->triMatCol[i] == spA->triMatRow[i] ) {
                nrmA += fabs(spA->triMatElem[i]);
            } else {
                nrmA += fabs(spA->triMatElem[i]) * 2;
            }
        }
    }
    
    return nrmA;
}

/** @brief Calculate norm of dense matrix A
 *  @param[in] A sdpDenseData pointer
 *  @param[in] type type of norm
 *  @return norm of A
 */
extern double dataMatNormDense( void *A, int type ) {
    
    sdpDenseData *dnsA = (sdpDenseData *)A;
    double nrmA = 0.0;
    int incx = 1;
    int i, j, colStart, colLen;
    
    if (type == 2) {
        for ( i = 0; i < dnsA->nSDPCol; ++i ) {
            colStart = i * (dnsA->nSDPCol * 2 - i + 1) / 2;
            colLen = dnsA->nSDPCol - i;
            nrmA += dnsA->dsMatElem[colStart] * dnsA->dsMatElem[colStart];
            for ( j = i + 1; j < colLen; ++j ) {
                nrmA += dnsA->dsMatElem[colStart + j] * dnsA->dsMatElem[colStart + j] * 2;
            }
        }
        nrmA = sqrt(nrmA);
    } else if (type == 1) {
        for ( i = 0; i < dnsA->nSDPCol; ++i ) {
            colStart = i * (dnsA->nSDPCol * 2 - i + 1) / 2;
            colLen = dnsA->nSDPCol - i - 1;
            nrmA += fabs(dnsA->dsMatElem[colStart]) + sum1(&colLen, dnsA->dsMatElem + colStart + 1, &incx) * 2;
        }
    }
    
    return nrmA;
}

/** @brief Calculate norm of rank one sparse matrix A
 *  @param[in] A sdpRankOneSparseData pointer
 *  @param[in] type type of norm
 *  @return norm of A
 */
extern double dataMatNormRankOneSparse( void *A, int type ) {
    
    sdpRankOneSparseData *spR1A = (sdpRankOneSparseData *)A;
    double nrmA = 0.0;
    int incx = 1;
    
    if (type == 2) {
        nrmA = nrm2(&spR1A->nSpR1FactorElem, spR1A->spR1MatElem, &incx);
        nrmA = nrmA * nrmA;
    } else if (type == 1) {
        nrmA = sum1(&spR1A->nSpR1FactorElem, spR1A->spR1MatElem, &incx);
        nrmA = nrmA * nrmA;
    }
    
    return nrmA;
}

/** @brief Calculate norm of rank one densematrix A
 *  @param[in] A sdpRankOneDenseData pointer
 *  @param[in] type type of norm
 *  @return norm of A
 */
extern double dataMatNormRankOneDense( void *A, int type ) {
    
    sdpRankOneDenseData *dnsR1A = (sdpRankOneDenseData *)A;
    double nrmA = 0.0;
    int incx = 1;
    
    if (type == 2) {
        nrmA = nrm2(&dnsR1A->nSDPCol, dnsR1A->r1MatFactor, &incx);
        nrmA = nrmA * nrmA;
    } else if (type == 1) {
        nrmA = sum1(&dnsR1A->nSDPCol, dnsR1A->r1MatFactor, &incx);
        nrmA = nrmA * nrmA;
    }
    
    return nrmA;
}

/** @brief Calculate number of nonzero elements in sparse matrix A
 *  @param[in] A sdpSparseData pointer
 *  @return number of nonzero elements in A
 */
extern int dataMatGetNnzSparse( void *A ) {
    
    sdpSparseData *spA = (sdpSparseData *)A;
    
    return spA->nTriMatElem;
}

/** @brief Calculate number of nonzero elements in dense matrix A
 *  @param[in] A sdpDenseData pointer
 *  @return number of nonzero elements in A
 */
extern int dataMatGetNnzDense( void *A ) {
    
    sdpDenseData *dnsA = (sdpDenseData *)A;
    
    return dnsA->nSDPCol * (dnsA->nSDPCol + 1) / 2;
}

/** @brief Calculate number of nonzero elements in rank one sparse matrix A
 *  @param[in] A sdpRankOneSparseData pointer
 *  @return number of nonzero elements in A
 */
extern int dataMatGetNnzRankOneSparse( void *A ) {
    
    sdpRankOneSparseData *spR1A = (sdpRankOneSparseData *)A;
    
    return spR1A->nSpR1FactorElem * (spR1A->nSpR1FactorElem + 1) / 2;
}

/** @brief Calculate number of nonzero elements in rank one dense matrix A
 *  @param[in] A sdpRankOneDenseData pointer
 *  @return number of nonzero elements in A
 */
extern int dataMatGetNnzRankOneDense( void *A ) {
    
    sdpRankOneDenseData *dnsR1A = (sdpRankOneDenseData *)A;
    
    return dnsR1A->nSDPCol * (dnsR1A->nSDPCol + 1) / 2;
}

/** @brief Construct full matrix from sparse matrix A
 *  @param[in] A sdpSparseData pointer
 *  @param[out] v store full matrix
 */
extern void dataMatDumpSparse( void *A, double *v ) {
    
    sdpSparseData *spA = (sdpSparseData *)A;
    HDSDP_ZERO(v, double, spA->nSDPCol * spA->nSDPCol);
    
    for ( int i = 0; i < spA->nTriMatElem; ++i ) {
        v[spA->triMatCol[i] * spA->nSDPCol + spA->triMatRow[i]] = spA->triMatElem[i];
        v[spA->triMatRow[i] * spA->nSDPCol + spA->triMatCol[i]] = spA->triMatElem[i];
    }
    
    return;
}

/** @brief Construct full matrix from dense matrix A
 *  @param[in] A sdpDenseData pointer
 *  @param[out] v store full matrix
 */
extern void dataMatDumpDense( void *A, double *v ) {
    
    sdpDenseData *dnsA = (sdpDenseData *)A;
    int nElem = dnsA->nSDPCol * (dnsA->nSDPCol + 1) / 2;
    HDSDP_MEMCPY(v, dnsA->dsMatElem, double, nElem);
    int i, j;
    
    for ( i = 0; i < dnsA->nSDPCol; ++i) {
        for ( j = i; j < dnsA->nSDPCol; ++j ) {
            v[i * dnsA->nSDPCol + j] = dnsA->dsMatElem[i * (dnsA->nSDPCol * 2 - i - 1) / 2 + j];
            v[j * dnsA->nSDPCol + i] = v[i * dnsA->nSDPCol + j];
        }
    }
    
    return;
}

/** @brief Construct full matrix from rank one sparse matrix A
 *  @param[in] A sdpRankOneSparseData pointer
 *  @param[out] v store full matrix
 */
extern void dataMatDumpRankOneSparse( void *A, double *v ) {
    
    sdpRankOneSparseData *spR1A = (sdpRankOneSparseData *)A;
    HDSDP_ZERO(v, double, spR1A->nSDPCol * spR1A->nSDPCol);
    
    int i, j;
    
    for ( i = 0; i < spR1A->nSpR1FactorElem; ++i ) { // colIdx
        for ( j = i; j < spR1A->nSpR1FactorElem; ++j ) { // rowIdx
            v[spR1A->spR1MatIdx[i] * spR1A->nSDPCol + spR1A->spR1MatIdx[j]] = spR1A->spR1FactorSign * spR1A->spR1MatElem[i] * spR1A->spR1MatElem[j];
            v[spR1A->spR1MatIdx[j] * spR1A->nSDPCol + spR1A->spR1MatIdx[i]] = v[spR1A->spR1MatIdx[i] * spR1A->nSDPCol + spR1A->spR1MatIdx[j]];
        }
    }
    
    return;
}

/** @brief Construct full matrix from rank one dense matrix A
 *  @param[in] A sdpRankOneDenseData pointer
 *  @param[out] v store full matrix
 */
extern void dataMatDumpRankOneDense( void *A, double *v ) {
    
    sdpRankOneDenseData *dnsR1A = (sdpRankOneDenseData *)A;
    
    int i, j;
    
    for ( i = 0; i < dnsR1A->nSDPCol; ++i ) { // colIdx
        for ( j = i; j < dnsR1A->nSDPCol; ++j ) { // rowIdx
            v[i * dnsR1A->nSDPCol + j] = dnsR1A->r1FactorSign * dnsR1A->r1MatFactor[i] * dnsR1A->r1MatFactor[j];
            v[j * dnsR1A->nSDPCol + i] = v[i * dnsR1A->nSDPCol + j];
        }
    }
    
    return;
}
