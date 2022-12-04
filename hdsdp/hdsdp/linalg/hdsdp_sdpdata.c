/** @file hdsdp\_sdpdata.h
 *  @brief Implenent the SDP Coeffcient data operations
 *  @date 11/24/2022
 */

#include "hdsdp_sdpdata.h"
#include "vec_opts.h"
#include "hdsdp_utils.h"

#include <math.h>

extern hdsdp_retcode sdpDataMatCreate( sdp_coeff **psdpCoeff ) {
    
    hdsdp_retcode retcode = HDSDP_RETCODE_OK;
    
    if ( !psdpCoeff ) {
        retcode = HDSDP_RETCODE_FAILED;
        goto exit_cleanup;
    }
    
    sdp_coeff *sdpCoeff = NULL;
    HDSDP_INIT(sdpCoeff, sdp_coeff, 1);
    HDSDP_MEMCHECK(sdpCoeff);
    
    *psdpCoeff = sdpCoeff;
    
exit_cleanup:
    
    return retcode;
}

#define SPARSE_THRESHOLD (0.3)

/** @brief Set SDP data matrix. This routine selects data type and assiciate structure with their operations
 * 
 *
 */
extern hdsdp_retcode sdpDataMatSetData( sdp_coeff *sdpCoeff, int nSDPCol, int dataMatNnz, int *dataMatIdx, double *dataMatElem ) {
    
    hdsdp_retcode retcode = HDSDP_RETCODE_OK;
    
    
    
    
    
    
exit_cleanup:
    
    return retcode;
}

/* Wrapper for the operations */
extern void sdpDataMataApB( sdp_coeff *A, double alpha, void *B ) {
    
    A->dataMataApB(A->dataMat, alpha, B);
    
    return;
}

/* TODO: Wrapper for the other operations */


/* TODO: The name dataMat may be changed for SDP context */
extern void dataMatScalZeroImpl( void *A, double alpha ) {
    
    return;
}

/** @brief A = alpha \* A for sparse matrix A
 *  @param[in] A sdpSparseData pointer
 *  @param[in] alpha Scale parameter
 */
extern void dataMatScalSparseImpl( void *A, double alpha ) {
    
    sdp_coeff_sparse *spA = (sdp_coeff_sparse *) A;
    int incx = 1;
    
    scal(&spA->nTriMatElem, &alpha, spA->triMatElem, &incx);
    
    return;
}

/** @brief A = alpha \* A for dense matrix A
 *  @param[in] A sdpDenseData pointer
 *  @param[in] alpha Scale parameter
 */
extern void dataMatScalDenseImpl( void *A, double alpha ) {
    
    sdp_coeff_dense *dsA = (sdp_coeff_dense *) A;
    int nElem = PACK_NNZ(dsA->nSDPCol);
    int incx = 1;
    
    scal(&nElem, &alpha, dsA->dsMatElem, &incx);
    
    return;
}

/** @brief A = alpha \* A for rank one sparse matrix A
 *  @param[in] A sdpRankOneSparseData pointer
 *  @param[in] alpha Scale parameter
 */
extern void dataMatScalRankOneSparseImpl( void *A, double alpha ) {
    
    /* Avoid alpha = 0.0 case */
    assert( fabs(alpha) >= 1e-08 );
    
    sdp_coeff_spr1 *spR1A = (sdp_coeff_spr1 *) A;
    
    spR1A->spR1FactorSign = spR1A->spR1FactorSign * alpha;
    
    return;
}

/** @brief A = alpha \* A for rank one dense matrix A
 *  @param[in] A sdpRankOneDenseData pointer
 *  @param[in] alpha Scale parameter
 */
extern void dataMatScalRankOneDenseImpl( void *A, double alpha ) {
    
    /* Avoid alpha = 0.0 case */
    assert( fabs(alpha) >= 1e-08 );
    
    sdp_coeff_dsr1 *dsR1A = (sdp_coeff_dsr1 *) A;
    dsR1A->r1FactorSign = dsR1A->r1FactorSign * alpha;
    
    return;
}

extern double dataMatNormZeroImpl( void *A, int type ) {
    
    return 0.0;
}

/** @brief Calculate norm of sparse matrix A
 *  @param[in] A sdpSparseData pointer
 *  @param[in] type Type of norm
 *  @return Norm of A
 */
extern double dataMatNormSparseImpl( void *A, int type ) {
    
    sdp_coeff_sparse *spA = (sdp_coeff_sparse *) A;
    double nrmA = 0.0;
    
    if ( type == FRO_NORM ) {
        for ( int i = 0; i < spA->nTriMatElem; ++i ) {
            nrmA += ( spA->triMatCol[i] == spA->triMatRow[i] ) ?
                      spA->triMatElem[i] * spA->triMatElem[i] :
                      spA->triMatElem[i] * spA->triMatElem[i] * 2.0;
        }
        nrmA = sqrt(nrmA);
    } else if ( type == ABS_NORM ) {
        for ( int i = 0; i < spA->nTriMatElem; ++i ) {
            nrmA += ( spA->triMatCol[i] == spA->triMatRow[i] ) ?
                      fabs(spA->triMatElem[i]) :
                      fabs(spA->triMatElem[i]) * 2.0;
        }
    }
    
    return nrmA;
}

/** @brief Calculate norm of dense matrix A
 *  @param[in] A sdpDenseData pointer
 *  @param[in] type Type of norm
 *  @return Norm of A
 */
extern double dataMatNormDenseImpl( void *A, int type ) {
    
    sdp_coeff_dense *dsA = (sdp_coeff_dense *) A;
    double nrmA = 0.0;
    int nCol = dsA->nSDPCol;
    int incx = 1;
    int colLen; ///< Exclude diagonal
    double colNrm; ///< Exclude diagonal
    double *p = dsA->dsMatElem;
    
    if ( type == FRO_NORM ) {
        for ( int i = 0; i < nCol; ++i ) {
            nrmA += p[0] * p[0];
            colLen = nCol - i - 1;
            colNrm = nrm2(&colLen, p + 1, &incx);
            nrmA += colNrm * colNrm * 2.0;
            p = p + colLen + 1;
        }
        nrmA = sqrt(nrmA);
    } else if ( type == ABS_NORM ) {
        for ( int i = 0; i < nCol; ++i ) {
            nrmA += fabs(p[0]);
            colLen = nCol - i - 1;
            nrmA += nrm1(&colLen, p + 1, &incx) * 2;
            p = p + colLen + 1;
        }
    }
    
    return nrmA;
}

/** @brief Calculate norm of rank one sparse matrix A
 *  @param[in] A sdpRankOneSparseData pointer
 *  @param[in] type Type of norm
 *  @return Norm of A
 */
extern double dataMatNormRankOneSparseImpl( void *A, int type ) {
    
    sdp_coeff_spr1 *spR1A = (sdp_coeff_spr1 *) A;
    double nrmA = 0.0;
    int incx = 1;
    
    if ( type == FRO_NORM ) {
        nrmA = nrm2(&spR1A->nSpR1FactorElem, spR1A->spR1MatElem, &incx);
        nrmA = nrmA * nrmA * fabs(spR1A->spR1FactorSign);
    } else if ( type == ABS_NORM ) {
        nrmA = nrm1(&spR1A->nSpR1FactorElem, spR1A->spR1MatElem, &incx);
        nrmA = nrmA * nrmA * fabs(spR1A->spR1FactorSign);
    }
    
    return nrmA;
}

/** @brief Calculate norm of rank one densematrix A
 *  @param[in] A sdpRankOneDenseData pointer
 *  @param[in] type Type of norm
 *  @return Norm of A
 */
extern double dataMatNormRankOneDenseImpl( void *A, int type ) {
    
    sdp_coeff_dsr1 *dsR1A = (sdp_coeff_dsr1 *) A;
    double nrmA = 0.0;
    int incx = 1;
    
    if ( type == FRO_NORM ) {
        nrmA = nrm2(&dsR1A->nSDPCol, dsR1A->r1MatFactor, &incx);
        nrmA = nrmA * nrmA * fabs(dsR1A->r1FactorSign);
    } else if ( type == ABS_NORM ) {
        nrmA = nrm1(&dsR1A->nSDPCol, dsR1A->r1MatFactor, &incx);
        nrmA = nrmA * nrmA * fabs(dsR1A->r1FactorSign);
    }
    
    return nrmA;
}

extern int dataMatGetNnzZeroImpl( void *A ) {
    
    return 0;
}

/** @brief Calculate number of nonzero elements in sparse matrix A
 *  @param[in] A sdpSparseData pointer
 *  @return Number of nonzero elements in A
 */
extern int dataMatGetNnzSparseImpl( void *A ) {
    
    sdp_coeff_sparse *spA = (sdp_coeff_sparse *) A;
    
    return spA->nTriMatElem;
}

/** @brief Calculate number of nonzero elements in dense matrix A
 *  @param[in] A sdpDenseData pointer
 *  @return Number of nonzero elements in A
 */
extern int dataMatGetNnzDenseImpl( void *A ) {
    
    sdp_coeff_dense *dsA = (sdp_coeff_dense *) A;
    
    return PACK_NNZ(dsA->nSDPCol);
}

/** @brief Calculate number of nonzero elements in rank one sparse matrix A
 *  @param[in] A sdpRankOneSparseData pointer
 *  @return Number of nonzero elements in A
 */
extern int dataMatGetNnzRankOneSparseImpl( void *A ) {
    
    sdp_coeff_spr1 *spR1A = (sdp_coeff_spr1 *) A;
    
    return PACK_NNZ(spR1A->nSpR1FactorElem);
}

/** @brief Calculate number of nonzero elements in rank one dense matrix A
 *  @param[in] A sdpRankOneDenseData pointer
 *  @return Number of nonzero elements in A
 */
extern int dataMatGetNnzRankOneDenseImpl( void *A ) {
    
    sdp_coeff_dsr1 *dsR1A = (sdp_coeff_dsr1 *) A;
    
    return PACK_NNZ(dsR1A->nSDPCol);
}

extern void dataMatDumpZeroImpl( void *A, double *v ) {
    
    sdp_coeff_zero *zeroA = (sdp_coeff_zero *) A;
    HDSDP_ZERO(v, double, zeroA->nSDPCol * zeroA->nSDPCol);
    
    return;
}

/** @brief Construct full matrix from sparse matrix A
 *  @param[in] A sdpSparseData pointer
 *  @param[out] v Store full matrix
 */
extern void dataMatDumpSparseImpl( void *A, double *v ) {
    
    sdp_coeff_sparse *spA = (sdp_coeff_sparse *) A;
    HDSDP_ZERO(v, double, spA->nSDPCol * spA->nSDPCol);
    int nCol = spA->nSDPCol;
    
    for ( int i = 0; i < spA->nTriMatElem; ++i ) {
        v[spA->triMatCol[i] * nCol + spA->triMatRow[i]] = spA->triMatElem[i];
        if ( spA->triMatCol[i] != spA->triMatRow[i] ) {
            v[spA->triMatRow[i] * nCol + spA->triMatCol[i]] = spA->triMatElem[i];
        }
    }
    
    return;
}

/** @brief Construct full matrix from dense matrix A
 *  @param[in] A sdpDenseData pointer
 *  @param[out] v Store full matrix
 */
extern void dataMatDumpDenseImpl( void *A, double *v ) {
    
    sdp_coeff_dense *dsA = (sdp_coeff_dense *) A;
    int nCol = dsA->nSDPCol;
    int dsIdx; ///< Index for (col, row) = (i, j) in dsMatElem
    
    for ( int i = 0, j; i < nCol; ++i) { // Column index
        dsIdx = i * (nCol * 2 - i + 1) / 2;
        v[i * nCol + i] = dsA->dsMatElem[dsIdx];
        for ( j = i + 1; j < nCol; ++j ) { // Row index
            dsIdx = i * (nCol * 2 - i - 1) / 2 + j;
            v[i * nCol + j] = dsA->dsMatElem[dsIdx];
            v[j * nCol + i] = v[i * nCol + j];
        }
    }
    
    return;
}

/** @brief Construct full matrix from rank one sparse matrix A
 *  @param[in] A sdpRankOneSparseData pointer
 *  @param[out] v Store full matrix
 */
extern void dataMatDumpRankOneSparseImpl( void *A, double *v ) {
    
    sdp_coeff_spr1 *spR1A = (sdp_coeff_spr1 *) A;
    HDSDP_ZERO(v, double, spR1A->nSDPCol * spR1A->nSDPCol);
    int nCol = spR1A->nSDPCol;
    int *spR1MatIdx = spR1A->spR1MatIdx;
    
    for ( int i = 0, j; i < spR1A->nSpR1FactorElem; ++i ) { // Column index
        v[spR1MatIdx[i] * nCol + spR1MatIdx[i]] = \
        spR1A->spR1FactorSign * spR1A->spR1MatElem[i] * spR1A->spR1MatElem[i];
        for ( j = i + 1; j < spR1A->nSpR1FactorElem; ++j ) { // Row index
            v[spR1MatIdx[i] * nCol + spR1MatIdx[j]] = \
            spR1A->spR1FactorSign * spR1A->spR1MatElem[i] * spR1A->spR1MatElem[j];
            v[spR1MatIdx[j] * nCol + spR1MatIdx[i]] = \
            v[spR1MatIdx[i] * nCol + spR1MatIdx[j]];
        }
    }
    
    return;
}

/** @brief Construct full matrix from rank one dense matrix A
 *  @param[in] A sdpRankOneDenseData pointer
 *  @param[out] v Store full matrix
 */
extern void dataMatDumpRankOneDenseImpl( void *A, double *v ) {
    
    sdp_coeff_dsr1 *dsR1A = (sdp_coeff_dsr1 *) A;
    
    for ( int i = 0, j; i < dsR1A->nSDPCol; ++i ) { // Column index
        v[i * dsR1A->nSDPCol + i] = \
        dsR1A->r1FactorSign * dsR1A->r1MatFactor[i] * dsR1A->r1MatFactor[i];
        for ( j = i + 1; j < dsR1A->nSDPCol; ++j ) { // Row index
            v[i * dsR1A->nSDPCol + j] = \
            dsR1A->r1FactorSign * dsR1A->r1MatFactor[i] * dsR1A->r1MatFactor[j];
            v[j * dsR1A->nSDPCol + i] = v[i * dsR1A->nSDPCol + j];
        }
    }
    
    return;
}
