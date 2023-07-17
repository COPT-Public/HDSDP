/** @file hdsdp\_sdpdata.h
 *  @brief Implenent the SDP Coeffcient data operations
 *  @date 11/24/2022
 */

#ifdef HEADERPATH
#include "interface/hdsdp_utils.h"
#include "interface/def_hdsdp_sdpdata.h"
#include "linalg/hdsdp_sdpdata.h"
#include "linalg/vec_opts.h"
#include "linalg/sparse_opts.h"
#include "linalg/dense_opts.h"
#include "linalg/r1_opts.h"
#include "external/hdsdp_cs.h"
#else
#include "hdsdp_utils.h"
#include "def_hdsdp_sdpdata.h"
#include "hdsdp_sdpdata.h"
#include "vec_opts.h"
#include "sparse_opts.h"
#include "dense_opts.h"
#include "r1_opts.h"
#include "hdsdp_cs.h"
#endif

#include <math.h>

static hdsdp_retcode dataMatCreateZeroImpl( void **pA, int nSDPCol, int dataMatNnz, int *dataMatIdx, double *dataMatElem ) {
    
    hdsdp_retcode retcode = HDSDP_RETCODE_OK;
    
    (void) dataMatNnz;
    (void) dataMatIdx;
    (void) dataMatElem;
    
    if ( !pA ) {
        retcode = HDSDP_RETCODE_FAILED;
        goto exit_cleanup;
    }
    
    sdp_coeff_zero *zero = NULL;
    HDSDP_INIT(zero, sdp_coeff_zero, 1);
    HDSDP_MEMCHECK(zero);
    
    zero->nSDPCol = nSDPCol;
    *pA = (void *) zero;
    
exit_cleanup:
    
    return retcode;
}

static void dataMatViewZeroImpl( void *A );
static void dataMatViewSparseImpl( void *A );
static void dataMatViewDenseImpl( void *A );
static void dataMatViewRankOneSparseImpl( void *A );
static void dataMatViewRankOneDenseImpl( void *A );

static hdsdp_retcode dataMatCreateSparseImpl( void **pA, int nSDPCol, int dataMatNnz, int *dataMatIdx, double *dataMatElem ) {
    
    hdsdp_retcode retcode = HDSDP_RETCODE_OK;
    
    if ( !pA ) {
        retcode = HDSDP_RETCODE_FAILED;
        goto exit_cleanup;
    }
    
    sdp_coeff_sparse *sparse = NULL;
    HDSDP_INIT(sparse, sdp_coeff_sparse, 1);
    HDSDP_MEMCHECK(sparse);
    
    sparse->nSDPCol = nSDPCol;
    sparse->nTriMatElem = dataMatNnz;
    
    HDSDP_INIT(sparse->triMatRow, int, dataMatNnz);
    HDSDP_INIT(sparse->triMatCol, int, dataMatNnz);
    HDSDP_INIT(sparse->triMatElem, double, dataMatNnz);
    
    HDSDP_MEMCHECK(sparse->triMatRow);
    HDSDP_MEMCHECK(sparse->triMatCol);
    HDSDP_MEMCHECK(sparse->triMatElem);
    
    if ( !HUtilCheckIfAscending(dataMatNnz, dataMatIdx) ) {
        HUtilSortDblByInt(dataMatElem, dataMatIdx, 0, dataMatNnz - 1);
    }
    
    tsp_decompress(sparse->nSDPCol, sparse->nTriMatElem, dataMatIdx, dataMatElem,
                   sparse->triMatRow, sparse->triMatCol, sparse->triMatElem);
    
    *pA = (void *) sparse;
    
exit_cleanup:
    
    return retcode;
}

static hdsdp_retcode dataMatCreateDenseImpl( void **pA, int nSDPCol, int dataMatNnz, int *dataMatIdx, double *dataMatElem ) {
    
    hdsdp_retcode retcode = HDSDP_RETCODE_OK;
    
    if ( !pA ) {
        retcode = HDSDP_RETCODE_FAILED;
        goto exit_cleanup;
    }
    
    sdp_coeff_dense *dense = NULL;
    HDSDP_INIT(dense, sdp_coeff_dense, 1);
    HDSDP_MEMCHECK(dense);
    
    dense->nSDPCol = nSDPCol;
    HDSDP_INIT(dense->dsMatElem, double, PACK_NNZ(nSDPCol));
    HDSDP_MEMCHECK(dense->dsMatElem);
    
    pds_decompress(dataMatNnz, dataMatIdx, dataMatElem, dense->dsMatElem);

    *pA = (void *) dense;
    
exit_cleanup:
               
    return retcode;
}

/* Rank-one matrices cannot be directly created */
static hdsdp_retcode dataMatCreateRankOneSparseImpl( void **pA, int nSDPCol, int dataMatNnz, int *dataMatIdx, double *dataMatElem ) {
    
    return HDSDP_RETCODE_FAILED;
}

static hdsdp_retcode dataMatCreateRankOneDenseImpl( void **pA, int nSDPCol, int dataMatNnz, int *dataMatIdx, double *dataMatElem ) {
    
    return HDSDP_RETCODE_FAILED;
}

static void dataMatScalZeroImpl( void *A, double alpha ) {
    
    return;
}

/** @brief A = alpha \* A for sparse matrix A
 *  @param[in] A sdpSparseData pointer
 *  @param[in] alpha Scale parameter
 */
static void dataMatScalSparseImpl( void *A, double alpha ) {
    
    sdp_coeff_sparse *spA = (sdp_coeff_sparse *) A;
    int incx = 1;
    
    scal(&spA->nTriMatElem, &alpha, spA->triMatElem, &incx);
    
    return;
}

/** @brief A = alpha \* A for dense matrix A
 *  @param[in] A sdpDenseData pointer
 *  @param[in] alpha Scale parameter
 */
static void dataMatScalDenseImpl( void *A, double alpha ) {
    
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
static void dataMatScalRankOneSparseImpl( void *A, double alpha ) {
    
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
static void dataMatScalRankOneDenseImpl( void *A, double alpha ) {
    
    /* Avoid alpha = 0.0 case */
    assert( fabs(alpha) >= 1e-08 );
    
    sdp_coeff_dsr1 *dsR1A = (sdp_coeff_dsr1 *) A;
    dsR1A->r1FactorSign = dsR1A->r1FactorSign * alpha;
    
    return;
}

static double dataMatNormZeroImpl( void *A, int type ) {
    
    return 0.0;
}

/** @brief Calculate norm of sparse matrix A
 *  @param[in] A sdpSparseData pointer
 *  @param[in] type Type of norm
 *  @return Norm of A
 */
static double dataMatNormSparseImpl( void *A, int type ) {
    
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
static double dataMatNormDenseImpl( void *A, int type ) {
    
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
static double dataMatNormRankOneSparseImpl( void *A, int type ) {
    
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
static double dataMatNormRankOneDenseImpl( void *A, int type ) {
    
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

static hdsdp_retcode dataMatBuildUpEigZeroImpl( void *A, int *rank, double *auxiFullMat, double **eVals, double **eVecs ) {
    
    *rank = 0;
    
    return HDSDP_RETCODE_OK;
}

static hdsdp_retcode dataMatBuildUpEigSparseImpl( void *A, int *rank, double *auxiFullMat, double **eVals, double **eVecs ) {
    
    hdsdp_retcode retcode = HDSDP_RETCODE_OK;
    
    if ( !eVals || !eVecs ) {
        retcode = HDSDP_RETCODE_FAILED;
        goto exit_cleanup;
    }
    
    sdp_coeff_sparse *sparse = (sdp_coeff_sparse *) A;
    
    double sgn = 0.0;
    HDSDP_ZERO(auxiFullMat, double, sparse->nSDPCol);
    int isRankOne = tsp_r1_extract(sparse->nSDPCol, sparse->nTriMatElem, sparse->triMatRow,
                                   sparse->triMatCol, sparse->triMatElem, &sgn, auxiFullMat);
    
    if ( !isRankOne ) {
        goto exit_cleanup;
    }
    
    *rank = 1;
    HDSDP_INIT(*eVals, double, 1);
    HDSDP_INIT(*eVecs, double, sparse->nSDPCol);
    HDSDP_MEMCHECK(*eVecs);
    (*eVals)[0] = sgn;
    HDSDP_MEMCPY(*eVecs, auxiFullMat, double, sparse->nSDPCol);
    
exit_cleanup:
    
    return retcode;
}

static hdsdp_retcode dataMatBuildUpEigDenseImpl( void *A, int *rank, double *auxiFullMat, double **eVals, double **eVecs ) {
    
    hdsdp_retcode retcode = HDSDP_RETCODE_OK;
    
    if ( !eVals || !eVecs ) {
        retcode = HDSDP_RETCODE_FAILED;
        goto exit_cleanup;
    }
    
    sdp_coeff_dense *dense = (sdp_coeff_dense *) A;
    
    double sgn = 0.0;
    int isRankOne = pds_r1_extract(dense->nSDPCol, dense->dsMatElem, &sgn, auxiFullMat);
    
    if ( !isRankOne ) {
        goto exit_cleanup;
    }
    
    *rank = 1;
    HDSDP_INIT(*eVals, double, 1);
    HDSDP_INIT(*eVecs, double, dense->nSDPCol);
    HDSDP_MEMCHECK(*eVecs);
    (*eVals)[0] = sgn;
    HDSDP_MEMCPY(*eVecs, auxiFullMat, double, dense->nSDPCol);
    
exit_cleanup:
    
    return retcode;
}

static hdsdp_retcode dataMatBuildUpRankOneSparseImpl( void *A, int *dummy1, double *dummy2, double **dummy3, double **dummy4 ) {
    
    return HDSDP_RETCODE_FAILED;
}

static hdsdp_retcode dataMatBuildUpRankOneDenseImpl( void *A, int *dummy1, double *dummy2, double **dummy3, double **dummy4 ) {
    
    return HDSDP_RETCODE_FAILED;
}

static int dataMatGetNnzZeroImpl( void *A ) {
    
    return 0;
}

/** @brief Calculate number of nonzero elements in sparse matrix A
 *  @param[in] A sdpSparseData pointer
 *  @return Number of nonzero elements in A
 */
static int dataMatGetNnzSparseImpl( void *A ) {
    
    sdp_coeff_sparse *spA = (sdp_coeff_sparse *) A;
    
    return spA->nTriMatElem;
}

/** @brief Calculate number of nonzero elements in dense matrix A
 *  @param[in] A sdpDenseData pointer
 *  @return Number of nonzero elements in A
 */
static int dataMatGetNnzDenseImpl( void *A ) {
    
    sdp_coeff_dense *dsA = (sdp_coeff_dense *) A;
    
    return PACK_NNZ(dsA->nSDPCol);
}

/** @brief Calculate number of nonzero elements in rank one sparse matrix A
 *  @param[in] A sdpRankOneSparseData pointer
 *  @return Number of nonzero elements in A
 */
static int dataMatGetNnzRankOneSparseImpl( void *A ) {
    
    sdp_coeff_spr1 *spR1A = (sdp_coeff_spr1 *) A;
    
    return PACK_NNZ(spR1A->nSpR1FactorElem);
}

/** @brief Calculate number of nonzero elements in rank one dense matrix A
 *  @param[in] A sdpRankOneDenseData pointer
 *  @return Number of nonzero elements in A
 */
static int dataMatGetNnzRankOneDenseImpl( void *A ) {
    
    sdp_coeff_dsr1 *dsR1A = (sdp_coeff_dsr1 *) A;
    
    return PACK_NNZ(dsR1A->nSDPCol);
}

static void dataMatDumpZeroImpl( void *A, double *v ) {
    
    sdp_coeff_zero *zeroA = (sdp_coeff_zero *) A;
    HDSDP_ZERO(v, double, zeroA->nSDPCol * zeroA->nSDPCol);
    
    return;
}

/** @brief Construct full matrix from sparse matrix A
 *  @param[in] A sdpSparseData pointer
 *  @param[out] v Store full matrix
 */
static void dataMatDumpSparseImpl( void *A, double *v ) {
    
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
static void dataMatDumpDenseImpl( void *A, double *v ) {
    
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
static void dataMatDumpRankOneSparseImpl( void *A, double *v ) {
    
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
static void dataMatDumpRankOneDenseImpl( void *A, double *v ) {
    
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

static void dataMatGetZeroSparsityImpl( void *A, int *spout ) {
    
    return;
}

static void dataMatGetSparseSparsityImpl( void *A, int *spout ) {
    
    sdp_coeff_sparse *spA = (sdp_coeff_sparse *) A;
    
    for ( int k = 0; k < spA->nTriMatElem; ++k ) {
        spout[PACK_IDX(spA->nSDPCol, spA->triMatRow[k], spA->triMatCol[k])] = 1;
    }

    return;
}

static void dataMatGetDenseSparsityImpl( void *A, int *spout ) {
    
    sdp_coeff_dense *dsA = (sdp_coeff_dense *) A;
    
    for ( int i = 0; i < PACK_NNZ(dsA->nSDPCol); ++i ) {
        spout[i] = 1;
    }
    
    /* This method should not be invoked */
    assert( 0 );
    return;
}

static void dataMatGetRankOneSparseSparsityImpl( void *A, int *spout ) {
    
    sdp_coeff_spr1 *spR1A = (sdp_coeff_spr1 *) A;
    
    for ( int i = 0, j; i < spR1A->nSpR1FactorElem; ++i ) {
        for ( j = 0; j < i + 1; ++j ) {
            spout[PACK_IDX(spR1A->nSDPCol, spR1A->spR1MatIdx[i], spR1A->spR1MatIdx[j])] = 1;
        }
    }
    
    return;
}

static void dataMatGetRankOneDenseSparsityImpl( void *A, int *spout ) {
    
    sdp_coeff_dsr1 *dsR1A = (sdp_coeff_dsr1 *) A;
    
    for ( int i = 0; i < PACK_NNZ(dsR1A->nSDPCol); ++i ) {
        spout[i] = 1;
    }
    
    /* This method should not be invoked */
    assert( 0 );
    return;
}

static void dataMatAddZeroToBufferImpl( void *A, double a, int *spmap, double *buffer ) {
    /* Let buffer <- buffer + a * A */
    
    return;
}

static void dataMatAddSparseToBufferImpl( void *A, double a, int *spmat, double *buffer ) {
    /* If spmat exists, then the buffer uses sparse data structure,
       Otherwise (if spmat is NULL), buffer is a dense array
     */
    
    sdp_coeff_sparse *spA = (sdp_coeff_sparse *) A;
    
    if ( !spmat ) {
        /* Simply add it */
        for ( int i = 0; i < spA->nTriMatElem; ++i ) {
            FULL_ENTRY(buffer, spA->nSDPCol, spA->triMatRow[i], spA->triMatCol[i]) += \
            a * spA->triMatElem[i];
        }
    } else {
        /* Add according to the mapping */
        for ( int i = 0; i < spA->nTriMatElem; ++i ) {
            buffer[spmat[PACK_IDX(spA->nSDPCol, spA->triMatRow[i], spA->triMatCol[i])]] += \
            a * spA->triMatElem[i];
        }
    }
    
    return;
}

static void dataMatAddDenseToBufferImpl( void *A, double a, int *spmat, double *buffer ) {
    /* Whenever a dense matrix exists, the dual matrix must be dense */
    sdp_coeff_dense *dsA = (sdp_coeff_dense *) A;
    
    if ( !spmat ) {
        /* We are adding a dense packed matrix to a dense full matrix */
        double *dPackElem = dsA->dsMatElem;
        double *dFullElem = buffer;
        for ( int i = 0; i < dsA->nSDPCol; ++i ) {
            for ( int j = 0; j < dsA->nSDPCol - i; ++j ) {
                dFullElem[j] += a * dPackElem[j];
            }
            dFullElem += dsA->nSDPCol + 1;
            dPackElem += (dsA->nSDPCol - i);
        }
        
    } else {
        /* This case should never happen */
        assert( 0 );
    }
    
    return;
}

static void dataMatAddRankOneSparseToBufferImpl( void *A, double a, int *spmat, double *buffer ) {
    
    sdp_coeff_spr1 *spR1A = (sdp_coeff_spr1 *) A;
    
    double aAsign = a * spR1A->spR1FactorSign;
    
    if ( !spmat ) {
        for ( int i = 0; i < spR1A->nSpR1FactorElem; ++i ) {
            for ( int j = 0; j < i + 1; ++j ) {
                FULL_ENTRY(buffer, spR1A->nSDPCol, spR1A->spR1MatIdx[i], spR1A->spR1MatIdx[j]) += \
                aAsign * spR1A->spR1MatElem[i] * spR1A->spR1MatElem[j];
            }
        }
    } else {
        for ( int i = 0; i < spR1A->nSpR1FactorElem; ++i ) {
            for ( int j = 0; j < i + 1; ++j ) {
                buffer[spmat[PACK_IDX(spR1A->nSDPCol, spR1A->spR1MatIdx[i], spR1A->spR1MatIdx[j])]] += \
                aAsign * spR1A->spR1MatElem[i] * spR1A->spR1MatElem[j];
            }
        }
    }
    
    return;
}

static void dataMatAddRankOneDenseToBufferImpl( void *A, double a, int *spmat, double *buffer ) {
    
    sdp_coeff_dsr1 *dsR1A = (sdp_coeff_dsr1 *) A;
    
    double aAsign = a * dsR1A->r1FactorSign;
    
    if ( !spmat ) {
        syr(&HCharConstantUploLow, &dsR1A->nSDPCol, &aAsign,
            dsR1A->r1MatFactor, &HIntConstantOne, buffer, &dsR1A->nSDPCol);
    } else {
        /* This case should never happen */
        assert( 0 );
    }
    
    return;
}

static void dataMatClearZeroImpl( void *A ) {
    
    if ( !A ) {
        return;
    }
    
    HDSDP_ZERO(A, sdp_coeff_zero, 1);
    
    return;
}

static void dataMatClearSparseImpl( void *A ) {
    
    if ( !A ) {
        return;
    }
    
    sdp_coeff_sparse *sparse = (sdp_coeff_sparse *) A;
    
    HDSDP_FREE(sparse->triMatRow);
    HDSDP_FREE(sparse->triMatCol);
    HDSDP_FREE(sparse->triMatElem);
    HDSDP_ZERO(A, sdp_coeff_sparse, 1);
    
    return;
}

static void dataMatClearDenseImpl( void *A ) {
    
    if ( !A ) {
        return;
    }
    
    sdp_coeff_dense *dense = (sdp_coeff_dense *) A;
    
    HDSDP_FREE(dense->dsMatElem);
    HDSDP_ZERO(A, sdp_coeff_dense, 1);
    
    return;
}

static void dataMatClearRankOneSparseImpl( void *A ) {
    
    if ( !A ) {
        return;
    }
    
    sdp_coeff_spr1 *spr1 = (sdp_coeff_spr1 *) A;
    
    HDSDP_FREE(spr1->spR1MatElem);
    HDSDP_FREE(spr1->spR1MatIdx);
    
    HDSDP_ZERO(A, sdp_coeff_spr1, 1);
    
    return;
}

static void dataMatClearRankOneDenseImpl( void *A ) {
    
    if ( !A ) {
        return;
    }
    
    sdp_coeff_dsr1 *dsr1 = (sdp_coeff_dsr1 *) A;
    
    HDSDP_FREE(dsr1->r1MatFactor);
    HDSDP_ZERO(A, sdp_coeff_dsr1, 1);
    
    return;
}

static void dataMatDestroyZeroImpl( void **pA ) {
    
    if ( !pA ) {
        return;
    }
    
    dataMatClearZeroImpl(*pA);
    HDSDP_FREE(*pA);
    
    return;
}

static void dataMatDestroySparseImpl( void **pA ) {
    
    if ( !pA ) {
        return;
    }
    
    dataMatClearSparseImpl(*pA);
    HDSDP_FREE(*pA);
    
    return;
}

static void dataMatDestroyDenseImpl( void **pA ) {
    
    if ( !pA ) {
        return;
    }
    
    dataMatClearDenseImpl(*pA);
    HDSDP_FREE(*pA);
    
    return;
}

static void dataMatDestroyRankOneSparseImpl( void **pA ) {
    
    if ( !pA ) {
        return;
    }
    
    dataMatClearRankOneSparseImpl(*pA);
    HDSDP_FREE(*pA);
    
    return;
}

static void dataMatDestroyRankOneDenseImpl( void **pA ) {
    
    if ( !pA ) {
        return;
    }
    
    dataMatClearRankOneDenseImpl(*pA);
    HDSDP_FREE(*pA);
    
    return;
}

static void dataMatViewZeroImpl( void *A ) {
    
    printf("Zero matrix of size %d L1 = [%5.3e] L2 = [%5.3e] \n",
           ((sdp_coeff_zero *) A)->nSDPCol, 0.0, 0.0);
    
    return;
}

static void dataMatViewSparseImpl( void *A ) {
    
    sdp_coeff_sparse *sparse = (sdp_coeff_sparse *) A;
    
    printf("Sparse matrix of size %d and %d nnzs L1 = [%5.3e] L2 = [%5.3e]. \n",
           sparse->nSDPCol, sparse->nTriMatElem,
           dataMatNormSparseImpl(sparse, ABS_NORM),
           dataMatNormSparseImpl(A, FRO_NORM));
    
#if 0
    dcs dcsMat;
    dcsMat.n = sparse->nSDPCol;
    dcsMat.m = sparse->nSDPCol;
    dcsMat.nz = sparse->nTriMatElem;
    dcsMat.nzmax = sparse->nTriMatElem;
    dcsMat.p = sparse->triMatCol;
    dcsMat.i = sparse->triMatRow;
    dcsMat.x = sparse->triMatElem;
    dcs_print(&dcsMat, 1);
#endif
    
    return;
}

static void dataMatViewDenseImpl( void *A ) {
    
    sdp_coeff_dense *dense = (sdp_coeff_dense *) A;
    printf("Dense matrix of size %d L1 = [%5.3e] L2 = [%5.3e]. \n",
           dense->nSDPCol, dataMatNormDenseImpl(dense, ABS_NORM),
           dataMatNormDenseImpl(dense, FRO_NORM));
    
    return;
}

static void dataMatViewRankOneSparseImpl( void *A ) {
    
    sdp_coeff_spr1 *spr1 = (sdp_coeff_spr1 *) A;
    printf("Sparse rank-one matrix of size %d L1 = [%5.3e] L2 = [%5.3e]. \n", spr1->nSDPCol,
           dataMatNormRankOneSparseImpl(spr1, ABS_NORM),
           dataMatNormRankOneSparseImpl(spr1, FRO_NORM));
    
    return;
}

static void dataMatViewRankOneDenseImpl( void *A ) {
    
    sdp_coeff_dsr1 *dsr1 = (sdp_coeff_dsr1 *) A;
    printf("Dense rank-one matrix of size %d L1 = [%5.3e] L2 = [%5.3e]. \n", dsr1->nSDPCol,
           dataMatNormRankOneDenseImpl(dsr1, ABS_NORM),
           dataMatNormRankOneDenseImpl(dsr1, FRO_NORM));
    
    return;
}

/* Scale a from A = sgn * a * a' and put numerical difficulties in sgn */
static void dataMatNormalizeRankOneSparseImpl( sdp_coeff_spr1 *A ) {
    
    int incx = 1;
    double aScal = nrm2(&A->nSpR1FactorElem, A->spR1MatElem, &incx);
    A->spR1FactorSign *= (aScal * aScal);
    rscl(&A->nSpR1FactorElem, &aScal, A->spR1MatElem, &incx);
    
    return;
}

static void dataMatNormalizeRankOneDenseImpl( sdp_coeff_dsr1 *A ) {
    
    int incx = 1;
    double aScal = nrm2(&A->nSDPCol, A->r1MatFactor, &incx);
    A->r1FactorSign *= (aScal * aScal);
    rscl(&A->nSDPCol, &aScal, A->r1MatFactor, &incx);
    
    return;
}
/*==========================================================================================*/
/* KKT operation: compute v = S^-1 a if A = sign * a * a' */
/*==========================================================================================*/
static void dataMatZeroKKT2SolveRankOneImpl( void *A, hdsdp_linsys *S, double *Sinv, double *sign, double *v ) {
    
    assert( 0 );
    return;
}

static void dataMatSparseKKT2SolveRankOneImpl( void *A, hdsdp_linsys *S, double *Sinv, double *sign, double *v ) {
    
    assert( 0 );
    return;
}

static void dataMatDenseKKT2SolveRankOneImpl( void *A, hdsdp_linsys *S, double *Sinv, double *sign, double *v ) {
    
    assert( 0 );
    return;
}

static void dataMatRankOneSparseKKT2SolveRankOneImpl( void *A, hdsdp_linsys *S, double *Sinv, double *sign, double *v ) {
    /* Implement v = S^-1 a.
       In HDSDP we have both the factorization of S and S^-1 in explicit form
       
       If a itself is sparse, then S^-1 a can be done by taking a sparse linear combination of columns of S^-1
       If a is dense, then we do direct solve S \ v
     
     */
    sdp_coeff_spr1 *spr1 = (sdp_coeff_spr1 *) A;
    
    if ( spr1->nSpR1FactorElem >= 0.3 * spr1->nSDPCol ) {
        /* If a is sparse but not that sparse */
        HFpLinsysSolve(S, 1, spr1->spR1MatFactor, v);
    } else {
        /* If a contains few nonzeros, we take linear combination of Sinv */
        HDSDP_ZERO(v, double, spr1->nSDPCol);
        for ( int i = 0; i < spr1->nSpR1FactorElem; ++i ) {
            double ai = spr1->spR1MatElem[i];
            axpy(&spr1->nSDPCol, &ai, Sinv + spr1->spR1MatIdx[i] * spr1->nSDPCol,
                 &HIntConstantOne, v, &HIntConstantOne);
        }
    }
    
    *sign = spr1->spR1FactorSign;
    
    return;
}

static void dataMatRankOneDenseKKT2SolveRankOneImpl( void *A, hdsdp_linsys *S, double *Sinv, double *sign, double *v ) {
    
    sdp_coeff_dsr1 *dsr1 = (sdp_coeff_dsr1 *) A;
    
    if ( S->LinType == HDSDP_LINSYS_SPARSE_DIRECT ) {
        HFpLinsysSolve(S, 1, dsr1->r1MatFactor, v);
    } else {
        fds_symv(dsr1->nSDPCol, 1.0, Sinv, dsr1->r1MatFactor, 0.0, v);
    }
    
    *sign = dsr1->r1FactorSign;
    
    return;
}
/*==========================================================================================*/
/* KKT operation: compute trace: trace(A * S^-1) = sign * a' * v */
/*==========================================================================================*/
static double dataMatZeroKKT2TraceASinvImpl( void *A, double *v ) {
    
    assert( 0 );
    return 0.0;
}

static double dataMatSparseKKT2TraceASinvImpl( void *A, double *v ) {
    
    assert( 0 );
    return 0.0;
}

static double dataMatDenseKKT2TraceASinvImpl( void *A, double *v ) {
    
    assert( 0 );
    return 0.0;
}

static double dataMatRankOneSparseKKT2TraceASinvImpl( void *A, double *v ) {
    
    sdp_coeff_spr1 *spr1 = (sdp_coeff_spr1 *) A;
    double asinv = 0.0;
    
    for ( int i = 0; i < spr1->nSpR1FactorElem; ++i ) {
        asinv += v[spr1->spR1MatIdx[i]] * spr1->spR1MatElem[i];
    }
    
    return spr1->spR1FactorSign * asinv;
}

static double dataMatRankOneDenseKKT2TraceASinvImpl( void *A, double *v ) {
    
    sdp_coeff_dsr1 *dsr1 = (sdp_coeff_dsr1 *) A;
    return dsr1->r1FactorSign * dot(&dsr1->nSDPCol, dsr1->r1MatFactor, &HIntConstantOne, v, &HIntConstantOne);
}

/* KKT operation: compute quadratic form v' * A * v */
static double dataMatZeroKKT2ComputeQuadForm( void *A, double *v, double *aux ) {
    
    return 0.0;
}

static double dataMatSparseKKT2ComputeQuadForm( void *A, double *v, double *aux ) {
    
    sdp_coeff_sparse *sparse = (sdp_coeff_sparse *) A;
    
    return tsp_quadform(sparse->nSDPCol, sparse->nTriMatElem, sparse->triMatRow,
                        sparse->triMatCol, sparse->triMatElem, v);
}

static double dataMatDenseKKT2ComputeQuadForm( void *A, double *v, double *aux ) {
    
    sdp_coeff_dense *dense = (sdp_coeff_dense *) A;
    
    return pds_quadform(dense->nSDPCol, dense->dsMatElem, v, aux);
}

static double dataMatRankOneSparseKKT2ComputeQuadForm( void *A, double *v, double *aux ) {
    
    sdp_coeff_spr1 *spr1 = (sdp_coeff_spr1 *) A;
    
    return spr1_quadform(spr1->nSDPCol, spr1->spR1FactorSign,
                         spr1->nSpR1FactorElem, spr1->spR1MatIdx, spr1->spR1MatElem, v);
}

static double dataMatRankOneDenseKKT2ComputeQuadForm( void *A, double *v, double *aux ) {
    
    sdp_coeff_dsr1 *dsr1 = (sdp_coeff_dsr1 *) A;
    
    return dsr1_quadform(dsr1->nSDPCol, dsr1->r1FactorSign, dsr1->r1MatFactor, v);
}
/*==========================================================================================*/
/* KKT operation: Compute S^-1 A * S^-1 and simultaneously compute trace(A * S^-1) */
/*==========================================================================================*/
static double dataMatZeroKKT3ComputeSinvASinvImpl( void *A, hdsdp_linsys *S, double *Sinv, double *aux, double *B ) {
    
    return 0.0;
}

static double dataMatSparseKKT3ComputeSinvASinvImpl( void *A, hdsdp_linsys *S, double *Sinv, double *aux, double *B ) {
    
    sdp_coeff_sparse *sparse = (sdp_coeff_sparse *) A;
    
    int iRow, iCol;
    
    double aVal = 0.0;
    double *SinvACol = NULL;
    double *SinvARow = NULL;
    double *SinvRow = NULL;
    double *SinvCol = NULL;
    double *SinvA = aux;
    double dTraceASinv = 0.0;
    
    HDSDP_ZERO(SinvA, double, sparse->nSDPCol * sparse->nSDPCol);
    
    /* First we set up S^-1 * A by combining columns of S^-1 */
    for ( int iElem = 0; iElem < sparse->nTriMatElem; ++iElem ) {
        
        iRow = sparse->triMatRow[iElem];
        iCol = sparse->triMatCol[iElem];
        aVal = sparse->triMatElem[iElem];
        SinvACol = SinvA + sparse->nSDPCol * iCol;
        SinvRow = Sinv + sparse->nSDPCol * iRow;
        axpy(&sparse->nSDPCol, &aVal, SinvRow, &HIntConstantOne, SinvACol, &HIntConstantOne);
        
        /* Map to the symmetric component of the sparse matrix */
        if ( iRow != iCol ) {
            SinvACol = SinvA + sparse->nSDPCol * iRow;
            SinvRow = Sinv + sparse->nSDPCol * iCol;
            axpy(&sparse->nSDPCol, &aVal, SinvRow, &HIntConstantOne, SinvACol, &HIntConstantOne);
        }
    }
    
    /* Then we do matrix-matrix multiplication B = SinvA * S^-1 = S^-1 * SinvA'
       Note that we know that the result of computation will be symmetric. Using level-2 or level-3 blas
       will make the method not symmetry aware. Thus we only resort to level-1 blas.
       Also the stride of ddot will be > 0 and we think it unavoidable now.
     
       When computing the final result, we simultaneously compute trace(A * S^-1),
       which is the sum of diagonal of SinvA
     */
    
    for ( iRow = 0; iRow < sparse->nSDPCol; ++iRow ) {
        SinvARow = SinvA + iRow;
        /* Reach the diagonal */
        dTraceASinv += SinvARow[iRow * sparse->nSDPCol];
        for ( iCol = 0; iCol <= iRow; ++iCol ) {
            SinvCol = Sinv + sparse->nSDPCol * iCol;
            FULL_ENTRY(B, sparse->nSDPCol, iRow, iCol) = \
            dot(&sparse->nSDPCol, SinvARow, &sparse->nSDPCol, SinvCol, &HIntConstantOne);
        }
    }
    
    return dTraceASinv;
}

static double dataMatDenseKKT3ComputeSinvASinvImpl( void *A, hdsdp_linsys *S, double *Sinv, double *aux, double *B ) {
    
    sdp_coeff_dense *dense = (sdp_coeff_dense *) A;
    
    /* Different from the implementation for the sparse matrix, in the dense matrix case we first compute
       A * S^-1 so that level-2 blas can be employed */
    
    double *SinvCol = NULL;
    double *SinvRow = NULL;
    double *ASinvCol = NULL;
    double *ASinvRow = NULL;
    double *ASinv = aux;
    double dTraceASinv = 0.0;
    
    HDSDP_ZERO(ASinv, double, dense->nSDPCol * dense->nSDPCol);
    
    /* First compute A * S^-1 column by column */
    for ( int iCol = 0; iCol < dense->nSDPCol; ++iCol ) {
        SinvCol = Sinv + dense->nSDPCol * iCol;
        ASinvCol = ASinv + dense->nSDPCol * iCol;
        pds_spmv(UPLOLOW, dense->nSDPCol, 1.0, dense->dsMatElem, SinvCol, 1, 0.0, ASinvCol, 1);
    }
    
    /* Then we are ready to compute S^-1 * ASinv.
       This part of the code involves two conversions.
       
       The loop order acturally targets the upper-triangular of the matrix
       and we map the element symmetrically
     */
    for ( int iCol = 0; iCol < dense->nSDPCol; ++iCol ) {
        ASinvCol = ASinv + dense->nSDPCol * iCol;
        dTraceASinv += ASinvCol[iCol];
        for ( int iRow = 0; iRow <= iCol; ++iRow ) {
            /* We note that S^-1 is symmetric. So the i-th row is exactly the i-th column
               and the stride becomes 1 */
            SinvCol = Sinv + dense->nSDPCol * iRow;
            FULL_ENTRY(B, dense->nSDPCol, iCol, iRow) = \
            dot(&dense->nSDPCol, ASinvCol, &HIntConstantOne, SinvCol, &HIntConstantOne);
        }
    }
    
    return dTraceASinv;
}

static double dataMatRankOneSparseKKT3ComputeSinvASinvImpl( void *A, hdsdp_linsys *S, double *Sinv, double *aux, double *B ) {
    
    sdp_coeff_spr1 *spr1 = (sdp_coeff_spr1 *) A;
    /*
       Generally this part of code will not be invoked.
       For rank-one matrices, we have S^-1 * A * S^-1 = sign * v * v'
       and trace(A * S^-1) = sign * a' * v
     
       To compute v we resort to the KKT M2 routine
     */
    
    double *dSinvAVecBuffer = aux;
    double dFactorSign = 0.0;
    double dTraceASinv = 0.0;
    
    HDSDP_ZERO(B, double, spr1->nSDPCol * spr1->nSDPCol);
    
    /* First we get v and it's stored in aux */
    dataMatRankOneSparseKKT2SolveRankOneImpl(A, S, Sinv, &dFactorSign, dSinvAVecBuffer);
    /* Then we compute sign * a' * v */
    dTraceASinv = dataMatRankOneSparseKKT2TraceASinvImpl(A, dSinvAVecBuffer);
    /* Finally we do rank-one update */
    pds_syr(UPLOLOW, spr1->nSDPCol, dFactorSign, dSinvAVecBuffer, 1, B, spr1->nSDPCol);
    
    return dTraceASinv;
}

static double dataMatRankOneDenseKKT3ComputeSinvASinvImpl( void *A, hdsdp_linsys *S, double *Sinv, double *aux, double *B ) {
    
    sdp_coeff_dsr1 *dsr1 = (sdp_coeff_dsr1 *) A;
    
    /* The same routine as the sparse case */
    double *dSinvAVecBuffer = aux;
    double dFactorSign = 0.0;
    double dTraceASinv = 0.0;
    
    HDSDP_ZERO(B, double, dsr1->nSDPCol * dsr1->nSDPCol);
    /* First we get v and it's stored in aux */
    dataMatRankOneDenseKKT2SolveRankOneImpl(A, S, Sinv, &dFactorSign, dSinvAVecBuffer);
    /* Then we compute sign * a' * v */
    dTraceASinv = dataMatRankOneDenseKKT2TraceASinvImpl(A, dSinvAVecBuffer);
    /* Finally we do rank-one update */
    pds_syr(UPLOLOW, dsr1->nSDPCol, dFactorSign, dSinvAVecBuffer, 1, B, dsr1->nSDPCol);
    
    return dTraceASinv;
}

/*==========================================================================================*/
/* KKT operation: Compute matrix buffer dot product */
/*==========================================================================================*/
static double dataMatZeroDotDenseKKT3Impl( void *A, double *B, double *aux ) {
    
    return 0.0;
}

static double dataMatSparseDotDenseKKT3Impl( void *A, double *B, double *aux ) {
    
    sdp_coeff_sparse *sparse = (sdp_coeff_sparse *) A;
    
    double dAdotB = 0.0;
    
    int iRow = 0;
    int iCol = 0;
    
    for ( int iElem = 0; iElem < sparse->nTriMatElem; ++iElem ) {
        
        iRow = sparse->triMatRow[iElem];
        iCol = sparse->triMatCol[iElem];
        
        if ( iRow == iCol ) {
            dAdotB += 0.5 * sparse->triMatElem[iElem] * FULL_ENTRY(B, sparse->nSDPCol, iRow, iCol);
        } else {
            dAdotB += sparse->triMatElem[iElem] * FULL_ENTRY(B, sparse->nSDPCol, iRow, iCol);
        }
    }
    
    return 2.0 * dAdotB;
}

static double dataMatDenseDotDenseKKT3Impl( void *A, double *B, double *aux ) {
    
    sdp_coeff_dense *dense = (sdp_coeff_dense *) A;
    /* Compute the dot product between a full buffer and a packed dense matrix */
    double dAdotB = 0.0;
    
    int iRow = 0;
    int iCol = 0;
    
    double *dACol = dense->dsMatElem;
    double *dBCol = B;
    
    for ( iCol = 0; iCol < dense->nSDPCol; ++iCol ) {
        /* The first element of each column is diagonal */
        dAdotB += 0.5 * dACol[0] * dBCol[0];
        for ( iRow = 1; iRow < dense->nSDPCol - iCol; ++iRow ) {
            dAdotB += dACol[iRow] * dBCol[iRow];
        }
        /* Move pointer to the start of each column */
        dACol += dense->nSDPCol - iCol;
        dBCol += dense->nSDPCol + 1;
    }
    
    return dAdotB;
}

static double dataMatRankOneSparseDotDenseKKT3Impl( void *A, double *B, double *aux ) {
    
    sdp_coeff_spr1 *spr1 = (sdp_coeff_spr1 *) A;
    
    double dAdotB = 0.0;
    
    int iRowElem = 0;
    int iColElem = 0;
    
    for ( iColElem = 0; iColElem < spr1->nSpR1FactorElem; ++iColElem ) {
        dAdotB += 0.5 * spr1->spR1MatElem[iColElem] * spr1->spR1MatElem[iColElem] * \
                    FULL_ENTRY(B, spr1->nSDPCol, iColElem, iColElem);
        for ( iRowElem = iColElem + 1; iRowElem < spr1->nSpR1FactorElem; ++iRowElem ) {
            dAdotB += spr1->spR1MatElem[iRowElem] * spr1->spR1MatElem[iColElem] * \
            FULL_ENTRY(B, spr1->nSDPCol, iRowElem, iColElem);
        }
    }
    
    return 2.0 * spr1->spR1FactorSign * dAdotB;
}

static double dataMatRankOneDenseDotDenseKKT3Impl( void *A, double *B, double *aux ) {
    
    sdp_coeff_dsr1 *dsr1 = (sdp_coeff_dsr1 *) A;
    fds_symv(dsr1->nSDPCol, 1.0, B, dsr1->r1MatFactor, 0.0, aux);
    
    return dsr1->r1FactorSign * dot(&dsr1->nSDPCol, &dsr1->r1FactorSign, &HIntConstantOne, aux, &HIntConstantOne);
}

/*==========================================================================================*/
/* KKT operation: Compute A * S^-1 and simultaneously give trace(A * S^-1)  */
/*==========================================================================================*/
static double dataMatZeroKKT4ComputeASinvImpl( void *A, hdsdp_linsys *S, double *Sinv, double *aux, double *B ) {
    
    assert( 0 );
    return 0.0;
}

static double dataMatSparseKKT4ComputeASinvImpl( void *A, double *Sinv, double *aux, double *B ) {
    
    sdp_coeff_sparse *sparse = (sdp_coeff_sparse *) A;
    /* Compute A * S^-1 for sparse A. We employ the same strategy as we did when setting up
       S^-1 * A * S^-1 in M3 strategy
     */
    
    int iRow, iCol;
    
    double aVal = 0.0;
    double *SinvACol = NULL;
    double *SinvRow = NULL;
    double *SinvA = B;
    double dTraceASinv = 0.0;
    
    HDSDP_ZERO(B, double, sparse->nSDPCol * sparse->nSDPCol);
    
    /* First we set up S^-1 * A by combining columns of S^-1 */
    for ( int iElem = 0; iElem < sparse->nTriMatElem; ++iElem ) {
        
        iRow = sparse->triMatRow[iElem];
        iCol = sparse->triMatCol[iElem];
        aVal = sparse->triMatElem[iElem];
        SinvACol = SinvA + sparse->nSDPCol * iCol;
        SinvRow = Sinv + sparse->nSDPCol * iRow;
        axpy(&sparse->nSDPCol, &aVal, SinvRow, &HIntConstantOne, SinvACol, &HIntConstantOne);
        
        /* Map to the symmetric component of the sparse matrix */
        if ( iRow != iCol ) {
            SinvACol = SinvA + sparse->nSDPCol * iRow;
            SinvRow = Sinv + sparse->nSDPCol * iCol;
            axpy(&sparse->nSDPCol, &aVal, SinvRow, &HIntConstantOne, SinvACol, &HIntConstantOne);
        }
    }
    
    /* Compute trace(B) */
    for ( iCol = 0; iCol < sparse->nSDPCol; ++iCol ) {
        dTraceASinv += B[iCol * sparse->nSDPCol + iCol];
    }
    
    return dTraceASinv;
}

static double dataMatDenseKKT4ComputeASinvImpl( void *A, hdsdp_linsys *S, double *Sinv, double *aux, double *B ) {
    
    sdp_coeff_dense *dense = (sdp_coeff_dense *) A;
    
    double *SinvCol = NULL;
    double *ASinvCol = NULL;
    double dTraceASinv = 0.0;
    
    HDSDP_ZERO(B, double, dense->nSDPCol * dense->nSDPCol);
    
    /* First compute A * S^-1 column by column */
    for ( int iCol = 0; iCol < dense->nSDPCol; ++iCol ) {
        SinvCol = Sinv + dense->nSDPCol * iCol;
        ASinvCol = B + dense->nSDPCol * iCol;
        pds_spmv(UPLOLOW, dense->nSDPCol, 1.0, dense->dsMatElem, SinvCol, 1, 0.0, ASinvCol, 1);
    }
    
    /* Compute trace(B) */
    for ( int iCol = 0; iCol < dense->nSDPCol; ++iCol ) {
        dTraceASinv += B[iCol * dense->nSDPCol + iCol];
    }
    
    return dTraceASinv;
}

static double dataMatRankOnesSparseKKT4ComputeASinvImpl( void *A, hdsdp_linsys *S, double *Sinv, double *aux, double *B ) {
    
    sdp_coeff_spr1 *spr1 = (sdp_coeff_spr1 *) A;
    
    /* Compute A * S^-1 = sign * a * a' * S^-1 = sign * a * (S^-1 * a)' */
    HDSDP_ZERO(B, double, spr1->nSDPCol * spr1->nSDPCol);
    
    /* Get buffers ready */
    double *dSinvAVec = aux;
    double dFactorSign = 0.0;
    double dTraceASinv = 0.0;
    double *SinvCol = NULL;
    double *ASinvRow = NULL;
    
    /* We employ two strategies to set up A * S^-1 */
    if ( spr1->nSpR1FactorElem >= sqrt(spr1->nSDPCol) ) {
        /* If A is rather dense, first we set up S^-1 a and get sign */
        dataMatRankOneSparseKKT2SolveRankOneImpl(A, S, Sinv, &dFactorSign, dSinvAVec);
        /* Then we do rank-one update */
        fds_ger(spr1->nSDPCol, spr1->nSDPCol, dFactorSign,
                spr1->spR1MatFactor, 1, dSinvAVec, 1, B, spr1->nSDPCol);
    } else {
        /* We directly compute A * S^-1 */
        for ( int iRowElem = 0; iRowElem < spr1->nSpR1FactorElem; ++iRowElem ) {
            for ( int iColElem = 0; iColElem < spr1->nSpR1FactorElem; ++iColElem ) {
                double dAElem = spr1->spR1FactorSign * spr1->spR1MatElem[iRowElem] * spr1->spR1MatElem[iColElem];
                SinvCol = Sinv + spr1->nSDPCol * spr1->spR1MatIdx[iColElem];
                ASinvRow = B + spr1->spR1MatIdx[iRowElem];
                axpy(&spr1->nSDPCol, &dAElem, SinvCol, &HIntConstantOne, ASinvRow, &HIntConstantOne);
            }
        }
    }
    
    /* Compute trace(B) */
    for ( int iCol = 0; iCol < spr1->nSDPCol; ++iCol ) {
        dTraceASinv += B[iCol * spr1->nSDPCol + iCol];
    }
    
    return dTraceASinv;
}

static double dataMatRankOnesDenseKKT4ComputeASinvImpl( void *A, hdsdp_linsys *S, double *Sinv, double *aux, double *B ) {
    
    sdp_coeff_dsr1 *dsr1 = (sdp_coeff_dsr1 *) A;
    
    HDSDP_ZERO(B, double, dsr1->nSDPCol * dsr1->nSDPCol);
    
    double *dSinvAVec = aux;
    double dFactorSign = 0.0;
    double dTraceASinv = 0.0;
    
    dataMatRankOneSparseKKT2SolveRankOneImpl(A, S, Sinv, &dFactorSign, dSinvAVec);
    fds_ger(dsr1->nSDPCol, dsr1->nSDPCol, dFactorSign,
            dsr1->r1MatFactor, 1, dSinvAVec, 1, B, dsr1->nSDPCol);
    
    for ( int iCol = 0; iCol < dsr1->nSDPCol; ++iCol ) {
        dTraceASinv += B[iCol * dsr1->nSDPCol + iCol];
    }
    
    return dTraceASinv;
}

/*==========================================================================================*/
/* KKT operation: Compute trace(A * S^-1 * ASinv)  */
/*==========================================================================================*/
static double dataMatZeroKKT4ADotSinvBImpl( void *A, double *Sinv, double *ASinv, double *aux ) {
    
    return 0.0;
}

static double dataMatSparseKKT4ADotSinvBImpl( void *A, double *Sinv, double *ASinv, double *aux ) {
    
    sdp_coeff_sparse *sparse = (sdp_coeff_sparse *) A;
    
    return 0.0;
}

static double dataMatDenseKKT4ADotSinvBImpl( void *A, double *Sinv, double *ASinv, double *aux ) {
    
    sdp_coeff_dense *dense = (sdp_coeff_dense *) A;
    
    return 0.0;
}

static double dataMatRankOnesSparseKKT4ADotSinvBImpl( void *A, double *Sinv, double *ASinv, double *aux ) {
    
    sdp_coeff_spr1 *spr1 = (sdp_coeff_spr1 *) A;
    
    return 0.0;
}

static double dataMatRankOnesDenseKKT4ADotSinvBImpl( void *A, double *Sinv, double *ASinv, double *aux ) {
    
    sdp_coeff_dsr1 *dsr1 = (sdp_coeff_dsr1 *) A;
    
    return 0.0;
}

/*==========================================================================================*/
/* KKT operation: template */
/*==========================================================================================*/
static void dataMatZeroKKTImpl( void *A ) {
    
    return;
}

static void dataMatSparseKKTImpl( void *A ) {
    
    sdp_coeff_sparse *sparse = (sdp_coeff_sparse *) A;
    
    return;
}

static void dataMatDenseKKTImpl( void *A ) {
    
    sdp_coeff_dense *dense = (sdp_coeff_dense *) A;
    
    return;
}

static void dataMatRankOnesSparseKKTImpl( void *A ) {
    
    sdp_coeff_spr1 *spr1 = (sdp_coeff_spr1 *) A;
    
    return;
}

static void dataMatRankOnesDenseKKTImpl( void *A ) {
    
    sdp_coeff_dsr1 *dsr1 = (sdp_coeff_dsr1 *) A;
    
    return;
}

static void sdpDataMatIChooseType( sdp_coeff *sdpCoeff, sdp_coeff_type dataType ) {
    
    sdpCoeff->dataType = dataType;
     
    switch ( dataType ) {
        case SDP_COEFF_ZERO:
            sdpCoeff->create = dataMatCreateZeroImpl;
            sdpCoeff->scal = dataMatScalZeroImpl;
            sdpCoeff->norm = dataMatNormZeroImpl;
            sdpCoeff->eig = dataMatBuildUpEigZeroImpl;
            sdpCoeff->getnnz = dataMatGetNnzZeroImpl;
            sdpCoeff->dump = dataMatDumpZeroImpl;
            sdpCoeff->getmatnz = dataMatGetZeroSparsityImpl;
            sdpCoeff->add2buffer = dataMatAddZeroToBufferImpl;
            sdpCoeff->destroy = dataMatDestroyZeroImpl;
            sdpCoeff->view = dataMatViewZeroImpl;
            sdpCoeff->kkt2r1solve = dataMatZeroKKT2SolveRankOneImpl;
            sdpCoeff->kkt2quadform = dataMatZeroKKT2ComputeQuadForm;
            sdpCoeff->kkt2r1asinv = dataMatZeroKKT2TraceASinvImpl;
            sdpCoeff->kkt3sinvAsinv = dataMatZeroKKT3ComputeSinvASinvImpl;
            sdpCoeff->kkt3AdotB = dataMatZeroDotDenseKKT3Impl;
            sdpCoeff->kkt4Asinv = dataMatZeroKKT4ComputeASinvImpl;
            sdpCoeff->kkt4AdotsinvB = NULL;
            break;
        case SDP_COEFF_SPARSE:
            sdpCoeff->create = dataMatCreateSparseImpl;
            sdpCoeff->scal = dataMatScalSparseImpl;
            sdpCoeff->norm = dataMatNormSparseImpl;
            sdpCoeff->eig = dataMatBuildUpEigSparseImpl;
            sdpCoeff->getnnz = dataMatGetNnzSparseImpl;
            sdpCoeff->dump = dataMatDumpSparseImpl;
            sdpCoeff->getmatnz = dataMatGetSparseSparsityImpl;
            sdpCoeff->add2buffer = dataMatAddSparseToBufferImpl;
            sdpCoeff->destroy = dataMatDestroySparseImpl;
            sdpCoeff->view = dataMatViewSparseImpl;
            sdpCoeff->kkt2r1solve = dataMatSparseKKT2SolveRankOneImpl;
            sdpCoeff->kkt2quadform = dataMatSparseKKT2ComputeQuadForm;
            sdpCoeff->kkt2r1asinv = dataMatSparseKKT2TraceASinvImpl;
            sdpCoeff->kkt3sinvAsinv = dataMatSparseKKT3ComputeSinvASinvImpl;
            sdpCoeff->kkt3AdotB = dataMatSparseDotDenseKKT3Impl;
            sdpCoeff->kkt4Asinv = dataMatSparseKKT4ComputeASinvImpl;
            sdpCoeff->kkt4AdotsinvB = NULL;
            break;
        case SDP_COEFF_DENSE:
            sdpCoeff->create = dataMatCreateDenseImpl;
            sdpCoeff->scal = dataMatScalDenseImpl;
            sdpCoeff->norm = dataMatNormDenseImpl;
            sdpCoeff->eig = dataMatBuildUpEigDenseImpl;
            sdpCoeff->getnnz = dataMatGetNnzDenseImpl;
            sdpCoeff->dump = dataMatDumpDenseImpl;
            sdpCoeff->getmatnz = dataMatGetDenseSparsityImpl;
            sdpCoeff->add2buffer = dataMatAddDenseToBufferImpl;
            sdpCoeff->destroy = dataMatDestroyDenseImpl;
            sdpCoeff->view = dataMatViewDenseImpl;
            sdpCoeff->kkt2r1solve = dataMatDenseKKT2SolveRankOneImpl;
            sdpCoeff->kkt2quadform = dataMatDenseKKT2ComputeQuadForm;
            sdpCoeff->kkt2r1asinv = dataMatDenseKKT2TraceASinvImpl;
            sdpCoeff->kkt3sinvAsinv = dataMatDenseKKT3ComputeSinvASinvImpl;
            sdpCoeff->kkt3AdotB = dataMatDenseDotDenseKKT3Impl;
            sdpCoeff->kkt4Asinv = dataMatDenseKKT4ComputeASinvImpl;
            sdpCoeff->kkt4AdotsinvB = NULL;
            break;
        case SDP_COEFF_SPR1:
            sdpCoeff->create = dataMatCreateRankOneSparseImpl;
            sdpCoeff->scal = dataMatScalRankOneSparseImpl;
            sdpCoeff->norm = dataMatNormRankOneSparseImpl;
            sdpCoeff->eig = dataMatBuildUpRankOneSparseImpl;
            sdpCoeff->getnnz = dataMatGetNnzRankOneSparseImpl;
            sdpCoeff->dump = dataMatDumpRankOneSparseImpl;
            sdpCoeff->getmatnz = dataMatGetRankOneSparseSparsityImpl;
            sdpCoeff->add2buffer = dataMatAddRankOneSparseToBufferImpl;
            sdpCoeff->destroy = dataMatDestroyRankOneSparseImpl;
            sdpCoeff->view = dataMatViewRankOneSparseImpl;
            sdpCoeff->kkt2r1solve = dataMatRankOneSparseKKT2SolveRankOneImpl;
            sdpCoeff->kkt2quadform = dataMatRankOneSparseKKT2ComputeQuadForm;
            sdpCoeff->kkt2r1asinv = dataMatRankOneSparseKKT2TraceASinvImpl;
            sdpCoeff->kkt3sinvAsinv = dataMatRankOneSparseKKT3ComputeSinvASinvImpl;
            sdpCoeff->kkt3AdotB = dataMatRankOneSparseDotDenseKKT3Impl;
            sdpCoeff->kkt4Asinv = dataMatRankOnesSparseKKT4ComputeASinvImpl;
            sdpCoeff->kkt4AdotsinvB = NULL;
            break;
        case SDP_COEFF_DSR1:
            sdpCoeff->create = dataMatCreateRankOneDenseImpl;
            sdpCoeff->scal = dataMatScalRankOneDenseImpl;
            sdpCoeff->norm = dataMatNormRankOneDenseImpl;
            sdpCoeff->eig = dataMatBuildUpRankOneDenseImpl;
            sdpCoeff->getnnz = dataMatGetNnzRankOneDenseImpl;
            sdpCoeff->dump = dataMatDumpRankOneDenseImpl;
            sdpCoeff->getmatnz = dataMatGetRankOneDenseSparsityImpl;
            sdpCoeff->add2buffer = dataMatAddRankOneDenseToBufferImpl;
            sdpCoeff->destroy = dataMatDestroyRankOneDenseImpl;
            sdpCoeff->view = dataMatViewRankOneDenseImpl;
            sdpCoeff->kkt2r1solve = dataMatRankOneDenseKKT2SolveRankOneImpl;
            sdpCoeff->kkt2quadform = dataMatRankOneDenseKKT2ComputeQuadForm;
            sdpCoeff->kkt2r1asinv = dataMatRankOneDenseKKT2TraceASinvImpl;
            sdpCoeff->kkt3sinvAsinv = dataMatRankOneDenseKKT3ComputeSinvASinvImpl;
            sdpCoeff->kkt3AdotB = dataMatRankOneDenseDotDenseKKT3Impl;
            sdpCoeff->kkt4Asinv = dataMatRankOnesDenseKKT4ComputeASinvImpl;
            sdpCoeff->kkt4AdotsinvB = NULL;
            break;
        default:
            assert( 0 );
            break;
    }
    
    return;
}

/* External methods for the SDP data */
extern hdsdp_retcode sdpDataMatCreate( sdp_coeff **psdpCoeff ) {
    
    hdsdp_retcode retcode = HDSDP_RETCODE_OK;
    
    if ( !psdpCoeff ) {
        retcode = HDSDP_RETCODE_FAILED;
        goto exit_cleanup;
    }
    
    sdp_coeff *sdpCoeff = NULL;
    HDSDP_INIT(sdpCoeff, sdp_coeff, 1);
    HDSDP_MEMCHECK(sdpCoeff);
    
    /* -1 means there is no eigen-decomposition available */
    sdpCoeff->eigRank = -1;
    *psdpCoeff = sdpCoeff;
    
exit_cleanup:
    
    return retcode;
}

/** @brief Set SDP data matrix. This routine selects data type and assiciate structure with their operations
 * 
 *
 */
extern hdsdp_retcode sdpDataMatSetData( sdp_coeff *sdpCoeff, int nSDPCol, int dataMatNnz, int *dataMatIdx, double *dataMatElem ) {
    
    hdsdp_retcode retcode = HDSDP_RETCODE_OK;
    sdpCoeff->nSDPCol = nSDPCol;
    
    /* Choose data matrix type */
    int nPack = PACK_NNZ(nSDPCol);
    
    /* At this stage, only sparse, zero and dense matrices are classified */
    if ( dataMatNnz == 0 ) {
        sdpDataMatIChooseType(sdpCoeff, SDP_COEFF_ZERO);
    } else if ( dataMatNnz > 0.3 * nPack ) {
        sdpDataMatIChooseType(sdpCoeff, SDP_COEFF_DENSE);
    } else {
        sdpDataMatIChooseType(sdpCoeff, SDP_COEFF_SPARSE);
    }
    
    /* Create data */
    HDSDP_CALL(sdpCoeff->create(&sdpCoeff->dataMat, nSDPCol,
                                dataMatNnz, dataMatIdx, dataMatElem));
    
exit_cleanup:
    
    return retcode;
}

extern int sdpDataMatGetRank( sdp_coeff *sdpCoeff ) {
    
    if ( sdpCoeff->dataType == SDP_COEFF_ZERO ) {
        return 0;
    } else if ( sdpCoeff->dataType == SDP_COEFF_DSR1 || sdpCoeff->dataType == SDP_COEFF_SPR1 ) {
        return 1;
    } else if ( sdpCoeff->eigRank != -1 ) {
        return sdpCoeff->eigRank;
    }
    
    return sdpCoeff->nSDPCol;
}

extern void sdpDataMatScal( sdp_coeff *sdpCoeff, double scal ) {
    
    sdpCoeff->scal(sdpCoeff->dataMat, scal);
    
    return;
}

extern double sdpDataMatNorm( sdp_coeff *sdpCoeff, int type ) {
    
    return sdpCoeff->norm(sdpCoeff->dataMat, type);
}

/* Used in presolving and low-rank structure detection */
extern hdsdp_retcode sdpDataMatBuildUpEigs( sdp_coeff *sdpCoeff, double *dAuxFullMatrix ) {
    
    hdsdp_retcode retcode = HDSDP_RETCODE_OK;
    
    HDSDP_ZERO(dAuxFullMatrix, double, sdpCoeff->nSDPCol);
    HDSDP_CALL(sdpCoeff->eig(sdpCoeff->dataMat, &sdpCoeff->eigRank, dAuxFullMatrix,
                             &sdpCoeff->eigVals, &sdpCoeff->eigVecs));
    
    if ( sdpCoeff->eigRank != 1 ) {
        goto exit_cleanup;
    }
    
    /* The matrix is rank-one */
    /* Count nnz in the decomposition */
    int nRankOneNz = 0;
    for ( int iRow = 0; iRow < sdpCoeff->nSDPCol; ++iRow ) {
        if ( fabs(sdpCoeff->eigVecs[iRow]) > 1e-10 ) {
            nRankOneNz += 1;
        }
    }
    
    /* We no longer need it */
    sdpCoeff->destroy(&sdpCoeff->dataMat);
    
    int useDense = 0;
    if ( nRankOneNz > 0.5 * sdpCoeff->nSDPCol ) {
        useDense = 1;
    }
    
    if ( useDense ) {
        
        sdp_coeff_dsr1 *dsr1 = NULL;
        
        HDSDP_INIT(dsr1, sdp_coeff_dsr1, 1);
        HDSDP_MEMCHECK(dsr1);
        
        dsr1->nSDPCol = sdpCoeff->nSDPCol;
        dsr1->r1FactorSign = sdpCoeff->eigVals[0];
        
        HDSDP_INIT(dsr1->r1MatFactor, double, sdpCoeff->nSDPCol);
        HDSDP_MEMCHECK(dsr1->r1MatFactor);
        HDSDP_MEMCPY(dsr1->r1MatFactor, sdpCoeff->eigVecs, double, sdpCoeff->nSDPCol);
        
        dataMatNormalizeRankOneDenseImpl(dsr1);
        sdpCoeff->dataMat = dsr1;
        sdpDataMatIChooseType(sdpCoeff, SDP_COEFF_DSR1);
        
    } else {
        
        sdp_coeff_spr1 *spr1 = NULL;
        
        HDSDP_INIT(spr1, sdp_coeff_spr1, 1);
        HDSDP_MEMCHECK(spr1);
        
        spr1->nSDPCol = sdpCoeff->nSDPCol;
        spr1->spR1FactorSign = sdpCoeff->eigVals[0];
        spr1->nSpR1FactorElem = nRankOneNz;
        
        HDSDP_INIT(spr1->spR1MatIdx, int, nRankOneNz);
        HDSDP_INIT(spr1->spR1MatElem, double, nRankOneNz);
        
        int iNz = 0;
        for ( int iRow = 0; iRow < sdpCoeff->nSDPCol; ++iRow ) {
            if ( fabs(sdpCoeff->eigVecs[iRow] ) > 1e-10 ) {
                spr1->spR1MatIdx[iNz] = iRow;
                spr1->spR1MatElem[iNz] = sdpCoeff->eigVecs[iRow];
                iNz += 1;
            }
        }
        
        assert( iNz == nRankOneNz );
        dataMatNormalizeRankOneSparseImpl(spr1);
        sdpCoeff->dataMat = spr1;
        sdpDataMatIChooseType(sdpCoeff, SDP_COEFF_SPR1);
    }
    
    /* Reset rank */
    sdpCoeff->eigRank = -1;
    HDSDP_FREE(sdpCoeff->eigVals);
    HDSDP_FREE(sdpCoeff->eigVecs);
    
exit_cleanup:
    
    return retcode;
}

extern inline int sdpDataMatGetNnz( sdp_coeff *sdpCoeff ) {
    
    return sdpCoeff->getnnz(sdpCoeff->dataMat);
}

extern inline void sdpDataMatDump( sdp_coeff *sdpCoeff, double *dFullMatrix ) {
    
    return sdpCoeff->dump(sdpCoeff->dataMat, dFullMatrix);
}

extern void sdpDataMatGetMatNz( sdp_coeff *sdpCoeff, int *iMatSpsPattern ) {
    
    sdpCoeff->getmatnz(sdpCoeff->dataMat, iMatSpsPattern);
    
    return;
}

extern void sdpDataMatAddToBuffer( sdp_coeff *sdpCoeff, double dElem, int *iMatSpsPattern, double *dBuffer ) {
    
    if ( dElem ) {
        sdpCoeff->add2buffer(sdpCoeff->dataMat, dElem, iMatSpsPattern, dBuffer);
    }
    
    return;
}

extern inline sdp_coeff_type sdpDataMatGetType( sdp_coeff *sdpCoeff ) {
    
    return sdpCoeff->dataType;
}

extern inline void sdpDataMatClear( sdp_coeff *sdpCoeff ) {
    
    if ( !sdpCoeff ) {
        return;
    }
    
    sdpCoeff->destroy(&sdpCoeff->dataMat);
    
    if ( sdpCoeff->eigRank != -1 ) {
        HDSDP_FREE(sdpCoeff->eigVals);
        HDSDP_FREE(sdpCoeff->eigVecs);
    }
    
    HDSDP_ZERO(sdpCoeff, sdp_coeff, 1);
    
    return;
}

extern inline void sdpDataMatDestroy( sdp_coeff **psdpCoeff ) {
    
    if ( !psdpCoeff ) {
        return;
    }
    
    sdpDataMatClear(*psdpCoeff);
    HDSDP_FREE(*psdpCoeff);
    
    return;
}

extern inline void sdpDataMatView( sdp_coeff *sdpCoeff ) {
    
    sdpCoeff->view(sdpCoeff->dataMat);
    
    return;
}

/* Data-dependent KKT operations */
extern void sdpDataMatKKT2SolveRankOne( sdp_coeff *sdpCoeff, hdsdp_linsys *dualFactor, double *dInvMatrix, double *dInvFactor, double *sign ) {
    
    sdpCoeff->kkt2r1solve(sdpCoeff->dataMat, dualFactor, dInvMatrix, dInvFactor, sign);
    return;
}

extern double sdpDataMatKKT2QuadForm( sdp_coeff *sdpCoeff, double *dQuadVector, double *dAuxiVec ) {
    
    return sdpCoeff->kkt2quadform(sdpCoeff->dataMat, dQuadVector, dAuxiVec);
}

extern double sdpDataMatKKT2TraceASinv( sdp_coeff *sdpCoeff, double *dSinvAVec ) {
    
    return sdpCoeff->kkt2r1asinv(sdpCoeff->dataMat, dSinvAVec);
}

extern double sdpDataMatKKT3ComputeSinvASinv( sdp_coeff *sdpCoeff, hdsdp_linsys *dualFactor, double *dInvMatrix,
                                             double *dAuxiMat, double *dSinvASinvBuffer ) {
    
    return sdpCoeff->kkt3sinvAsinv(sdpCoeff->dataMat, dualFactor, dInvMatrix, dAuxiMat, dSinvASinvBuffer);
}

extern double sdpDataMatKKT3TraceABuffer( sdp_coeff *sdpCoeff, double *B, double *aux ) {
    
    return sdpCoeff->kkt3AdotB(sdpCoeff->dataMat, B, aux);
}
