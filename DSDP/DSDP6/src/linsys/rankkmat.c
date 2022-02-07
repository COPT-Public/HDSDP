#include <stdio.h>
#include "rankonemat.h"
#include "rankkmat.h"
#include "dsdplapack.h"
#include "dsdpdata.h"

static char etype[] = "Rank-k matrix computation";

extern DSDP_INT rkMatInit( rkMat *R ) {
    // Initialize rank k matrix
    R->rank = 0;
    R->dim  = 0;
    R->isdata = FALSE;
    R->data = NULL;
    R->mattype = DSDP_UNKNOWN;
    R->origdata = NULL;
    return DSDP_RETCODE_OK;
}

extern DSDP_INT rkMatAllocIter( rkMat *R, DSDP_INT n ) {
    
    assert( n > 0 );
    R->dim = n;
    R->isdata = FALSE;
    R->rank = n;
    
    R->data = (r1Mat **) calloc(n, sizeof(r1Mat *));
    if (!R->data) return DSDP_RETCODE_FAILED;
    r1Mat *r1data = NULL;
    
    for (DSDP_INT i = 0; i < n; ++i) {
        r1data = (r1Mat *) calloc(1, sizeof(r1Mat));
        R->data[i] = r1data;
        r1MatInit(r1data);
        r1MatAlloc(r1data, n);
    }
    
    return DSDP_RETCODE_OK;
}

extern DSDP_INT rkMatAllocAndSetData( rkMat *R, DSDP_INT n, DSDP_INT rank,
                                     double *eigvals, double *eigvecs ) {
    // Allocate memory for the r1Mat array
    DSDP_INT retcode = DSDP_RETCODE_OK;
    
    assert( R->dim == 0 && rank <= n );
    R->dim  = n;
    R->data = (r1Mat **) calloc(rank, sizeof(r1Mat *));
    R->rank = rank;
    R->isdata = TRUE;
    
    r1Mat *r1data = NULL;
    
    for (DSDP_INT i = 0; i < rank; ++i) {
        r1data = (r1Mat *) calloc(1, sizeof(r1Mat));
        R->data[i] = r1data;
        r1MatInit(r1data);
        r1MatAlloc(r1data, n);
        retcode = r1MatSetData(r1data, eigvals[i], &eigvecs[i * n]);
    }
    
    return retcode;
}

extern DSDP_INT rkMatStoreOriginalData( rkMat *R, DSDP_INT mattype, void *data ) {
    
    assert( data || mattype == MAT_TYPE_RANKK );
    R->origdata = data;
    R->mattype = mattype;
    
    return DSDP_RETCODE_OK;
}

extern DSDP_INT rkMatAllocAndSelectData( rkMat *R, DSDP_INT n, DSDP_INT rank, double thresh,
                                         double *eigvals, double *eigvecs ) {
    
    DSDP_INT retcode = DSDP_RETCODE_OK;
    R->dim  = n;
    R->data = (r1Mat **) calloc(rank, sizeof(r1Mat *));
    R->rank = rank;
    R->isdata = TRUE;
    r1Mat *r1data = NULL;
    DSDP_INT counter = 0;
    
    for (DSDP_INT i = 0; i < n; ++i) {
        if (fabs(eigvals[i]) > thresh) {
            r1data = (r1Mat *) calloc(1, sizeof(r1Mat));
            R->data[counter] = r1data;
            r1MatInit(r1data);
            r1MatAlloc(r1data, n);
            retcode = r1MatSetData(r1data, eigvals[i], &eigvecs[i * n]);
            counter += 1;
            if (counter == rank) {
                break;
            }
        }
    }
    
    return retcode;
}

extern DSDP_INT rkMatrkTrace( rkMat *R1, rkMat *R2, double *trace ) {
    
    // Compute the inner product between two rank-k matrices
    /*
      trace(R1 * R2) = trace( \sum_i \sum_j c_i * d_j a_i a_i' * b_j * b_j')
                     = \sum_i \sum_j c_i * d_j trace( a_i a_i' * b_j * b_j')
                     = \sum_i \sum_j c_i * d_j (a_i' * b_j)^2
     
     Implemented by calling r1Matr1Trace
     
     When this routine is called, R1 is an iterator (SinvASinv) and R2 is data (A)
    */
    
    DSDP_INT retcode = DSDP_RETCODE_OK;
    double res = 0.0, tmp = 0.0;
    DSDP_INT n = R1->dim, i, j, useRk = TRUE;
    if (R2->mattype == MAT_TYPE_SPARSE) {
        spsMat *spsdata = R2->origdata;
        if (spsdata->nnz < 2 * n) {
            useRk = FALSE;
            for (i = 0; i < R2->rank; ++i) {
                spsMatxTAx(spsdata, R2->data[i]->x, &tmp);
                res += tmp;
            }
        }
    }
    
    if (useRk) {
        for (i = 0; i < R1->rank; ++i) {
            for (j = 0; j < R2->rank; ++j) {
                retcode = r1Matr1Trace(R1->data[i], R2->data[j], &tmp);
                res += tmp;
            }
        }
    }
    
    *trace = res;
    return retcode;
}

extern DSDP_INT rkMatdenseTrace( rkMat *R, dsMat *A, double *trace ) {
    
    // Compute the innter product between rank-k and dense matrices
    /*
     trace(R1 * D) = trace( \sum_i d_i * a_i * a_i' * D )
                   = \sum_i trace( a_i * a_i' * D )
                   = \sum_i d_i * (a_i' * D * a_i)
     
     Implemented by calling r1MatdenseTrace
     
   */
    
    DSDP_INT retcode = DSDP_RETCODE_OK;
    assert( R->dim == A->dim );
    DSDP_INT rank = R->rank;
    double res = 0.0, tmp;
    
    for (DSDP_INT i = 0; i < rank; ++i) {
        retcode = r1MatdenseTrace(R->data[i], A, &tmp); checkCode;
        res += tmp;
    }
    
    *trace = res;
    return retcode;
}

extern DSDP_INT rkMatspsTrace( rkMat *R, spsMat *A, double *trace ) {
    
    // Compute the inner product between rank-1 and sparse matrix
    /*
     Implemented by calling r1MatspsTrace
    */
    
    DSDP_INT retcode = DSDP_RETCODE_OK;
    assert( R->dim == A->dim );
    double res = 0.0, tmp;
    DSDP_INT rank = R->rank;
    
    for (DSDP_INT i = 0; i < rank; ++i) {
        retcode = r1MatspsTrace(R->data[i], A, &tmp); checkCode;
        res += tmp;
    }
    
    *trace = res;
    return retcode;
}

extern DSDP_INT rkMatdiagTrace( rkMat *R, double diag, double *trace ) {
    // Compute trace( \sum_i d_i * a_i * a_i' diag ) = \sum_i diag * d_i * norm(a_i)^2
    DSDP_INT retcode = DSDP_RETCODE_OK;
    
    if (diag == 0.0) {
        *trace = 0.0;
        return retcode;
    }
    
    double res = 0.0, tmp;
    DSDP_INT rank = R->rank;
    
    for (DSDP_INT i = 0; i < rank; ++i) {
        retcode = r1MatdiagTrace(R->data[i], diag, &tmp); checkCode;
        res += tmp;
    }
    
    *trace = res;
    return retcode;
}

extern DSDP_INT rkMatCountNnz( rkMat *R ) {
    
    for (DSDP_INT i = 0; i < R->rank; ++i) {
        r1MatCountNnz(R->data[i]);
    }
    
    return DSDP_RETCODE_OK;
}

extern DSDP_INT rkMatFree( rkMat *R ) {
    
    if (R) {
        if (R->isdata) {
            for (DSDP_INT i = 0; i < R->rank; ++i) {
                r1MatFree(R->data[i]);
                DSDP_FREE(R->data[i]);
            }
        } else {
            for (DSDP_INT i = 0; i < R->dim; ++i) {
                r1MatFree(R->data[i]);
                DSDP_FREE(R->data[i]);
            }
        }
    }
    
    R->dim    = 0;
    R->rank   = 0;
    R->isdata = FALSE;
    
    if (R->mattype == MAT_TYPE_DENSE) {
        denseMatFree(R->origdata);
        DSDP_FREE(R->origdata);
    }
    
    if (R->mattype == MAT_TYPE_SPARSE) {
        spsMatFree(R->origdata);
        DSDP_FREE(R->origdata);
    }
        
    return DSDP_RETCODE_OK;
}

extern DSDP_INT rkMatFnorm( rkMat *R, double *fnrm ) {
    /* Compute the Frobenius norm of rank k matrix
            
      ||\sum_i d_i a_i a_i||F^2 = trace ((\sum_i d_i a_i a_i) (\sum_i d_i a_i a_i))
                                = d_i^2 norm(a_i)^4
     */
    
    assert( R->dim );
    
    double res = 0.0, tmp;
    DSDP_INT rank = R->rank;
    
    for (DSDP_INT i = 0; i < rank; ++i) {
        r1MatFnorm(R->data[i], &tmp);
        res += tmp * tmp;
    }
    
    *fnrm = sqrt(res);
    return DSDP_RETCODE_OK;
}

extern DSDP_INT rkMatRscale( rkMat *R, double r ) {
    
    DSDP_INT retcode = DSDP_RETCODE_OK;
    DSDP_INT rank = R->rank;
    for (DSDP_INT i = 0; i < rank; ++i) {
        r1MatRscale(R->data[i], r);
    }
    
    switch (R->mattype) {
        case MAT_TYPE_SPARSE:
            spsMatRscale(R->origdata, r);
            break;
        case MAT_TYPE_DENSE:
            denseMatRscale(R->origdata, r);
            break;
        case MAT_TYPE_RANKK:
            break;
        case MAT_TYPE_ZERO:
            break;
        default:
            error(etype, "Invalid Matrix type. \n");
            break;
    }
    
    return retcode;
}

extern DSDP_INT rkMatisRank1( rkMat *R, DSDP_INT *isRank1 ) {
    
    if (R->rank == 1) {
        *isRank1 = TRUE;
    } else {
        *isRank1 = FALSE;
    }
    
    return DSDP_RETCODE_OK;
}

extern DSDP_INT rkMatView( rkMat *R ) {
    
    assert( R->dim );
    
    printf("Matrix View: \n");
    printf("Matrix rank "ID" : \n", R->rank);
    for (DSDP_INT i = 0; i < R->rank; ++i) {
        r1MatView(R->data[i]);
    }
    
    return DSDP_RETCODE_OK;
}
