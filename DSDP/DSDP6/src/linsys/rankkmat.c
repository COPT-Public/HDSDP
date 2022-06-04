#include <stdio.h>
#include "rankonemat.h"
#include "rankkmat.h"
#include "dsdplapack.h"
#include "dsdpdata.h"

static char etype[] = "Rank-k matrix computation";
static char uplolow = DSDP_MAT_LOW;
static DSDP_INT one = 1;

extern DSDP_INT rkMatInit( rkMat *R ) {
    // Initialize rank k matrix
    R->rank = 0; R->dim  = 0; R->isdata = FALSE; R->data = NULL;
    return DSDP_RETCODE_OK;
}

extern DSDP_INT rkMatAllocIter( rkMat *R, DSDP_INT n ) {
    
    assert( n > 0 );
    R->dim = n; R->isdata = FALSE; R->rank = n;
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
    R->dim  = n; R->data = (r1Mat **) calloc(rank, sizeof(r1Mat *));
    R->rank = rank; R->isdata = TRUE;
    r1Mat *r1data = NULL;
    for (DSDP_INT i = 0; i < rank; ++i) {
        r1data = (r1Mat *) calloc(1, sizeof(r1Mat));
        R->data[i] = r1data; r1MatInit(r1data); r1MatAlloc(r1data, n);
        retcode = r1MatSetData(r1data, eigvals[i], &eigvecs[i * n]);
    }
    
    return retcode;
}

extern DSDP_INT rkMatAllocAndSelectData( rkMat *R, DSDP_INT n, DSDP_INT rank, double thresh,
                                         double *eigvals, double *eigvecs ) {
    
    DSDP_INT retcode = DSDP_RETCODE_OK;
    R->dim  = n; R->data = (r1Mat **) calloc(rank, sizeof(r1Mat *));
    R->rank = rank; R->isdata = TRUE; r1Mat *r1data = NULL;
    DSDP_INT counter = 0;
    
    for (DSDP_INT i = 0; i < n; ++i) {
        if (fabs(eigvals[i]) > thresh) {
            r1data = (r1Mat *) calloc(1, sizeof(r1Mat));
            R->data[counter] = r1data;
            r1MatInit(r1data); r1MatAlloc(r1data, n);
            retcode = r1MatSetData(r1data, eigvals[i], &eigvecs[i * n]);
            counter += 1; if (counter == rank) { break; }
        }
    }
    
    return retcode;
}

extern void rkMatrkTrace( rkMat *R1, rkMat *R2, double *trace ) {
    
    // Compute the inner product between two rank-k matrices
    /*
      trace(R1 * R2) = trace( \sum_i \sum_j c_i * d_j a_i a_i' * b_j * b_j')
                     = \sum_i \sum_j c_i * d_j trace( a_i a_i' * b_j * b_j')
                     = \sum_i \sum_j c_i * d_j (a_i' * b_j)^2
     
     Implemented by calling r1Matr1Trace
     
     When this routine is called, R1 is an iterator (SinvASinv) and R2 is data (A)
    */
    
    double res = 0.0; DSDP_INT i, j;
    for (i = 0; i < R1->rank; ++i) {
        for (j = 0; j < R2->rank; ++j) {
            res += r1Matr1Trace(R1->data[i], R2->data[j]);
        }
    }
    *trace = res;
}

extern void rkMatdenseUpdate( dsMat *dAMat, rkMat *rkBMat ) {
    // Compute rank-1 update A = A + \sum_i alpha_i * b_i * b_i'
    DSDP_INT dim = rkBMat->dim; r1Mat *r1data = NULL;
    for (DSDP_INT i = 0; i < rkBMat->rank; ++i) {
        r1data = rkBMat->data[i];
        if (r1data->nnz > dim / 4) {
            dspr(&uplolow, &dim, &r1data->sign, r1data->x, &one, dAMat->array);
        } else {
            double sign = r1data->sign, *rx = r1data->x;
            DSDP_INT *nzIdx, j, k;
            nzIdx = r1data->nzIdx;
            for (j = 0; j < r1data->nnz; ++j) {
                for (k = 0; k <= j; ++k) {
                    packIdx(dAMat->array, dim, nzIdx[i], nzIdx[j]) += \
                    sign * rx[nzIdx[i]] * rx[nzIdx[j]];
                }
            }
        }
    }
}

extern void rkMatdenseTrace( rkMat *R, dsMat *A, double *trace ) {
    // Compute the innter product between rank-k and dense matrices
    double res = 0.0, tmp;
    for (DSDP_INT i = 0; i < R->rank; ++i) {
        r1MatdenseTrace(R->data[i], A, &tmp); res += tmp;
    }
    *trace = res;
}

extern void rkMatdiagTrace( rkMat *R, double diag, double *trace ) {
    // Compute trace( \sum_i d_i * a_i * a_i' diag ) = \sum_i diag * d_i * norm(a_i)^2
    if (diag == 0.0) {
        *trace = 0.0;
    }
    double res = 0.0, tmp; DSDP_INT rank = R->rank;
    for (DSDP_INT i = 0; i < rank; ++i) {
        r1MatdiagTrace(R->data[i], diag, &tmp); res += tmp;
    }
    *trace = res;
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
                r1MatFree(R->data[i]); DSDP_FREE(R->data[i]);
            }
        } else {
            for (DSDP_INT i = 0; i < R->dim; ++i) {
                r1MatFree(R->data[i]); DSDP_FREE(R->data[i]);
            }
        }
        R->dim = 0; R->rank = 0; R->isdata = FALSE;
    }
    
    return DSDP_RETCODE_OK;
}

extern DSDP_INT rkMatFnorm( rkMat *R, double *fnrm ) {
    /* Compute the Frobenius norm of rank k matrix */
    double res = 0.0, tmp;
    DSDP_INT rank = R->rank;
    
    for (DSDP_INT i = 0; i < rank; ++i) {
        r1MatFnorm(R->data[i], &tmp); res += tmp * tmp;
    }
    
    *fnrm = sqrt(res);
    return DSDP_RETCODE_OK;
}

extern DSDP_INT rkMatScale( rkMat *R, double a ) {
    
    for (DSDP_INT i = 0; i < R->rank; ++i) {
        r1MatScale(R->data[i], a);
    }
    return DSDP_RETCODE_OK;
}

extern DSDP_INT rkMatRscale( rkMat *R, double r ) {
    
    for (DSDP_INT i = 0; i < R->rank; ++i) {
        r1MatRscale(R->data[i], r);
    }
    
    return DSDP_RETCODE_OK;
}

extern DSDP_INT rkMatisRank1( rkMat *R, DSDP_INT *isRank1 ) {
    
    if (R->rank == 1) {
        *isRank1 = TRUE;
    } else {
        *isRank1 = FALSE;
    }
    
    return DSDP_RETCODE_OK;
}

extern DSDP_INT rkMatGetRank( rkMat *R ) {
    return R->rank;
}

extern r1Mat *rkMatGetBase( rkMat *R, DSDP_INT i) {
    assert( i < R->rank );
    return R->data[i];
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
