#include <stdio.h>
#include "rankonemat.h"
#include "rankkmat.h"
#include "dsdplapack.h"
#include "dsdpdata.h"

static char etype[] = "Rank-k matrix computation";
static char uplolow = DSDP_MAT_LOW;
static DSDP_INT one = 1;

extern void rkMatInit( rkMat *R ) {
    // Initialize rank k matrix
    R->rank = 0; R->dim  = 0; R->isdata = FALSE; R->data = NULL;
}

extern DSDP_INT rkMatAllocIter( rkMat *R, DSDP_INT n ) {
    
    DSDP_INT retcode = DSDP_RETCODE_OK;
    r1Mat *r1data = NULL;
    R->dim = n; R->isdata = FALSE; R->rank = n;
    R->data = (r1Mat **) calloc(n, sizeof(r1Mat *));
    if (!R->data) { retcode = DSDP_RETCODE_FAILED; return retcode; }
    for (DSDP_INT i = 0; i < n; ++i) {
        r1data = (r1Mat *) calloc(1, sizeof(r1Mat));
        R->data[i] = r1data; r1MatInit(r1data);
        if (r1MatAlloc(r1data, n) != DSDP_RETCODE_OK) {
            retcode = DSDP_RETCODE_FAILED; return retcode;
        }
    }
    return retcode;
}

extern DSDP_INT rkMatAllocAndSetData( rkMat *R, DSDP_INT n, DSDP_INT rank,
                                     double *eigvals, double *eigvecs ) {
    // Allocate memory for the r1Mat array
    DSDP_INT retcode = DSDP_RETCODE_OK;
    
    R->dim  = n; R->data = (r1Mat **) calloc(rank, sizeof(r1Mat *));
    R->rank = rank; R->isdata = TRUE;
    r1Mat *r1data = NULL;
    for (DSDP_INT i = 0; i < rank; ++i) {
        r1data = (r1Mat *) calloc(1, sizeof(r1Mat));
        R->data[i] = r1data; r1MatInit(r1data);
        retcode = r1MatAlloc(r1data, n);
        if (retcode != DSDP_RETCODE_OK) { return retcode; }
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
            r1MatInit(r1data);
            retcode = r1MatAlloc(r1data, n);
            if (retcode != DSDP_RETCODE_OK) { return retcode; }
            retcode = r1MatSetData(r1data, eigvals[i], &eigvecs[i * n]);
            if (retcode != DSDP_RETCODE_OK) { return retcode; }
            counter += 1; if (counter == rank) { break; }
        }
    }
    
    return retcode;
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
    if (diag == 0.0) { *trace = 0.0; }
    double res = 0.0, tmp; DSDP_INT rank = R->rank;
    for (DSDP_INT i = 0; i < rank; ++i) {
        r1MatdiagTrace(R->data[i], diag, &tmp); res += tmp;
    }
    *trace = res;
}

extern void rkMatCountNnz( rkMat *R ) {
    
    for (DSDP_INT i = 0; i < R->rank; ++i) {
        r1MatCountNnz(R->data[i]);
    }
}

extern void rkMatFree( rkMat *R ) {
    
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
}

extern void rkMatFnorm( rkMat *R, double *fnrm ) {
    // Frobenius norm of rank k matrix
    double res = 0.0, tmp;
    DSDP_INT rank = R->rank;
    
    for (DSDP_INT i = 0; i < rank; ++i) {
        r1MatFnorm(R->data[i], &tmp); res += tmp * tmp;
    }
    
    *fnrm = sqrt(res);
}

extern void rkMatScale( rkMat *R, double a ) {
    
    for (DSDP_INT i = 0; i < R->rank; ++i) { r1MatScale(R->data[i], a); }
}

extern void rkMatRscale( rkMat *R, double r ) {
    
    for (DSDP_INT i = 0; i < R->rank; ++i) { r1MatRscale(R->data[i], r); }
}

extern void rkMatisRank1( rkMat *R, DSDP_INT *isRank1 ) {
    
    *isRank1 = (R->rank == 1) ? TRUE : FALSE;
}

extern DSDP_INT rkMatGetRank( rkMat *R ) {
    
    return R->rank;
}

extern r1Mat *rkMatGetBase( rkMat *R, DSDP_INT i) {
    
    return R->data[i];
}

extern void rkMatCheckSparsity( rkMat *R, DSDP_INT *isdense, double thresh ) {
    
    for (DSDP_INT i = 0; i < R->rank; ++i) {
        r1MatCheckSparsity(R->data[i], isdense, thresh);
        if (isdense) { break; }
    }
}

extern void rkMatGetSymbolic( rkMat *R, DSDP_INT *hash, DSDP_INT *firstNnz, DSDP_INT *nnzs ) {
    
    for (DSDP_INT i = 0; i < R->rank; ++i) {
        r1MatGetSymbolic(R->data[i], hash, firstNnz, nnzs);
    }
}

extern DSDP_INT rkMatIsConstant( rkMat *R ) {
    
    if (R->rank != 1) {
        return FALSE;
    } else {
        return r1MatIsConstant(R->data[0]);
    }
}

extern void rkMatView( rkMat *R ) {
    
    printf("Matrix View: \n");
    printf("Matrix rank "ID" : \n", R->rank);
    for (DSDP_INT i = 0; i < R->rank; ++i) {
        r1MatView(R->data[i]);
    }
}
