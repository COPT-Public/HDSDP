/** @file hdsdp\_file\_io.c
 *  @author Wenzhi Gao
 * HDSDP input and output for SDPA format
 *
 */

#include "interface/hdsdp_utils.h"
#include "interface/hdsdp_file_io.h"
#include "interface/hdsdp_user_data.h"

#include "external/hdsdp_cs.h"

#include <stdio.h>
#include <math.h>

#ifndef SDPA_BUFFERSIZE
#define SDPA_BUFFERSIZE 1024
#define read_line(file, target, line) fgets(target, SDPA_BUFFERSIZE, file); \
                                      ++line;
#endif


/** @function HReadSDPA
 * Read SDPA data from file
 *
 * In a word we insert data into CSparse triplet format and then compress it
 */
hdsdp_retcode HReadSDPA( char *fname, int *pnConstrs, int *pnBlks, int **pblkDims, double **prowRHS,
                         int ***pconeMatBeg, int ***pconeMatIdx, double ***pconeMatElem, int *pnCols, int *pnLPCols,
                         int **pLpMatBeg, int **pLpMatIdx, double **pLpMatElem, int *pnElems ) {
    
    hdsdp_retcode retcode = HDSDP_RETCODE_OK;
    FILE *file;
    
    char fileLine[SDPA_BUFFERSIZE] = "*";
    
    int whichLine = 0;
    int nConstrs = 0;
    int nCols = 0;
    int nLPCols = 0;
    int nBlks = 0;
    int nElems = 0;
    int iBlk = 0;
    int iCon = 0;
    int blkDim = 0;
    
    /* Allocated memory */
    int *blkDims = NULL;
    double *rowRHS = NULL;
    int **coneMatBeg = NULL;
    int **coneMatIdx = NULL;
    int *LpMatBeg = NULL;
    int *LpMatIdx = NULL;
    double *LpMatElem = NULL;
    double **coneMatElem = NULL;
    
    /* Auxiliary memory */
    dcs *LpConeData = NULL;
    dcs **pSDPConeData = NULL;
    
    /* Open SDPA file for reading */
    file = fopen(fname, "r");
    
    if ( !file ) {
        HDSDP_ERROR_TRACE;
        retcode = HDSDP_RETCODE_FAILED;
        goto exit_cleanup;
    }
    
    /* Jump through comments */
    while ( !feof(file) && (fileLine[0] == '*' ||
                            fileLine[0] == '"') ) {
        read_line(file, fileLine, whichLine);
    }
    
    /* Get number of constraints */
    if ( sscanf(fileLine, "%d", &nConstrs) != 1 ) {
        HDSDP_ERROR_TRACE;
        retcode = HDSDP_RETCODE_FAILED;
        goto exit_cleanup;
    }
    
    read_line(file, fileLine, whichLine);
    
    /* Get number of blocks */
    if ( sscanf(fileLine, "%d", &nBlks) != 1 ) {
        HDSDP_ERROR_TRACE;
        retcode = HDSDP_RETCODE_FAILED;
        goto exit_cleanup;
    }
    
    /* Get all the dimensions */
    ++whichLine;
    HDSDP_INIT(blkDims, int, nBlks);
    HDSDP_MEMCHECK(blkDims);
    
    for ( iBlk = 0; iBlk < nBlks - 1; ++iBlk ) {
        if ( fscanf(file, "{") == 1 || fscanf(file, "(") == 1 ||
             fscanf(file, "'") == 1 ) {
            --iBlk;
        } else if ( fscanf(file, "%d", &blkDim) == 1 ) {
            if ( blkDim <= 0 ) {
                /* Reason: only one diagonal block is supported and
                           it must appear at the end
                 */
                HDSDP_ERROR_TRACE;
                retcode = HDSDP_RETCODE_FAILED;
                goto exit_cleanup;
            }
            nCols += blkDim * blkDim;
            blkDims[iBlk] = blkDim;
            
        } else {
            HDSDP_ERROR_TRACE;
            retcode = HDSDP_RETCODE_FAILED;
            goto exit_cleanup;
        }
    }
    
    /* Special treatment for LP */
    if ( fscanf(file, "%d", &blkDim) == 1 ) {
        if ( blkDim < 0 ) {
            --nBlks;
            nLPCols = -blkDim;
        } else {
            blkDims[iBlk] = blkDim;
            nCols += blkDim * blkDim;
        }
    }
    
    /* Read dual objective / primal RHS */
    read_line(file, fileLine, whichLine);
    HDSDP_INIT(rowRHS, double, nConstrs);
    HDSDP_MEMCHECK(rowRHS);
    
    int iCol;
    int iRow;
    double dElem;
    char charTmp;
    
    for ( iRow = 0; iRow < nConstrs; ++iRow ) {
        
        if ( fscanf(file, ",") == 1 ) {
            --iRow;
            continue;
        }
        
        while ( fscanf(file, "%lg", &dElem) != 1 ) {
            fscanf(file, "%c", &charTmp);
            if ( charTmp == '\n' ) {
                HDSDP_ERROR_TRACE;
                retcode = HDSDP_RETCODE_FAILED;
                goto exit_cleanup;
            }
        }
        
        rowRHS[iRow] = dElem;
    }
    
    /* Reset */
    fgets(fileLine, SDPA_BUFFERSIZE, file);
    int lineNum = whichLine;
    fseek(file, 0, SEEK_SET);
    whichLine = 0;
    
    for ( int iLine = 0; iLine < lineNum; ++iLine ) {
        charTmp = '*';
        while ( charTmp != '\n' ) {
            fscanf(file, "%c", &charTmp);
            ++whichLine;
        }
    }
    
    /* Read conic data via CSparse */
    int isMemOK = 1;
    
    if ( nLPCols > 0 ) {
        LpConeData = dcs_spalloc(nLPCols, nConstrs + 1, nConstrs, 1, 1);
        HDSDP_MEMCHECK(LpConeData);
    }
    
    HDSDP_INIT(pSDPConeData, dcs *, nBlks);
    HDSDP_MEMCHECK(pSDPConeData);
    
    for ( iBlk = 0; iBlk < nBlks; ++iBlk ) {
        pSDPConeData[iBlk] = dcs_spalloc(PACK_NNZ(blkDims[iBlk]), nConstrs + 1, 100, 1, 1);
        HDSDP_MEMCHECK(pSDPConeData[iBlk]);
    }
    
    int LpConeID = -1;
    if ( nLPCols > 0 ) {
        LpConeID = nBlks;
    }
    
    while( !feof(file) ) {
        fileLine[0] = '\0';
        read_line(file, fileLine, whichLine);
        
        if ( sscanf(fileLine, "%d %d %d %d %lg", &iCon, &iBlk, &iRow, &iCol, &dElem) != 5 ) {
            
            if ( feof(file) ) {
                break;
            } else {
                HDSDP_ERROR_TRACE;
                retcode = HDSDP_RETCODE_FAILED;
                goto exit_cleanup;
            }
            
        } else {
        
            /* To 0-based indexing */
            iBlk -= 1;
            iRow -= 1;
            iCol -= 1;
            
            if ( fabs(dElem) < 1e-10 ) {
                printf("[Warning] Entry smaller than 1e-10 is ignored. \n");
                continue;
            }
            
            if ( iBlk == LpConeID ) {
                isMemOK = dcs_entry(LpConeData, iRow, iCol, dElem);
            } else {
                
                if ( iRow > iCol ) {
                    int iTmp = iRow;
                    iRow = iCol;
                    iCol = iTmp;
                }
                isMemOK = dcs_entry(pSDPConeData[iBlk], PACK_IDX(blkDims[iBlk], iCol, iRow), iCon, dElem);
            }
            
            nElems += 1;
        }
    }
    
    if ( !isMemOK ) {
        HDSDP_ERROR_TRACE;
        retcode = HDSDP_RETCODE_MEMORY;
        return retcode;
    }
    
    /* Compress data and complete conversion */
    dcs *cscMat = NULL;
    for ( iBlk = 0; iBlk < nBlks; ++iBlk ) {
        cscMat = dcs_compress(pSDPConeData[iBlk]);
        HDSDP_MEMCHECK(cscMat);
        dcs_spfree(pSDPConeData[iBlk]);
        pSDPConeData[iBlk] = cscMat;
    }
    
    if ( nLPCols > 0 ) {
        cscMat = dcs_compress(LpConeData);
        dcs_spfree(LpConeData);
        LpConeData = cscMat;
    }
    
    /* Copy the data out */
    HDSDP_INIT(coneMatBeg, int *, nBlks);
    HDSDP_INIT(coneMatIdx, int *, nBlks);
    HDSDP_INIT(coneMatElem, double *, nBlks);
    
    HDSDP_MEMCHECK(coneMatBeg);
    HDSDP_MEMCHECK(coneMatIdx);
    HDSDP_MEMCHECK(coneMatElem);
    
    for ( iBlk = 0; iBlk < nBlks; ++iBlk ) {
        
        cscMat = pSDPConeData[iBlk];
        HDSDP_INIT(coneMatBeg[iBlk], int, cscMat->n + 1);
        HDSDP_INIT(coneMatIdx[iBlk], int, cscMat->p[cscMat->n]);
        HDSDP_INIT(coneMatElem[iBlk], double, cscMat->p[cscMat->n]);
        
        HDSDP_MEMCHECK(coneMatBeg[iBlk]);
        HDSDP_MEMCHECK(coneMatIdx[iBlk]);
        HDSDP_MEMCHECK(coneMatElem[iBlk]);
        
        HDSDP_MEMCPY(coneMatBeg[iBlk], cscMat->p, int, cscMat->n + 1);
        HDSDP_MEMCPY(coneMatIdx[iBlk], cscMat->i, int, cscMat->p[cscMat->n]);
        HDSDP_MEMCPY(coneMatElem[iBlk], cscMat->x, double, cscMat->p[cscMat->n]);
        
    }
    
    if ( nLPCols > 0 ) {
        
        HDSDP_INIT(LpMatBeg, int, LpConeData->n + 1);
        HDSDP_INIT(LpMatIdx, int, LpConeData->p[LpConeData->n]);
        HDSDP_INIT(LpMatElem, double, LpConeData->p[LpConeData->n]);
        
        HDSDP_MEMCHECK(LpMatBeg);
        HDSDP_MEMCHECK(LpMatIdx);
        HDSDP_MEMCHECK(LpMatElem);
        
        HDSDP_MEMCPY(LpMatBeg, LpConeData->p, int, LpConeData->n + 1);
        HDSDP_MEMCPY(LpMatIdx, LpConeData->i, int, LpConeData->p[cscMat->n]);
        HDSDP_MEMCPY(LpMatElem, LpConeData->x, double, LpConeData->p[cscMat->n]);
    }
     
    /* Get the outputs */
    *pnConstrs = nConstrs;
    *pnBlks = nBlks;
    *pblkDims = blkDims;
    *prowRHS = rowRHS;
    *pconeMatBeg = coneMatBeg;
    *pconeMatIdx = coneMatIdx;
    *pconeMatElem = coneMatElem;
    *pnCols = nCols;
    *pnLPCols = nLPCols;
    *pLpMatBeg = LpMatBeg;
    *pLpMatIdx = LpMatIdx;
    *pLpMatElem = LpMatElem;
    *pnElems = nElems;
    
exit_cleanup:
    
    /* Auxiliary memory, always freed */
    dcs_spfree(LpConeData);
    
    if ( pSDPConeData ) {
        for ( iBlk = 0; iBlk < nBlks; ++iBlk ) {
            dcs_spfree(pSDPConeData[iBlk]);
        }
        HDSDP_FREE(pSDPConeData);
    }
    
    if ( retcode != HDSDP_RETCODE_OK ) {
        /* Free all the internal memories if failed */
        HDSDP_FREE(blkDims);
        HDSDP_FREE(rowRHS);
        
        HDSDP_FREE(LpMatBeg);
        HDSDP_FREE(LpMatIdx);
        HDSDP_FREE(LpMatElem);
        
        if ( coneMatBeg ) {
            for ( iBlk = 0; iBlk < nBlks; ++iBlk ) {
                HDSDP_FREE(coneMatBeg[iBlk]);
            }
            HDSDP_FREE(coneMatBeg);
        }
        
        if ( coneMatIdx ) {
            for ( iBlk = 0; iBlk < nBlks; ++iBlk ) {
                HDSDP_FREE(coneMatIdx[iBlk]);
            }
            HDSDP_FREE(coneMatIdx);
        }
        
        if ( coneMatElem ) {
            for ( iBlk = 0; iBlk < nBlks; ++iBlk ) {
                HDSDP_FREE(coneMatElem[iBlk]);
            }
            HDSDP_FREE(coneMatElem);
        }
        
    }
    
    return retcode;
}
