/** @file hdsdp\_solver.c
 *  @brief Implement the conversion interface between COPT and HDSDP
 *
 * The routines implements the following operations
 *
 * - Converting of COPT data to HDSDP format
 * - Invoking the HDSDP solver to solve the problem
 * - Pass solution statistics back
 *
 */


#include "lp/cone_solver.h"

#include "base/copt_prob.h"
#include "base/copt_qmatrix.h"

#include "tool/tool_define.h"
#include "tool/tool_math.h"
#include "tool/tool_memory.h"
#include "tool/tool_message.h"
#include "tool/tool_matconv.h"
#include "tool/tool_index.h"
#include "tool/tool_sort.h"

#include "dsdphsd.h"
#include "dsdpcs.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>


static int ss_HDSDPDataInsert(int iRow, int iCol, double iElem, dcs *sdpBlkData) {
    
    if (dcs_entry(sdpBlkData, iRow, iCol, iElem) != 1)
    {
        dcs_spfree(sdpBlkData);
        return (-1);
    }
    
    return 0;
}

static void ss_HDSDPSol2COPT(double* coptSol, double *hdsdpSol, int dim)
{
    
    
    return;
}

#ifndef HDSDP_CALL
#define HDSDP_CALL(func)                                                                                               \
  do                                                                                                                   \
  {                                                                                                                    \
    if ((dsdpcode = (func)) != DSDP_RETCODE_OK)                                                                        \
    {                                                                                                                  \
      retcode = RETCODE_ERROR;                                                                                         \
      goto exit_cleanup;                                                                                               \
    }                                                                                                                  \
  } while (0)

#define CS_CALL(func)                                                                                               \
  do                                                                                                                   \
  {                                                                                                                    \
    if ((func) != 0)                                                                        \
    {                                                                                                                  \
      retcode = RETCODE_ERROR;                                                                                         \
      goto exit_cleanup;                                                                                               \
    }                                                                                                                  \
  } while (0)
#endif

ss_retcode ss_Prob_SolveHDSDP(copt_prob* prob)
{
    ss_retcode retcode = RETCODE_OK;
    HSDSolver *hdsdp = NULL;
    DSDP_INT dsdpcode = DSDP_RETCODE_OK;
    dcs *sdpBlkData = NULL;
    dcs *cscSDPBlkData = NULL;
    
    int nDualCol = 0;  /* Number of dual columns */
    int nDualRow = 0;  /* Number of dual rows */
    int nDualElem = 0; /* Number of nonzeros in A */

    int nSlackCol = 0; /* Number of slack columns */
    int nObjElem = 0;  /* Number of nonzeros in linear objective */

    /* 'b' */
    double* rowRHS = NULL;

    /* 'A' */
    int* rowMatBeg = NULL;
    int* rowMatIdx = NULL;
    double* rowMatElem = NULL;

    /* 'A_ij' and 'C_j' */
    ss_index** psdColIdx = NULL;
    ss_index** symMatIdx = NULL;
    int matLen = 0;
    int* matPos = NULL;
    double* matElem = NULL;
    double *lpObj = NULL;
    
    /* Count number of slack columns */
    for (int i = 0; i < prob->nPSDConstr; ++i)
    {
      ss_psdmatrix* psdconstr = prob->psdconstr[i];
      int hasLower = psdconstr->dRowLower > -prob->param.dInfBound;
      int hasUpper = psdconstr->dRowUpper < +prob->param.dInfBound;
      if ((hasLower && !hasUpper) || (!hasLower && hasUpper))
      {
        nSlackCol++;
      }
    }
    for (int i = 0; i < prob->nRow; ++i)
    {
      int hasLower = prob->rowLower[i] > -prob->param.dInfBound;
      int hasUpper = prob->rowUpper[i] < +prob->param.dInfBound;
      if ((hasLower && !hasUpper) || (!hasLower && hasUpper))
      {
        nSlackCol++;
      }
    }

    /* Count number of nonzeros in linear objective */
    for (int i = 0; i < prob->nCol; ++i)
    {
      if (fabs(prob->colCost[i]) > prob->param.dMatrixTol)
      {
        ++nObjElem;
      }
    }

    /* Count number of matrix elements */
    if (prob->attr.iHasPSDObj)
    {
      ss_psdmatrix* psdobj = prob->psdobj;
      for (int i = 0; i < psdobj->nPSDMatCnt; ++i)
      {
        ss_qmatrix* symmat = prob->symmat[psdobj->symMatIdx[i]];
        matLen += symmat->nQMatCnt;
      }
    }
    for (int i = 0; i < prob->nPSDConstr; ++i)
    {
      ss_psdmatrix* psdconstr = prob->psdconstr[i];
      for (int j = 0; j < psdconstr->nPSDMatCnt; ++j)
      {
        ss_qmatrix* symmat = prob->symmat[psdconstr->symMatIdx[j]];
        matLen += symmat->nQMatCnt;
      }
    }

    /* Dimension statistics */
    nDualCol = prob->nCol + nSlackCol;
    nDualRow = prob->nPSDConstr + prob->nRow;
    nDualElem += nObjElem;
    for (int i = 0; i < prob->nPSDConstr; ++i)
    {
      ss_psdmatrix* psdconstr = prob->psdconstr[i];
      nDualElem += psdconstr->nRowMatCnt;
    }
    nDualElem += prob->nElem;
    nDualElem += nSlackCol;

    /* Allocate memory */
    SS_INIT(rowRHS, nDualRow);

    SS_INIT(rowMatBeg, nDualRow + 2);
    SS_INIT(rowMatIdx, nDualElem);
    SS_INIT(rowMatElem, nDualElem);
    SS_INIT_ZERO(lpObj, nDualCol);

    SS_INIT_ZERO(psdColIdx, prob->nPSD);
    SS_INIT_ZERO(symMatIdx, prob->nPSD);
    SS_INIT(matPos, matLen);
    SS_INIT(matElem, matLen);

    /**
     * Build 'b'
     */
    for (int i = 0; i < prob->nPSDConstr; ++i)
    {
      ss_psdmatrix* psdconstr = prob->psdconstr[i];
      int hasLower = psdconstr->dRowLower > -prob->param.dInfBound;
      int hasUpper = psdconstr->dRowUpper < +prob->param.dInfBound;
      if (hasLower && !hasUpper)
      {
        rowRHS[i] = psdconstr->dRowLower;
      }
      else
      {
        rowRHS[i] = psdconstr->dRowUpper;
      }
    }
    for (int i = 0; i < prob->nRow; ++i)
    {
      int hasLower = prob->rowLower[i] > -prob->param.dInfBound;
      int hasUpper = prob->rowUpper[i] < +prob->param.dInfBound;
      if (hasLower && !hasUpper)
      {
        rowRHS[prob->nPSDConstr + i] = prob->rowLower[i];
      }
      else
      {
        rowRHS[prob->nPSDConstr + i] = prob->rowUpper[i];
      }
    }

    /**
     * Build 'c' + 'A'
     */
    int iRowIdx = 0;
    int iElemIdx = 0;
    int iSlackIdx = 0;

    /* The linear objective */
    rowMatBeg[iRowIdx] = 0;
    iRowIdx++;

    if (prob->attr.iObjSense == COPT_MAXIMIZE)
    {
      for (int i = 0; i < prob->nCol; ++i)
      {
        /* Nonzeros: linear objective */
        if (fabs(prob->colCost[i]) > prob->param.dMatrixTol)
        {
          rowMatIdx[iElemIdx] = i;
          rowMatElem[iElemIdx] = -prob->colCost[i];
          iElemIdx++;
        }
      }
    }
    else
    {
      for (int i = 0; i < prob->nCol; ++i)
      {
        /* Nonzeros: linear objective */
        if (fabs(prob->colCost[i]) > prob->param.dMatrixTol)
        {
          rowMatIdx[iElemIdx] = i;
          rowMatElem[iElemIdx] = prob->colCost[i];
          iElemIdx++;
        }
      }
    }

    /* The PSD part */
    for (int i = 0; i < prob->nPSDConstr; ++i)
    {
      ss_psdmatrix* psdconstr = prob->psdconstr[i];
      int hasLower = psdconstr->dRowLower > -prob->param.dInfBound;
      int hasUpper = psdconstr->dRowUpper < +prob->param.dInfBound;

      /* Begin pointer */
      rowMatBeg[iRowIdx] = iElemIdx;
      iRowIdx++;

      /* Nonzeros: linear part */
      for (int j = 0; j < psdconstr->nRowMatCnt; ++j)
      {
        rowMatIdx[iElemIdx] = psdconstr->rowMatIdx[j];
        rowMatElem[iElemIdx] = psdconstr->rowMatElem[j];
        iElemIdx++;
      }
      /* Nonzeros: slack column */
      if (hasLower && !hasUpper)
      {
        rowMatIdx[iElemIdx] = prob->nCol + iSlackIdx;
        rowMatElem[iElemIdx] = -1.0;
        iSlackIdx++;
        iElemIdx++;
      }
      else if (!hasLower && hasUpper)
      {
        rowMatIdx[iElemIdx] = prob->nCol + iSlackIdx;
        rowMatElem[iElemIdx] = 1.0;
        iSlackIdx++;
        iElemIdx++;
      }
    }

    /* The linear part */
    for (int i = 0; i < prob->nRow; ++i)
    {
      int hasLower = prob->rowLower[i] > -prob->param.dInfBound;
      int hasUpper = prob->rowUpper[i] < +prob->param.dInfBound;

      /* Begin pointer */
      rowMatBeg[iRowIdx] = iElemIdx;
      iRowIdx++;

      /* Nonzeros: linear part */
      int iRowBeg = prob->rowMatBeg[i];
      int iRowEnd = prob->rowMatCnt[i] + iRowBeg;
      for (int j = iRowBeg; j < iRowEnd; ++j)
      {
        rowMatIdx[iElemIdx] = prob->rowMatIdx[j];
        rowMatElem[iElemIdx] = prob->rowMatElem[j];
        iElemIdx++;
      }
      /* Nonzeros: slack column */
      if (hasLower && !hasUpper)
      {
        rowMatIdx[iElemIdx] = prob->nCol + iSlackIdx;
        rowMatElem[iElemIdx] = -1.0;
        iSlackIdx++;
        iElemIdx++;
      }
      else if (!hasLower && hasUpper)
      {
        rowMatIdx[iElemIdx] = prob->nCol + iSlackIdx;
        rowMatElem[iElemIdx] = 1.0;
        iSlackIdx++;
        iElemIdx++;
      }
    }
    rowMatBeg[iRowIdx] = iElemIdx;
    SS_ASSERT(iSlackIdx == nSlackCol);

    /**
     * Build 'psdColIdx' and 'psdMatIdx' in dual problem
     */

    /* Create index structures */
    for (int i = 0; i < prob->nPSD; ++i)
    {
      SS_CALL(ss_Index_Create(&psdColIdx[i], 0));
      SS_CALL(ss_Index_Create(&symMatIdx[i], 0));
    }

    /* The PSD objective */
    if (prob->attr.iHasPSDObj)
    {
      ss_psdmatrix* psdobj = prob->psdobj;
      for (int i = 0; i < psdobj->nPSDMatCnt; ++i)
      {
        int iPsd = psdobj->psdColIdx[i];
        SS_CALL(ss_Index_Add(psdColIdx[iPsd], 0));
        SS_CALL(ss_Index_Add(symMatIdx[iPsd], psdobj->symMatIdx[i]));
      }
    }

    /* The PSD constraints */
    for (int i = 0; i < prob->nPSDConstr; ++i)
    {
      ss_psdmatrix* psdconstr = prob->psdconstr[i];
      for (int j = 0; j < psdconstr->nPSDMatCnt; ++j)
      {
        int iPsd = psdconstr->psdColIdx[j];
        SS_CALL(ss_Index_Add(psdColIdx[iPsd], i + 1));
        SS_CALL(ss_Index_Add(symMatIdx[iPsd], psdconstr->symMatIdx[j]));
      }
    }
    
    /* Create HDSDP solver */
    DSDP_INT nNzs = nDualElem + matLen;
    DSDP_INT nBlocks = prob->nPSD;
    DSDP_INT nSDPVars = 0;
    DSDP_INT nConstrs = prob->nPSDConstr;
    DSDP_INT nLpVars = prob->nCol;
    
    for (int i = 0; i < nBlocks; ++i)
    {
        nSDPVars += prob->psdDim[i];
    }
    
    HDSDP_CALL(HDSDPCreate(&hdsdp, NULL));
    HDSDP_CALL(DSDPSetDim(hdsdp, nSDPVars, nBlocks, nConstrs, nLpVars, &nNzs));
    
    /* Set Conic data */
    for (int i = 0; i < nBlocks; ++i)
    {
        int psdDim = prob->psdDim[i];
        int flatPsdDim = (int) (psdDim * (psdDim + 1) / 2);
        int nIndex = ss_Index_Get_Num(psdColIdx[i]);
        int* colIndex = ss_Index_Get_Index(psdColIdx[i]);
        int* matIndex = ss_Index_Get_Index(symMatIdx[i]);
        
        /* CSparse entry routine builds the block data */
        sdpBlkData = dcs_spalloc(flatPsdDim, nConstrs + 1, 1000, TRUE, TRUE);
        
        if (sdpBlkData == NULL)
        {
            retcode = RETCODE_ERROR;
            goto exit_cleanup;
        }
        
        for (int j = 0; j < nIndex; ++j)
        {
            ss_qmatrix* symmat = prob->symmat[matIndex[j]];
            
            /* C_j appears in the last column in HDSDP structure */
            if (colIndex[j] == 0)
            {
                /* Flip the sign for 'C_j' when maximization */
                if (prob->attr.iObjSense == COPT_MAXIMIZE)
                {
                    for (int k = 0; k < symmat->nQMatCnt; ++k)
                    {
                      int iRow = symmat->qMatRow[k];
                      int iCol = symmat->qMatCol[k];
                      int iFlat = ((2 * psdDim - iCol - 1) * iCol / 2) + iRow;
                      CS_CALL(ss_HDSDPDataInsert(iFlat, nConstrs,
                                                 -symmat->qMatElem[k], sdpBlkData));
                    }
                }
                else
                {
                    for (int k = 0; k < symmat->nQMatCnt; ++k)
                    {
                      int iRow = symmat->qMatRow[k];
                      int iCol = symmat->qMatCol[k];
                      int iFlat = ((2 * psdDim - iCol - 1) * iCol / 2) + iRow;
                      CS_CALL(ss_HDSDPDataInsert(iFlat, nConstrs,
                                                 symmat->qMatElem[k], sdpBlkData));
                    }
                }
            }
            else
            {
                /* Constraint data */
                for (int k = 0; k < symmat->nQMatCnt; ++k)
                {
                  int iRow = symmat->qMatRow[k];
                  int iCol = symmat->qMatCol[k];
                  int iFlat = ((2 * psdDim - iCol - 1) * iCol / 2) + iRow;
                  CS_CALL(ss_HDSDPDataInsert(iFlat, colIndex[j] - 1,
                                             symmat->qMatElem[k], sdpBlkData));
                }
            }
        }

        /* Compress data and set it */
        cscSDPBlkData = dcs_compress(sdpBlkData);
        dcs_spfree(sdpBlkData); sdpBlkData = NULL;
        HDSDP_CALL(DSDPSetSDPConeData(hdsdp, i, psdDim, cscSDPBlkData->p,
	    			      cscSDPBlkData->i, cscSDPBlkData->x));
    	/* Free the current block */
        dcs_spfree(cscSDPBlkData); cscSDPBlkData = NULL;
    }
    
    /* LP Cone data */
    /* Build c and skip the first column */
    for (int i = 0; i < rowMatBeg[1]; ++i)
    {
        lpObj[rowMatIdx[i]] = rowMatElem[i];
    }
    
    for (int i = 1; i < nDualRow + 2; ++i)
    {
        rowMatBeg[i] -= nObjElem;
    }
    
    if (nDualCol > 0) 
    {
	HDSDP_CALL(DSDPSetLPData(hdsdp, nDualCol, rowMatBeg + 1,
                             rowMatIdx + nObjElem, rowMatElem + nObjElem,
                             lpObj));
    }
    
    /* Dual objective data */
    HDSDP_CALL(DSDPSetObj(hdsdp, rowRHS));
    
    /* Call HDSDP optimization routine */
    HDSDP_CALL(DSDPOptimize(hdsdp));
    
    /* TODO: Status conversion and solution extraction */
    prob->attr.iHasLpSol = 0;
    
    
    /* Get DIMACS errors */
    double dimacs[6] = {0.0};
    
    /* TODO: Get DIMACS Errors */
    
    /* Pass back dimacs numbers */
    prob->attr.dDimacPInf = dimacs[0];
    prob->attr.dDimacDInf = dimacs[1];
    prob->attr.dDimacMinPEV = dimacs[2];
    prob->attr.dDimacMinDEV = dimacs[3];
    prob->attr.dDimacRelGap = dimacs[4];
    prob->attr.dDimacRelXDotS = dimacs[5];
    
    /* Get Objective */
    
    
exit_cleanup:

    if (sdpBlkData) {
        dcs_spfree(sdpBlkData);
        sdpBlkData = NULL;
    }
    
    if (cscSDPBlkData) {
        dcs_spfree(cscSDPBlkData);
        cscSDPBlkData = NULL;
    }
    
    HDSDPDestroy(hdsdp);
    
    /* Free memory */
    SS_FREE(rowRHS);
    SS_FREE(lpObj);

    SS_FREE(rowMatBeg);
    SS_FREE(rowMatIdx);
    SS_FREE(rowMatElem);

    for (int i = 0; i < prob->nPSD; ++i)
    {
      ss_Index_Delete(&psdColIdx[i]);
    }
    SS_FREE(psdColIdx);
    for (int i = 0; i < prob->nPSD; ++i)
    {
      ss_Index_Delete(&symMatIdx[i]);
    }
    SS_FREE(symMatIdx);
    SS_FREE(matPos);
    SS_FREE(matElem);
    
    return retcode;
}

