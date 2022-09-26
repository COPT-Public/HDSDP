#include "dsdpreadsdpa.h"
#include "dsdplapack.h"
#include "dsdpcs.h"
#include "dsdplog.h"
#include "dsdpdata.h"

// A simple SDPA reader for HDSDP
static char etype[] = "SDPA File Reader";

#define BFSIZE 1024

static DSDP_INT DSDPPrepareSDPData( char     *filename,    // 'xxx.dat-s'
                                    double   **dObjVec,    // DSDP dual objective
                                    DSDP_INT ***sdpAp,     // SDP Ap data
                                    DSDP_INT ***sdpAi,     // SDP Ai data
                                    double   ***sdpAx,     // SDP Ax data
                                    DSDP_INT **lpAp,       // LP  Ap data
                                    DSDP_INT **lpAi,       // LP  Ai data
                                    double   **lpAx,       // LP  Ax data
                                    double   **lpObj,      // LP  coefficient data
                                    DSDP_INT *nBlock,      // Number of SDP blocks
                                    DSDP_INT **nBlockVars, // Number of variables in each block
                                    DSDP_INT *nLPVars,     // Number of LP variables
                                    DSDP_INT *nVars,       // Total number of variables
                                    DSDP_INT *nConstr,     // Total number of constraints
                                    DSDP_INT *nNzs ) {     // Number of nonzero elements
    
    // Adapted from readsdpa.c of DSDP5.8
    DSDP_INT retcode = DSDP_RETCODE_OK;
    FILE *file; char chartmp, thisline[BFSIZE] = "*";
    DSDP_INT i, j, k, ngot, blockid, constrid, m, n, line = 0, tline = 0, lpidx = -1, memerr = FALSE;
    DSDP_INT nblock = 0, *blocksizes = NULL, nvars = 0, nlpvars = 0, nconstr = 0, nnz = 0, lpexist;
    double *dObj = NULL, *c = NULL, val = 0.0; dcs **sdpAs = NULL, *lpA = NULL;
    file = fopen(filename, "r");
    printf("| Reading data from %s \n", filename);
    if (!file) { printf("| Failed to open file. \n"); fatal_error_msg(etype); }
    // Jump through comments
    while(!feof(file) && (thisline[0] == '*' || thisline[0] == '"')) {
        fgets(thisline, BFSIZE, file); ++line;
    }
    // Read the number of constraints
    if (sscanf(thisline, ID, &nconstr) < 1) {
        printf("| Failed to read number of constraints from line "ID".\n", line); fatal_error_msg(etype);
    }
    // Read the number of blocks
    fgets(thisline, BFSIZE, file); ++line;
    if (sscanf(thisline, ID, &nblock) != 1) {
        printf("| Failed to read number of blocks from line "ID".\n", line); fatal_error_msg(etype);
    }
    // Allocate blocksize vector and read block sizes
    DSDP_INT *tmpblocksize = (DSDP_INT *) calloc(nblock, sizeof(DSDP_INT)); ++line;
    if (!tmpblocksize) { printf("| Failed to allocate buffer for SDPA reader. \n"); fatal_error_msg(etype); }
    
    for (i = 0; i < nblock; ++i) {
        if (fscanf(file, "{") == 1 || fscanf(file, "(") == 1 || fscanf(file, ",") == 1) {
            --i;
        } else if (fscanf(file, ID, &n) == 1) {
            if (n > 0) {
                nvars += n; tmpblocksize[i] = n;
            } else if (n < 0) {
                if (lpidx != -1) {
                    printf("| Invalid LP data. Only one diagonal block is supported. \n");
                    fatal_error_msg(etype);
                }
                lpidx = i; nlpvars = -n; tmpblocksize[i] = -n;
            } else {
                printf("| Empty block detected. \n"); fatal_error_msg(etype);
            }
        } else {
            printf("| Failed to read blocksize vector. \n"); fatal_error_msg(etype);
        }
    }
    
    // Is there a diagonal block ?
    lpexist = (lpidx == -1) ? 0 : 1;
    if (lpexist) {
        nblock -= 1;
        blocksizes = (DSDP_INT *) calloc(nblock, sizeof(DSDP_INT));
        memcpy(blocksizes, tmpblocksize, sizeof(DSDP_INT) * lpidx);
        for (i = lpidx; i < nblock; ++i) {
            blocksizes[i] = tmpblocksize[i + 1];
        }
        DSDP_FREE(tmpblocksize);
    } else {
        blocksizes = tmpblocksize;
    }
    
    // Read objective
    fgets(thisline, BFSIZE, file); ++line;
    dObj = (double *) calloc(nconstr, sizeof(double));
    if (!dObj) {
        printf("| Failed to allocate space for dual objective vector. \n");
        fatal_error_msg(etype);
    }
    
    for (i = 0; i < nconstr; ++i) {
        if (fscanf(file, ",") == 1) { --i; continue; }
        while (fscanf(file, "%lg", &val) != 1) {
            fscanf(file, "%c", &chartmp);
            if (chartmp == '\n') {
                printf("| Failed to read objective. \n");
                fatal_error_msg(etype);
            }
        }
        dObj[i] = val;
    }
    
    // Read data
    sdpAs = (dcs **) calloc(nblock, sizeof(dcs *));
    n = nconstr;
    for (i = 0; i < nblock; ++i) {
        m = nsym(blocksizes[i]);
        sdpAs[i] = dcs_spalloc(m, n + 1, 1000, TRUE, TRUE); // Triplet matrix with value
        if (!sdpAs[i]) { memerr = TRUE; }
    }
    
    if (memerr) {
        printf("| Failed to allocate space for SDP dcs entry matrix. \n");
        fatal_error_msg(etype);
    }
    
    // Allocate LP data
    if (lpexist) {
        lpA = dcs_spalloc(nlpvars, nconstr, 1000, TRUE, TRUE);
        if (!lpA) {
            printf("| Failed to allocate space for LP dcs entry matrix. \n");
            fatal_error_msg(etype);
        }
        c = (double *) calloc(nlpvars, sizeof(double));
        if (!c) {
            printf("| Failed to allocate space for LP objective vector. \n");
            fatal_error_msg(etype);
        }
    }
    
    fgets(thisline, BFSIZE, file); tline = line; fseek(file, 0, SEEK_SET); line = 0;
    for (i = 0; i < tline; ++i) {
        chartmp = '*'; while (chartmp != '\n') { fscanf(file, "%c", &chartmp); } ++line;
    }
    
    if (lpexist) {
        while (!feof(file)) {
            thisline[0] = '\0'; blockid = constrid = -1; i = j = -1; val = 0.0;
            fgets(thisline, BFSIZE, file); ++line;
            ngot = sscanf(thisline, "%d %d %d %d %lg", &constrid, &blockid, &i, &j, &val);
            
            if (i > j) {
                printf("| Warning: non-upper triangular element detected. \n");
                k = i; i = j; j = k;
            }
            
            if (ngot != 5) {
                if (feof(file)) {
                    break;
                } else {
                    if (ngot > 0) {
                        printf("| Line %d is invalid \n", line); fatal_error_msg(etype);
                    }
                }
            } else if (val != 0.0) {
                // Insert data
                if (blockid == lpidx + 1) {
                    if (constrid == 0) {
                        c[i - 1] = -val;
                    } else {
                        dcs_entry(lpA, i - 1, constrid - 1, val);
                    }
                } else {
                    if (blockid > lpidx + 1) { blockid -= 1; }
                    n = blocksizes[blockid - 1];
                    i = (DSDP_INT) ((2 * n - i) * (i - 1) / 2) + (j - 1); // Transform upper into lower
                    if (constrid == 0) {
                        dcs_entry(sdpAs[blockid - 1], i, nconstr, -val);
                    } else {
                        dcs_entry(sdpAs[blockid - 1], i, constrid - 1, val);
                    }
                }
                nnz += 1;
            }
        }
    } else {
        while (!feof(file)) {
            thisline[0] = '\0'; blockid = constrid = -1; i = j = -1; val = 0.0;
            fgets(thisline, BFSIZE, file); ++line;
            ngot = sscanf(thisline, "%d %d %d %d %lg", &constrid, &blockid, &i, &j, &val);
            if (i > j) {
                printf("| Warning: non-upper triangular element detected. \n");
                k = i; i = j; j = k;
            }
            
            if (ngot != 5) {
                if (feof(file)) {
                    break;
                } else {
                    if (ngot > 0) {
                        printf("| Line %d is invalid \n", line); fatal_error_msg(etype);
                    }
                }
            } else if (val != 0.0) {
                // Insert data
                n = blocksizes[blockid - 1];
                i = (DSDP_INT) ((2 * n - i) * (i - 1) / 2) + (j - 1); // Transform upper into lower
                if (constrid == 0) {
                    dcs_entry(sdpAs[blockid - 1], i, nconstr, -val);
                } else {
                    dcs_entry(sdpAs[blockid - 1], i, constrid - 1, val);
                }
                nnz += 1;
            }
        }
    }
    
    // Compress data
    *sdpAp = (DSDP_INT **) calloc(nblock, sizeof(DSDP_INT *));
    *sdpAi = (DSDP_INT **) calloc(nblock, sizeof(DSDP_INT *));
    *sdpAx = (double   **) calloc(nblock, sizeof(double *));
    
    if (!sdpAp || !sdpAi || !sdpAx) {
        printf("| Failed to allocate space for SDP input dcs array \n");
        fatal_error_msg(etype);
    }
    
    dcs *tmp = NULL;
    
    for (i = 0; i < nblock; ++i) {
        tmp = dcs_compress(sdpAs[i]); dcs_spfree(sdpAs[i]); sdpAs[i] = tmp;
        (*sdpAp)[i] = (DSDP_INT *) calloc(tmp->n + 1, sizeof(DSDP_INT));
        (*sdpAi)[i] = (DSDP_INT *) calloc(tmp->p[tmp->n], sizeof(DSDP_INT));
        (*sdpAx)[i] = (double   *) calloc(tmp->p[tmp->n], sizeof(double));
        if (!(*sdpAp)[i] || !(*sdpAp)[i] || !(*sdpAp)[i]) {
            printf("| Failed to allocate space for SDP input dcs matrix \n");
            fatal_error_msg(etype);
        }
        memcpy((*sdpAp)[i], tmp->p, sizeof(DSDP_INT) * (tmp->n + 1));
        memcpy((*sdpAi)[i], tmp->i, sizeof(DSDP_INT) * tmp->p[tmp->n]);
        memcpy((*sdpAx)[i], tmp->x, sizeof(double)   * tmp->p[tmp->n]);
        dcs_spfree(tmp);
    }
    
    if (lpexist) {
        tmp = dcs_compress(lpA); dcs_spfree(lpA);
        *lpAp  = (DSDP_INT *) calloc(tmp->n + 1, sizeof(DSDP_INT));
        *lpAi  = (DSDP_INT *) calloc(tmp->p[tmp->n], sizeof(DSDP_INT));
        *lpAx  = (double   *) calloc(tmp->p[tmp->n], sizeof(double));
        *lpObj = (double   *) calloc(tmp->m, sizeof(double));
        if (!lpAp || !lpAi || !lpAx || !lpObj) {
            printf("| Failed to allocate space for LP input data. \n");
            fatal_error_msg(etype);
        }
        memcpy(*lpAp, tmp->p, sizeof(DSDP_INT) * (tmp->n + 1));
        memcpy(*lpAi, tmp->i, sizeof(DSDP_INT) * tmp->p[tmp->n]);
        memcpy(*lpAx, tmp->x, sizeof(double) * tmp->p[tmp->n]);
        memcpy(*lpObj, c, sizeof(double) * tmp->m);
        dcs_spfree(tmp);
    }
    
    *dObjVec = (double *) calloc(nconstr, sizeof(double));
    if (!dObjVec) {
        printf("| Failed to allocate dual objective input data. \n");
        fatal_error_msg(etype);
    }
    memcpy(*dObjVec, dObj, sizeof(double) * nconstr);
    *nBlock = nblock; *nBlockVars = (DSDP_INT *) calloc(nblock, sizeof(DSDP_INT));
    memcpy(*nBlockVars, blocksizes, sizeof(DSDP_INT) * nblock);
    *nVars = nvars; *nConstr = nconstr; *nNzs = nnz; *nLPVars = nlpvars;
    DSDP_FREE(blocksizes); DSDP_FREE(dObj); DSDP_FREE(c); DSDP_FREE(sdpAs);
    return retcode;
}

static DSDP_INT extractSDPAfname( char *path, char *file ) {
    // Extract filename
    DSDP_INT retcode = DSDP_RETCODE_OK;
    char sdpasuffix[] = ".dat-s", *where; strcpy(file, path);
    where = strstr(file, sdpasuffix);
    if (!where) {
        printf("| Invalid SDPA file name. \n");
        retcode = DSDP_RETCODE_FAILED; return retcode;
    }
    memset(where, 0, sizeof(char)); return retcode;
}

extern DSDP_INT DSDPAnalyzeSDPA(int argc, char *argv[]) {
    // Analyze an SDPA file and export the dual symbolic structure to .csv file
    DSDP_INT retcode = DSDP_RETCODE_OK;
    FILE *file; char filename[100], thisline[100], export[100];
    if (argc < 2) {
        DSDPPrintVersion(); return retcode;
    } else {
        strncpy(thisline, argv[1], 90); file = fopen(thisline, "r");
    }
    retcode = extractSDPAfname(thisline, export);
    
    if (retcode != DSDP_RETCODE_OK) {
        printf("| Failed to extract SDPA filename. \n"); return retcode;
    }
    
    HSDSolver *analyzer = NULL, **panalyzer = &analyzer;
    
    retcode = DSDPCreate(panalyzer, export);
    if (retcode != DSDP_RETCODE_OK) {
        printf("| Failed to invoke solver. \n"); return retcode;
    }
    
    strncpy(filename, argv[1], 90);
    
    // Prepare conic data
    DSDP_INT **coneAp = NULL, **coneAi = NULL, *lpAp = NULL, *lpAi = NULL;
    double   **coneAx = NULL, *lpAx = NULL, *lpObj = NULL, *dObj = NULL;
    DSDP_INT *blocksizes = NULL, nConstrs = 0, nBlocks = 0, nSDPVars = 0, nLPVars = 0, nNz = 0, i;
    
    DSDPPrintVersion();
    printf("| Running SDPA structure analysis. \n");
    retcode = DSDPPrepareSDPData(filename, &dObj, &coneAp, &coneAi, &coneAx, &lpAp, &lpAi, &lpAx, &lpObj,
                                 &nBlocks, &blocksizes, &nLPVars, &nSDPVars, &nConstrs, &nNz);
    
    if (retcode != DSDP_RETCODE_OK) {
        printf("| Failed to prepare SDP data. \n");
        goto exit_cleanup;
    }
    
    retcode = DSDPSetDim(analyzer, nSDPVars, nBlocks, nConstrs, nLPVars, &nNz);
    for (i = 0; i < nBlocks; ++i) {
        retcode = DSDPSetSDPConeData(analyzer, i, blocksizes[i],
                                     coneAp[i], coneAi[i], coneAx[i]);
        if (retcode != DSDP_RETCODE_OK) {
            break;
        }
    }
    
    if (retcode != DSDP_RETCODE_OK) {
        printf("| Failed to set conic data. \n");
        goto exit_cleanup;
    }
    
    if (nLPVars > 0) {
        retcode = DSDPSetLPData(analyzer, nLPVars, lpAp, lpAi, lpAx, lpObj);
        if (retcode != DSDP_RETCODE_OK) {
            printf("| Failed to set LP data. \n");
            goto exit_cleanup;
        }
    }
    
    retcode = DSDPSetObj(analyzer, dObj);
    if (retcode != DSDP_RETCODE_OK) {
        printf("| Failed to set dual objective. \n");
        goto exit_cleanup;
    }
    
    retcode = DSDPExport(analyzer, DSDP_EXPORT_DSYMBOLIC, export);
    if (retcode != DSDP_RETCODE_OK) {
        printf("| Failed to export analysis results. \n");
        goto exit_cleanup;
    }
    showBeautifulDashlines();
    
exit_cleanup:
    DSDPDestroy(analyzer);
    for (DSDP_INT i = 0; i < nBlocks; ++i) {
        DSDP_FREE(coneAi[i]); DSDP_FREE(coneAp[i]); DSDP_FREE(coneAx[i]);
    }
    DSDP_FREE(dObj); DSDP_FREE(blocksizes);
    return retcode;
}

extern DSDP_INT DSDPSolveSDPA(int argc, char *argv[]) {
    
    DSDP_INT retcode = DSDP_RETCODE_OK;
    FILE *file; char filename[100], thisline[100], maxtime[100], export[100];
    double tmax = 0.0;
    if (argc < 2) {
        DSDPPrintVersion(); return retcode;
    } else {
        strncpy(thisline, argv[1], 90); file = fopen(thisline, "r");
        if (argc == 4 && strcmp(argv[2], "-timelimit") == 0) {
            strncpy(maxtime, argv[3], 20); tmax = atof(maxtime);
        }
    }
    retcode = extractSDPAfname(thisline, export);
    
    if (!retcode) {
        printf("| Failed to extract SDPA filename. Using 'model' by default. \n");
        char name[] = "model"; memset(export, 0, 100 * sizeof(char));
        strcat(export, name);
    }
    
    if (tmax <= 0.0) { tmax = 15000.0; }
    
    HSDSolver *hsdSolver = NULL, **phsdSolver = &hsdSolver;
    
    retcode = DSDPCreate(phsdSolver, NULL);
    if (retcode != DSDP_RETCODE_OK) {
        printf("| Failed to invoke solver. \n");
        fatal_error;
    }
    
    strncpy(filename, argv[1], 90);
    
    // Prepare conic data
    DSDP_INT **coneAp = NULL, **coneAi = NULL, *lpAp = NULL, *lpAi = NULL;
    double   **coneAx = NULL, *lpAx = NULL, *lpObj = NULL, *dObj = NULL;
    DSDP_INT *blocksizes = NULL, nConstrs = 0, nBlocks = 0, nSDPVars = 0, nLPVars = 0, nNz = 0, i;
    
    DSDPPrintVersion();
    double start = my_clock();
    retcode = DSDPPrepareSDPData(filename, &dObj, &coneAp, &coneAi, &coneAx, &lpAp, &lpAi, &lpAx, &lpObj,
                                 &nBlocks, &blocksizes, &nLPVars, &nSDPVars, &nConstrs, &nNz);
    if (retcode != DSDP_RETCODE_OK) {
        printf("| Failed to prepare SDP data. \n");
        goto exit_cleanup;
    }
    
    retcode = DSDPSetDim(hsdSolver, nSDPVars, nBlocks, nConstrs, nLPVars, &nNz);
    
    for (i = 0; i < nBlocks; ++i) {
        retcode = DSDPSetSDPConeData(hsdSolver, i, blocksizes[i],
                                     coneAp[i], coneAi[i], coneAx[i]);
        if (retcode != DSDP_RETCODE_OK) {
            break;
        }
    }
    
    if (retcode != DSDP_RETCODE_OK) {
        printf("| Failed to set conic data. \n");
        goto exit_cleanup;
    }
    
    if (nLPVars > 0) {
        retcode = DSDPSetLPData(hsdSolver, nLPVars, lpAp, lpAi, lpAx, lpObj);
        if (retcode != DSDP_RETCODE_OK) {
            printf("| Failed to set LP data. \n");
            goto exit_cleanup;
        }
    }
    
    printf("| Data read into solver."
           " Elapsed Time: %3.3f seconds. \n", my_clock() - start);
    
    retcode = DSDPSetObj(hsdSolver, dObj);
    if (retcode != DSDP_RETCODE_OK) {
        printf("| Failed to set dual objective. \n");
        goto exit_cleanup;
    }
    
    DSDPSetDblParam(hsdSolver, DBL_PARAM_TIMELIMIT, tmax);
    retcode = DSDPOptimize(hsdSolver);
    
    if (retcode != DSDP_RETCODE_OK) {
        printf("| Optimization failed. \n");
    }
    
exit_cleanup:
    DSDPDestroy(hsdSolver);
    for (i = 0; i < nBlocks; ++i) {
        DSDP_FREE(coneAi[i]); DSDP_FREE(coneAp[i]); DSDP_FREE(coneAx[i]);
    }
    DSDP_FREE(dObj); DSDP_FREE(blocksizes);
    return retcode;
}
