#include <string.h>
#include "dsdpreadsdpa.h"
#include "dsdphsd.h"
#include "dsdpsolver.h"
#include "dsdplog.h"
// A simple SDPA reader for HDSDP

static char etype[] = "SDPA File Reader";

#define BFSIZE 4096 // 4 * 1024

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
    
    FILE *file;
    char chartmp, thisline[BFSIZE] = "*";
    DSDP_INT i, j, ngot, blockid, constrid, m, n, line = 0, tline = 0, lpidx = -1;
    DSDP_INT nblock = 0, nvars = 0, nlpvars = 0, nconstr = 0, nnz = 0, lpexist;
    DSDP_INT *blocksizes = NULL;
    double *dObj = NULL, val = 0.0;
    cs **sdpAs = NULL, *lpA = NULL;
    double *c = NULL;
    
    file = fopen(filename, "r");
    
    printf("| Reading data from %s \n", filename);
    if (!file) {
        error_clean(etype, "Unable to open file. \n");
    }
    
    // Jump through comments
    while(!feof(file) && (thisline[0] == '*' || thisline[0] == '"')) {
        fgets(thisline, BFSIZE, file); ++line;
    }
    
    // Read nConstrs
    if (sscanf(thisline, ID, &nconstr) < 1) {
        printf("[%s]: Failed to read number of constraints "
               "from line "ID".\n", etype, line);
        error_clean(etype, "Failed to read SDPA data \n");
    }
    
    // Read nBlocks
    fgets(thisline, BFSIZE, file); ++line;
    if (sscanf(thisline, ID, &nblock) != 1) {
        printf("[%s]: Failed to read number of blocks "
               "from line "ID".\n", etype, line);
        error_clean(etype, "Failed to read SDPA data \n");
    }
    
    // Allocate blocksize vector and read block sizes
    DSDP_INT *tmpblocksize = (DSDP_INT *) calloc(nblock, sizeof(DSDP_INT)); ++line;
    
    for (i = 0; i < nblock; ++i) {
        if (fscanf(file, "{") == 1 ||
            fscanf(file, "(") == 1 ||
            fscanf(file, ",") == 1) {
            --i;
        } else if (fscanf(file, ID, &n) == 1) {
            if (n > 0) {
                nvars += n;
                tmpblocksize[i] = n;
            } else if (n < 0) {
                if (lpidx != -1) {
                    error_clean(etype, "Invalid LP data. \n");
                }
                lpidx = i;
                nlpvars = -n;
//                nvars -= n;
                tmpblocksize[i] = -n;
            } else {
                error_clean(etype, "Empty block detected. \n");
            }
        } else {
            error_clean(etype, "Failed to read blocksize vector. \n");
        }
    }
//    lpidx = -1;
//    nlpvars = 0;
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
    
    for (i = 0; i < nconstr; ++i) {
        if (fscanf(file, ",") == 1) {
            --i;
            continue;
        }
        while (fscanf(file, "%lg", &val) != 1) {
            fscanf(file, "%c", &chartmp);
            if (chartmp == '\n') {
                error_clean(etype, "Failed to read objective. \n");
            }
        }
        dObj[i] = val;
    }
    
    // Read data
    sdpAs = (cs **) calloc(nblock, sizeof(cs *));
    n = nconstr;
    
    for (i = 0; i < nblock; ++i) {
        m = nsym(blocksizes[i]);
        // Triplet matrix with value
        sdpAs[i] = cs_di_spalloc(m, n + 1, 1000, TRUE, TRUE);
    }
    
    if (lpexist) {
        lpA = cs_di_spalloc(nlpvars, nconstr, 1000, TRUE, TRUE);
        c = (double *) calloc(nlpvars, sizeof(double));
    }
    
    fgets(thisline, BFSIZE, file);
    tline = line;
    fseek(file, 0, SEEK_SET);
    line = 0;
    
    for (i = 0; i < tline; ++i) {
        chartmp = '*';
        while (chartmp != '\n') {
            fscanf(file, "%c", &chartmp);
        }
        ++line;
    }
    
    if (lpexist) {
        while (!feof(file)) {
            thisline[0] = '\0'; blockid = constrid = -1;
            i = j = -1; val = 0.0;
            fgets(thisline, BFSIZE, file); ++line;
            ngot = sscanf(thisline, "%d %d %d %d %lg",
                          &constrid, &blockid, &i, &j, &val);
            if (ngot != 5) {
                if (feof(file)) {
                    break;
                } else {
                    printf("Line %d is invalid \n", line);
                    error_clean(etype, "Failed to extract data. \n");
                }
            } else if (val != 0.0) {
                // Insert data
                if (blockid == lpidx + 1) {
                    assert( i == j );
                    if (constrid == 0) {
                        c[i - 1] = -val;
                    } else {
                        cs_di_entry(lpA, i - 1, constrid - 1, val);
                    }
                } else {
                    if (blockid > lpidx + 1) {
                        blockid -= 1;
                    }
                    n = blocksizes[blockid - 1];
                    i = (DSDP_INT) ((2 * n - i) * (i - 1) / 2) + (j - 1); // Transform upper into lower
                    if (constrid == 0) {
                        cs_di_entry(sdpAs[blockid - 1], i, nconstr, -val);
                    } else {
                        cs_di_entry(sdpAs[blockid - 1], i, constrid - 1, val);
                    }
                }
                nnz += 1;
            }
        }
    } else {
        while (!feof(file)) {
            thisline[0] = '\0'; blockid = constrid = -1;
            i = j = -1; val = 0.0;
            fgets(thisline, BFSIZE, file); ++line;
            ngot = sscanf(thisline, "%d %d %d %d %lg",
                          &constrid, &blockid, &i, &j, &val);
            if (ngot != 5) {
                if (feof(file)) {
                    break;
                } else {
                    if (ngot > 0) {
                        printf("Line %d is invalid \n", line);
                        error_clean(etype, "Failed to extract data. \n");
                    }
                }
            } else if (val != 0.0) {
                // Insert data
                n = blocksizes[blockid - 1];
                i = (DSDP_INT) ((2 * n - i) * (i - 1) / 2) + (j - 1); // Transform upper into lower
                if (constrid == 0) {
                    cs_di_entry(sdpAs[blockid - 1], i, nconstr, -val);
                } else {
                    cs_di_entry(sdpAs[blockid - 1], i, constrid - 1, val);
                }
                nnz += 1;
            }
        }
    }
    
    // Compress data
    *sdpAp = (DSDP_INT **) calloc(nblock, sizeof(DSDP_INT *));
    *sdpAi = (DSDP_INT **) calloc(nblock, sizeof(DSDP_INT *));
    *sdpAx = (double   **) calloc(nblock, sizeof(double *));
    
    cs *tmp = NULL;
    
    for (i = 0; i < nblock; ++i) {
        tmp = cs_compress(sdpAs[i]); cs_spfree(sdpAs[i]); sdpAs[i] = tmp;
        (*sdpAp)[i] = (DSDP_INT *) calloc(tmp->n + 1, sizeof(DSDP_INT));
        (*sdpAi)[i] = (DSDP_INT *) calloc(tmp->p[tmp->n], sizeof(DSDP_INT));
        (*sdpAx)[i] = (double   *) calloc(tmp->p[tmp->n], sizeof(double));
        memcpy((*sdpAp)[i], tmp->p, sizeof(DSDP_INT) * (tmp->n + 1));
        memcpy((*sdpAi)[i], tmp->i, sizeof(DSDP_INT) * tmp->p[tmp->n]);
        memcpy((*sdpAx)[i], tmp->x, sizeof(double)   * tmp->p[tmp->n]);
        cs_spfree(tmp);
    }
    
    if (lpexist) {
        tmp = cs_compress(lpA); cs_spfree(lpA);
        *lpAp = (DSDP_INT *) calloc(tmp->n + 1, sizeof(DSDP_INT));
        *lpAi = (DSDP_INT *) calloc(tmp->p[tmp->n], sizeof(DSDP_INT));
        *lpAx = (double   *) calloc(tmp->p[tmp->n], sizeof(double));
        *lpObj = (double   *) calloc(tmp->m, sizeof(double));
        memcpy(*lpAp, tmp->p, sizeof(DSDP_INT) * (tmp->n + 1));
        memcpy(*lpAi, tmp->i, sizeof(DSDP_INT) * tmp->p[tmp->n]);
        memcpy(*lpAx, tmp->x, sizeof(double) * tmp->p[tmp->n]);
        memcpy(*lpObj, c, sizeof(double) * tmp->m);
        cs_spfree(tmp);
    }
    
    *dObjVec = (double *) calloc(nconstr, sizeof(double));
    memcpy(*dObjVec, dObj, sizeof(double) * nconstr);
    
    *nBlock = nblock; 
    *nBlockVars = (DSDP_INT *) calloc(nblock, sizeof(DSDP_INT));
    memcpy(*nBlockVars, blocksizes, sizeof(DSDP_INT) * nblock);
    *nVars = nvars; *nConstr = nconstr; *nNzs = nnz; *nLPVars = nlpvars;
    
    DSDP_FREE(blocksizes);
    DSDP_FREE(dObj); DSDP_FREE(c);
    DSDP_FREE(sdpAs);
    return retcode;
    
exit_cleanup:
    
    DSDP_FREE(blocksizes);
    DSDP_FREE(dObj);
    
    for (i = 0; i < nblock; ++i) {
        cs_spfree(sdpAs[i]);
    }
    DSDP_FREE(sdpAs);
    return retcode;
}

static DSDP_INT extractSDPAfname( char *path, char *file ) {
    // Extract filename
    DSDP_INT retcode = DSDP_RETCODE_OK;
    char sdpasuffix[] = ".dat-s", *where;
    strcpy(file, path);
    where = strstr(file, sdpasuffix);
    if (!where) {
        printf("| Fatal error. Invalid SDPA file name. \n");
        retcode = DSDP_RETCODE_FAILED; return retcode;
    }
    memset(where, 0, sizeof(char));
    return retcode;
}

extern DSDP_INT DSDPAnalyzeSDPA(int argc, char *argv[]) {
    DSDP_INT retcode = DSDP_RETCODE_OK;
    
    FILE *file;
    char filename[100], thisline[100], export[100];
    if (argc < 2) {
        DSDPPrintVersion(); return retcode;
    } else {
        strncpy(thisline, argv[1], 90); file = fopen(thisline, "r");
    }
    
    retcode = extractSDPAfname(thisline, export); checkCode;
    
    Solver *hsdSolver = NULL;
    Solver **phsdSolver = &hsdSolver;
    
    retcode = DSDPCreate(phsdSolver, export);
    
    strncpy(filename, argv[1], 90);
    
    DSDP_INT **coneAp = NULL;
    DSDP_INT **coneAi = NULL;
    double   **coneAx = NULL;
    
    DSDP_INT *lpAp  = NULL;
    DSDP_INT *lpAi  = NULL;
    double   *lpAx  = NULL;
    double   *lpObj = NULL;
    
    double   *dObj    = NULL;
    DSDP_INT *blocksizes = NULL;
    DSDP_INT nConstrs = 0, nBlocks = 0, nSDPVars = 0, nLPVars = 0, nNz = 0;
    
    DSDPPrintVersion();
    printf("| Running SDPA structure analysis. \n");
    DSDPPrepareSDPData(filename, &dObj,
                       &coneAp, &coneAi, &coneAx,
                       &lpAp, &lpAi, &lpAx, &lpObj,
                       &nBlocks, &blocksizes, &nLPVars,
                       &nSDPVars, &nConstrs, &nNz);
    retcode = DSDPSetDim(hsdSolver, nSDPVars, nBlocks, nConstrs, nLPVars, &nNz);
    for (DSDP_INT i = 0; i < nBlocks; ++i) {
        retcode = DSDPSetSDPConeData(hsdSolver, i, blocksizes[i], NULL,
                                     coneAp[i], coneAi[i], coneAx[i]);
    }
    
    if (nLPVars > 0) {
        retcode = DSDPSetLPData(hsdSolver, nLPVars, lpAp, lpAi, lpAx, lpObj);
    }
    
    retcode = DSDPSetObj(hsdSolver, dObj);
    retcode = DSDPExport(hsdSolver, DSDP_EXPORT_DSYMBOLIC, export);
    dsdpshowdash();
    
exit_cleanup:
    DSDPDestroy(hsdSolver);
    // Free data
    for (DSDP_INT i = 0; i < nBlocks; ++i) {
        DSDP_FREE(coneAi[i]); DSDP_FREE(coneAp[i]); DSDP_FREE(coneAx[i]);
    }
    DSDP_FREE(dObj); DSDP_FREE(blocksizes);
    
    return retcode;
}

extern DSDP_INT DSDPSolveSDPA(int argc, char *argv[]) {
    
    DSDP_INT retcode = DSDP_RETCODE_OK;
    
    FILE *file;
    char filename[100], thisline[100], maxtime[100];
    double tmax = 0.0;
    if (argc < 2) {
        DSDPPrintVersion(); return retcode;
    } else {
        strncpy(thisline, argv[1], 90); file = fopen(thisline, "r");
        if (argc == 4 && strcmp(argv[2], "-timelimit") == 0) {
            strncpy(maxtime, argv[3], 20); tmax = atof(maxtime);
        }
    }
    
    if (tmax <= 0.0) { tmax = 15000.0; }
    
    Solver *hsdSolver = NULL;
    Solver **phsdSolver = &hsdSolver;
    
    retcode = DSDPCreate(phsdSolver, NULL);
    
    strncpy(filename, argv[1], 90);
    
    DSDP_INT **coneAp = NULL;
    DSDP_INT **coneAi = NULL;
    double   **coneAx = NULL;
    
    DSDP_INT *lpAp  = NULL;
    DSDP_INT *lpAi  = NULL;
    double   *lpAx  = NULL;
    double   *lpObj = NULL;
    
    double   *dObj    = NULL;
    DSDP_INT *blocksizes = NULL;
    DSDP_INT nConstrs = 0, nBlocks = 0, nSDPVars = 0, nLPVars = 0, nNz = 0;
    
    DSDPPrintVersion();
    double start = my_clock();
    
    DSDPPrepareSDPData(filename, &dObj,
                       &coneAp, &coneAi, &coneAx,
                       &lpAp, &lpAi, &lpAx, &lpObj,
                       &nBlocks, &blocksizes, &nLPVars,
                       &nSDPVars, &nConstrs, &nNz);
    
    retcode = DSDPSetDim(hsdSolver, nSDPVars, nBlocks, nConstrs, nLPVars, &nNz);
    
    for (DSDP_INT i = 0; i < nBlocks; ++i) {
        retcode = DSDPSetSDPConeData(hsdSolver, i, blocksizes[i], NULL,
                                     coneAp[i], coneAi[i], coneAx[i]);
    }
    
    if (nLPVars > 0) {
        retcode = DSDPSetLPData(hsdSolver, nLPVars, lpAp, lpAi, lpAx, lpObj);
    }
    
    printf("| Data read into solver."
           " Elapsed Time: %3.3f seconds. \n", my_clock() - start);
    
    retcode = DSDPSetObj(hsdSolver, dObj);
    DSDPSetDblParam(hsdSolver, DBL_PARAM_TIMELIMIT, tmax);
    retcode = DSDPOptimize(hsdSolver);
    
    // mwTrace("End Profiling. \n");
    
exit_cleanup:
    
    retcode = DSDPDestroy(hsdSolver);
    
    // Free data
    for (DSDP_INT i = 0; i < nBlocks; ++i) {
        DSDP_FREE(coneAi[i]);
        DSDP_FREE(coneAp[i]);
        DSDP_FREE(coneAx[i]);
    }
    
    DSDP_FREE(dObj);
    DSDP_FREE(blocksizes);
    
    return retcode;
}
