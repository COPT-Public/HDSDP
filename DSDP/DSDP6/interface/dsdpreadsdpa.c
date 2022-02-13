#include "dsdpreadsdpa.h"
#include "dsdphsd.h"
#include "dsdpsolver.h"
// A simple SDPA reader for HDSDP

static char etype[] = "SDPA File Reader";

#define BFSIZE 4096 // 4 * 1024

static DSDP_INT DSDPPrepareSDPData( char     *filename,    // 'xxx.dat-s'
                                    double   **dObjVec,    // SDP dual objective
                                    DSDP_INT ***sdpAp,     // SDP Ap data
                                    DSDP_INT ***sdpAi,     // SDP Ai data
                                    double   ***sdpAx,     // SDP Ax data
                                    DSDP_INT *nBlock,      // Number of SDP blocks
                                    DSDP_INT **nBlockVars, // Number of variables in each block
                                    DSDP_INT *nVars,       // Total number of variables
                                    DSDP_INT *nConstr,     // Total number of constraints
                                    DSDP_INT *nNzs ) {     // Number of nonzero elements
    
    // Adapted from readsdpa.c of DSDP5.8
    DSDP_INT retcode = DSDP_RETCODE_OK;
    
    FILE *file;
    char chartmp, thisline[BFSIZE] = "*";
    DSDP_INT i, j, ngot, blockid, constrid, m, n, line = 0, tline = 0;
    DSDP_INT nblock = 0, nvars = 0, nconstr = 0, nnz = 0;
    DSDP_INT *blocksizes = NULL;
    double *dObj = NULL, val = 0.0;
    cs **sdpAs = NULL;
    
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
    blocksizes = (DSDP_INT *) calloc(nblock, sizeof(DSDP_INT)); ++line;
    for (i = 0; i < nblock; ++i) {
        if (fscanf(file, "{") == 1 ||
            fscanf(file, "(") == 1 ||
            fscanf(file, ",") == 1) {
            --i;
        } else if (fscanf(file, ID, &n) == 1) {
            if (n > 0) {
                nvars += n;
                blocksizes[i] = n;
            } else if (n < 0) {
                nvars -= n;
                blocksizes[i] = - n;
            } else {
                error_clean(etype, "Empty block detected. \n");
            }
        } else {
            error_clean(etype, "Failed to read blocksize vector. \n");
        }
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
    
    while (!feof(file)) {
        thisline[0] = '\0';
        blockid = constrid = -1;
        i = j = -1;
        val = 0.0;
        
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
        } else if (val != 0) {
            // Insert data
            n = blocksizes[blockid - 1];
            // Transform upper into lower
            i = (DSDP_INT) ((2 * n - i) * (i - 1) / 2) + (j - 1);
            
            if (constrid == 0) {
                cs_di_entry(sdpAs[blockid - 1], i, nconstr, -val);
            } else {
                cs_di_entry(sdpAs[blockid - 1], i, constrid - 1, val);
            }
            nnz += 1;
        }
    }
    
    // Compress data
    *sdpAp = (DSDP_INT **) calloc(nblock, sizeof(DSDP_INT *));
    *sdpAi = (DSDP_INT **) calloc(nblock, sizeof(DSDP_INT *));
    *sdpAx = (double   **) calloc(nblock, sizeof(DSDP_INT *));
    
    cs *tmp = NULL;
    
    for (i = 0; i < nblock; ++i) {
        tmp = cs_compress(sdpAs[i]);
        cs_spfree(sdpAs[i]);
        
        sdpAs[i] = tmp;
        (*sdpAp)[i] = (DSDP_INT *) calloc(tmp->n + 1, sizeof(DSDP_INT));
        (*sdpAi)[i] = (DSDP_INT *) calloc(tmp->p[tmp->n], sizeof(DSDP_INT));
        (*sdpAx)[i] = (double   *) calloc(tmp->p[tmp->n], sizeof(double));
        
        memcpy((*sdpAp)[i], tmp->p, sizeof(DSDP_INT) * (tmp->n + 1));
        memcpy((*sdpAi)[i], tmp->i, sizeof(DSDP_INT) * tmp->p[tmp->n]);
        memcpy((*sdpAx)[i], tmp->x, sizeof(double)   * tmp->p[tmp->n]);
        
        cs_spfree(tmp);
    }
    
    *dObjVec = (double *) calloc(nconstr, sizeof(double));
    memcpy(*dObjVec, dObj, sizeof(double) * nconstr);
    
    *nBlock = nblock;
    *nBlockVars = (DSDP_INT *) calloc(nblock, sizeof(DSDP_INT));
    memcpy(*nBlockVars, blocksizes, sizeof(DSDP_INT) * nblock);
    *nVars = nvars;
    *nConstr = nconstr;
    *nNzs = nnz;
    
    DSDP_FREE(blocksizes);
    DSDP_FREE(dObj);
    
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


extern DSDP_INT DSDPSolveSDPA(int argc, char *argv[]) {
    
    DSDP_INT retcode = DSDP_RETCODE_OK;
    
    FILE *file;
    char filename[100], thisline[100];
    
    if (argc < 2) {
        DSDPPrintVersion();
        return retcode;
    } else {
        strncpy(thisline, argv[1], 90);
        file = fopen(thisline, "r");
    }
    
    Solver *hsdSolver = NULL;
    Solver **phsdSolver = &hsdSolver;
    retcode = DSDPCreate(phsdSolver);
    
    strncpy(filename, argv[1], 90);
    
    DSDP_INT **coneAp = NULL;
    DSDP_INT **coneAi = NULL;
    double   **coneAx = NULL;
    double   *dObj    = NULL;
    DSDP_INT *blocksizes = NULL;
    DSDP_INT nConstrs = 0, nBlocks = 0, nSDPVars = 0, nNz = 0;
    
    DSDPPrintVersion();
    clock_t start = clock();
    DSDPPrepareSDPData(filename, &dObj, &coneAp, &coneAi, &coneAx,
                       &nBlocks, &blocksizes, &nSDPVars, &nConstrs, &nNz);
    retcode = DSDPSetDim(hsdSolver, nSDPVars, nBlocks, nConstrs, 0, &nNz);
    
    for (DSDP_INT i = 0; i < nBlocks; ++i) {
        retcode = DSDPSetSDPConeData(hsdSolver, i, blocksizes[i], NULL,
                                     coneAp[i], coneAi[i], coneAx[i]);
    }
    printf("| Data read into solver."
           " Elapsed Time: %10.3e. \n", (double) (clock() - start) / CLOCKS_PER_SEC);
    
    retcode = DSDPSetObj(hsdSolver, dObj);
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
