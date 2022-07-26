#include <stdio.h>
#include "dsdpoutput.h"

#define SUFFIX_DUALSYM "-dsym.csv"
#define SUFFIX_SOL     ".csv"

static void dumpdMatSizes( HSDSolver *dsdpSolver, FILE *output ) {
    // Export block statistics
    for (DSDP_INT i = 0; i < dsdpSolver->nBlock; ++i) {
        fprintf(output, "%d,%d,%d\n", i + 1, dsdpSolver->S[i]->dim, !dsdpSolver->S[i]->nominalsps);
    }
    fprintf(output, "\n");
}

static DSDP_INT dumpDualBlock( DSDP_INT blockid, DSDP_INT n, DSDP_INT *dSym, FILE *output ) {
    // Export a block
    DSDP_INT retcode = DSDP_RETCODE_OK;
    if (!dSym || !output) { retcode = DSDP_RETCODE_FAILED; return retcode; }
    DSDP_INT i, j;
    for (i = 0; i < n; ++i) {
        for (j = 0; j <= i; ++j) {
            if (packIdx(dSym, n, i, j)) {
                fprintf(output, "%d,%d,%d\n",
                        blockid + 1, i + 1, j + 1);
            }
        }
    }
    return retcode;
}

/* Implement the output interface for HDSDP */
extern DSDP_INT dumpDualSol( HSDSolver *dsdpSolver, char *fname ) {
    // Dump the dual solution y
    DSDP_INT retcode = DSDP_RETCODE_OK;
    if (dsdpSolver->insStatus != DSDP_STATUS_SOLVED) {
        printf("| The instance is not solved. \n");
        retcode = DSDP_RETCODE_FAILED; return retcode;
    }
    char filename[100] = "", suffix[] = SUFFIX_SOL;
    strcat(filename, fname); strcat(filename, suffix);
    FILE *fp; fp = fopen(filename, "w+");
    
    if (!fp) {
        printf("| Failed to create symbolic file. \n");
        retcode = DSDP_RETCODE_FAILED; return retcode;
    }
    // Dual solution y
    for (DSDP_INT i = 0; i < dsdpSolver->y->dim; ++i) {
        fprintf(fp, "%10.10e,", dsdpSolver->y->x[i]);
    }
    retcode = (fclose(fp) == 0) ? DSDP_RETCODE_OK : DSDP_RETCODE_FAILED;
    if (retcode == DSDP_RETCODE_OK) {
        printf("| Successfully exported dual solution y to %s. \n", filename);
    } else {
        printf("| Error in export. \n");
    }
    return retcode;
}

extern DSDP_INT dumpDualSymbolic( HSDSolver *dsdpSolver, char *fname ) {
    // Dump the dual symbolic structure
    DSDP_INT retcode = DSDP_RETCODE_OK;
    char filename[100] = "", filename2[100] = "", suffix[] = SUFFIX_DUALSYM;
    strcat(filename, fname); strcat(filename, suffix);
    strcat(filename2, fname); strcat(filename2, "-blk"); strcat(filename2, suffix);
    FILE *fp, *fp2; fp = fopen(filename, "w+"); fp2 = fopen(filename2, "w+");
    if (!fp || !fp2) { printf("| Failed to create symbolic file. \n"); }
    
    // Dual symbolic
    fprintf(fp, "%s,%s,%s\n", "block", "col", "row");
    for (DSDP_INT i = 0; i < dsdpSolver->nBlock; ++i) {
        if (dsdpSolver->S[i]->nominalsps) {
            fprintf(fp, "%d,%d,%d\n", i + 1, -1, -1); continue;
        }
        retcode = dumpDualBlock(i, dsdpSolver->sdpData[i]->dimS, dsdpSolver->symS[i], fp);
        if (retcode != DSDP_RETCODE_OK) {
            printf("| Export symbolic structure of block %d fails. \n", i + 1);
            break;
        }
    }
    // Block sizes
    fprintf(fp2, "%s,%s,%s\n", "block", "dim", "sparse");
    dumpdMatSizes(dsdpSolver, fp2);
    retcode = (fclose(fp) == 0 && fclose(fp2) == 0) ? DSDP_RETCODE_OK : DSDP_RETCODE_FAILED;
    if (retcode == DSDP_RETCODE_OK) {
        printf("| Successfully exported symbolic structure to %s. \n", filename);
        printf("| Successfully exported block structure to %s. \n", filename2);
    } else {
        printf("| Error in export. \n");
    }
    return retcode;
}
