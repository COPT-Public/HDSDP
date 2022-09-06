#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "rot3solve.h"

#define check_code(code) if (retcode != DSDP_RETCODE_OK) \
                            { printf("Optimization failed on line %d", __LINE__); \
                              goto cleanup;}

#ifdef ROT_35
static DSDP_INT Asp[37] = {
    0, 1, 2, 4, 5, 6, 7, 9, 11, 12, 14, 16, 18,
    19, 20, 21, 22, 24, 26, 27, 29, 32, 34,
    36, 38, 39, 41, 43, 45, 47, 49, 51, 52,
    53, 54, 55, 55
};

static DSDP_INT Asi[55] = {
    54, 53, 51, 52, 50, 49, 48, 44, 47, 43,
    46, 42, 39, 45, 38, 41, 37, 40, 36, 35,
    34, 33, 26, 32, 25, 31, 24, 18, 30, 17,
    23, 29, 16, 22, 15, 28, 14, 21, 13, 9,
    27, 8, 20, 7, 19, 6, 12, 5, 11, 4, 10,
    3, 2, 1, 0
};

static double Asx[55] = {
    -1.00, -1.00, -1.00, -1.00, -1.00, -1.00, -1.00, -1.00,
    -1.00, -1.00, -1.00, -1.00, -1.00, -1.00, -1.00, -1.00,
    -1.00, -1.00, -1.00, -1.00, -1.00, -1.00, -1.00, -1.00,
    -1.00, -1.00, -1.00, -1.00, -1.00, -1.00, -1.00, -1.00,
    -1.00, -1.00, -1.00, -1.00, -1.00, -1.00, -1.00, -1.00,
    -1.00, -1.00, -1.00, -1.00, -1.00, -1.00, -1.00, -1.00,
    -1.00, -1.00, -1.00, -1.00, -1.00, -1.00, -1.00
};

static DSDP_INT Alp[36] = {
    0,
    1, 1,
    2, 2,
    3, 3, 3, 3, 3,
    4, 4,
    5, 5, 5,
    6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6,
    7, 7,
    8, 8, 8,
    9, 9, 9, 9,
    10
};

DSDP_INT Ali[MATRIX_DIM] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
double   Alx[MATRIX_DIM] = {-1.0, -2.0, -1.0, -2.0, -2.0, -1.0, -2.0, -2.0, -2.0, -1.0};
double lpObj[1] = {-1.0};
double btmp[NUM_CONSTR] = {0.0};
double sol[NUM_CONSTR]  = {0.0};
#else

static DSDP_INT Asp[17] = {
    0, 1, 2, 4, 5, 6, 7, 9, 11, 12, 14, 16, 18,
    19, 20, 21, 21
};

static DSDP_INT Asi[21] = {
    20,
    19, 17,
    18,
    16,
    15, 14,
    13, 10,
    12, 9,
    8,
    11, 5,
    7, 4,
    6, 3,
    2,
    1,
    0
};

static double Asx[21] = {
    -1.00, -1.00, -1.00, -1.00, -1.00, -1.00, -1.00, -1.00,
    -1.00, -1.00, -1.00, -1.00, -1.00, -1.00, -1.00, -1.00,
    -1.00, -1.00, -1.00, -1.00, -1.00
};

static DSDP_INT Alp[16] = {
    0,
    1, 1,
    2, 2,
    3, 3, 3, 3, 3,
    4, 4,
    5, 5, 5,
    6
};

DSDP_INT Ali[MATRIX_DIM] = {0, 0, 0, 0, 0, 0};
double   Alx[MATRIX_DIM] = {-1.0, -2.0, -1.0, -2.0, -2.0, -1.0};
double lpObj[1] = {-1.0};
double btmp[NUM_CONSTR] = {0.0};
double sol[NUM_CONSTR]  = {0.0};
#endif

static DSDP_INT extractName( char *path, char *file ) {
    
    DSDP_INT retcode = DSDP_RETCODE_OK;
    char csvsuffix[] = ".csv", *where; strcpy(file, path);
    where = strstr(file, csvsuffix);
    if (!where) {
        printf("| Invalid SDPA file name. \n");
        retcode = DSDP_RETCODE_FAILED; return retcode;
    }
    memset(where, 0, sizeof(char)); return retcode;
}

static DSDP_INT readBvec( char *fname, double *b ) {
    
    FILE *fp = fopen(fname, "r");
    
    if (!fp) {
        printf("| Failed to open %s. \n", fname);
        return DSDP_RETCODE_FAILED;
    }
    
    char line[2048]; fgets(line, 2048, fp);
    const char *tok; DSDP_INT i = 0;
    
    for (tok = strtok( line, ","); tok != NULL; tok = strtok(NULL, ",")) {
        if ( i >= NUM_CONSTR ) {
            printf("More than %d entries found", NUM_CONSTR);
            return DSDP_RETCODE_FAILED;
        }
        b[i] = atof(tok); ++i;
    }
    
    return DSDP_RETCODE_OK;
}

static void printDual( DSDP_INT n, double *y ) {
    
    for (DSDP_INT i = 0; i < n; ++i) {
        printf("y[%2d] = %+10.6e, ", i + 1, y[i]);
        if ( i % 5 == 4) {
            printf("\n");
        }
    }
    
    return;
}

static DSDP_INT exportDual( DSDP_INT n, double *y, char *fname ) {
    
    char filename[100] = "", csv[] = "-out.csv", tmp[] = "rot_y-out";
    
    if (!fname) {
        strcat(filename, tmp);
    } else {
        strcat(filename, fname);
    }
    strcat(filename, csv);
    
    FILE *fp = fopen(filename, "w+");
    
    if (!fp) {
        printf("| Failed to open %s to write. \n", filename);
        return DSDP_RETCODE_FAILED;
    }
    
    for (DSDP_INT i = 0; i < n - 1; ++i) {
        fprintf(fp, "%10.10e,", y[i]);
    }
    
    fprintf(fp, "%10.10e", y[n - 1]); fclose(fp);
    printf("\nExported solution to %s. \n", filename);
    return DSDP_RETCODE_OK;
}

extern DSDP_INT solveRot( double *b, char *fname ) {
    
    DSDP_INT retcode = DSDP_RETCODE_OK;
    HSDSolver *dsdp = NULL;
    
    retcode = DSDPCreate(&dsdp, NULL);
    check_code(retcode);
    
    retcode = DSDPSetDim(dsdp, MATRIX_DIM, 1, NUM_CONSTR, NUM_LPVAR, NULL);
    check_code(retcode);
    
    retcode = DSDPSetSDPConeData(dsdp, 0, MATRIX_DIM, Asp, Asi, Asx);
    check_code(retcode);
    
    retcode = DSDPSetLPData(dsdp, NUM_LPVAR, Alp, Ali, Alx, lpObj);
    check_code(retcode);
    
    retcode = DSDPSetObj(dsdp, b);
    check_code(retcode);
     
    DSDPSetDblParam(dsdp, DBL_PARAM_REL_OPTTOL, 1e-12);
    DSDPSetDblParam(dsdp, DBL_PARAM_ABS_OPTTOL, 1e-12);
    DSDPSetDblParam(dsdp, DBL_PARAM_PRLX_PENTALTY, 1e+04);
    
    retcode = DSDPOptimize(dsdp);
    check_code(retcode);
    
    retcode = DSDPGetDual(dsdp, sol, NULL);
    exportDual(NUM_CONSTR, sol, fname);
    printDual(NUM_CONSTR, sol);
    check_code(retcode);
    
cleanup:
    DSDPDestroy(dsdp); DSDP_FREE(dsdp);
    return retcode;
}

extern DSDP_INT solveRotfromFile( int argc, char **argv ) {
    
    DSDP_INT retcode = DSDP_RETCODE_OK;
    char thisline[100], export[100];
    
    if (argc < 2) {
        DSDPPrintVersion(); return retcode;
    }
    
    strncpy(thisline, argv[1], 90);
    retcode = extractName(thisline, export);
    
    if (retcode != DSDP_RETCODE_OK) {
        printf("| Failed to extract .csv filename. \n");
        return retcode;
    }
    
    retcode = readBvec(thisline, btmp);
    
    if (retcode != DSDP_RETCODE_OK) {
        return retcode;
    }
    
    retcode = solveRot(btmp, export);
    
    return retcode;
}
