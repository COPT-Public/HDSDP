/** @file hdsdp\_file\_io.h
 * HDSDP input and output for SDPA format
 */
#ifndef hdsdp_file_reader_h
#define hdsdp_file_reader_h

hdsdp_retcode HReadSDPA( char *fname, int *pnConstrs, int *pnBlks, int **pblkDims, double **prowRHS,
                         int ***pconeMatBeg, int ***pconeMatIdx, double ***pconeMatElem, int *pnCols, int *pnLPCols,
                         int **pLpMatBeg, int **pLpMatIdx, double **pLpMatElem, int *pnElems );

#endif /* hdsdp_file_reader_h */
