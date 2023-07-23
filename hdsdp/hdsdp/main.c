#include <stdio.h>
#include <stdlib.h>

int test_file_io( char *fname );
int test_solver( char *fname );

int main(int argc, const char * argv[]) {
 
    char *fname = "/Users/gaowenzhi/Desktop/gwz/DSDP/DSDP6/matlab/benchmark/sdplib/buck3.dat-s";
    return test_solver(fname);
}
