#include <stdio.h>
#include <stdlib.h>

int test_file_io( char *fname );
int test_solver( char *fname );

int main(int argc, const char * argv[]) {
 
    if ( argc > 1 ) {
        return test_solver(argv[1]);
    } else {
        char *fname = "/Users/gaowenzhi/Desktop/gwz/benchmark/sdplib/hamming_8_3_4.dat-s";
        return test_solver(fname);
    }
}
