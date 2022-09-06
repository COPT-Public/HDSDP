cmake_minimum_required(VERSION 3.21)

set(CMAKE_OSX_ARCHITECTURES x86_64)
project(hdsdplinsys C)
set(CMAKE_C_STANDARD 11)

# Configure MKL paths
set(ENV{MKL_LIB_PATH} /opt/intel/oneapi/mkl/2021.2.0/lib)
set(ENV{MKL_OMP_PATH} /opt/intel/oneapi/compiler/2021.2.0/mac/compiler/lib)
set(ENV{HDSDP_HOME} /Users/gaowenzhi/Desktop/gwz/DSDP/DSDP6)

# Include project directory
include_directories($ENV{HDSDP_HOME}/externals)
include_directories($ENV{HDSDP_HOME}/include)
include_directories($ENV{HDSDP_HOME}/interface)
include_directories($ENV{HDSDP_HOME}/src/linsys)

add_library(hdsdplinsys SHARED
        $ENV{HDSDP_HOME}/externals/config.c
        $ENV{HDSDP_HOME}/externals/config.h
        $ENV{HDSDP_HOME}/externals/cs.c
        $ENV{HDSDP_HOME}/externals/cs.h
        $ENV{HDSDP_HOME}/include/dsdplapack.h
        $ENV{HDSDP_HOME}/interface/dsdphsd.h
        $ENV{HDSDP_HOME}/src/linsys/denseopts.c
        $ENV{HDSDP_HOME}/src/linsys/denseopts.h
        $ENV{HDSDP_HOME}/src/linsys/adpcg.c
        $ENV{HDSDP_HOME}/src/linsys/adpcg.h
        $ENV{HDSDP_HOME}/src/linsys/dsdpfeast.h
        $ENV{HDSDP_HOME}/src/linsys/lanczos.c
        $ENV{HDSDP_HOME}/src/linsys/lanczos.h
        $ENV{HDSDP_HOME}/src/linsys/pardiso.h
        $ENV{HDSDP_HOME}/src/linsys/rankkopts.c
        $ENV{HDSDP_HOME}/src/linsys/rankkopts.h
        $ENV{HDSDP_HOME}/src/linsys/rank1opts.c
        $ENV{HDSDP_HOME}/src/linsys/rank1opts.h
        $ENV{HDSDP_HOME}/src/linsys/schurmat.c
        $ENV{HDSDP_HOME}/src/linsys/schurmat.h
        $ENV{HDSDP_HOME}/src/linsys/sparseopts.c
        $ENV{HDSDP_HOME}/src/linsys/sparseopts.h
        $ENV{HDSDP_HOME}/src/linsys/speigs.c
        $ENV{HDSDP_HOME}/src/linsys/speigs.h
        $ENV{HDSDP_HOME}/src/linsys/structs.h
        $ENV{HDSDP_HOME}/src/linsys/vec.c
        $ENV{HDSDP_HOME}/src/linsys/vec.h
        $ENV{HDSDP_HOME}/src/linsys/linsysconfig.c
        $ENV{HDSDP_HOME}/src/linsys/linsysconfig.h)

target_link_libraries(hdsdplinsys $ENV{MKL_LIB_PATH}/libmkl_core.a)
target_link_libraries(hdsdplinsys $ENV{MKL_LIB_PATH}/libmkl_intel_lp64.a)
target_link_libraries(hdsdplinsys $ENV{MKL_LIB_PATH}/libmkl_intel_thread.a)
target_link_libraries(hdsdplinsys $ENV{MKL_OMP_PATH}/libiomp5.dylib)
