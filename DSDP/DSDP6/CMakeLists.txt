cmake_minimum_required(VERSION 3.21)

set(CMAKE_OSX_ARCHITECTURES x86_64)
project(hdsdp C)
set(CMAKE_C_STANDARD 11)

# Configure MKL paths
set(ENV{MKL_LIB_PATH} /opt/intel/oneapi/mkl/2021.2.0/lib)
set(ENV{MKL_OMP_PATH} /opt/intel/oneapi/compiler/2021.2.0/mac/compiler/lib)

# Include project directory
include_directories(externals)
include_directories(include)
include_directories(interface)
include_directories(src/data)
include_directories(src/dsalg)
include_directories(src/linsys)
include_directories(src/presolve)

add_executable(hdsdp
        externals/config.c
        externals/config.h
        externals/cs.c
        externals/cs.h
        include/dsdplapack.h
        interface/dsdphsd.c
        interface/dsdphsd.h
        interface/dsdplog.c
        interface/dsdplog.h
        interface/dsdpoutput.c
        interface/dsdpoutput.h
        interface/dsdpparam.c
        interface/dsdpparam.h
        interface/dsdpreadsdpa.c
        interface/dsdpreadsdpa.h
        interface/dsdpsolver.h
        interface/dsdpstats.c
        interface/dsdpstats.h
        interface/main.c
        src/data/dsdpdata.c
        src/data/dsdpdata.h
        src/data/dsdpsort.c
        src/data/dsdpsort.h
        src/dsalg/dsdpcorrector.c
        src/dsalg/dsdpcorrector.h
        src/dsalg/dsdpdinfeas.c
        src/dsalg/dsdpdinfeas.h
        src/dsalg/dsdpinitializer.c
        src/dsalg/dsdpinitializer.h
        src/dsalg/dsdppfeas.c
        src/dsalg/dsdppfeas.h
        src/dsalg/dsdppfeascheck.c
        src/dsalg/dsdppfeascheck.h
        src/dsalg/dsdpproxmeasure.c
        src/dsalg/dsdpproxmeasure.h
        src/dsalg/dsdppsol.c
        src/dsalg/dsdppsol.h
        src/dsalg/dsdputils.c
        src/dsalg/dsdputils.h
        src/dsalg/heurpool.c
        src/dsalg/heurpool.h
        src/dsalg/schurmatops.c
        src/dsalg/schurmatops.h
        src/dsalg/stepdirection.c
        src/dsalg/stepdirection.h
        src/dsalg/stepheur.c
        src/dsalg/stepheur.h
        src/dsalg/symschur.c
        src/dsalg/symschur.h
        src/linsys/denseopts.c
        src/linsys/denseopts.h
        src/linsys/adpcg.c
        src/linsys/adpcg.h
        src/linsys/dsdpfeast.h
        src/linsys/lanczos.c
        src/linsys/lanczos.h
        src/linsys/pardiso.h
        src/linsys/rankkopts.c
        src/linsys/rankkopts.h
        src/linsys/rank1opts.c
        src/linsys/rank1opts.h
        src/linsys/schurmat.c
        src/linsys/schurmat.h
        src/linsys/sparseopts.c
        src/linsys/sparseopts.h
        src/linsys/speigs.c
        src/linsys/speigs.h
        src/linsys/structs.h
        src/linsys/vec.c
        src/linsys/vec.h
        src/linsys/linsysconfig.h
        src/linsys/linsysconfig.c
        src/presolve/dsdppresolve.c
        src/presolve/dsdppresolve.h)

add_library(lhdsdp SHARED
        externals/config.c
        externals/config.h
        externals/cs.c
        externals/cs.h
        include/dsdplapack.h
        interface/dsdphsd.c
        interface/dsdphsd.h
        interface/dsdplog.c
        interface/dsdplog.h
        interface/dsdpoutput.c
        interface/dsdpoutput.h
        interface/dsdpparam.c
        interface/dsdpparam.h
        interface/dsdpreadsdpa.c
        interface/dsdpreadsdpa.h
        interface/dsdpsolver.h
        interface/dsdpstats.c
        interface/dsdpstats.h
        interface/main.c
        src/data/dsdpdata.c
        src/data/dsdpdata.h
        src/data/dsdpsort.c
        src/data/dsdpsort.h
        src/dsalg/dsdpcorrector.c
        src/dsalg/dsdpcorrector.h
        src/dsalg/dsdpdinfeas.c
        src/dsalg/dsdpdinfeas.h
        src/dsalg/dsdpinitializer.c
        src/dsalg/dsdpinitializer.h
        src/dsalg/dsdppfeas.c
        src/dsalg/dsdppfeas.h
        src/dsalg/dsdppfeascheck.c
        src/dsalg/dsdppfeascheck.h
        src/dsalg/dsdpproxmeasure.c
        src/dsalg/dsdpproxmeasure.h
        src/dsalg/dsdppsol.c
        src/dsalg/dsdppsol.h
        src/dsalg/dsdputils.c
        src/dsalg/dsdputils.h
        src/dsalg/heurpool.c
        src/dsalg/heurpool.h
        src/dsalg/schurmatops.c
        src/dsalg/schurmatops.h
        src/dsalg/stepdirection.c
        src/dsalg/stepdirection.h
        src/dsalg/stepheur.c
        src/dsalg/stepheur.h
        src/dsalg/symschur.c
        src/dsalg/symschur.h
        src/linsys/denseopts.c
        src/linsys/denseopts.h
        src/linsys/adpcg.c
        src/linsys/adpcg.h
        src/linsys/dsdpfeast.h
        src/linsys/lanczos.c
        src/linsys/lanczos.h
        src/linsys/pardiso.h
        src/linsys/rankkopts.c
        src/linsys/rankkopts.h
        src/linsys/rank1opts.c
        src/linsys/rank1opts.h
        src/linsys/schurmat.c
        src/linsys/schurmat.h
        src/linsys/sparseopts.c
        src/linsys/sparseopts.h
        src/linsys/speigs.c
        src/linsys/speigs.h
        src/linsys/structs.h
        src/linsys/vec.c
        src/linsys/vec.h
        src/linsys/linsysconfig.h
        src/linsys/linsysconfig.c
        src/presolve/dsdppresolve.c
        src/presolve/dsdppresolve.h)

target_link_libraries(hdsdp $ENV{MKL_LIB_PATH}/libmkl_core.a)
target_link_libraries(hdsdp $ENV{MKL_LIB_PATH}/libmkl_intel_lp64.a)
target_link_libraries(hdsdp $ENV{MKL_LIB_PATH}/libmkl_intel_thread.a)
target_link_libraries(hdsdp $ENV{MKL_OMP_PATH}/libiomp5.dylib)

target_link_libraries(lhdsdp $ENV{MKL_LIB_PATH}/libmkl_core.a)
target_link_libraries(lhdsdp $ENV{MKL_LIB_PATH}/libmkl_intel_lp64.a)
target_link_libraries(lhdsdp $ENV{MKL_LIB_PATH}/libmkl_intel_thread.a)
target_link_libraries(lhdsdp $ENV{MKL_OMP_PATH}/libiomp5.dylib)

