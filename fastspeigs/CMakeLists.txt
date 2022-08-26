cmake_minimum_required(VERSION 3.21)

set(CMAKE_OSX_ARCHITECTURES x86_64)
project(SPEIGS C)
set(CMAKE_C_STANDARD 11)

# Configure MKL paths
set(ENV{MKL_LIB_PATH} /opt/intel/oneapi/mkl/2021.2.0/lib)
set(ENV{MKL_OMP_PATH} /opt/intel/oneapi/compiler/2021.2.0/mac/compiler/lib)

# Include project directory
include_directories(src)

add_library(speigs
        src/spinfo.h
        src/speigs.c
        src/speigs.h)

target_link_libraries(speigs $ENV{MKL_LIB_PATH}/libmkl_core.a)
target_link_libraries(speigs $ENV{MKL_LIB_PATH}/libmkl_intel_lp64.a)
target_link_libraries(speigs $ENV{MKL_LIB_PATH}/libmkl_intel_thread.a)
target_link_libraries(speigs $ENV{MKL_OMP_PATH}/libiomp5.dylib)

add_executable(speigs_exp
        src/spinfo.h
        src/speigs.c
        src/speigs.h
        src/example.c
        src/example.h)

target_link_libraries(speigs_exp $ENV{MKL_LIB_PATH}/libmkl_core.a)
target_link_libraries(speigs_exp $ENV{MKL_LIB_PATH}/libmkl_intel_lp64.a)
target_link_libraries(speigs_exp $ENV{MKL_LIB_PATH}/libmkl_intel_thread.a)
target_link_libraries(speigs_exp $ENV{MKL_OMP_PATH}/libiomp5.dylib)