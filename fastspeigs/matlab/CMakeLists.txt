cmake_minimum_required(VERSION 3.21)

set(CMAKE_OSX_ARCHITECTURES x86_64)
project(SPEIGS C)
set(CMAKE_C_STANDARD 11)

find_package(Matlab)
add_compile_definitions(MATLAB_MEX_FILE)


# Configure MKL paths
set(ENV{MKL_LIB_PATH} /opt/intel/oneapi/mkl/2021.2.0/lib)
set(ENV{MKL_OMP_PATH} /opt/intel/oneapi/compiler/2021.2.0/mac/compiler/lib)

# Include project directory
include_directories(../src)
include_directories(externals)

matlab_add_mex(NAME speigs  SRC ../src/spinfo.h ../src/speigs.c ../src/speigs.h mex_speigs.c externals/matrix.h externals/mex.h externals/tmwtypes.h)
target_link_libraries(speigs $ENV{MKL_LIB_PATH}/libmkl_core.a)
target_link_libraries(speigs $ENV{MKL_LIB_PATH}/libmkl_intel_ilp64.a)
target_link_libraries(speigs $ENV{MKL_LIB_PATH}/libmkl_intel_thread.a)
target_link_libraries(speigs $ENV{MKL_OMP_PATH}/libiomp5.dylib)