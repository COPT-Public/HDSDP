cmake_minimum_required(VERSION 3.21)

set(CMAKE_OSX_ARCHITECTURES x86_64)
project(POTLP C)
set(CMAKE_C_STANDARD 11)

find_package(Matlab)
add_compile_definitions(MATLAB_MEX_FILE)

# Configure MKL paths
set(ENV{MKL_LIB_PATH} /opt/intel/oneapi/mkl/2021.2.0/lib)
set(ENV{MKL_OMP_PATH} /opt/intel/oneapi/compiler/2021.2.0/mac/compiler/lib)

# Include project directory
include_directories(.)
include_directories(../interface)
include_directories(../linalg)
include_directories(../include)
include_directories(../potreduce)

matlab_add_mex(NAME potlp SRC ./potlp.c ./matrix.h ./tmwtypes.h ./mex.h
			  ../interface/lp_solver.c ../interface/lp_solver.h ../interface/pot_param.h  ../interface/pot_solver.c ../interface/pot_solver.h
			  ../linalg/vec_mat.c ../linalg/vec_mat.h
			  ../include/pot_def.h 
			  ../potreduce/pot_constr_mat.c ../potreduce/pot_constr_mat.h ../potreduce/pot_lanczos.c ../potreduce/pot_lanczos.h
			  ../potreduce/pot_objfunc.c ../potreduce/pot_objfunc.h ../potreduce/pot_structs.h ../potreduce/pot_utils.c ../potreduce/pot_utils.h
			  ../potreduce/pot_vector.c ../potreduce/pot_vector.h)

target_link_libraries(potlp $ENV{MKL_LIB_PATH}/libmkl_core.a)
target_link_libraries(potlp $ENV{MKL_LIB_PATH}/libmkl_intel_lp64.a)
target_link_libraries(potlp $ENV{MKL_LIB_PATH}/libmkl_sequential.a)