# SPEIGS
#### An efficient preprocessor for VERY SParse EIgen-decomposition problem.

SPEIGS implements an efficient eigen-decomposition pre-processor for *Extremely* sparse matrices (arising from SDPs). Given a real symmetric matrix, it computes the full eigen-decomposition $A = V\Lambda V^T$ . It internally uses Lapack *deyevr* as the subroutine and detects the possible special structures within the matrix to factorize. It is **recommended** to try the package if $A$ satisfies either of the conditions below

1. $A$ is EXTREMELY sparse. e.g.., there may exist empty row/columns
2. $A$ is possibly low-rank. e.g. Some matrices may be approximately rank-one.

If neither of the two cases is satisfied, we recommend the use of other state-of-the-art sparse eigen-decomposition libraries such as ARPACK and MKL FEAST. 

#### Origin

SPEIGS originates from **DSDP5.8** (https://www.mcs.anl.gov/hs/software/DSDP/), a semi-definite programming solver by Steve Benson and is formalized as a library in its successor **HDSDP** (https://github.com/COPT-Public/HDSDP).

#### Current release

The current version of SPEIGS is 1.0.0 and can be called from C and MATLAB Mex interface.

#### Installation

SPEIGS is built via CMAKE system and is linked to MKL library for the implementation of *dsyevr*. The user can also switch to other implementations by modifying the paths and linked libraries from the **CMakeLists.txt**.

```cmake
# Option 1. Link with MKL
set(ENV{MKL_LIB_PATH} YOUR_MKL_PATH)
set(ENV{MKL_OMP_PATH} YOUR_OMP_PATH)
target_link_libraries(speigs $ENV{MKL_LIB_PATH}/libmkl_core.a)
target_link_libraries(speigs $ENV{MKL_LIB_PATH}/libmkl_intel_lp64.a)
target_link_libraries(speigs $ENV{MKL_LIB_PATH}/libmkl_intel_thread.a)
target_link_libraries(speigs $ENV{MKL_OMP_PATH}/libiomp5.dylib)

# Option 2. Modify the paths to link with other implementations
# set(ENV{LAPACK_BLAS_PATH} YOUR_LAPACK_BLAS_PATH)
# target_link_libraries(speigs $ENV{LAPACK_BLAS_PATH}/liblapack.a)
# target_link_libraries(speigs $ENV{LAPACK_BLAS_PATH}/liblas.a)
```

After configuring the installation paths.  Users can execute

```bash
mkdir build
cmake ..
make
```

in the command line and build the SPEIG library.

#### Documentation

The interface of SPDEIGS is well-documented using **doxygen** system and the users can run

```bash
cd doc
doxygen .
```

in the command line to generate HTML or LaTex documents to the interface.

#### Examples

The examples for SPEIGS are available at 

```Ebash
src/example.h
src/example.c
matlab/mex_speigs.c
```

and in a word, SPEIGS runs in a two-phase fashion

- An analysis phase that detects special structure within the matrix 
- A factorization phase that extracts the decomposition exploiting the structures from analysis phase 

The users can flexibly decide whether to use SPEIGS to factorize based on the result from the analysis phase.

#### Use in MATLAB

SPEIGS is callable from MATLAB by the MEX interface. Users can either build the mexfile using **CMakeLists.txt** from `matlab`  directory or download the pre-built mex files from https://github.com/leavesgrp/SPEIGS/releases/tag/v1.0.0 and run

```
install_mex
```

in Matlab. On successful installation, SPEIGS can be uses as `eig` function by

```
[V, e] = speigs(tril(A), opts);
```

and `test_speigs.m` , `help speigs`  would provide help on how to use the routine.

#### Performance

Since SPEIGS serves as a pre-processor that targets special structures of matrices. If there does exist structures to exploit, SPEIG might be 1000x faster than conventional eigen solvers. SPEIG would be less efficient without structure and users can decide after the analysis phase.

#### Maintainer

SPEIGS is a by-product of the **HDSDP** solver, which is maintained by Wenzhi Gao from Shanghai University of Finance and Economics.

