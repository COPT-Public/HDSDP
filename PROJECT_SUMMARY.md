# HDSDP project summary for COIN-OR

- **Purpose:** Open-source software for sparse semidefinite programming,
  comprising a C library and `sdpasolve` command-line solver.
- **Use for operations research:** Solves SDPs supplied in SDPA `.dat-s`
  format and linear/conic inputs in MPS `.mps` format; includes example
  instances and a user manual.
- **Maturity:** The project includes a C library, command-line solver,
  examples, and a user manual. The bundled SDP smoke test checks one example;
  additional solver features are under development.
- **Build environment verified:** macOS on Apple silicon, Apple Clang,
  CMake 4.4.3, BLAS/LAPACK from Accelerate. Optional Intel MKL Pardiso is
  disabled by default and is not covered by the smoke test.
- **Build and installation checks:** CTest runs `sdpasolve` on `truss1.dat-s`,
  requiring exit code zero and a primal-dual optimal status message. The
  installed executable was also checked with the same script.
- **Dependencies:** External BLAS and LAPACK; optional Intel MKL. Bundled
  SuiteSparse CSparse/CXSparse and LDL, OSQP QDLDL, and the separate SPEIGS
  component are listed in `THIRD_PARTY_NOTICES`.
- **Project page:** https://github.com/COPT-Public/HDSDP
- **Project manager:** Wenzhi Gao, `gwz@stanford.edu`.
