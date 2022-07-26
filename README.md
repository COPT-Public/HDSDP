## HDSDP 

#### A Homogeneous Dual-scaling Interior Point Solver for Sparse Semi-definite Programming

HDSDP implements a dual-scaling interior point method for sparse Semi-definite programming problems (SDPs). The solver interits several features from  **DSDP5.8** [1] and presents some new aspects for the dual-scaling algorithm. HDSDP is written in ANSI C and is maintained by Cardinal Operations COPT development team. More features are still under active development.

#### Solver overview

HDSDP solves standard form semi-definite programs expressed by

$$\min_{\mathcal{A} \mathbf{X} = \mathbf{b}, \mathbf{X}\in \mathbb{S}_+^n} \left\langle \mathbf{C}, \mathbf{X} \right\rangle  \\
  \\ 
$$

where we do linear optimization subject to affine constraints over the cone of positive-semidefinite matrices.

#### Current release

HDSDP is now under active development and releases a pre-built binary (v0.9.2) which reads SDPA (.dat-s) format files, solves the SDP problem, and dumps the solution. Users testing solvers in Linux can download the binary from https://github.com/COPT-Public/HDSDP/releases/download/v0.9.2/hdsdp_20220714_v_0_9_2-alpha.out.gz .

In the version to come, HDSDP will support a complete user interface in C and support for MATLAB is also under development.

#### Getting started

After downloading the binary from the release site, you could execute 

```bash
gunzip hdsdp_20220714_v_0_9_2-alpha.out.gz
```

in the command line to extract the binary and then execute 

```
chmod +x hdsdp_20220714_v_0_9_2-alpha.out
```

to allow it to run on your device. If successful, running the binary directly 

```
./hdsdp_20220714_v_0_9_2-alpha.out
```

would output a header

```
--------------------------------------------------------------------------------------------------
| Checking running environment. 
| Number of threads available: xx.  Optimizing over xx threads. 
--------------------------------------------------------------------------------------------------
| Homogeneous Dual Scaling Interior Point Solver. Version 0.9.2 (Build date 7.14 2022)                                   
--------------------------------------------------------------------------------------------------
```

in the console. Then by running

```
./hdsdp_20220714_v_0_9_2-alpha.out SDPAFILE.dat-s
```

we can solve SDPs represented in standard SDPA format.  HDSDP has an internal time limit of 15000 seconds or you could optionally set a time limit by 

```
./hdsdp_20220714_v_0_9_2-alpha.out SDPAFILE.dat-s -timelimit 15000
```

If everything goes well, we would see logs like below.

```
--------------------------------------------------------------------------------------------------
| Homogeneous Dual Scaling Interior Point Solver. Version 0.9.2 (Build date 7.14 2022)                                   
--------------------------------------------------------------------------------------------------
| Reading data from theta1.dat-s 
--------------------------------------------------------------------------------------------------
| nSDPBlock: 1 | nConstrs: 104 | LP. Dim: 0 | SDP. Dim: 50 | Nonzeros: 1428 
| Data read into solver. Elapsed Time: 0.047 seconds. 
--------------------------------------------------------------------------------------------------
| Start presolving 
| ...
| Presolve Ends. Elapsed Time: 0.000 seconds 
--------------------------------------------------------------------------------------------------
.
.
.
Phase A Log: 'P': Primal solution found. '*': Phase A ends. 'F': Error happens. 'M': Max iteration
--------------------------------------------------------------------------------------------------
| Iter |         pObj |         dObj |      dInf |      k/t |       mu |   step |    Pnorm |   E |
--------------------------------------------------------------------------------------------------
|    1 |  1.00000e+05 |  0.00000e+00 |  7.07e+05 | 0.00e+00 | 7.75e+02 |   0.00 |  1.0e+20 |     |
|    2 |  1.61239e+05 | -1.99998e+05 |  0.00e+00 | 0.00e+00 | 7.36e+02 |   1.00 |  1.1e+01 |   * |
--------------------------------------------------------------------------------------------------
Phase A Log Ends. 

--------------------------------------------------------------------------------------------------
| DSDP Phase A ends with status: DSDP_PRIMAL_DUAL_FEASIBLE                                         
| Elapsed Time: 0.001 seconds                                                                   
--------------------------------------------------------------------------------------------------
| DSDP Phase A certificates primal-dual feasibility                                                
| Primal relaxation penalty is set to  2.000e+06 
| Perturbing dual iterations by  0.000e+00 
| DSDP Phase B starts. Restarting dual-scaling                                                     
| Heuristic start: mu:  2.800e+02 pObj:  1.612e+05  dObj: -2.000e+05                              
--------------------------------------------------------------------------------------------------

Phase B Log: 'P': Primal solution found. '*': Optimal. 'F': Error happens. 'M': Max iteration
--------------------------------------------------------------------------------------------------
| Iter |              pObj |              dObj |       pInf |       mu |   step |    Pnorm |   E |
--------------------------------------------------------------------------------------------------
|    1 |  1.6123930939e+05 | -1.9999819855e+05 |  2.000e+00 | 2.80e+02 |   1.00 |  1.1e+01 |     |
|    2 |  1.6123930939e+05 | -9.0898247344e+03 |  2.000e+00 | 1.73e+02 |   0.15 |  4.5e+01 |     |
|    3 |  3.5989199238e+04 | -4.3785145278e+02 |  8.693e-05 | 8.33e+00 |   0.23 |  3.0e+01 |   P |
...
|   13 | -2.2999855599e+01 | -2.3000001332e+01 |  2.877e-13 | 1.13e-07 |   0.19 |  3.5e+00 |   P |
|   14 | -2.2999971556e+01 | -2.3000000127e+01 |  5.666e-14 | 2.31e-08 |   0.12 |  1.9e+01 |   * |
--------------------------------------------------------------------------------------------------
Phase B Log Ends. 

--------------------------------------------------------------------------------------------------
| DSDP Phase B ends with status: DSDP_OPTIMAL                                                      
| Elapsed Time: 0.043 seconds                                                                   
--------------------------------------------------------------------------------------------------
| DSDP Ends.                                                                                        
--------------------------------------------------------------------------------------------------
| Primal solution is extracted.                                                                    
| Final pObj: -2.29999e+01   dObj: -2.30000e+01 
--------------------------------------------------------------------------------------------------
| DSDP Time Summary: 
--------------------------------------------------------------------------------------------------
|           Event |    Time(s) | 
--------------------------------------------------------------------------------------------------
|        Presolve |      0.000 | 
|  Phase A (iter) |      0.001 | (2) 
|  Phase B (iter) |      0.043 | (14) 
|           Get X |      0.010 | 
|       Postsolve |      0.000 | 
|             All |      0.055 | (16) 
--------------------------------------------------------------------------------------------------
| DIMACS Error:
--------------------------------------------------------------------------------------------------
| Err1: 4.67e-09 | Err2: 0.0e+00 | Err3: 0.0e+00 | Err4: 0.0e+00 | Err5: 2.6e-06  | Err6: 2.6e-06  
--------------------------------------------------------------------------------------------------
```

#### Contributing

HDSDP is still in its preliminary release and will start accepting pull requests in a future release.

#### Developers

HDSDP is developed by 

- Wenzhi Gao

  Cardinal Operations (Intern)

- Dongdong Ge

  Cardinal Operations

- Yinyu Ye

  The original author of DSDP5.8

and is now maintained by the Cardinal Operations COPT development team.

#### License 

HDSDP is licensed under the MIT License. Check license.md for more details

#### Acknowledgement

We thank all the users using and making suggestings during the development of DSDP/HDSDP solver. Specially we than Dr. Qi Huangfu from COPT development team for his constructive suggestions during the solver development. We sincerely respect the developers of DSDP for their precious suggestions [1] and their invaluable efforts getting DSDP through all along the way. It is the efficient and elegant implementation from DSDP5.8 that guides HDSDP to where it is.

Interested users can refer to the following literature and links for the theory and elegant implementation behind **DSDP (v5.8)**. 

#### References

[1] DSDP5.8 https://www.mcs.anl.gov/hs/software/DSDP/.

[2] Benson, Steven J., and Yinyu Ye. "Algorithm 875: DSDP5—software for semidefinite programming." *ACM Transactions on Mathematical Software (TOMS)* 34.3 (2008): 1-20.

[3] Benson, Steven J., Yinyu Ye, and Xiong Zhang. "Solving large-scale sparse semidefinite programs for combinatorial optimization." *SIAM Journal on Optimization* 10.2 (2000): 443-461.
