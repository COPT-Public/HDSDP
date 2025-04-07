## HDSDP 

#### A Homogeneous Dual-scaling Interior Point Solver for Sparse Semi-definite Programming

HDSDP implements a dual-scaling interior point method for sparse Semi-definite programming problems (SDPs). The solver inherits several features from  **DSDP5.8** [1] and presents some new aspects for the dual-scaling algorithm. HDSDP is written in ANSI C and is maintained by Cardinal Operations COPT development team. More features are still under active development.

#### Solver overview

HDSDP solves standard-form semi-definite programs expressed by

$$\min_{\mathcal{A} \mathbf{X} = \mathbf{b}, \mathbf{X}\in \mathbb{S}_+^n} \left\langle \mathbf{C}, \mathbf{X} \right\rangle  \\
  \\ 
$$

where we do linear optimization subject to affine constraints over the cone of positive-semidefinite matrices.

#### License 

HDSDP is licensed under the MIT License. Check license.md for more details

#### Cite

```
@article{gao2025algorithm,
  title={Algorithm xxxx: HDSDP: Software for Semidefinite Programming},
  author={Gao, Wenzhi and Ge, Dongdong and Ye, Yinyu},
  journal={ACM Transactions on Mathematical Software},
  year={2025},
  publisher={ACM New York, NY}
}
```

#### References

[1] DSDP5.8 https://www.mcs.anl.gov/hs/software/DSDP/.

[2] Benson, Steven J., and Yinyu Ye. "Algorithm 875: DSDP5—software for semidefinite programming." *ACM Transactions on Mathematical Software (TOMS)* 34.3 (2008): 1-20.

[3] Benson, Steven J., Yinyu Ye, and Xiong Zhang. "Solving large-scale sparse semidefinite programs for combinatorial optimization." *SIAM Journal on Optimization* 10.2 (2000): 443-461.
