### DSDP 调参记录

1. cnpil 系列问题的 C 矩阵为 0，为原始问题的可行性问题，在得到第一个原始问题可行解之后可以直接停止

3. thetaG11/thetaG51 问题需要从大的 Infeasibility 开始初始化 (kappa 与 mu 取大, init beta 1e+03)

3. buck 问题需要 Phase A 的步长较小，否则会由于步长过小出现问题（现在可以加入 Phase A 专用 Corrector 步解决）

4. rose 问题后期会出现由于数值问题导致的 Primal 超过 Dual 的情况，现在通过降低 gap 的精度解决

4. biggs 问题的 Dual infeasibility 非常难以消除，后续需要特殊处理

7. lzc 问题可以在 450s 解决，但目前为 inaccurate，可以考虑在完成 Primal relaxation 后增加精度

8. Phase A 的 proximity norm 可能出现 nan，可能数值问题

9. m 巨大的问题需要加入 CG Reuse 参数的启发式调整（重复使用 15 次左右，重复使用完后先回到 Diagonal scaling 再考虑 Cholesky 分解）

   ldc1024 hamming_8_3_4 theta12 等



在 DSDP5.8 的 DSDPSetup 中可以获取其参数启发式选取策略





#### 需要解决的问题

1. Debug Schur 系统
2. HSD 内加入 Bound cone 还需要额外考虑



### Mittelmann Hans 调参（无 Corrector）

| 问题名称      | $\kappa$ | $\tau$ | $\beta$ | $\mu$ | $\rho$ | $\alpha_A$ | CG Reuse | 提升                     | Iteration  | 备注                                                         |
| ------------- | -------- | ------ | ------- | ----- | ------ | ---------- | -------- | ------------------------ | :--------: | ------------------------------------------------------------ |
| 1dc.1024      | 1.0      | 1.0    | 1.0     | 1e+10 | 5.0    | 0.75       | 6        | 可以求解                 | 50(2550s)  | 1. m 巨大  2. 降低 CG 精度/增加迭代次数 3.提高利用次数       |
| 1et.1024      | 1.0      | 1.0    | 1.0     | 1e+10 | 5.0    | 0.75       | 10       | 可以求解                 |    315s    | 1. m 较大 2. MKL 函数                                        |
| 1tc.1024      | 1.0      | 1.0    | 1.0     | 1e+10 | 5.0    | 0.75       | 10       | 195s                     |  58(195s)  | 1. m 较大 2. 迭代次数较多，比原版略慢，需加入 Corrector      |
| 1zc.1024      | 1.0      | 1.0    | 1.0     | 1e+10 | 5.0    | 0.75       | 10       | 提速(1400s 到450s)       |  44(450s)  | 1. m 巨大  2. 降低 CG 精度/增加迭代次数 3.提高利用次数       |
| hamming_8_3_4 | 1.0      | 1.0    | 1.0     | 1e+10 | 5.0    | 0.75       | 10       | 提速(243s 到 25s)        |  12(25s)   | 1. m 巨大 2.  CG 条件数很好，一般不需要 Cholesky 分解 3. MKL 函数失效，需要调用 Lanczos |
| hamming_9_5_6 |          |        |         |       |        |            |          |                          |            | m 超出内存范围，无法求解                                     |
| buck3         | 1.0      | 1.0    | 1.0     | 1e+10 | 4.0    | 0.5        | 10       | 可以更精确求解           | 158（15s） | Phase A 步长需要较小（或加入 Corrector 步）                  |
| buck4         | 1.0      | 1.0    | 1.0     | 1e+10 | 4.0    | 0.5        | 10       | 可以更精确求解           |  185(79s)  | 1. 对初始化极为敏感 2.Phase A 步长需要较小（或加入 Corrector 步） |
| buck5         |          |        |         |       |        |            |          | 调参中                   |            | Phase A 步长需要较小（或加入 Corrector 步）                  |
| trto3         | 1.0      | 1.0    | 1.0     | 1e+10 | 4.0    | 0.5        | 10       | 可以更精确求解（7s)      |     92     | 对 Phase A 步长较敏感                                        |
| trto4         | 1.0      | 1.0    | 1.0     | 1e+10 | 4.0    | 0.5        | 10       | 可以更精确求解           |    139     | 对 Phase A 步长较敏感                                        |
| trto5         |          |        |         |       |        |            |          |                          |            |                                                              |
| inc_600       | 1e+10    | 1.0    | 1.0     | 1e+10 | 3.0    | 0.75       | 10       | 可以低精度求解 (1e-02)   |  86(72s)   | 1. 选择 mu 时 MKL 函数经常出现数值错误，需要使用 Lanczos 2. 在 Pardiso 未修复时需要在 Potential recduction 与选择 mu 时非常保守 |
| inc_1200      | 1e+10    | 1.0    | 1.0     | 1e+10 | 3.0    | 0.75       | 10       | 可以较低精度求解 (1e-02) |  72(446s)  | 1. 与 inc_600 一样，需要检查 selectMu 函数及 Potential reduction 策略 |
| theta12       | 1e+10    | 1.0    | 1.0     | 1e+10 | 4.0    | 0.75       | 10       | 844s                     |  29(844s)  | 1. C 矩阵为常数矩阵，初始化需要较大残差 2. m 巨大，需要考虑特殊的 CG 策略 |
| theta102      |          |        |         |       |        |            |          |                          |            | m 超出内存范围，无法求解                                     |
| theta123      |          |        |         |       |        |            |          |                          |            | m 超出内存范围，无法求解                                     |
| rose13        | 1e+10    | 1.0    | 1.0     | 1e+10 | 4.0    | 0.75       | 10       | 可以更精确求解           |  115(46s)  |                                                              |
| rose15        |          |        |         |       |        |            |          |                          |            |                                                              |
| G40_mb        | 1e+10    | 1.0    | 1.0     | 1e+10 | 4.0    | 0.75       | 10       | 调参中                   | 202(202s)  | G 系列问题求解难度较低，但目前迭代次数较多，需要加入 Corrector 步 |
| G48_mb        | 1e+10    | 1.0    | 1.0     | 1e+10 | 4.0    | 0.75       | 10       | 可以更精确求解           |  32(264s)  |                                                              |
| G40mc         | 1e+10    | 1.0    | 1.0     | 1e+10 | 4.0    | 0.75       | 10       | 调参中                   |  71 (70s)  |                                                              |
| G48mc         | 1e+10    | 1.0    | 1.0     | 1e+10 | 4.0    | 0.75       | 10       |                          | 16(10.4s)  |                                                              |
| G55mc         |          |        |         |       |        |            |          |                          |            |                                                              |
| G59mc         |          |        |         |       |        |            |          |                          |            |                                                              |
| torusg3-15    | 1e+10    | 1.0    | 1.0     | 1e+10 | 5.0    | 0.75       | 10       | 186s                     |  26(186s)  | 1. 需要 scaling 求解                                         |
| hand          |          |        |         |       |        |            |          |                          |            |                                                              |
| foot          |          |        |         |       |        |            |          |                          |            |                                                              |





- Alh_1.r20  

- BH2.r14    

- Bex2_1_5   

- Bst_jcbpaf2

- CH2_1.r14  

  

- H30_.r16   

- NH2-.r14   

- NH3.r16    

- NH4.r18   

- biggs      

- broyden25  

  

- cancer     

- checker    

- chs_5000   

- cnhil10    

- cnhil8     

- cphil10    

- cphil12    

- dia_patch  

- e*quad*2   

- e*stable*2 

- ice_2.0     

- mater-4    

- mater-5    

- mater-6    

- neosfbr25  

- neosfbr30e8

- neu1       

- neu1g      

- neu2       

- neu2c      

- neu2g      

- neu3       

- neu3g      

- nonc_500   

- p_auss2    

- prob_1_2_0 

- prob_1_2_1 

- prob_2_4_0 

- prob_2_4_1 

- rabmo      

- reimer5    

- r1_2000    

- ros_500    

- ros_2000        

- sensor_1000

- sensor_500 

- shmup3     

- shmup4     

- shmup5     

- spar060    

- swissroll  

- taha1a     

- taha1b     

- taha1c     

- 

- t_texture  

- 

- vibra3     

- vibra4     

- vibra5     

- yalsdp     