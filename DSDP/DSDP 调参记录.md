### DSDP 调参记录

1. cnpil 系列问题的 C 矩阵为 0，为原始问题的可行性问题，在得到第一个原始问题可行解之后可以之间停止

2. theta 系列问题容易在 Phase A 找到可行解后停留很久，可以考虑关闭 Attempt 

   theta12 问题的约束数极多，最好考虑调整 CG 的精度与迭代次数，让 Schur 系统的求解集中在 CG 中

3. thetaG11/thetaG51 问题需要从大的 Infeasibility 开始初始化 (kappa 与 mu 取大, init beta 1e+03)

3. buck 问题需要 Phase A 的步长较小，否则会由于步长过小出现问题（现在可以加入 Phase A 专用 Corrector 步解决）

4. rose 问题后期会出现由于数值问题导致的 Primal 超过 Dual 的情况，现在通过降低 gap 的精度解决

4. biggs 问题的 Dual infeasibility 非常难以消除，后续需要特殊处理

5. lzc 问题可以在 450s 解决，但目前为 inaccurate，可以考虑在完成 Primal relaxation 后增加精度

6. Phase A upperbound 的公式可能有误（导致 Primal 下降超过 Dual），需要检查修正

7. Phase A 的 proximity norm 可能出现 nan，可能是 Schur 矩阵计算有误（M5 Technique 的数值问题）



在 DSDP5.8 的 DSDPSetup 中可以获取其参数启发式选取策略



#### 需要解决的问题

1. Debug Schur 系统
2. Bound Cone 的构建方法