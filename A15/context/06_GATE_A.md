# Gate A：材料—电流—观测机制的冻结验证

**裁决：PARTIAL PASS；不得进入 Gate B、reference posterior 或 Flow。** 局部测量信息谱与留出任务剩余谱分离，有限幅度下也有少数 full-wave witness；但五个预登记 validation regime 都无法满足“每个 seed 至少 16/32 场景”的必要条件。锁定 test 没有打开。此裁决依据 [冻结协议](../EXPERIMENT_PROTOCOL.md)、[逐场景原始数据](../results/gate_a/) 与 [可核验的提前停止证明](../results/gate_a_futility_certificate.json)，不是事后改阈值。

## 离散物理和切线先过关

复用 A13 的 corrected-self-cell scalar Helmholtz 后端。实际离散方程是 J=χE、E=Einc+GD J，因此 M=I−Dχ GD，材料切线满足 M δJ=D_E δχ。接收切线为 GS δJ。实参数 VJP 取复内积实部。八场景、每场景八随机方向、四个读出与七个步长的正式差分，最优误差逐项最大约 2.03×10⁻¹⁰；伴随配对最大 2.53×10⁻¹²。self-cell 积分、零 contrast、Born 极限、多次散射、线性残差、32/48 网格与 Haar 谱不变性见 [实现验收](IMPLEMENTATION_VALIDATION.md)。这些只验证离散实现。

## 预定主设置的两个谱

32×32、观测照明 0/π、32 receiver、训练观测 RMS 的 5% 实分量噪声，留出照明 π/2。96 个 validation 场景（11/29/47 各 32）中，top measurement-information eigenvalue λ 的场景中位数为 2.91×10⁴；top task-residual eigenvalue ηJ、ηY 中位数分别为 0.00591、0.00297。top LIS 与 top J-task / Y-task 的第一个 principal angle 中位数约 89.97° / 89.87°。[测量谱](../figures/fig02_measurement_spectrum.png)、[top-3 任务谱](../figures/fig03_task_spectra.png)、[首个预定场景完整任务谱](../figures/fig03b_full_task_spectrum.png) 与 [方向图](../figures/fig04_eigenvectors_and_angles.png) 已目视审查。经验 covariance 的 ridge 让 J 谱有小尾部；它不应被解释为独立可实现的材料自由度。

局部 Gaussian 近似的总任务方差在训练 task RMS 单位下，J 中位数约 0.0420，Y 中位数约 0.00743。因而“两个谱不同”在代数意义上成立，但 generalized eigenproblem 本来就惩罚测量强方向；正交性不能替代任务效果量或 posterior calibration。

## 有限幅度 full-wave 对象对

所有方向先满足 vᵀ(I+AᵀA)v=1，再按冻结的七个幅度构造 χ± 并分别真实求解。witness 同时要求局部 measurement fraction≤0.25、两观测分布 KL≤1，以及归一化留出任务差异≥0.25。主设置 J 的三 seed 计数为 **4/32、4/32、1/32**；新 receiver 任务 Y 为 **0/32、0/32、0/32**。J 在 KL≤1 下的最佳任务差异中位数约 0.158，Y 约 0.0647；J 的所选有限差分 secant 误差中位数约 0.00121。当前失败主要是可达到的任务效果与跨场景频率不足，而非这些幅度上普遍丧失线性化。

在 KL 预算 0.1、0.3、1 的场景平均包络里，task-aware J 为 0.107，best-of-32 random 为 0.0050；Y 为 0.0553 对 0.0051。按完整 scene 配对 bootstrap、对 J/Y 作 Bonferroni 后的 97.5% 区间，J 的差值为 [0.0956,0.1087]，Y 为 [0.0437,0.0567]。相对排名很强，但绝对效果与 prevalence 不达门槛；不能用相对倍数掩盖这一点。[有限幅度 trade-off](../figures/fig05_measurement_vs_task.png) 及 [affine 对象对](../figures/fig06_explicit_material_pair.png) 有代表性，但该对象对仅属于经验 Gaussian surrogate 椭球，不自动属于原始 blob prior。

## 原始 prior 支持审计

先在原 blob 生成器的合法幅度、中心、角度内扫描，同时固定宽度和轴比：96 场景中 J/Y witness 均为 0。再开放所有六个连续参数（幅度、两中心、角度、宽度、轴比），保持 blob 数、合法范围及尺寸—符号规则：J 为 **0/32、1/32、0/32**，Y 为 **0/96**。唯一合法 J witness 的 scene 410004 在幅度 1 下有 KL=0.816、J 差异=0.255、measurement fraction=0.227；同 scene 的 affine 扫描未找到 witness。完整合法扫描 96 场景的最佳 J 差异中位数仅 0.0625，只有 7 个场景超过 affine 扫描的最佳值。两扫描的搜索空间不同，此计数不是算法公平排名。这说明机制在原生成器支持内确实可能出现，但目前只得到稀有、刚过阈值的实例。[支持审计](../figures/fig06d_full_prior_support_audit.png)、[合法对象对](../figures/fig06e_full_prior_pair.png) 均已打开核对。Top-3 方向与固定幅度网格不穷尽完整非 Gaussian posterior。

## 五种 regime 与提前停止

| Regime | seed 11 J/Y witness | seed 29 J/Y | seed 47 J/Y | Gate 作用 |
|---|---:|---:|---:|---|
| main，32 receiver、5% noise | 4/0 | 4/0 | 1/0 | 三 seed 完整 |
| sparse8，8 receiver、5% noise | 3/0 | 4/0 | 3/0 | 三 seed 完整；低观测边界 |
| weak，0.25χ | 5/0 | 未完整 | 未运行 | seed 11 已排除 |
| strong，2χ | 5/1 | 未完整 | 未运行 | seed 11 已排除 |
| low_noise，1% noise | 1/0 | 未完整 | 未运行 | seed 11 已排除 |

每格是 32 scene 中的计数；“未完整”保留原始已完成场景，但不当作 32 场景估计。每一个 regime 的 seed 11 都低于必要的 16/32；后续 seed 不可能救活该 regime。弱/强/低噪声按 [提前停止记录](../results/gate_a_early_stop.json) 中止，节省本地求解成本，而不改协议或把缺失记为阴性。[Regime 图](../figures/fig_gate_a_regime_comparison.png) 只绘制完整 seed 块。[完整分析 JSON](../results/gate_a_validation_summary.json) 逐任务列出随机、POD、LIS 的配对比较与 bootstrap；因 prevalence 已失败，无 regime 被选入锁定 test。

## 边界、反例和解释

训练 prior 中满足 minχ>−0.95 的比例，三 seed 分别约 48.95%、48.41%、49.37%；它只是条件子集，不能替代原 signed-contrast 主 prior。3π/2 留出照明的 24 场景检查中 J 差异中位数约 0.150、最大 0.234，没有救活门槛。探索性的 Born-to-full receiver 差与最佳 J 效果 Spearman ρ≈0.35，提示多次散射可能改变强度，但不能从相关性推因果或事后重选 regime。

预定的 sign-size 关联反转 OOD 独立 96 场景已完成，J witness 三 seed 为 **3/32、10/32、10/32**，Y 仍 0/96；最多的两组仍低于 16/32。OOD 的 J 最佳有限任务效果场景中位数分别为 0.139、0.146、0.162，未随 witness 计数同步普遍增加。与 ID 是不同随机 scene，不能按 scene 配对或解释为单一几何符号因素的因果效应；这里也没有对 OOD 原 prior 支持的参数扫描。[OOD 边界图](../figures/fig_gate_a_ood_boundary.png) 与 [原始汇总](../results/gate_a_ood_diagnostic.json) 是诊断，不参与验证选择。

主动反例已出现：LIS/task 方向确实分离，但任务效果不够普遍；Gaussian surrogate 的可行对象对通常比合法 blob 扰动强，也会漏掉极少数合法方向。测量拟合好只说明当前 acquisition 难分辨，不证明未来任务有重要不确定性。当前数据不支持“应使用 Hybrid posterior coordinates”的主 claim；后续 Gate 均为 **NOT RUN**，参见 [Gate B](GATE_B.md)、[reference](POSTERIOR_REFERENCE.md)、[Flow](FLOW_RESULTS.md)。
