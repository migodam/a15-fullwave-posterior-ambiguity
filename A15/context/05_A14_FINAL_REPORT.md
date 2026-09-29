# A14：Task-Relevant Posterior Coordinates 数值研究最终报告

## 1. Executive verdict

**NO-GO，针对当前预定 acquisition、signed-blob prior 和“常见且任务显著的剩余不确定性应由 Hybrid/Flow 保留”这一主线。** Gate A 为 PARTIAL PASS：离散切线正确，局部测量信息谱与任务剩余谱不同，有限幅度全波求解中也确有少量 measurement-weak/task-strong 对象；但五个预登记 regime 的必要跨场景门槛全部失败，原生成器合法连续参数扫描在 96 个场景仅找到 1 个刚过阈值的内部电流 witness。按冻结协议，Gate B、reference posterior 和 Flow 均未执行；因此本轮没有表示或后验生成的成功证据。锁定 test 未打开。

这是否定这套预先定义的广泛机制主张与后续投入依据，**不是**全 Maxwell 或所有照明/目标任务都不存在这种方向的不可能性定理。[Gate A 完整证据](GATE_A.md)、[实现验收](IMPLEMENTATION_VALIDATION.md)、[图像逐张检查](VISUAL_REVIEW.md) 与原始 JSON 是本报告依据。

项目入口为 [README](../README.md)；[最终文件清单](../results/final_manifest.json) 列出代码、数据、报告和图片的 SHA256，冻结协议 SHA256 为 4da141608105c7ab511b22e06762d741f559c2da809edb106edb1705fc41e4ec。

## 2. A13 failure carried forward

[A13 最终报告](../../FlowMultiDecode_Local/reports/SHARED_PHYSICAL_FLOW_FEASIBILITY.md) 对任意压缩 material latent 加独立 J/E/χ/ρ heads 的同样本物理一致性作出 NO-GO。A13 的 Toy 条件后验通过了自己的门槛，但物理 heads 对 J−χE、E−Einc−GDJ 等方程残差不达标，也没有稳定的多 mode 独立预算收益。A14 没有复用旧 encoder、heads 或 posterior checkpoint；canonical state 是真实 χ。所有 J/E/y 均由同一 χ 的复制且校验过的 scalar Helmholtz 求解器产生。这解决读出构造的一致性定义，却**不自动证明**任务后验表示有效。

## 3. Material-current tangent derivation verification

A13 离散 solver 的符号是 J=Dχ(Einc+GDJ)，故 M=I−DχGD；在固定入射波下对实材料求导，MδJ=D_Eδχ，δE=GDδJ，δy=GSδJ。A14 使用同一 LU 解各照明和切线 RHS，复数观测拆成 Re/Im，实参数伴随取复内积实部。代码见 [物理接口](../src/physics.py) 与 [复用后端](../src/physics_backend.py)。

正式预飞行的 8 scenes×8 directions×4 readouts×7 差分步长全部满足相邻两个步长误差低于 10⁻⁵；各项最优误差最大约 2.03×10⁻¹⁰。JVP/VJP 配对最大相对误差 2.53×10⁻¹²，线性性 4.48×10⁻¹⁵。self-cell 积分误差约 4.41×10⁻¹³，零 contrast 电流为零，弱 Born 极限与多次散射差异符合预期；三 scene 的 32/48 网格接收响应相对差为 0.064–0.086%。这只支持当前离散 2D scalar model 的实现正确性，不能上升为连续 3D Maxwell 证明。

## 4. Measurement information spectrum

用 4096 个独立 training-prior 材料估计 C0，加冻结的相对 ridge；训练 128 个物理场景固定观测与任务尺度。主设置 32×32、k=2π、两观测照明 0/π、32 receivers，实观测分量独立 Gaussian noise，σ 为训练观测 RMS 的 5%。定义 A=Γn⁻¹ᐟ² JF C0¹ᐟ²；96 个 validation 场景的 top λ 中位数约 2.91×10⁴。[测量信息谱图](../figures/fig02_measurement_spectrum.png) 显示前部到数值尾部，但小尾部不能被称作精确 null。该谱刻画当前观测对训练 prior 二阶尺度中各方向的压缩，不是 raw Jacobian 排名，也不是非 Gaussian posterior 的解析真值。

## 5. Task ambiguity spectrum

对新照明 π/2 的内部 current Jnew 与半接收器间距偏移阵列的新响应 ynew，分别计算 BQᵀWQBQv=η(I+AᵀA)v。主设置 top ηJ、ηY 场景中位数分别约 0.00591、0.00297；top LIS 与 top J/Y task 子空间第一 principal angle 中位数约 89.97°/89.87°。局部 Gaussian 剩余 task 方差的迹，在训练任务 RMS 单位下，中位数 J 约 0.0420、Y 约 0.00743。[task 谱](../figures/fig03_task_spectra.png)、[完整首场景谱](../figures/fig03b_full_task_spectrum.png) 和 [方向/夹角图](../figures/fig04_eigenvectors_and_angles.png) 都已目视检查。

两个谱在此离散 setting 确实不同，但 generalized eigensystem 本身会偏向弱观测方向；高维中的方向夹角也不是足够的实用效果。经验 covariance ridge 的小 J 尾部更不能当作原 blob prior 的高维独立概率质量。是否值得保留这些方向要看有限幅度对象对和任务 posterior，而不能只看图 2–4。

新照明 internal E 的切线已纳入预飞行差分与伴随检查；由于 J/Y 两项主任务已在冻结 Gate A 失败，没有追加独立 E 谱的全场景 sweep，不能替 E 作机制裁决。

## 6. Finite-amplitude validation

每方向以 vᵀ(I+AᵀA)v=1 归一化，幅度按冻结网格 0.03125–2 双侧扫描；对 χ± 独立全波求解。测量分布同协方差时，KL=0.5‖Γn⁻¹ᐟ²(F(χ+)−F(χ−))‖²。witness 必须同时有 measurement fraction≤0.25、真实 KL≤1、留出任务差异≥0.25 个训练 prior 任务 RMS。

主设置三 seed 32 scene 的 J witness 为 4、4、1；Y 为 0、0、0。J 在 KL≤1 下的最佳差异场景中位数约 0.158，Y 约 0.0647。所选 J secant 误差中位数约 0.00121，说明此失败主要不是有限幅度立刻推翻局部线性化。task-aware 相比 POD/LIS/best-of-32 random 的三 KL 预算平均包络有显著配对优势，但 J 绝对平均效果约 0.107、Y 约 0.055，远低于单场景 0.25 的筛选效果阈值。最佳随机方向的搜索成本已包含 32 次方向尝试。不能以很大的相对倍数替代跨场景有效性。

8 receiver 边界亦完整三 seed：J witness 为 3、4、3，Y 仍全零。弱散射、强散射、1% 噪声各自的 seed 11 完整 32 场景计数分别是 J/Y=5/0、5/1、1/0，都严格小于所需 16/32。因此任一后续 seed 的结果都无法使这些 regime 通过全 seed Gate。其余场景按预定必要条件的精确 futility 中止，已完成的部分 seed 文件保留且未混入完整估计。[regime 图](../figures/fig_gate_a_regime_comparison.png) 清楚区分数据覆盖。锁定 test 未打开。

## 7. Representation audit

**未运行 Gate B。** 项目实现 raw/KL、正交可逆 Haar、POD、LIS、task-only、Hybrid 与先验剩余量接口，并完成 Haar 旋转的 full-rank 谱不变性检查（相对误差约 8.43×10⁻¹⁶）。这只证明坐标实现基本正确，不是 rank 16/32/64/128 的 forward distortion、task covariance 或 posterior variance 对比。没有可报告的 Hybrid 优于 POD、Wavelet-only、LIS-only、raw/full 的实验证据。确定性 LIS 截断删掉的方差不能拿来代表保留补空间先验的经典 LIS。[Gate B 状态](GATE_B.md)。

## 8. Explicit ambiguity pairs and prior support

[Fig 6 affine pair](../figures/fig06_explicit_material_pair.png) 对同一 validation 场景展示 χ±、留出 J 与观测幅值：真实复数观测 KL=0.326，J 差异=0.276；它是经验 Gaussian covariance 椭球的代理对，不保证属于原生成器。固定宽度/轴比的合法 blob 参数扫描在 96 场景 0 个 witness。开放所有六个连续 blob 参数但保留合法区间、个数和符号—尺寸规则后，J 为 0/32、1/32、0/32，Y 为 0/96。唯一合法 witness 是 [Fig 6e](../figures/fig06e_full_prior_pair.png) 中的 scene 410004：KL=0.816、J=0.255、measurement fraction=0.227。它证明原 prior 支持内有条件地存在这种对，但不是稳定可复现的占优 regime。完整合法扫描的最佳 J 差异中位数 0.0625；只有 7/96 场景优于 affine 扫描的最佳值。两者搜索空间不等，不能当公平算法赛。

## 9. Reference posterior

**未运行。** 16 维受控物理 prior、三次独立 tempered SMC、ESS/祖先退化/接受率、粒子加倍、129²/257² 切片 quadrature 全部受 Gate A/B 阻断。源码原型和单元测试不等于 posterior truth；没有原始 1024D 或受控 16D 的 posterior reference。[阶段记录](POSTERIOR_REFERENCE.md)。

## 10. Flow experiments

**未训练、未评价。** F0 MAP/Laplace 与 F1–F4 conditional Flow 比较只有受阻断的实现准备，没有真实 posterior energy score、coverage、mode probability、NFE 或摊销交叉点。本轮不能声称 Flow 必需，也不能声称不必需；只可说当前 Gate A 失败使训练它缺乏预注册的依据。[阶段记录](FLOW_RESULTS.md)。

## 11. Task posterior fidelity

**未验证。** 虽然从每个给定 χ 通过真实 solver 取得的 J/E/y 天然属于同一物理 realization，本轮没有经过诊断的 posterior samples，因此无法评估 p(Jnew,Enew,ynew|y) 的覆盖率、方差、交叉协方差或分支概率。局部 Gaussian 的 η 与迹是 mechanism 诊断，不能替代这些分布指标，也不能继承 A13 Toy 的 posterior success。

## 12. Computational cost

全部任务留在当前 Mac，复杂线性系统用 CPU float64/complex128；MPS 仅用于未训练小网络的 FP32 数值一致性检查。逐场景程序保存协议 hash、split/scene seed、谱、方向、七幅度全波扫描和耗时。已完成场景的 LU/RHS 与累计 scene 秒数由 [成本 JSON](../results/computational_costs.json) 按源码调用结构重建；它是下界，排除有意中断场景、训练尺度生成与一次性诊断，且并发 scene 秒数不能当作实际经过时间。预飞行内存峰值约 0.94 GB，两个后期并发 Gate A 进程各约 0.24–0.32 GB RSS。没有 SMC/Flow 训练或推断时间可以报告。

具体完成量为 Gate A **459 个场景、至少 283,203 次 LU 和 2,731,050 个线性 RHS、累计 15,264 scene 秒**；固定参数与完整合法参数支持审计各 96 场景，分别至少 6,268/5,546 次 LU。12 个训练尺度块另至少 1,536 次 LU，但其耗时未单独仪表化。按这些已记录块相加，至少约 296,553 次 LU、2,784,082 个 RHS；该合计仍遗漏中断时的部分工作和预飞行等一次性计算。最终项目目录约 0.1 GB；20 项单元检查通过。

## 13. Failure cases and boundaries

第一，signed-contrast 原 prior 的约半数训练样本有 minχ≤−0.95，不能把所有场景解释为被动介质。第二，经验 Gaussian 白化只复现二阶尺度；affine proxy 的许多方向离开原 blob 生成器支持。第三，唯一真实合法 witness 刚越过 0.25 门槛，未说明其 posterior mass 或出现概率。第四，扫描只覆盖 top 3 generalized directions 和预设七个幅度，不是对全部非 Gaussian posterior 的穷举。第五，3π/2 新照明 24 个 validation 场景的 J 效果中位数约 0.150、最大 0.234；没有在该边界救活主张。第六，Born-to-full 差与 J 效果的探索性 Spearman ρ≈0.35 可能指向多次散射条件，但散点分散、witness 稀少，不能推因果。第七，符号—尺寸关联反转 OOD 的独立 96 场景中 J witness 为 3/32、10/32、10/32，Y 仍 0/96；两组虽比 ID 多，但仍低于 16/32。OOD 的每 seed 最佳 J 效果中位数为 0.139、0.146、0.162，并未同步普遍增强。ID/OOD 不同随机 scene，且 OOD 不是原 prior 的确认集；[描述性图](../figures/fig_gate_a_ood_boundary.png) 不作因果或 posterior claim。

## 14. Literature and novelty implications

[LIS 原始工作](https://arxiv.org/abs/1403.4680) 明确保留补空间的 marginal prior；[非 Gaussian prior 正规化 LIS](https://arxiv.org/abs/2202.00074) 进一步限制仅凭协方差白化的解释。[Goal-oriented Bayesian reduction](https://arxiv.org/abs/1607.01881) 已研究目标量 posterior covariance/mean 的低秩近似；[null-space FM 近邻](https://arxiv.org/abs/2604.22005) 使“Flow 只放 null”不能算主要新意。故理论上有价值的候选只能是明确 full-wave 材料—电流—接收器链中、原 prior 支持内、有限幅度与 task posterior 都成立的谱错位，而本轮未满足后两项的可重复门槛。[有界文献审计](LITERATURE_AUDIT.md) 因 ScholarQA 查询遇到 429，并非穷尽 prior art。

## 15. Connection to N-Port / overlap / A9–A12

A14 的 D_E 注入、M⁻¹ 电流传递、GS 测量过滤，与 A9–A12 的 material/current tangent 和 measurement filtering 在算子结构上直接相接；本轮又添加先验、噪声和 downstream task，问 posterior 残留方向。这个连接不是 A9–A12 已证明 A14 posterior claim 的证据。本轮没有调用 N-Port、Maxwell 后端或远程机器；详细边界见 [接口报告](CONNECTION_TO_OVERLAP_AND_TRANSFER.md)。

## 16. Final recommendation and Q1–Q7

| 问题 | 本轮直接回答 |
|---|---|
| Q1：测量信息谱与任务剩余谱不同？ | **局部 Gaussian/离散切线意义上是。** 方向和特征值不同，但这不足以保证原 prior 下任务显著。 |
| Q2：有限幅度 full-wave 仍有意义？ | **稀有条件下有，未跨场景复现。** affine J 9/96，合法 blob J 1/96，Y 0/96；五 regime 必要门槛全失败。 |
| Q3：Hybrid 是否优于 POD/Wavelet/LIS/raw？ | **未验证。** Gate B 被冻结规则关闭。 |
| Q4：需要 Flow，还是 Laplace/SMC 已足够？ | **未验证。** 没有可信 reference，也没有 Flow 比较。 |
| Q5：若需 Flow 应承担哪部分不确定性？ | **无数据裁决。** 理论候选是测量更新后仍影响留出任务的相关条件剩余，而非独立补一个 null marginal。 |
| Q6：联合新 J/E/y 后验是否正确？ | **未验证。** 单一 χ 物理读出保证样本内方程成立，不保证分布、coverage 或 mode 概率正确。 |
| Q7：等级？ | **NO-GO，针对本次预定义的广泛 A14 主 claim 与后续 Hybrid/Flow 投入。** 保留稀有合法 witness 作为重新定义物理 regime 的线索，不后验改 Gate。 |

下一轮若继续，应先用独立验证协议研究原 blob 参数流形中的稀有 witness、特定多次散射场景及新 acquisition/task 的效果边界，再决定是否建立 reference posterior。此轮既没有证明 Flow 的价值，也没有用大模型去补救失败的物理门槛。
