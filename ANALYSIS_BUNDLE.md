# A15 analysis bundle

Public snapshot, 2026-09-29. Read [START_HERE](https://raw.githubusercontent.com/migodam/a15-fullwave-posterior-ambiguity/main/START_HERE.md). The following are existing reports, not new experiments.


---

Source: [A15_FULLWAVE_POSTERIOR_REPORT.md](https://raw.githubusercontent.com/migodam/a15-fullwave-posterior-ambiguity/main/A15/reports/A15_FULLWAVE_POSTERIOR_REPORT.md)

# A15 被动双谐振体：全波后验任务歧义与采集设计

## 1. Executive verdict

**NO-GO：本轮不支持进入 Flow 或按已验证 bright–dark population mechanism 撰文。**

这里的 NO-GO 是对当前实现、冻结 family 与证据强度的研究决策。**后验参考不足不等于物理不存在。** 本轮数值实施已执行至预算/验收停止条件；确认测试和条件性 Flow 没有执行，不计作完成的正面实验。

| 项目 | 裁决 |
|---|---|
| 物理求解 | 当前圆柱标量全波数值验收通过；原VIE网格不足 |
| Physics mechanism | 当前两模态定量解释未通过，存在不可激发暗方向反例 |
| Posterior risk | 参考未收敛，population claim 未决；Gate A未建立 |
| Acquisition design | 有一致的未认证诊断线索，确认结论未决 |
| Gate B | FAIL：未满足通过条件，含参考未决 |
| Flow | **FLOW NOT RUN — GATE CLOSED** |
| TAP-paper potential | NO-GO for current standalone mechanism claim；文献核验也不完整 |

## 2. 模型、prior 与冻结边界

两个被动介质圆，r=.1、εbase=12、k0=8.625。θ=(s,Δ,γ) 为独立uniform，范围 [.02,.16]×[−3,3]×[.001,.08]。Q=ln固定目标圆盘平均|E|²/入射平均|Einc|²，是场能增强代理量，可为负。εimag按γωref/ω跨频变化，非负；所有readout源自同一个θ。

主观测为主频两个复通道，σ=0.004051675402/实分量；邻频/近场/额外同频均新增两个复数值。prior盒在所选family内保留全部质量，没有validation后缩窄。确认盒相对于六模板均匀混合才是1/6，而非普适物理总体概率。

协议 SHA256 `eddc3a97a3c39f0bdead44956a3bffadaaa103d2a8bb2568e0d1b9b9658cb38b`。源原文和A14不修改。每次compact/阶段恢复完整重读prompt的规则已写入AGENTS.md与恢复记录。

## 3. 物理与数值验收

复用Green/self-cell后发现32–128网格不能满足低噪声精度，因此在确认前冻结全圆柱多次散射后端。order20对24的24场景最大白化输出误差4.04e−10、Q误差4.19e−12、闭合面功率平衡≤5.62e−15。零contrast、Born、互易性、单圆柱、active/full一致性、几何面积、JVP/VJP及相邻差分步长核验通过。详见 PHYSICS_VALIDATION.md 与44项测试。

这些是限定数值证据，不是普遍Maxwell定理。等面积圆盘self-cell不是精确方形cell积分，VIE模态仅作有偏诊断。

## 4. Discovery 与反例

六模板、每模板24状态×41频率；前三候选再做9³粗积分和独立screening查询。未把这些查询计入prevalence。选择标准偏向未消除任务风险，而非最大场强。

实测图显示冻结k0更接近目标内部能量谷；没有可认证的两暗/亮共振峰指认。81点细化仅在冻结后作诊断，不改变频率。此实施偏离原定“先细化再冻结”的时间顺序已保留；不把它追认为选型证据。

## 5. 后验参考

33³与41³耦合表共314,574次全波求解；另33³去cross-block反事实107,811次。表全部复制回本地，节点、正积分权重、三频y/Q/吸收/电流能量、残差、配置和源码hash齐全。连续truth独立于积分节点，3×128 validation queries。

主噪声基础acquisition仅94/384通过全部33/41差异阈值；最大均值差0.11836、端点差0.15384。4个独立SMC运行共32,768次求解；两个query正权自适应cubature追加19,387次。仍没有达到足以认证全体查询的参考精度。不能靠有效节点数、温度达到1或少数好看query宣告收敛。

## 6. Population risk（仅诊断）

41³主噪声表：R_Q=0.01292496，posterior-mean MSE=0.01345738，90% coverage=86.20%，ambiguous=0.00%（0/384）。33³为1/384。**这些未通过积分验收，不能当真实总体结论或阴性证明。** 5%边界出现更高表上比例，但不能以此修改主噪声。

整个query做5,000次paired bootstrap。没有删除未收敛query后重定义population，没有把节点、粒子、像素或接收器当独立重复。锁定3×512 test未开启。

## 7. 等预算采集

| 方案 | 未认证 R_Q | 表上风险下降 |
|---|---:|---:|
| base1 | 0.0209962 | -62.45% |
| base | 0.0129250 | 0.00% |
| same_angles | 0.0098180 | 24.04% |
| optimized_angles | 0.0057502 | 55.51% |
| neighbor | 0.0038403 | 70.29% |
| near | 0.0006875 | 94.68% |
| random | 0.0062992 | 51.26% |
| base8 | 0.0056772 | 56.08% |

表上near最好，neighbor其次，optimized同频优于固定同频。neighbor/near相对optimized多降约14.78/39.17个百分点；独立seed方向一致，但CI不含数值积分偏差，不能裁决确认性Holm PASS。

事前排序来自全波局部切线proxy，不能归功于尚未通过的coupled-mode公式。邻频设备/源成本与近场布设成本另列，不把通道数公平称为所有工程成本公平。

## 8. 模态与耦合归因

两resolvent态在16/16场景外部误差超过噪声；response-POD在14/16超过。快照电流能量高保留不保证likelihood/task准确。更严格hidden-fiber margin没有建立，任意两极点拟合也未穷尽。

在当前反射对称geometry和入射下，奇对称电流模式不可被θ方向注入，receiver和ROI一阶task也可能同时不敏感。暗并不等于后验任务歧义。Fig11四项与余项记录支持这个反例。

去cross-block表上risk大降、neighbor对优化同频优势消失，与事前方向一致；由于反事实积分未独立收敛，不能形成归因证书。大间隙实物对照显示反馈减弱但gap1.2仍非可忽略，不能将其冒称uncoupled真实装置。

## 9. 显式材料与场

Fig5–7展示seed11、主噪声validation内最大表上方差query。只称低/高Q节点，不称两物理分支。10/90%节点ROI比约1.71，但基础观测KL约4.23：它不是一个强不可区分对象对。图上幅值相近不能替代含相位、噪声白化的距离。

## 10. Flow

**FLOW NOT RUN — GATE CLOSED。** 不写FLOW NOT NEEDED，因为Laplace/tiny mixture的正式充分性没有验证；也不写Flow needed。NF/FM、50k/200k训练、分支概率和摊销交叉点都未运行。

## 11. Q1–Q10

1. **Q1：非小prior mass外等价/内不同family？** 未建立；当前候选没有合格确认性证据，参考不足。
2. **Q2：risk多大？** 表上约0.012925（ln-energy²），未认证；不能当最终Bayes risk。
3. **Q3：何种来源？** 有receiver filtering和反馈诊断线索；不可激发暗模是反例；不能把普通有限观测与bright–dark反馈完全区分。
4. **Q4：两模态能定量解释？** 当前两种投影不能在噪声预算内解释；不排除其他更合适的模型，但本轮不继续无限调参。
5. **Q5：哪个新增测量最有效？** 未认证排序near > neighbor > optimized同频；真实population优劣未确认。
6. **Q6：事前预测成功？** 局部全波proxy预测排序一致；two-mode理论预测成功未建立。
7. **Q7：够独立TAP mechanism paper吗？** 不够；缺可信population/reference及定量解释，也缺完整近邻文献核验。
8. **Q8：A9/A12关系？** material→current transfer→receiver filtering是自然数学接口，目前只宜discussion/methodological connection。
9. **Q9：需要Flow吗？** 未知；缺通过gate的后验复杂性证据。
10. **Q10：Flow状态？** **FLOW NOT RUN — GATE CLOSED**。未验证简单近似，不能误写FLOW NOT NEEDED。

## 12. 计算、复现与剩余缺口

主运行合计约49.4万次forward calls（另有少量pilot/benchmark及worker探索，精确全生命周期峰值/调用数未从第一步统一采集，详见 COST_AND_PROVENANCE.md）。CPU complex128用于表；CUDA/GMRES benchmark保留。所有数据在本地可重读，远端只使用隔离目录，不改全局环境。

table可按块恢复，协议/源码不符则拒绝恢复。test/Flow有上游gate guard。已有测试仅支持实现合同，不代表科研gate通过。

尚未完成且不计成功：锁定测试、全population可信参考、完整uncoupled加密、posterior physical branch识别、精确hidden-fiber margin、通过的两态模型、Holm确认比较、Laplace/mixture充分性、Flow、TAP novelty审查。没有用剩余预算重新寻找positive family。

## 13. 最终建议

保留本轮作为可复算的NO-GO/参考失败与对称性反例档案。不要升级模型大小、启动Flow或扩大positive叙事。如果后续重新立项，首先改善全支持后验积分并为reference单列预算；之后再判断更换保持/打破对称的物理family是否值得。新研究需新冻结协议，本轮阈值与test状态保持不变。


---

Source: [POSTERIOR_REFERENCE.md](https://raw.githubusercontent.com/migodam/a15-fullwave-posterior-ambiguity/main/A15/reports/POSTERIOR_REFERENCE.md)

# A15 后验参考验收（参考不足）

**状态：未通过。锁定测试未开启，Flow gate 关闭。**

两个全波表均完整生成：33³=35,937、41³=68,921 个合法材料节点，每节点三个频率，合计 314,574 次正演。没有删除失败节点重新归一化。源代码、协议及每个表的原始文件 SHA256 保存在各表 status.json；本地保留完整表。

独立连续 truth 的 validation 为三 seed ×128 queries。主噪声、基础 acquisition 的严格逐 query 比较只有 94/384 通过全部容差；最大 Q 均值差 0.11836、最大区间端点差 0.15384。新增 acquisition 同样未收敛。详见 `results/reference_validation_gate.json`。有效积分节点数不是 ESS 收敛证明。

已追加两个预定 query 的独立 tempered SMC（各两个独立 256 粒子运行，共 32,768 次全波求解），采用正确 conditional ESS、自适应温度、pCN 和祖先诊断。温度到达 1 不等于参考收敛；重复运行的方差差异仍不可忽略。

正权自适应 cubature 保留全部 prior 支持，对两个相同 query 分别使用约 9.7k 全波求解。每个叶 cell 使用 3³/5³ Gauss–Legendre 规则，所有叶节点及权重保留。没有只截取 MAP 周围区域。预算停止时嵌套 raw-moment 相对误差估计仍约 4.3% 和 8.6%；区间与方差也不能仅用这些估计认证。完整轨迹与节点位于 `adaptive_reference.json` / `adaptive_q*.npz`。

这是积分分辨率和有限预算的失败，**不是物理不存在的证据**。表上 population risk、prevalence、采集排序只作 validation 诊断，不能使用它们宣告 Gate A/B PASS。禁止开启 locked test，禁止训练 Flow。其余预算用于规定的 coupling 消融和数值/机制边界检查，不修改 prior 或阈值。


---

Source: [GATE_A.md](https://raw.githubusercontent.com/migodam/a15-fullwave-posterior-ambiguity/main/A15/reports/GATE_A.md)

# Gate A：未建立（参考不足）

**未通过执行 gate；物理存在性结论未决。** locked test 未开启。

主噪声 validation 41³ 的未认证诊断：R_Q=0.01292496，90%覆盖率=86.1979%，posterior-mean MSE=0.01345738；ambiguous query 为0/384，三个 seed 均为0。33³ 为1/384。它们既未接近10%门槛，也未通过 reference 收敛；不能把0/384当作真实总体零概率。Wilson 和 bootstrap 区间仅包括有限查询抽样误差，不包括积分偏差。

1/2/5%噪声和1/2/4/8通道的全部 validation 诊断均保存。5%基础 acquisition 的未认证比例上升，说明门槛结果对噪声敏感；不得据此把主噪声改成5%。

参考失败明细见 POSTERIOR_REFERENCE.md。没有只保留通过积分的94个 query 重算总体，也没有把它们当新 population。后续 positive claim 均阻断。


---

Source: [GATE_B.md](https://raw.githubusercontent.com/migodam/a15-fullwave-posterior-ambiguity/main/A15/reports/GATE_B.md)

# Gate B：FAIL（未满足通过条件；含未决参考）

当前证据不满足 bright–dark 定量解释与 population/acquisition 确认的联合要求。这是执行 gate 的 FAIL，不是“所有被动物理系统不存在任务歧义”的定理。

16个预定校准场景，质量归一化圆柱 resolvent 的两模态 Galerkin 投影，外部误差为11.78–43.49个实分量噪声 RMS，0/16小于噪声；Q误差2.07–4.82。固定空间 n64 VIE response-POD 两态投影只有2/16外部误差小于噪声，Q误差0.0277–0.3477，即使 snapshot current energy 保留很高也不够。后者受VIE离散偏差限制，不能作精密解释证书。其他两极点拟合未穷尽，不能从这两种失败推出任何两态模型都不可能。

关键对称性反例：几何关于 y→−y 对称，照明 exp(ikx) 是偶函数，三个参数都保持该对称。奇对称电流方向不能被该参数族激发；它们在 y=0 的基础接收器上不可见，固定圆形 ROI 的一阶能量变化也抵消。因此“很暗”可同时意味着 injection≈0、task coupling≈0。Fig11 将四项分开，而不是以最大 resolvent norm 当证据；余项大小在 modal_audit.json。

删除跨体耦合的33³表完整生成，保留各体self interaction。表上基础风险从0.013041降到0.00030439，邻频相对优化同频的优势消失；符合事前方向，但该反事实的积分参考没有独立加密认证，仅能作为诊断。它不是另一真实装置。物理大间隙0.4/0.8/1.2对照也计算，反馈减弱但在有限间距并未可忽略：平均|Q−Q_uncoupled|仍约2.16/1.45/1.18。不得把gap=1.2称为无耦合真值。

精确 hidden-fiber margin 条件未建立；代表性 posterior 10/90% Q 节点的测量 KL≈4.23，不是几乎不可区分对象对。它们来自未收敛表，展示目的，不构成分支或 prevalence 证据。


---

Source: [PHYSICS_VALIDATION.md](https://raw.githubusercontent.com/migodam/a15-fullwave-posterior-ambiguity/main/A15/reports/PHYSICS_VALIDATION.md)

# A15 物理与实现验收

**圆柱全波后端：限定数值检查通过；VIE n≤128：不能承担本轮低噪声 posterior truth。**

采用 exp(-iωt)、外向 Hankel-1、ε=εreal+iγωref/ω，γ≥0。单位 E_z 入射，归一化入射强度 |E0|²/(2Z0)，取 Z0=1；吸收截面除以该入射强度。Q 是固定目标圆盘 mean |E|² 的自然对数相对增益，可为负，不是完整色散储能或温升。

A14 复数 Green kernel 与 corrected equal-area-disk self-cell 原样隔离复用；高层实材料限制未沿用。圆—格点 cell overlap 使用连续面积积分，总面积误差≤6.94e-18。active/full 系统 J/E 相对差≤9.3e-16。零 contrast、Born 极限、自单元数值积分、独立圆柱解析散射与互易性已检查。

32/48/64/96/128 单圆柱外部误差依次为 6.07%、2.67%、1.43%、0.732%、0.416%；不能对该 VIE 网格宣称满足 1% 噪声容差。因此在正式确认前冻结解析多圆柱展开后端：解所有截断圆柱谐波的多次散射耦合，不是 two-mode toy。order20 vs24 的24场景最大噪声白化 RMS 为4.04e-10、最大 Q 差4.19e-12；被动能量平衡相对误差≤5.62e-15。

8 scenes×8 参数方向的中心差分、相邻步长及实参数 adjoint 配对通过；全部链式导数包含移动几何、材料本构、坐标质量归一化、RHS 与读出。几何 Jacobian 的参考不是只改变 ε 的导数。44 项测试通过，记录 `results/all_tests.log`；其中10项是物理/协议合同检查，34项是后验数值工具检查。

局限：圆柱解析解只适用于当前均匀圆形 E_z 标量模型；不验证3D vector Maxwell。n64 VIE 模态图保留明确的空间离散偏差，不能用其误差大小替代圆柱 reference。互易性检查与单圆柱解析核查提供不同于复用 Green 代码的验证路线。


---

Source: [ACQUISITION_COMPARISON.md](https://raw.githubusercontent.com/migodam/a15-fullwave-posterior-ambiguity/main/A15/reports/ACQUISITION_COMPARISON.md)

# 等预算采集比较（未认证 validation）

每个新增方案恰好增加两个复数值；所有通道使用同一绝对实/虚噪声 σ。共同 truth 与基础噪声，384独立 query，5,000次完整 query 配对 bootstrap。粒子、积分节点、通道不作为重复试验。

| 方案 | 表上 R_Q | 相对基础风险下降 |
|---|---:|---:|
| base1 | 0.0209962 | -62.45% |
| base | 0.0129250 | 0.00% |
| same_angles | 0.0098180 | 24.04% |
| optimized_angles | 0.0057502 | 55.51% |
| neighbor | 0.0038403 | 70.29% |
| near | 0.0006875 | 94.68% |
| random | 0.0062992 | 51.26% |
| base8 | 0.0056772 | 56.08% |

邻频相对优化同频多降14.78个百分点，近场多降39.17个百分点；三个 seed 方向一致。但严格 reference 未通过，因此这些不是确认性效应；不宣告 Holm 检验通过。原始配对95%与保守97.5%区间保存在 population_diagnostics.json，均未计积分偏差。

discovery 事前 proxy 排序为 near→neighbor+→optimized angle，与表上排序一致。这只支持一个待确认的全波局部设计线索；不是两模态理论预测成功。邻频需要额外频点和相应源/设备成本，近场需要布设接近材料的传感器；相同复通道预算不等于相同工程成本。固定近场点在整个 prior 支持外部。


---

Source: [VISUAL_REVIEW.md](https://raw.githubusercontent.com/migodam/a15-fullwave-posterior-ambiguity/main/A15/reports/VISUAL_REVIEW.md)

# A15 visual review（逐图实际打开检查）

2026-09-26，Codex逐张打开11幅PNG；改动后重新打开关键图。下述是实际视觉观察，不是仅文件存在检查。所有后验/风险图均有PROVISIONAL标记，不提供未验证分支图。

| 图 | 实际观察 | 解释/下一步 |
|---|---|---|
| 1 geometry | 两圆几何、固定目标ROI及两基础接收器可辨；近场位于目标外；初版inset遮住右基础接收器，已移到中央空白区并重看 | 不随样本移动ROI，不把近场位置算成额外隐藏观测 |
| 2 spectrum | 冻结k0处目标Q靠近一个深谷；不是内部能量高峰。初版说明文字落在第二幅图内，已移到图外 | 不能用“强暗模增强”描述选择结果；81点细化不回改k0 |
| 3 modes | 两个小奇异方向具不同反射对称结构，放大系数并没有独立证明任务相关性 | 结合Fig11注入/可见性；n64有离散偏差，非完整模态证书 |
| 4 parameter map | gap强烈改变Q；gamma方向变化很弱；外部幅值与Q有明显关联而不是大面积横向分离 | 粗9³图只说明选型局限，不数object pairs；参数轴用真实GL坐标 |
| 5 theta posterior | gap集中在窄带，gamma近乎宽分布；节点分辨率造成可见条纹/空bin | 有效积分节点数不是收敛证据，不把条纹解释成物理多分支 |
| 6 Q posterior | 单个宽而起伏的直方图，不足以断言多峰；Q为负 | 标为未收敛后验；不将log能量的CV当物理能量CV |
| 7 explicit pair | 两场图由同一冻结family精确solver生成；目标内部较暗。10/90%节点能量比约1.71，但测量KL约4.23 | 不是强不可区分对象对，也不是两个已识别分支；幅值曲线近似重叠不等于复测量等价 |
| 8 risk distribution | 主2%噪声下384个query的表上SD均未跨阈值；多方法是同一query重复评估 | 不把方法数增加样本量；不能从未收敛表判定真实总体零发生率 |
| 9 acquisition | near最低、neighbor次之；5%边界的不确定性明显更大 | 不把主噪声换成5%追逐正面结果；不以bootstrap掩盖积分偏差 |
| 10 coupling | 移除cross-block后risk大降；图中两组同为33³规则 | 反事实不是完整物理装置，且未独立加密；不能给正式因果证书 |
| 11 four factors | 最暗方向常同时注入与任务耦合接近机器零；相邻可注入方向有更大泄漏 | 反射对称性提供可解释反例；必须保留注入和余项，不只画1/σ |

绘图实现审查发现并修复：JSON列表/字典schema不兼容导致5–8图失败；代表性query曾混用seed/noise；噪声未加回观测；field pair曾用默认半径而非冻结family；模态σ(M)曾误标resolvent奇异值；uncoupled与coupled图曾用不同积分阶。修复后加上query均值/方差回归测试及field/table Q一致性断言。首轮失败manifest保存在results/figures_first_run_failures.json。没有将这些作图错误生成的图用于科学裁决。

图片仍是诊断图而非论文定稿。预算内没有建立物理分支，故没有伪造两分支/Flow样本图。


---

Source: [COST_AND_PROVENANCE.md](https://raw.githubusercontent.com/migodam/a15-fullwave-posterior-ambiguity/main/A15/reports/COST_AND_PROVENANCE.md)

# A15 cost / provenance

主档在本地A15；XINAN只使用独立D:\AI\A15_BRIGHT_DARK_FULLWAVE。直连LAN，ProxyCommand/ProxyJump显式none；未改全局环境、原项目或A14。无需CUDA新安装。

| 已持久化工作 | 全波forward次数 | wall time |
|---|---:|---:|
| coupled33³（三频） | 107,811 | 194.67s |
| coupled41³（三频） | 206,763 | 385.45s |
| cross-block removed33³（三频） | 107,811 | 97.27s |
| SMC两query×两独立运行 | 32,768 | 139.03s |
| 正权adaptive两query | 19,387 | 5.28s |
| discovery主扫描+切线基点 | 约7,380 | discovery/summary.json |
| 三个候选9³×三频+screening | 6,945 | 原始screening表和日志 |
| 校准1024+设计切线基点 | 1,120 | freeze.log |
| 验证query coupled/uncoupled | 2,304 | 各summary.json |
| large-gap基点/切线 | 384 | large_gap.json |

计入物理验证和模态检查，主运行账目重建约494,424次；此数不含少量benchmark、checkpoint测试、绘图复算及已中止worker的探索性调用。它不是从每次启动起完整捕获的全生命周期精确计数。接近50万预算后停止扩展积分/搜索；未宣称精确global-call预算证书。SMC连续随机参数缓存命中少，而规则积分会复用单圆柱系数，wall time差异不代表SMC遗漏了物理求解。

table33/41曾并发运行，不能把两wall time相加当总elapsed。原始PID28244/12796已正常退出；第一个尝试的Windows detached进程在SSH结束后退出，未完成节点，记录保留。随后维持SSH父进程到结束；uncoupled也正常结束，所有表完整复制回本地并逐文件验证hash。

## 后端成本边界

50×50 pilot：CPU LU约62.7μs，GMRES302.6μs，CUDA单个853.4μs，batch32约71.1μs/系统。最终82×82矩阵核复核：CPU130.3μs，GMRES323.0μs，CUDA单个772.8μs，batch32为82.5μs/系统。**最终规模batch GPU矩阵核比CPU快，但该数排除了完整varying-θ装配和持续数据搬移**；不能宣称CPU在所有batch实现中绝对最快。本轮选择已验证的CPU装配/表路径，完整缓存场景约227μs；完整表实际吞吐由三频装配与I/O决定。

这是执行计划的一项局限：没有实现/验证全表GPU流水线，也没有从第一步记录每个阶段真实CPU/GPU峰值内存；峰值记为未采集，不伪造为0。磁盘与版本hash在RUN_MANIFEST.json。MPS未训练，CUDA只做benchmark。

## 验证与worker

44项测试通过，包括正权积分、CESS、SMC、噪声约定、互易性、单圆柱、self-cell、protocol、恢复、篡改拒绝和图中同一query一致性。test/Flow guard拒绝上游gate关闭时执行。

DeepSeek Flash承担普通posterior工具和画图初稿；Codex纠正CESS、核验公式、修复绘图实际数据错误并独立裁决。物理QA worker因范围扩张和family混用停止；其未核验结论未纳入报告。日志、parent reviews及输出路径保留。理论理解与最终判据未委托worker。

全部理论原件及A14来源与初始SHA256重新核对。运行生成器、源文件、配置、原始表文件、门控、报告和图的清单在RUN_MANIFEST.json；恢复时须先完整阅读A15_prompt.md，摘要不能替代原文。


---

Source: [LITERATURE_AUDIT.md](https://raw.githubusercontent.com/migodam/a15-fullwave-posterior-ambiguity/main/A15/reports/LITERATURE_AUDIT.md)

# 文献核验状态：不完整

使用 scholarqa-research 技能与 scholarqa-cli 进行了有界检索。Semantic Scholar snippet 与 paper 端点均返回 HTTP429；原始错误保留在 research/delegated/literature/retrieval.json。规划阶段另两次429也没有提供可引用证据。

本轮没有取得可验证的近邻论文证据，不能断言首创、没有相关工作或满足TAP novelty。A15.md 的 content-reference 占位符未转换成伪引用。有限预算留给数值证伪；文献审查不足单独列为 TAP 评估限制。用户原理论是待检验输入，不是外部学术证据。
