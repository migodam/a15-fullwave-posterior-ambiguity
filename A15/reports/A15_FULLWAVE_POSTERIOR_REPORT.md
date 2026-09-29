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
