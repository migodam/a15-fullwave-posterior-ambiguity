# A15 — Full-Wave Bright–Dark Task Ambiguity and Acquisition Design

你现在负责执行下一轮 full-wave 研究。

这不是继续 A14 的 Hybrid latent competition，也不是立即恢复 Flow Matching。

本轮核心问题只有一个：

\[
\boxed{
\textbf{
Can a passive full-wave electromagnetic system possess a
non-negligible-prior-mass family of objects that are nearly
indistinguishable to the current external measurement but have
substantially different internal task values?
}}
\]

进一步：

\[
\boxed{
\textbf{
Can the theory predict a new measurement that specifically breaks
this hidden family more efficiently than simply adding more similar
measurements?
}}
\]

只有这个物理问题被正面建立之后，才允许讨论 Flow。

---

# 0. 必须先读的文件

开始代码前先找到并完整阅读：

```text
A15.md
```

这是本轮理论研究结果。

同时读取：

```text
A14/A14.md
A14/Q14.md
A14/reports/A14_FINAL_REPORT.md
A14/reports/GATE_A.md
A14/src/physics.py
A14/src/physics_backend.py
```

以及 A9/A12 中与：

```text
material-current tangent
current transfer
receiver filtering
resolvent / hidden current modes
```

直接相关的最终理论文件。

不要修改已有 A14 数据和报告。

---

# 1. 新目录

本地创建：

```text
Flow Matching/A15/
```

建议结构：

```text
Flow Matching/A15/
    context/
    configs/
    src/
    experiments/
    results/
    figures/
    reports/
    tests/
```

复制只读理论 snapshot 到：

```text
Flow Matching/A15/context/
```

至少包括：

```text
A15.md
A14_FINAL_REPORT.md
```

保存 SHA256。

---

# 2. 远程算力

本轮允许并鼓励使用局域网直通：

```text
ssh xinan
```

调用 XINAN 上已有：

```text
RTX 4060
i9-14900 class CPU
CUDA environment
```

禁止走代理。

不要修改 XINAN 全局环境或已有项目。

远端新建独立工作目录，例如：

```text
D:\AI\A15_BRIGHT_DARK_FULLWAVE
```

具体 Windows/WSL/Python 环境先检测已有配置，再使用最稳定的现有环境。

本地 repo 保持 source of truth：

```text
Flow Matching/A15/
```

大型 raw tables/checkpoints 可以保留远端，但必须：

- 有 manifest；
- 有 config hash；
- 有 seed；
- 有代码 commit/hash；
- 有生成脚本；
- 将 summary JSON/CSV/figures/reports 同步回本地。

---

# 3. GPU 的第一用途不是 Flow

RTX 4060 首先用于：

```text
batched / accelerated full-wave solves
frequency sweeps
parameter sweeps
posterior-table generation
modal/resolvent analysis
```

如果现有 CPU complex128 solver 比 GPU 更可靠或实际更快：

可以继续用 CPU。

先 benchmark：

```text
CPU dense
GPU dense/batched
existing iterative implementation
```

然后选 wall-clock 最优且数值可靠的方案。

不要为了“必须用 GPU”而重写一套不可靠 solver。

---

# 4. 物理模型

第一阶段允许继续使用小型 2D full-wave scalar Helmholtz / volume-integral formulation 做机制证伪。

但必须从 A14 的 real signed contrast 升级成一个**明确被动的谐振散射系统**。

优先候选：

```text
two nearby passive dielectric resonators
```

例如两个圆形或椭圆介质体形成 bright/dark hybridized modes。

必须有明确：

```text
time convention
complex permittivity convention
passivity sign
input-power normalization
energy / absorption sanity checks
```

禁止把 gain medium 或非被动 contrast 当作 positive result。

---

# 5. 三个核心物理参数

冻结一个三维参数 family：

\[
\boxed{
\theta=(s,\Delta,\gamma)
}
\]

其中建议：

### \(s\)

resonator gap / separation / geometry coupling parameter，

主要控制 bright–dark coupling。

### \(\Delta\)

dark-resonator detuning parameter。

可以通过：

```text
radius
aspect ratio
real permittivity
```

其中一个实现。

第一版只选一种，不要同时让多个参数表达同一 detuning。

### \(\gamma\)

dark resonator loss parameter，例如：

```text
Im epsilon
loss tangent
conductive loss
```

必须满足被动性。

---

# 6. Discovery 与 confirmatory prior 必须分开

允许先在一个较宽 discovery box 中搜索：

```text
bright/dark resonant regime
```

这只是选 canonical geometry 和中心频率。

Discovery 数据不得作为最终 prevalence 证据。

找到候选后：

冻结：

```text
geometry family
theta ranges
prior
omega0
noise
measurement channels
ROI
task
success thresholds
```

写入：

```text
EXPERIMENT_PROTOCOL.md
configs/protocol.json
```

然后重新生成独立 validation / locked-test 数据。

不要在看到 validation 结果后继续移动 prior box。

---

# 7. Prior

主实验优先使用独立 bounded continuous prior：

\[
s\sim U[s_-,s_+],
\]

\[
\Delta\sim U[\Delta_-,\Delta_+],
\]

\[
\gamma\sim U[\gamma_-,\gamma_+].
\]

如果物理需要非均匀 prior，可以使用，但必须在协议冻结。

报告：

```text
prior width
physical interpretation
fraction of broad discovery prior occupied by frozen regime
```

不要只报告“在精选 resonance box 内 prevalence 很高”。

还要报告该 box 相对于更广 physical prior 的质量。

---

# 8. 找 bright/dark mode 的标准

不能只看一张 field plot。

至少使用两种独立证据。

建议：

## A. Full-wave current / field mode decomposition

检查响应是否主要由两个 resonant modes 主导。

可以使用：

```text
resolvent SVD
current covariance/SVD
rational pole fit
quasi-normal mode fit
```

选择最稳定的一种或两种。

## B. Radiation contrast

识别：

```text
bright mode:
large external coupling

dark mode:
large internal energy/current
small external coupling
```

## C. Two-mode reconstruction

尝试拟合：

\[
a_b=
\frac{\kappa s}
{b+g^2/d}.
\]

不要求符号完全采用理论 toy 的记号，

但需要一个可识别的 two-mode reduced model。

---

# 9. Two-mode explanation 的 Gate

为了把结果解释成 bright–dark mechanism，而不只是黑箱歧义，

要求在 confirmatory regime：

external response reduced-model error 小于：

```text
measurement noise scale
```

并且 preferably：

\[
E_{\rm CM,F}
<
0.25\times
\text{hidden-fiber measurement margin}.
\]

ROI task reduced-model error 也应显著小于：

```text
high-Q vs low-Q branch separation
```

如果 reduced model 不够准：

物理歧义本身仍可继续研究，

但不能宣称由当前 coupled-mode mechanism 已解释。

---

# 10. 主 measurement

主 measurement 使用：

\[
\boxed{\text{single-frequency external complex response at }\omega_0}
\]

不要一开始使用完整 360° 无限数据。

选择一个现实但有限的 acquisition。

优先：

```text
2 complex external channels
```

例如：

```text
forward/back scattering
```

或两个冻结角度的 complex receivers。

同时做 acquisition-budget boundary：

```text
1 complex channel
2 complex channels
4 complex channels
8 complex channels
```

但提前冻结 channel positions。

---

# 11. 主 task

第一主任务只定义一个：

\[
\boxed{
Q(\theta)
=
\log
\left[
\frac{
\int_{\Omega_{\rm ROI}} |E(\mathbf r)|^2\,d\mathbf r
}{
\int_{\Omega_{\rm ROI}} |E_{\rm inc}(\mathbf r)|^2\,d\mathbf r
}
\right]
}
\]

ROI 在协议中提前冻结。

优先：

```text
dark resonator volume
```

或：

```text
dark resonator + fixed gap ROI
```

不要按每个 sample 的热点位置动态选 ROI。

---

# 12. Secondary tasks

可以记录但不作为主 gate：

```text
ROI current energy
peak-smoothed field
total stored-energy proxy
local absorption if material model supports it
```

如果研究 absorption：

必须使用真正 lossy constitutive model。

禁止把：

```text
|J|^2
```

直接叫 heating。

---

# 13. Frequency sweep

对于每个 parameter state，除：

\[
\omega_0
\]

外至少计算：

\[
\omega_0-\Delta\omega,
\qquad
\omega_0+\Delta\omega.
\]

Discovery 阶段允许较密 sweep：

```text
21–81 frequency points
```

以识别 resonance/poles。

Confirmatory posterior table 不需要每个点这么密。

至少保存：

```text
omega0
omega0 - delta_omega
omega0 + delta_omega
```

并可另外保存少量 intermediate frequencies。

---

# 14. Posterior reference：优先利用只有 3 个参数这件事

本轮不要先用 neural posterior sampler。

因为：

\[
\dim\theta=3.
\]

这是优势。

建立高精度 reference posterior。

优先方案：

### structured Gauss–Legendre / adaptive tensor integration

例如先：

```text
33 × 33 × 33
```

必要时升级：

```text
41 × 41 × 41
```

或：

```text
49 × 49 × 49
```

每个节点预计算：

```text
y(omega0)
y(omega0-dw)
y(omega0+dw)
Q
secondary readouts
modal diagnostics
```

一个 table 可以服务很多 synthetic observations。

如果 tensor table 成本太高，

使用：

```text
adaptive sparse grid
Sobol prior table
tempered SMC
```

但必须给 convergence check。

---

# 15. 大规模计算预算

XINAN 可以承担大型 table。

初始目标：

```text
3e4 – 1e5 unique parameter states
```

和：

```text
3–9 relevant frequencies/state
```

根据 benchmark 自动调整。

总 full-wave solve 数允许到：

```text
1e5 – 5e5
```

量级，

只要：

```text
wall-clock
numerical residual
memory
checkpointing
```

合理。

不要为了追求数字大而浪费求解。

---

# 16. Posterior likelihood

主 measurement：

\[
y=F_{\omega_0}(\theta)+n,
\]

\[
n\sim\mathcal N_{\mathbb C}(0,\Gamma_n)
\]

或等价 Re/Im Gaussian 表示。

明确 noise normalization。

至少测试：

```text
1%
2%
5%
```

但选择一个 pre-registered primary noise level。

---

# 17. 不再数 object-pair witness

本轮的主量是：

\[
\boxed{
p(Q|y)
}
\]

和：

\[
\boxed{
\operatorname{Var}(Q|y).
}
\]

对每一个 synthetic observation：

保存：

```text
posterior mean Q
posterior std Q
credible interval
posterior branch weights if multimodal
posterior theta covariance
posterior predictive external response
```

---

# 18. Population-level risk

从 prior predictive 生成：

```text
512–2048
```

个 independent truth/noise queries。

对每一个求 reference posterior。

定义：

\[
R_Q
=
\mathbb E_Y
\operatorname{Var}(Q|Y).
\]

同时报告：

\[
P_Y
\left[
\operatorname{SD}(Q|Y)>\tau_Q
\right].
\]

阈值：

\[
\tau_Q
\]

必须在 locked test 前冻结。

可额外报告：

```text
posterior CV
95% interval width
high/low branch probability
```

---

# 19. 真正的 prior-mass Gate

建议 primary Gate：

### Gate A1 — task uncertainty

至少一个有实际意义的 query subset 满足：

\[
P_Y[
SD(Q|Y)>\tau_Q
]
\ge 0.10
\]

作为最低有质量门槛。

如果 discovery 结果支持，可以冻结更严格：

```text
0.15 or 0.20
```

但不能看 test 后改。

### Gate A2 — effect size

在 ambiguous query 中，

posterior 90% credible interval 对应至少：

```text
1.5× ROI-energy ratio
```

preferably：

```text
2×
```

或等价 log-energy separation。

### Gate A3 — physical support

全部 posterior mass 来自合法：

```text
passive
geometry-valid
solver-stable
```

对象。

---

# 20. Explicit fiber visualization

从 representative posterior 中选取：

```text
high-Q branch
low-Q branch
```

但这只是展示，不是统计单位。

画：

```text
geometry
material
field
current
external measurement
frequency spectrum
```

并显示：

\[
F_{\omega_0}(\theta_a)
\approx
F_{\omega_0}(\theta_b)
\]

而：

\[
Q(\theta_a)\ne Q(\theta_b).
\]

---

# 21. 核心消歧实验

在相同 additional measurement budget 下比较三种方案。

## Baseline A — more same-frequency external angles

在：

\[
\omega_0
\]

增加额外 receivers。

## Targeted B — neighboring frequency

用相同数量的 complex channels，

但放在：

\[
\omega_0+\Delta\omega
\]

或最优的一个预冻结邻频。

## Targeted C — near-field channel

如果模型允许，

使用一个或等预算 near-field channel。

---

# 22. measurement budget 必须公平

例如增加：

```text
2 complex values
```

则三种方法都只能增加：

```text
2 complex values
```

noise level 相同。

不能让 multi-frequency 方法偷偷获得更多 real measurements。

---

# 23. acquisition-design 主指标

对新增 measurement \(Z\)，

计算：

\[
R_Q^{\rm new}
=
\mathbb E_{Y,Z}
\operatorname{Var}(Q|Y,Z).
\]

报告：

\[
\boxed{
\text{risk reduction}
=
1-
\frac{R_Q^{\rm new}}{R_Q^{\rm base}}
}
\]

比较：

```text
same-frequency extra angles
neighbor frequency
near-field
random additional channel
```

如果理论正确，

targeted measurement 应显著优于 equal-budget generic channel。

---

# 24. 如果可能，计算 goal-oriented information gain

作为 secondary metric：

\[
I(Q;Z|Y).
\]

不要只计算：

\[
I(\theta;Z|Y).
\]

因为最终问题是 task risk，而不是 parameter entropy。

---

# 25. coupling ablation

必须有一个控制实验区分：

```text
普通有限采样欠定
```

和：

```text
bright-dark internal feedback mechanism
```

优先做：

### Ablation 1

在离散 Green interaction 中去除两个 resonator 之间的 cross-body coupling block，

保留各自 self interaction。

或者：

### Ablation 2

把 gap 增大到 coupling 可忽略但其他设置尽量保持。

比较：

```text
posterior task risk
modal structure
targeted-frequency advantage
```

如果 ambiguity 完全不受 coupling 影响，

bright–dark explanation 很可能不成立。

---

# 26. modal sufficient-condition audit

对于 high-risk scenes / posterior regions，

检查：

\[
M^{-1}C_sh
\]

的 modal decomposition。

记录：

\[
|u_k^\ast C_sh|
\]

材料 injection；

\[
\sigma_k
\]

resolvent amplification；

\[
\|Sv_k\|
\]

measurement visibility；

\[
\|L_Qv_k\|
\]

task visibility。

寻找是否存在明确：

\[
\boxed{
\text{injectable}
+
\text{amplified}
+
\text{measurement-dark}
+
\text{task-bright}
}
\]

mode。

---

# 27. 不允许只用最大 resolvent norm

禁止把：

\[
\|M^{-1}\|\gg1
\]

直接解释为 hidden ambiguity。

必须同时报告：

```text
injection coupling
receiver leakage
task coupling
residual modes
```

因为 resonance 也可能让 measurement 更容易看到对象。

---

# 28. Gate B — Maxwell mechanism verdict

定义：

### PASS

存在有非小 prior-predictive mass 的 full-wave task ambiguity，

且：

```text
bright-dark / resolvent mechanism
```

能够解释它；

并且 predicted targeted measurement 显著降低 task posterior risk。

### PARTIAL

真实 ambiguity 与 posterior risk 存在，

但 coupled-mode explanation 不够准确；

或：

theory explains it，但 prevalence 很窄。

### FAIL

只有 isolated object pairs / local tangent mismatch，

没有显著 posterior mass；

或 targeted acquisition 没有理论预言优势。

---

# 29. Flow Gate

只有 Gate B PASS 后才讨论 Flow。

先评价 posterior complexity。

对 representative queries：

比较：

```text
Laplace
1–4 component Gaussian mixture
direct quadrature/reference
tempered SMC
```

检查：

```text
number of physical posterior branches
branch probabilities
non-Gaussian skew/kurtosis
curved fibers
task proper scores
```

---

# 30. Flow 不必要的判据

如果：

```text
Laplace / tiny Gaussian mixture
```

已经能够以低成本得到正确：

```text
Q posterior mean
Q variance
coverage
branch probability
```

则明确结论：

\[
\boxed{\text{Flow is unnecessary for this canonical problem.}}
\]

不要为了项目名字强行训练。

---

# 31. 若 Flow Gate 通过

如果 posterior：

```text
multimodal
strongly curved
non-Gaussian
```

且 reference sampling 对大量 queries 昂贵，

才训练 Conditional Flow Matching。

Flow state 第一选择：

\[
\boxed{
\theta=(s,\Delta,\gamma)
}
\]

的真实 physical parameters。

不是：

```text
pixel chi
CNN latent
arbitrary learned embedding
```

---

# 32. parameter transform

如果 prior 是独立 bounded uniform，

优先先做 exact CDF transform：

\[
u_i=
\Phi^{-1}
\left[
F_{\theta_i}(\theta_i)
\right].
\]

使：

\[
u\sim N(0,I).
\]

Flow 在：

\[
u
\]

上运行。

decode：

\[
u\rightarrow\theta\rightarrow\chi\rightarrow\text{solver}.
\]

因此所有 samples 自动在 physical prior family 内。

---

# 33. Flow baselines

若真的进入 Flow：

比较：

```text
Laplace
Gaussian mixture
conditional normalizing flow
conditional FM
```

训练数据可以直接来自 simulator joint samples：

\[
(\theta,y).
\]

reference posterior 主要用于独立评价，

不是必须用 posterior particles 训练 FM。

---

# 34. Flow metrics

必须评价 posterior，而不是只评价 sample aesthetics：

```text
energy score
coverage
branch probability
Wasserstein / sliced Wasserstein
task posterior mean/std
posterior predictive check
```

尤其：

\[
p(Q|y)
\]

必须和 reference 匹配。

---

# 35. xinan 上允许的大训练

只有进入 Flow Gate 后，

可以使用 RTX 4060 训练：

```text
1–10M parameter conditional FM
```

数据规模：

```text
50k smoke
200k main
500k optional
```

根据已有 full-wave table 与 simulator 成本决定。

不要一开始就生成 500k。

---

# 36. Flow 的真正成功条件

不是：

```text
samples look diverse
```

而是：

\[
\boxed{
\text{在相同 task-posterior accuracy 下，
amortized query cost 明显低于 reference inference}.
}
\]

计算：

\[
C_{\rm offline}
+
N_{\rm query}C_{\rm FM}
\]

和：

\[
N_{\rm query}C_{\rm reference}.
\]

求 break-even：

\[
N_{\rm query}^\star.
\]

---

# 37. 主要 figures

至少准备：

### Fig 1
canonical two-resonator geometry + measurement + ROI。

### Fig 2
frequency spectrum，标出 bright/dark resonances。

### Fig 3
two dominant modes：
field/current + radiation visibility。

### Fig 4
external response vs internal ROI energy parameter map。

### Fig 5
同一/similar single-frequency measurement 下 posterior in \((s,\Delta,\gamma)\)。

### Fig 6
posterior \(p(Q|y)\)。

### Fig 7
两个 physically distinct posterior branches：
similar external measurement / different internal energy。

### Fig 8
prior-predictive distribution of posterior task risk。

### Fig 9
same-frequency additional channels vs neighboring frequency vs near-field 的 task-risk reduction。

### Fig 10
coupling ablation。

### Fig 11
material injection / resolvent amplification / receiver darkness / task brightness modal decomposition。

只有进入 Flow 后再增加 posterior-sampler figures。

---

# 38. Final report

生成：

```text
Flow Matching/A15/reports/A15_FULLWAVE_POSTERIOR_REPORT.md
```

必须回答：

## Q1
full-wave Maxwell/Helmholtz 中是否真的存在具有非小 prior mass 的 external-equivalent / internal-different family？

## Q2
posterior internal-task risk 有多大？

## Q3
该现象主要来自：
- ordinary undermeasurement，
- resonance，
- bright-dark coupling，
- receiver filtering，
- 还是它们组合？

## Q4
A15 的 coupled-mode mechanism 能否定量解释 full-wave data？

## Q5
哪一种额外 measurement 最有效降低 task risk？

## Q6
理论是否能在实验前预测这一 measurement？

## Q7
这种 ambiguity 的 population/prevalence 是否足以支撑独立 TAP mechanism paper？

## Q8
这条结果与 A9/A12 material-current transfer 是：
- natural extension，
- independent paper，
- 还是只适合 discussion？

## Q9
posterior 是否复杂到需要 Flow？

## Q10
若不需要 Flow，明确写：
```text
FLOW NOT NEEDED
```

若需要，再报告 FM 结果与 amortization break-even。

---

# 39. 最终 verdict

必须从以下选择：

```text
STRONG GO
GO
GO WITH REFRAMING
NO-GO
```

并分别给：

```text
physics mechanism verdict
posterior-risk verdict
acquisition-design verdict
Flow verdict
TAP-paper verdict
```

---

# 40. 研究纪律

不要因为解析 coupled-mode toy 是正面结果，就主动寻找只支持它的 full-wave scene。

主动找：

```text
no dark mode
bright leakage
extra radiation channel
posterior unimodal
targeted frequency no advantage
strong resonance but no hidden task risk
```

这些反例。

如果 full-wave 不支持 toy mechanism：

明确杀掉。

如果支持：

再扩大 prior 和 geometry family。

不要在第一轮成功后立刻宣称普遍 Maxwell 定理。

---

# 41. 最核心目标

本轮不是：

\[
\text{train a better Flow}
\]

而是：

\[
\boxed{
\textbf{
Show whether external electromagnetic data can determine an effective
scattering combination while leaving a non-negligible-prior-mass
internal-energy ambiguity, and whether a theoretically targeted
measurement can specifically remove that ambiguity.
}
}
\]

如果这条成立，

再问：

\[
\boxed{
\text{Is the remaining posterior sufficiently complex and repeatedly
queried that Flow Matching is actually useful?}
}
\]

这个顺序不能反。