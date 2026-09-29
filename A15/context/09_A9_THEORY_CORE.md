# A9_THEORY_CORE
## 材料切向激励、SOM/TSOM 模态失配与降阶逆更新保真
版本：2026-09-16。范围：以固定频率、固定几何的三维矢量 Maxwell **离散系统**为主；第4A节另给一个精确连续 Maxwell 构造及独立证明。多频可按下述“非共用 G”规则扩展。

本文的 `proved` 指在明确列出的数学假设下，本文给出了证明；不表示该命题首次出现于文献。连续 Maxwell、硬件噪声和具体离散化的一致性，不由有限矩阵证明自动获得。文献定位见配套 `A9_THEORY_CLAIM_LEDGER.md`。

**结论摘要。** Current/data observability 的排序一般不等于物理材料切向的测量能量排序。正确权重是 \(\tau_i=\sigma_i^2q_i\)，其中 \(q_i\) 是真实材料坐标、实际照明和全波重散射共同产生的切向激励量。该权重对固定 SOM 基底的能量选择精确最优，但不是对正则化更新的普适最优权重。真正的更新保证还需要实际降阶动力学误差、残差误差和逆问题稳定裕量。少量随机探针可以给出有置信水平的选择/拒绝与误差证书；不能保证任意问题在 64 个探针内获得明确排序。

---

# 0. 约定、边界与来源

除明确标出的连续 Maxwell 第4A节外，全文默认有限维实 Hilbert 空间。符号 \(M^*\) 表示当前坐标中正确的伴随：实正交坐标下为转置；只有在真正复线性且使用 Hermitian metric 时才为共轭转置。范数 \(\|\cdot\|_2\) 为算子范数，\(\|\cdot\|_F\) 为 Frobenius 范数。未标下标的向量范数是相应正交坐标的 Euclidean 范数。

材料参数、允许切空间、材料 metric、白化、照明权重和当前投影基底，在一次 GN/TR 线性化内固定。所有“最优”均相对于这些固定选择。改变 metric 或 tangent covariance 会改变 relevance，不能跨定义比较。

源码来源边界：本次连接可读取的旧 DDA 文件为
`research/a4_reliability_v1/src/dda.py`，其中明确使用 CM 极化率与辐射反作用修正。没有取得 A8 的精确极化率实现和旧 `eta_mix` 定义。因此下文给出必须由 A8 适配器实现的 \(a(\chi),Da(\chi)\) 合同，并给已核对的旧代码公式作实例；不宣称已审计 A8。旧代码把模型称作独立小型向量 DDA，而非连续模型证书。

---

# 1. 物理材料、电流、metric 与最终 \(K_p\)

## 定义 1.1：Maxwell contrast-source 系统

采用 \(e^{-i\omega t}\)，非磁性背景 \(\epsilon_b,\mu_b\) 已知。令
\[
 \chi=(\epsilon-\epsilon_b)/\epsilon_b,\qquad
 j=\chi E,\qquad J_{\rm pol}=-i\omega\epsilon_b\,j.
\]
\(j\) 是归一化极化电流变量，不是独立于材料的自由电流。对 \(N\) 个矢量单元，\(j_p,E_p\in\mathbb C^{3N}\)。离散 self term、单元体积和 \(k_b^2\) 因子必须一次性固定在 \(G_D,G_S\) 或电流变量中，不能重复计入。
\[
 E_p=e_p^{\rm inc}+G_Dj_p,\quad
 j_p=X(\chi)E_p,\quad
 L j_p=X e_p^{\rm inc},\quad L=I-XG_D.                 \tag{1.1}
\]
各向同性时 \(X=\operatorname{diag}(\chi_\ell)\otimes I_3\)。各向异性时是局部材料张量块。假设当前 \(L\) 可逆。测量为 \(f_p=G_S^{\rm phys}j_p\)，残差符号固定为 \(r_p=f_p-y_p\)。

### 引理 1.2：物理材料导数
固定几何、背景、频率、入射源和接收算子，对任意物理材料扰动 \(\delta\chi\)，
\[
 D j_p[\delta\chi]=L^{-1}\big(DX[\delta\chi]E_p\big).   \tag{1.2}
\]
**证明。** 对 \(Lj_p=Xe_p^{\rm inc}\) 求导：
\(L\,Dj_p-DX[\delta\chi]G_Dj_p=DX[\delta\chi]e_p^{\rm inc}\)。移项即得。无需 \(X\) 可逆，因此适用于 \(\chi=0\)。若几何、频率或入射源也是未知量，必须另加其导数项，不能继续使用只有 (1.2) 的公式。□

各向同性的 \(D_{E_p}\delta\chi:=DX[\delta\chi]E_p\) 是从 \(N\) 个材料自由度到 \(3N\) 个电流自由度的注入，逐单元为 \(\delta\chi_\ell E_{p,\ell}\)。不能误写成具有 \(3N\) 个独立材料分量的任意对角扰动。

## 定义 1.3：实化和噪声白化

\[
 \mathfrak r(v)=\binom{\Re v}{\Im v},\qquad
 \mathfrak R(M)=
 \begin{pmatrix}\Re M&-\Im M\\\Im M&\Re M\end{pmatrix}. \tag{1.3}
\]
有 \(\mathfrak r(Mv)=\mathfrak R(M)\mathfrak r(v)\)。
真实材料坐标取 \(x=(\Re\chi,\Im\chi)\in\mathbb R^{n_x}\)，或指定物理参数的实向量。若只允许实部变化，则只有相应列。

**重要区别：** 一个 \(D_c\in\mathbb C^{3N\times d}\) 若作用于实 \(z\in\mathbb R^d\)，其真实 Jacobian 是
\[
 \mathfrak V(D_c)=\binom{\Re D_c}{\Im D_c}\in\mathbb R^{6N\times d},
\]
不是 \(\mathfrak R(D_c)\)，后者错误地允许 \(2d\) 个输入方向。

令真实噪声 covariance 为 \(\Gamma_y\succ0\)，取 \(W_y\Gamma_y W_y^*=I\)。对 proper complex noise、Hermitian covariance \(\Gamma_c=\mathbb E[nn^H]\)，实化 covariance 为 \(\tfrac12\mathfrak R(\Gamma_c)\)。使用单位实噪声 covariance 的约定时，白化含 \(\sqrt2\) 因子。也可选另一目标函数归一化，但所有 \(\sigma,\tau,\mu,\Lambda\) 必须随之统一缩放。Improper noise 必须用完整实 covariance，不能只用 \(\Gamma_c\)。

若各照明有相同接收几何和相同空间白化，以下使用一个共同 \(G_S\)。标量照明权重 \(\sqrt{w_p}\) 吸收到 \(K_p,B_p,r_p\)，只吸收一次。可分离的 illumination covariance 可先在线性照明组合中白化；完全非可分离的联合白化一般破坏共同 SOM 基底，须改用不同 \(G_{S,p}\) 或总体堆叠算子，不得套用第 2–5 节的共同-\(G_S\) 结论。

## 定义 1.4：物理 metric 与 admissible tangent

物理材料 metric \(M_x\succ0\) 定义 \(\|\delta x\|_{M_x}^2=\delta x^*M_x\delta x\)。例如分片常数介电常数的 \(L^2\) metric 使用体积权重；不同单位的物理参数必须先声明单位/尺度。不能把参数化列长度或参数数目当作物理 metric。

电流 metric 为 \(M_c\succ0\)。电流密度的离散 \(L^2\) metric 通常用 \(\operatorname{diag}(v_\ell)\otimes I_3\)；若变量是单元偶极矩 \(p_\ell=v_\ell j_\ell\)，相应 metric 用 \(\operatorname{diag}(1/v_\ell)\otimes I_3\)。实化后各实虚块使用相同权重。取 \(C_c^*C_c=M_{c,R}\)。

给定物理 tangent 生成列 \(T_0\in\mathbb R^{n_x\times k}\)，先施加线性化等式约束/冻结 active face，再作
\[
 M_x^{1/2}T_0=U_d S_d W_d^*,\qquad
 Q_T=M_x^{-1/2}U_d=T_0W_dS_d^{-1}.                  \tag{1.4}
\]
仅保留数值非零奇异值。于是
\[
 Q_T^*M_xQ_T=I_d,\quad {\rm ran}\,Q_T={\rm ran}\,T_0,
 \quad \delta x=Q_Tz,\ \|\delta x\|_{M_x}=\|z\|.       \tag{1.5}
\]
这删除冗余列和真正物理零方向；**不能因为某物理方向的 \(Jz=0\) 就把它从 admissible tangent 删除**。

若 metric 半正定，先明确商空间和被商掉的零方向，再在商空间正交化；否则输出 `metric_singular`。若 \(d=0\)，输出 `empty_tangent`，所有能量为零，排序无意义。

不等式约束的边界产生 tangent cone，通常不是线性空间。对整个 cone 不存在既覆盖全部方向又具有对称各向同性探针的正交基。可选择冻结 active-face 的线性切空间，或在 cone 的 span 上做诊断并明确它不等于所有可行方向。不得默默把 cone 当线性空间。

### 引理 1.5：metric-dependent adjoint
\(D:X\to C\) 在原始物理坐标中的伴随是
\[
 D^\dagger=M_x^{-1}D^*M_c.                           \tag{1.6}
\]
正交化后为普通转置。**证明：** 直接代入
\(\langle Dx,c\rangle_{M_c}=\langle x,D^\dagger c\rangle_{M_x}\)。□

## 定义 1.6：最终统一算子

令
\[
 A=C_c\mathfrak R(XG_D)C_c^{-1},\quad R=(I-A)^{-1},
\]
\[
 B_pz=C_c\,\mathfrak r\!\left(DX[Q_Tz]E_p\right),\quad
 K_p=RB_p\in\mathbb R^{n_c\times d},                 \tag{1.7}
\]
\[
 G_S=W_y\mathfrak R(G_S^{\rm phys})C_c^{-1},\quad
 J_{T,p}=G_SK_p\in\mathbb R^{m\times d}.              \tag{1.8}
\]
其中 \(DX[Q_Tz]\) 包含实坐标到复材料的重组和参数链式导数。以下省略坐标变换符号，把正交坐标中的相容因子仍写作 \(X,G_D\)，使 \(A=XG_D\)。\(X\) 在变换后不一定保持原始对角外观。

给定 \(z\)，一个 probe 对每个照明只需要构造 \(B_pz\)，解一次 \( (I-A)u_p=B_pz\)，无需形成完整 \(K_p\)。同一次 probe 在所有照明使用相同 \(z\)。

## 1.7 DDA 极化率适配器

归一化偶极矩系统为
\[
 X_a=\operatorname{diag}(a_\ell(\chi_\ell))\otimes I_3,\quad
 (I-X_aG_{\rm off})p_p=X_ae_p^{\rm inc},\quad
 E_p^{\rm exc}=e_p^{\rm inc}+G_{\rm off}p_p.           \tag{1.9}
\]
切向为
\[
 Dp_p[\delta\chi]=(I-X_aG_{\rm off})^{-1}
                   DX_a[\delta\chi]E_p^{\rm exc}.    \tag{1.10}
\]
必须使用真正 \(Da(\chi)[\delta\chi]\)，不是 \(\delta a=\delta\chi\)。再用偶极矩 metric 转到 (1.7)。

固定单元体积 \(v\)、真空背景下，已核对旧实现的实例为
\[
 a_0=\frac{3v\chi}{\chi+3},\quad
 a=\frac{a_0}{1-i c_k a_0},\quad c_k=\frac{k^3}{6\pi},
\]
\[
 Da[\delta\chi]=
 \frac{9v}{(\chi+3)^2(1-i c_k a_0)^2}\,\delta\chi.    \tag{1.11}
\]
**证明：** \(a_0'=9v/(\chi+3)^2\)，且
\(d[a_0/(1-ic_ka_0)]/da_0=(1-ic_ka_0)^{-2}\)。□

若使用逆极化率形式，
\(D(a^{-1})=-a^{-1}(Da)a^{-1}\)；在 \(a=0\) 处不能这样计算，应采用 (1.9)–(1.10)。填充因子乘极化率、填充因子进入有效介电常数、动态去极化修正，是不同模型，必须对**实际代码公式**求导。体积、频率也变化时，还需它们的导数。

数值实现建议将 (1.11) 合并成
\[
 a=\frac{3v\chi}{\chi+3-i3c_kv\chi},\qquad
 a'=\frac{9v}{(\chi+3-i3c_kv\chi)^2}.
\]
它避免先计算 \(a_0\) 时在 \(\chi=-3\) 产生人为中间极点；只有合并分母为零才是真正模型奇点。

\(p/v\) 与平均极化电流密度的联系必须和 self-field convention 一致；DDA exciting field 不应误当连续单元内局部场。\(a\) 是由材料、形状、体积及 self correction 产生的离散散射响应，不是物理材料坐标的替代品。

## 1.8 复 SOM 的相位平面对照

也可用复 \(G_S=U\Sigma V^H\) 和 \(K_{c,p}\in\mathbb C^{n\times d}\) 计算每个**复 SOM 模态组**的能量；这对应真实电流中的二维相位平面，预算一个复模态通常保留两个实状态维度。重复奇异值须按组处理。

若 \(D_c\) 作用于实参数，\(H_c=D_cD_c^H\) 之外还需要 pseudo-covariance \(P_c=D_cD_c^T\)：
\[
 \mathfrak V(D_c)\mathfrak V(D_c)^T
 =\frac12\begin{pmatrix}
 \Re(H_c+P_c)&-\Im(H_c-P_c)\\
 \Im(H_c+P_c)& \Re(H_c-P_c)
 \end{pmatrix}.                                    \tag{1.12}
\]
**证明：** 写 \(D_c=A+iB\)，展开 \(AA^T,AB^T,BB^T\)。□
只有适当圆对称的成对材料列（例如 \([D,iD]\)）使 \(P_c=0\) 时，才能进一步简化为 \(\tfrac12\mathfrak R(H_c)\)。因此单独复 Hermitian covariance 的谱，不能一般冒充实 GN Jacobian 的谱。

---

# 2. 多照明：横向 covariance 与纵向 Jacobian

令 \(G_S=U\Sigma V^*\) 为 thin SVD，\(s=\operatorname{rank}G_S\)，
\(U\in\mathbb R^{m\times s},V\in\mathbb R^{n_c\times s}\)，
\(\Sigma=\operatorname{diag}(\sigma_1,\ldots,\sigma_s)>0\)，按降序排列。
\[
 \mathcal J=[G_SK_1,\ldots,G_SK_P]\in\mathbb R^{m\times Pd},\quad
 J_T=\begin{bmatrix}G_SK_1\\\vdots\\G_SK_P\end{bmatrix}
       \in\mathbb R^{Pm\times d}.                   \tag{2.1}
\]
\[
 C_T=\sum_{p=1}^P V^*K_pK_p^*V\in\mathbb R^{s\times s}. \tag{2.2}
\]
这里 \(C_T\) 是可观测电流模态坐标中的**聚合切向 covariance**，不是直接的材料 covariance。若任务指定材料 covariance \(\Omega\)，将 \(K_p\) 替换为 \(K_p\Omega^{1/2}\)；本文默认 \(\Omega=I_d\)。

## 定理 2.1：两种对象的精确恒等式
\[
 \mathcal J\mathcal J^*=U\Sigma C_T\Sigma U^*,
 \qquad (J_TJ_T^*)_{pq}=G_SK_pK_q^*G_S^*.           \tag{2.3}
\]
且 \(J_T^*J_T=\sum_pK_p^*G_S^*G_SK_p\)，
\(\|\mathcal J\|_F=\|J_T\|_F\)。
**证明。** 分块乘法。横向乘法只保留 \(p=q\) 的和；纵向乘法有全部跨照明块。Frobenius 范数均为 \(\sum_p\|G_SK_p\|_F^2\) 的平方根。□

## 定理 2.2：\(C_T=cI_s\) 的正确等价范围
设 \(c>0\)。则：
(a) 横向 \(\mathcal J\) 的左奇异子空间为 \(U\)，非零奇异值为 \(\sqrt c\,\sigma_i\)；
(b) 对所有共同秩 \(r\) 测量正交投影 \(\Pi\)，最小化
\(\sum_p\|(I-\Pi)G_SK_p\|_F^2\) 的最优解为 top-\(r\) 的 \(U\) 子空间（退化谱时不唯一）；
(c) 在固定 SOM 模态子集内，top-\(r\) 的 \(\sigma_i\) 最优；
(d) (b) 对 stacked 数据意味着使用 \(I_P\otimes\Pi\)，通常总秩 \(Pr\)，不是任意秩 \(r\) stacked projector。
**证明。** (a) 来自 (2.3)；(b) 将误差写成
\({\rm tr}(H)-{\rm tr}(\Pi H)\)，\(H=cU\Sigma^2U^*\)，用正交特征展开或 Ky Fan 极值原理；(c)(d) 是其限制与分块形式。□

这些结论**不能推出完整 \(J_T\) 的谱相同**，更不能推出相同材料右奇异向量、相同 Hessian 或相同 GN step。

## 定理 2.3：跨照明的可分离条件
若存在 \(H_P=(h_{pq})\succeq0\)，使
\[
 V^*K_pK_q^*V=h_{pq}I_s\quad\forall p,q,            \tag{2.4}
\]
则
\[
 J_TJ_T^*=H_P\otimes(U\Sigma^2U^*),\qquad
 C_T=({\rm tr}H_P)I_s.                              \tag{2.5}
\]
非零平方奇异值为所有正的 \(\lambda_a(H_P)\sigma_i^2\)。
特别地：
- \(\operatorname{rank}H_P=1,\ {\rm tr}H_P=c\)：stacked 与 horizontal 非零奇异值相同；
- \(H_P=(c/P)I_P\)：stacked 有 \(P\) 份 \(\sqrt{c/P}\sigma_i\)，不是一份 \(\sqrt c\,\sigma_i\)；
- 多个 \(\lambda_a>0\) 时，乘积可能交错，不能把共同 SOM top-\(r\) 等同于任意 stacked top-\(r\)。

**证明。** 将 (2.4) 代入 (2.3)，再将 \(H_P\) 和 \(U\Sigma^2U^*\) 分别谱分解。□

一个更弱但有用的条件是，对 \(i\ne j\)，
\(v_i^*K_pK_q^*v_j=0\) 对所有 \(p,q\) 成立。此时 stacked covariance 按 SOM mode 分成不混合的 \(P\) 维 sector，但每个 sector 的谱由整块跨照明 covariance 决定；\(q_i\) 只给这块的迹。\(K_p=a_pK_0,\ V^*K_0K_0^*V=c_0I\) 是 rank-one 可分离条件的简单充分构造。

## 反例 2.4：最小非平凡多模态反例
\[
 G_S=\operatorname{diag}(2,1),\quad d=1,\quad
 K_1=(1,0)^T,\quad K_2=(0,1)^T.
\]
于是 \(C_T=I_2\)，\(\mathcal J=\operatorname{diag}(2,1)\)，而
\(J_T=(2,0,0,1)^T\)。横向奇异值为 \(2,1\)，纵向只有一个 \(\sqrt5\)。
跨照明块 \(K_1K_2^*\ne hI\)，明确破坏 (2.4)。□

---

# 3. 固定 SOM 基底的精确切向能量定理

对 \(S\subset\{1,\ldots,s\}\)，定义
\(P_S=V_SV_S^*,\Pi_S=U_SU_S^*\)，有 \(G_SP_S=\Pi_SG_S\)。

## 定理 3.1：无需对角假设的能量恒等式
\[
 \sum_p\|G_S(I-P_S)K_p\|_F^2
 =\|(I-\Pi_S)\mathcal J\|_F^2
 =\sum_{i\notin S}\sigma_i^2(C_T)_{ii}.              \tag{3.1}
\]
stacked 的对应误差是
\(\|J_T-(I_P\otimes\Pi_S)J_T\|_F^2\)，等于同一个右端。
**证明。** 左乘 \(U^*\) 不改变可观测部分的 Frobenius 范数。第 \(i\) 行的平方范数为 \(\sigma_i^2\sum_p\|v_i^*K_p\|^2\)。不同输出行的平方范数相加，故所有 \(C_T\) 非对角元不进入迹；没有使用它们小或为零的假设。\(G_S\) 零空间部分本来测不到，不影响等式。□

定义
\[
 q_i=(C_T)_{ii}=\sum_p\|K_p^*v_i\|^2,\qquad
 \tau_i=\sigma_i^2q_i.                              \tag{3.2}
\]

## 推论 3.2：固定 cardinality 的最优选择
在固定 \(V\)、固定 \(|S|=r\) 下，选最大的 \(r\) 个 \(\tau_i\)，精确最小化 (3.1)。
**证明。** 总能量固定，最大化 \(\sum_{i\in S}\tau_i\)。任何未选权重大于已选权重的交换都严格改善；交换直至 top-\(r\)。边界相等时有多个最优集合。□

这不是任意旋转子空间下的最优解；后者取
\(H=\Sigma C_T\Sigma\) 的 top eigenspace。也不是 regularized GN step 的最优解。一个反例是
\(J=\operatorname{diag}(10,1)\)、\(r=(0,1)^T\)、\(\Lambda=\lambda I,\lambda>0\)。
最大能量行是第一行，但保留第一行得到零 step；保留第二行得到完整真实 step。relevance 还受 residual 和 regularizer 影响。

---

# 4. 排序反转：代数、物理与正面充分条件

## 定理 4.1：零激励和严格反转
(a) \(q_i=0\iff K_p^*v_i=0\ \forall p\)。这与 \(\sigma_i>0\) 相容。
(b) 若 \(\sigma_i>\sigma_j>0,\ q_i,q_j>0\)，则
\[
 \tau_i<\tau_j\iff
 \frac{q_j}{q_i}>\left(\frac{\sigma_i}{\sigma_j}\right)^2. \tag{4.1}
\]
若 \(q_i=0<q_j\)，反转自动成立；若 \(q_j=0\)，不能有此严格反转。
**证明。** (a) 是非负平方和为零的充要条件；(b) 对 \(\sigma_i^2q_i<\sigma_j^2q_j\) 除正数。□

## 定理 4.2：top-\(r\) 保持条件
令 \(S_\sigma=\{1,\ldots,r\}\)。其也是 tangent-energy 最优集合的充要条件为
\[
 \min_{i\in S_\sigma}\tau_i\ge\max_{j\notin S_\sigma}\tau_j. \tag{4.2}
\]
严格不等式给唯一跨边界集合。容易核验的充分条件包括：
\[
 0<q_L\le q_i\le q_U,\quad
 \sigma_r^2q_L>\sigma_{r+1}^2q_U;                    \tag{4.3}
\]
或
\[
 \|C_T-cI\|_2\le \rho c,\quad 0\le\rho<1,\quad
 \frac{\sigma_r^2}{\sigma_{r+1}^2}>\frac{1+\rho}{1-\rho}.
                                                               \tag{4.4}
\]
**证明。** (4.2) 用交换论证；(4.3) 逐项夹逼；(4.4) 用
\(|e_i^*(C_T-cI)e_i|\le\rho c\)。□

## 定理 4.3：注入 tight frame 与传播稳定性

在完整正交电流空间令
\[
 S_B=\sum_pB_pB_p^*,\quad L=I-A,\quad
 S_K=\sum_pK_pK_p^*=L^{-1}S_BL^{-*}.                 \tag{4.5}
\]
若 \(S_B=cLL^*\)，则 \(S_K=cI\)，从而 \(C_T=cI_s\)。
更实用地，若
\[
 c(1-\rho)I\preceq S_B\preceq c(1+\rho)I,\quad
 \|A\|_2\le a<1,
\]
则
\[
 \frac{c(1-\rho)}{(1+a)^2}\le q_i
 \le \frac{c(1+\rho)}{(1-a)^2}.                       \tag{4.6}
\]
**证明。** 第一结论直接相乘。第二结论先做 congruence，得到
\(c(1-\rho)RR^*\preceq S_K\preceq c(1+\rho)RR^*\)。
\(\sigma_{\min}(L)\ge1-a,\ \sigma_{\max}(L)\le1+a\)，再对单位 \(v_i\) 取 Rayleigh quotient。□

也可用已认证的 \(\sigma_{\min}(L),\sigma_{\max}(L)\) 替代 \(1\mp a\)，不要求 Neumann 收缩。完整 \(S_B\) 正定要求足够照明和材料自由度；不满足时不能假装有 \(q_L>0\)。

**Maxwell 解释。** 各向同性介电常数变化仅在局部总场方向注入电流；完整的电流模式并非都能被一个照明的材料扰动激励。多极化、足够独立的局部照明可形成近 tight frame；材料支持约束、场节点和极化缺口会破坏它。\(S_B=cLL^*\) 是严格的“传播预补偿激励”条件，不是声称任意天线阵都能实现的设计定理。

另一个精确式是
\[
 q_i=\sum_p\|B_p^*R^*v_i\|^2.                        \tag{4.7}
\]
它把 relevance 解释为：反传并重散射后的 SOM 电流模态，与局部总场加权的 admissible material tangent 的耦合。只看 \(\sigma_i\)，漏掉了注入与重散射两者。

## 定理 4.4：任意强代数反转
对 \(t\to\infty\)，取
\[
 G_t=\operatorname{diag}(1,t^{-1}),\quad
 K_t=\operatorname{diag}(t^{-2},1).
\]
则 \(\sigma_1/\sigma_2=t\to\infty\)，
\(q=(t^{-4},1)\)，\(\tau=(t^{-4},t^{-2})\)，
\(\tau_1/\tau_2=t^{-2}\to0\)。同时 \(\|G_t\|_2=\|K_t\|_2=1\)。
**证明。** 直接计算。无需靠不稳定 state resolvent 或无界接收增益。□

## 定理 4.5：被动、稳定矢量 Maxwell 点偶极模型中的反转族

考虑固定被动各向同性点偶极，\(\chi_0=1+0.1i\)，使用 (1.11) 且选 \(v,k\) 避免分母为零，\(a'(\chi_0)\ne0\)。无偶极间耦合，\(G_{\rm off}=0\)，其 (1.9) state resolvent 恒为 \(I\)。物理 metric 固定。

在 \(z\) 方向测量该偶极辐射的两个横向极化泛函，单位分析向量为
\(e_x\) 和 \(\cos\theta\,e_x+\sin\theta\,e_y\)，独立白化接收噪声作为测量模型假设。共同几何因子 \(g\ne0\) 可含传播和物理常数：
\[
 G_\theta=g
 \begin{pmatrix}1&0&0\\\cos\theta&\sin\theta&0\end{pmatrix},
 \quad 0<\theta<\pi/2.
\]
其两个正奇异值
\[
 \sigma_\pm=|g|\sqrt{1\pm\cos\theta},\quad
 v_+=(\cos(\theta/2),\sin(\theta/2),0)^T,\quad
 v_-=(-\sin(\theta/2),\cos(\theta/2),0)^T.
\]
令 \(\rho=\sigma_-/\sigma_+=\tan(\theta/2)\)，选择单位振幅、沿 \(z\) 传播的入射平面波，横向极化
\[
 E^{\rm inc}=\frac{\rho^2v_++v_-}{\sqrt{1+\rho^4}}.
\]
固定损耗虚部，仅改变材料实部的单位物理 tangent，则
\(K=a'(\chi_0)E^{\rm inc}\)，固定 metric 因子在比值中抵消。因此
\[
 \frac{q_+}{q_-}=\rho^4,\quad
 \frac{\tau_+}{\tau_-}=\rho^2\to0,\quad
 \frac{\sigma_+}{\sigma_-}=\rho^{-1}\to\infty.         \tag{4.8}
\]
**证明。** \(G_\theta^*G_\theta\) 的两个横向特征向量和特征值直接给上式；对 \(K\) 投影即得。源满足横向性，材料为被动耗散材料，state equation 不靠近多重散射共振。□

这是 Maxwell 偶极离散模型中**复 SOM 模态组**的物理实现；实化时按对应相位平面计秩，不对任意退化实奇异向量宣称 basis-independent 排序。

**适用边界。** 当 \(\theta\to0\)，接收问题退化，弱模态的绝对信息也趋零。该定理不是固定噪声下弱模态始终可恢复的声明。独立接收噪声和理想分析泛函也不自动包含真实硬件耦合。连续扩展物体的结论不能由本点偶极证明自动获得；第4A节另给精确连续球体证明。固定有限孔径、接收间距及统一观测condition number约束下的无界反转并未据此证明。对任一固定 \(\theta>0\) 的严格反转，若某体积分离散/小粒子近似对算子及材料导数的扰动足够小，严格不等式由连续性保留；近似误差的实际上界需另证或经网格/独立求解器验证。

若另强制 \(\sigma_{\min}^{+}(G)\ge s_0>0,\|G\|\le M\)，则 \(\sigma_i/\sigma_j\le M/s_0\)，任意强奇异值比本来就不可能。被动性和 state stability 不提供所有 \(q_i\) 的正下界；禁止把它们误当“不会反转”的理由。

---

# 4A. 连续、被动、非共振 Maxwell 的精确反转构造

本节独立证明连续 Maxwell 特例，不把点偶极结果自动外推。观测是给定方向的 far-field amplitude，允许两个近共线的理想接收极化泛函；材料 tangent 为球内均匀介电常数实部这一物理方向。

## 引理 4A.1：精确体积算子在球内 L2 空间有界
令 \(D\) 为固定有界球，\(k>0\)，\(g_k(r)=e^{ik|r|}/(4\pi|r|)\)。
在 \(\mathcal H=L^2(D;\mathbb C^3)\) 上把电流零延拓到全空间，定义包含正确 distributional self term 的
\[
 \mathcal G_D j=\left.(k^2I+\nabla\nabla)(g_k*j)\right|_D .
\]
它是有界算子。
**证明。** 分为 \(\nabla\nabla g_0*\) 和
\(k^2g_k I+\nabla\nabla(g_k-g_0)\) 两部分。第一部分的 Fourier multiplier 为
\(-\xi\xi^T/|\xi|^2\)，L2 operator norm不超过1；self term包括在分布导数中。第二部分在原点至多 \(O(1/|r|)\)，在有界差域 \(D-D\) 可积，满足有限 Schur bound。零延拓及限制不增加L2范数，故结论成立。□

## 定理 4A.2：固定被动球和固定state margin下的任意强反转
球体体积为 \(v_D\)，固定材料contrast为 \(\chi_0=i\eta,\eta>0\)，取
\[
 |\chi_0|\|\mathcal G_D\|\le1/10.
\]
使用精确Maxwell方程
\[
 j=\chi_0(E^{inc}+\mathcal G_Dj),\qquad
 R_{\chi_0}=(I-\chi_0\mathcal G_D)^{-1}.
\]
材料单位tangent为球内均匀 \(\delta\chi=1/\sqrt{v_D}\)，仅改变实部，保持正损耗；其材料L2范数为1。入射是沿z传播的单位振幅横向平面波。

令
\[
 w_x(r)=v_D^{-1/2}e^{ikz}e_x,\qquad
 w_y(r)=v_D^{-1/2}e^{ikz}e_y .
\]
沿 \(+\hat z\) 的两个far-field接收极化为
\(e_x\)、\(\cos\theta e_x+\sin\theta e_y\)，\(0<\theta<\pi/2\)。
精确current-to-data operator为
\[
 \mathcal G_{S,\theta}j
 =\gamma\sqrt{v_D}\begin{pmatrix}
 \langle w_x,j\rangle\\
 \cos\theta\langle w_x,j\rangle+\sin\theta\langle w_y,j\rangle
 \end{pmatrix},\qquad \gamma=k^2/(4\pi).              \tag{4A.1}
\]
内积约定 \(\langle w,j\rangle=\int_D\overline w\cdot j\)。相同白化尺度不改变下述比值。这里的两个SOM模态采用复表示；在实化代码中，每个复模态对应完整二维实相位平面，以下q和tau是该平面的总能量，不能仅取其中一条实轴。假设两个接收通道的噪声在这组读数坐标中等强、独立；非对角硬件噪声须重新白化再分析。

定义
\[
 w_+=\cos(\theta/2)w_x+\sin(\theta/2)w_y,\quad
 w_-=-\sin(\theta/2)w_x+\cos(\theta/2)w_y,\quad
 \rho=\tan(\theta/2).
\]
令 \(e_\pm\) 为对应横向极化，选择入射极化
\(E=(\rho^2e_++e_-)/\sqrt{1+\rho^4}\)。则两个正SOM奇异值和材料切向能量满足
\[
 \frac{\sigma_+}{\sigma_-}=\rho^{-1},\qquad
 \frac{q_+}{q_-}=\rho^4,\qquad
 \frac{\tau_+}{\tau_-}=\rho^2.                       \tag{4A.2}
\]
当theta趋零得到任意强反转，同时 \(\|R_{\chi_0}\|\le10/9\)，与theta无关。此处current、材料、观测都是连续Maxwell对象，不是Born/DDA近似。

**证明。**
1. 引理4A.1与Neumann级数给精确resolvent存在，范数≤10/9。材料虚部为正，被动耗散。
2. 均匀材料时 \(j(\chi)=\chi R_\chi E^{inc}\)。对单位材料tangent求导：
\[
 K(E)=R_{\chi_0}^2(E_xw_x+E_yw_y).                  \tag{4A.3}
\]
这是精确导数，保留全部多重散射项。
3. 球、背景和z轴对绕z旋转及包含z轴的镜面反射不变。因此
\(M_{ab}=\langle w_a,R_{\chi_0}^2w_b\rangle\)，\(a,b\in\{x,y\}\)，
与全部二维旋转及反射交换，故 \(M=hI_2\)。
与旋转交换只允许 \(aI+bJ_{90^\circ}\)，反射再强制b=0。
4. \(\|R_{\chi_0}^2-I\|\le(1-1/10)^{-2}-1=19/81\)，
故 \(|h-1|\le19/81\)、\(|h|\ge62/81>0\)。非零材料灵敏度有显式保证。
5. (4A.1)在正交 \(w_x,w_y\) 平面的矩阵为
\(\gamma\sqrt{v_D}\begin{pmatrix}1&0\\\cos\theta&\sin\theta\end{pmatrix}\)。
SVD给 \(\sigma_\pm=|\gamma|\sqrt{v_D}\sqrt{1\pm\cos\theta}\)。
步骤3给 \(\langle w_\pm,K(E)\rangle=hE_\pm\)，所以q比为rho^4，再乘sigma平方比即得tau比。□

**适用边界。** 任意强反转在一个明确的连续、被动、非共振类中成立，而不是只在任意矩阵中成立。它仍不声称固定接收condition number和固定绝对信号裕量下可有无界反转；不覆盖任意非球形/非均匀材料/完整材料tangent/全部照明；理想far-field泛函也不自动包含真实硬件互耦。无界比值来自接收极化趋共线及材料激励的协同变化，而不是state共振。有限theta时的绝对能量及噪声条件决定工程价值。

---

# 5. 非对角 covariance、weighted mixing 和子空间旋转

本节的角度结论针对聚合horizontal measurement covariance，不是完整stacked Jacobian的奇异子空间；后者需要第2节的跨照明条件或跨照明块分析。

令
\[
 C_T=D_T+E_T,\quad D_T=\operatorname{diag}(q_i),\quad
 H=\Sigma C_T\Sigma=H_0+F,\quad
 H_0=\operatorname{diag}(\tau_i),\quad F=\Sigma E_T\Sigma.
\]
\(E_T\) 和 \(F\) 为 Hermitian/实对称，且对角为零。

## 定理 5.1：真正的 mixing 大小
\[
 \|F\|_F^2=\sum_{i\ne j}\sigma_i^2\sigma_j^2|(C_T)_{ij}|^2,
 \quad \|F\|_\nu\le\sigma_1^2\|E_T\|_\nu,\quad \nu=2,F.
                                                               \tag{5.1}
\]
在所有 \(\sigma_i>0\) 时，也有
\(\|F\|_\nu\ge\sigma_s^2\|E_T\|_\nu\)。
**证明。** 逐元素计算 Frobenius 等式；上界由子乘性，下界对
\(E_T=\Sigma^{-1}F\Sigma^{-1}\) 使用上界。□

为避免猜测旧 `eta_mix`，本文重新声明：
\[
 \eta_{\rm raw}=\|E_T\|_F/\|C_T\|_F,\quad
 \eta_{\rm weighted}=\|F\|_F/\|H\|_F .
\]
非零分母时
\[
 \kappa_\Sigma^{-2}\eta_{\rm raw}
 \le\eta_{\rm weighted}
 \le\kappa_\Sigma^2\eta_{\rm raw}.                  \tag{5.2}
\]
故 raw mixing 小不等于 weighted mixing 小。若分母为零，返回 `zero_tangent_energy`，不要伪造数值比例。

## 定理 5.2：带 gap 的 principal-angle 界
将 \(\tau_i\) 重排降序，目标为 top-\(r\) 的坐标子空间 \(S_0\)，
\[
 g=\tau_{(r)}-\tau_{(r+1)},\quad e=\|F\|_2.
\]
若 \(g>0,\ e<g/2\)，令 \(S_H\) 为 \(H\) 的 top-\(r\) invariant subspace，则
\[
 |\lambda_i(H)-\tau_{(i)}|\le e,\qquad
 \|\sin\Theta(S_0,S_H)\|_2\le\frac{e}{g-e},
\]
\[
 \|\sin\Theta(S_0,S_H)\|_F\le\frac{\|F\|_F}{g-e}.      \tag{5.3}
\]
**证明。** Weyl 界来自对任意单位向量的 Rayleigh quotient 扰动不超过 \(e\)，再用 min–max。对 \(H\) 目标特征向量矩阵 \(Y\)，投影至 \(H_0\) 的非目标空间 \(X_\perp\)，得到 Sylvester 方程
\[
 H_{0,\perp}(X_\perp^*Y)-(X_\perp^*Y)\Lambda_H
 =-X_\perp^*FY.
\]
左谱至多 \(\tau_{(r+1)}\)，右谱至少 \(\tau_{(r)}-e\)，间距至少 \(g-e\)。两个实谱区间分离的 Sylvester 逆在算子及 Frobenius 范数下受 \(1/(g-e)\) 控制（也可由移位后的指数积分表示直接证明）。且
\(\|X_\perp^*Y\|=\|\sin\Theta\|\)，得到 (5.3)。□

**直接输出量**
\[
 {\tt mixing\_score}=e/\max_i\tau_i,\quad
 {\tt gap\_score}=g/(2e),\quad
 {\tt angle\_upper}=\arcsin\!\left(\min\{1,e/(g-e)\}\right)
\]
最后一式只有 \(g>2e\) 时标记 `certified`。\(e=0,g>0\) 时 angle=0、gap score=∞。\(g\le0\) 输出 `gap_unresolved`；\(0<g\le2e\) 输出 `angle_uncertified`，不能因 mixing 数字小而宣称角度小。

重复谱在目标组内部允许；跨保留边界存在重复则不具有唯一 top-\(r\) 子空间，须把整组纳入或仅报告能量。复 SOM 的相位平面对也应成组。

反例：\(H_0=I_2,F=\begin{pmatrix}0&\rho\\\rho&0\end{pmatrix}\)，任意小 \(\rho>0\) 的主向量都是 \(45^\circ\) 旋转。没有 gap，不能对单向量保证小角度。

## 5.3 只有 probes 时的 mixing 证书
用同一批已计算的 \(a_{p,k}=\Sigma V^*K_pz_k\)，对每个实向量 \(w\) 估计
\[
 h(w)=\sum_p\|K_p^*V\Sigma w\|^2=w^*Hw.
\]
第 8 节给每个 PSD quadratic trace 的置信区间。利用
\(H_{ij}=[h(e_i+e_j)-h(e_i-e_j)]/4\)，得到每个非对角元素区间及半宽 \(r_{ij}\)，从而
\[
 \|F-\widehat F\|_2\le
 \kappa_{\rm off}:=\sqrt{2\sum_{i<j}r_{ij}^2}.
\]
令 \(e_U=\|\widehat F\|_2+\kappa_{\rm off}\)，以及已认证 top subset 的
\(g_L=\min_{i\in S}L_{\tau_i}-\max_{j\notin S}U_{\tau_j}\)。
若 \(g_L>2e_U\)，(5.3) 可用 \(e_U,g_L\) 给带概率角度证书。所有二次型查询必须一起计入 union bound。这需要 \(O(s^2)\) 的小代数统计，但不增加 PDE tangent solves。否则仅输出估计值和 `mixing_uncertified`。对真实代码默认的非实复实现，应先实化或分别用 \(e_i\pm e_j,e_i\pm i e_j\) 恢复实虚交叉项。


# 6. 真正 TSOM 子空间、安全条件与泄漏机制

## 定义 6.1：TSOM 是子空间之和，不是投影乘积
在相容 metric 中计算 \(G_S\) 的 retained right basis \(V_S^+\) 和 \(G_D\) 的 retained right basis \(V_D^+\)。若为 \(G_D\) 选独立场 metric \(M_E=C_E^*C_E\)，应使用
\(G_D=C_E\mathfrak R(G_D^{\rm phys})C_c^{-1}\)、
\(X=C_c\mathfrak R(X^{\rm phys})C_E^{-1}\)；
这样 \(A=XG_D\) 与第 1 节一致。必须固定场 metric，否则 \(G_D\) 的 TSOM 排序本身会改变。

\[
 P_S=V_S^+(V_S^+)^*,\quad
 B_D=\operatorname{orth}((I-P_S)V_D^+),\quad
 P=P_S+B_DB_D^*,\quad Q=I-P.                         \tag{6.1}
\]
`orth` 是带秩揭示的正交化，丢弃零/冗余列。由于 \(P_SB_D=0\)，
\(P^*=P,P^2=P\)，且
\({\rm ran}P={\rm span}(V_S^+,V_D^+)\)。
一般 \(P_SP_D\) 不是正交投影，也不是交集投影。原始 SOM/TSOM 可将材料与 current 作为独立优化变量，中间迭代未必满足精确 state equation；本文 K 是 state-consistent material manifold 的 tangent，不是原算法任意 joint update 的方向。原始 TSOM 的实际优化方式不一定是本文分析的 Galerkin reduced solver；(6.1) 是分析其 retained current span 的对象，不能声称本定理直接替代原 TSOM 算法的收敛分析。

## 引理 6.2：精确 block-resolvent 恒等式
\(K_p=(I-A)^{-1}B_p\)。若 \(I_Q-QAQ\) 在 \({\rm ran}Q\) 可逆，定义
\(R_Q=(I_Q-QAQ)^{-1}\)，则
\[
 QK_p=R_Q\big(QB_p+QAP\,PK_p\big).                  \tag{6.2}
\]
这里 \(I_Q\) 是 discarded 子空间的恒等映射，不是在整个空间对奇异矩阵求逆。
**证明。** 对 \((I-A)K_p=B_p\) 左乘 \(Q\)，用 \(I=P+Q\)，移项并在 \(Q\) 空间求逆。□

完整 \(I-A\) 可逆不意味着该子块或 \(I_P-PAP\) 可逆。例如
\(I-A=\begin{pmatrix}0&1\\1&0\end{pmatrix}\)、\(P=\operatorname{diag}(1,0)\)，
完整最小奇异值为 1，但两个对角子块都奇异。

## 定理 6.3：direct injection + state-transfer 泄漏界
记 \(\gamma_Q=\|R_Q\|_2\)，\(\ell=\|QAP\|_2\)。对 \(\nu=2,F\)，
\[
 \boxed{\ \|QK_p\|_\nu
 \le\gamma_Q\big(\|QB_p\|_\nu+\ell\|PK_p\|_\nu\big).\ } \tag{6.3}
\]
若只知道完整 \(\gamma_L=\|(I-A)^{-1}\|_2\)，可再用
\(\|PK_p\|_\nu\le\gamma_L\|B_p\|_\nu\)。
**证明。** 对 (6.2) 使用三角不等式和相容范数的子乘性。□

特别地，\(QB_p=0\ \forall p\) 且 \(QAP=0\) 时，\(QK_p=0\)，无 material tangent 泄漏。近似满足且 \(\gamma_Q\) 有界时，泄漏受 (6.3) 控制。这是充分条件，不是必要条件；不同项可以相消。

如果还存在 \(R_P=(I_P-PAP)^{-1}\)，令
\(F_Q=R_QQAP R_PPAQ\)。若 \(\|F_Q\|_2=f_Q<1\)，则有无需完整 \(K_p\) 的界：
\[
 QK_p=(I_Q-F_Q)^{-1}R_Q
       \big(QB_p+QAP R_PPB_p\big),                  \tag{6.4}
\]
\[
 \|QK_p\|_\nu\le
 \frac{\gamma_Q}{1-f_Q}
 \big(\|QB_p\|_\nu+\|QAP R_P\|_2\|PB_p\|_\nu\big).
                                                               \tag{6.5}
\]
**证明。** \(P\) 块给 \(PK_p=R_P(PB_p+PAQ QK_p)\)，代回 (6.2)。
对 \(I-F_Q\) 用 Neumann 范数界。□
若 feedback 条件失败但完整系统稳定，(6.5) 仅无证书，不能据此宣布物理系统失稳。

## 定理 6.4：TSOM 的单向保护与缺失保护

令 \(\gamma_D\) 为 \(G_D\) 在 \(V_D^+\) 之外的最大奇异值；若取完全部 current directions，尾项为零。因为 \({\rm ran}Q\perp V_D^+\)，
\[
 \|G_DQ\|_2\le\gamma_D,\quad
 \|QAQ\|_2,\|PAQ\|_2\le\|X\|_2\gamma_D.             \tag{6.6}
\]
因此 \(a_Q=\|X\|\gamma_D<1\) 时，
\[
 \gamma_Q\le(1-a_Q)^{-1}.                           \tag{6.7}
\]
但 (6.6) **不控制** \(QXG_DP\)：它保护“discarded current 作为传播输入”这一侧，不自动保护“retained current 传播并经材料乘法后落到 discarded output”这一侧。

设
\(G_D=U_D^+\Sigma_D^+(V_D^+)^*+G_{D,\rm tail}\)，
\(\|G_{D,\rm tail}\|\le\gamma_D\)，则
\[
 \|QAP\|\le
 \|QXU_D^+\|\,\sigma_{D,1}+\|X\|\gamma_D.            \tag{6.8}
\]
将 (6.7)–(6.8) 代入 (6.3)/(6.5)，得到 TSOM-specific 可计算安全条件。
**证明。** (6.6) 用 SVD 尾子空间；(6.8) 将上述 SVD 分解代入 \(QXG_DP\)，再对子乘积取范数。□

### 推论 6.5：normal propagation + 近均匀材料的正面条件
若在所用相容坐标中 \(G_D\) normal，并且 \(V_D^+\) 包含完整的目标谱组，使它是 invariant spectral subspace，且 \(X=xI\)，则
\[
 \|QAP\|\le |x|\gamma_D.                            \tag{6.9}
\]
若 \(X=xI+\Delta X\)，则
\[
 \|QAP\|\le |x|\gamma_D+\|\Delta X\|\|G_DP\|.          \tag{6.10}
\]
**证明。** normal 情况下 \(G_DV_D^+\subset{\rm ran}V_D^+\subset{\rm ran}P\)，其输出被 \(Q\) 消去。剩余输入/输出位于尾谱，范数不超过 \(\gamma_D\)。对 \(\Delta X\) 单独用范数界。□
一般非均匀矢量 Maxwell 算子不必 normal；此推论是有物理含义的安全机制，不是无条件 Maxwell 属性。

由于 \(P\) 还包含 \(V_S^+\)，另有
\(\|G_SQ\|\le\gamma_S\)，其中 \(\gamma_S\) 为接收尾奇异值。
这可以减小直接测量泄漏，但不能删除第 9 节的反馈项。

## 6.6 三种机制的严格分离
在原始相容物理坐标写
\[
 B_p=B_p^{\rm inc}+B_p^{\rm loc},\quad
 B_p^{\rm inc}z=DX[Q_Tz]e_p^{\rm inc},\quad
 B_p^{\rm loc}z=DX[Q_Tz]G_Dj_p,
\]
\[
 K_p=B_p^{\rm inc}+B_p^{\rm loc}+(R-I)(B_p^{\rm inc}+B_p^{\rm loc}). \tag{6.11}
\]
必要时对各项施加第 1 节的同一正交坐标变换。
第一项是直接材料切向注入；第二项是基准总场中已有散射造成的局部场再分布；第三项是扰动电流产生后经历的多重散射传播。它们可能相干增强或相消，范数可上界相加，但能量不能无交叉项相加。

\([X,G_D]\ne0\) 既不是 failure 的充分条件，也不是必要条件：
可有非交换但所有注入和传播完全在 \(P\) 内；也可 \(X=I\)、交换子为零，但非 normal \(G_D\) 将 retained right mode 送入 discarded output。

## 6.7 最小安全/失败实例

统一 \(n=3,G_S=\operatorname{diag}(3,2,1)\)，SOM 只保留 \(e_1\)。

(a) 安全：\(G_D=\operatorname{diag}(0.1,0.3,0.01)\)，\(X=I\)，\(V_D^+=e_2\)，所以 \(P=\operatorname{diag}(1,1,0)\)。取 \(B=e_1\)，则 \(K=(10/9)e_1\)，\(QK=0\)。

(b) 直接注入 failure：保持 (a)，令 \(B=e_3\)。则 \(K=(100/99)e_3\)，全部 tangent 落在 discarded 子空间；state 仍稳定。

(c) 传播 transfer failure：\(G_D=g e_3e_2^*,X=I,B=e_2\)，\(g=1/2\)。dominant right singular vector 为 \(e_2\)，TSOM 仍保留 \(e_1,e_2\)，但
\(K=e_2+\tfrac12e_3\)，\(QB=0,QK\ne0\)，\([X,G_D]=0\)。
这个例子是精确线性代数机制例，不冒充互易自由空间 Maxwell discretization。

(d) 对称因子例：\(G_D=\operatorname{diag}(0.1,0.6,0.02)\)，
\[
 X=\begin{pmatrix}1&0&0\\0&1&1/2\\0&1/2&1\end{pmatrix},\quad B=e_2.
\]
\(G_D,X\) 都对称，\(X\succ0\)，\(P=\operatorname{diag}(1,1,0)\)。
\[
 K_2=0.98/0.389,\qquad K_3=0.3/0.389\ne0.
\]
这说明即使采用对称传播和材料因子，仍可因材料混合出现 transfer；仍须区分“对称矩阵”与“已实现的 Maxwell Green 矩阵”。

---

# 7. 同一正则器下的真实 inverse-update 扰动

下文 \(J=J_T\) 是真实实化、白化材料 tangent Jacobian，所有目标都在同一 \(z\) 坐标。考虑
\[
 \min_z\ \tfrac12\|Jz+r\|^2+
 R_{\rm reg}(z),\quad
 R_{\rm reg}(z)=\tfrac12z^*\Lambda z+g^*z+c,\quad\Lambda\succeq0.
\]
\[
 H=J^*J+\Lambda,\quad h=J^*r+g,\quad
 H\succeq\mu I,\quad\mu>0,\quad z^\star=-H^{-1}h.     \tag{7.1}
\]
近似系统为 \(\widetilde J,\widetilde r\)，但 \(\Lambda,g\) 不变。令
\[
 \epsilon=\|J-\widetilde J\|_2,\quad
 \delta=\|r-\widetilde r\|,\quad
 \beta=2\|J\|_2\epsilon+\epsilon^2.                   \tag{7.2}
\]

## 定理 7.1：候选 step bound 成立的精确条件
若 \(\beta<\mu\)，则近似解唯一，且
\[
 \boxed{\ 
 \|\widetilde z^\star-z^\star\|
 \le
 \frac{\epsilon\|r\|+(\|J\|+\epsilon)\delta
       +\beta\|z^\star\|}{\mu-\beta}.
 \ }                                                \tag{7.3}
\]
**证明。** 设 \(\Delta J=\widetilde J-J\)，
\[
 \Delta H=J^*\Delta J+\Delta J^*J+\Delta J^*\Delta J,\quad
 \|\Delta H\|\le\beta.
\]
故 \(\widetilde H\succeq(\mu-\beta)I\)。同时
\[
 \Delta h=\Delta J^*r+\widetilde J^*(\widetilde r-r),\quad
 \|\Delta h\|\le\epsilon\|r\|+(\|J\|+\epsilon)\delta.
\]
两套 normal equations 相减：
\[
 \widetilde H(\widetilde z^\star-z^\star)
 =-\Delta h-\Delta H z^\star.
\]
乘逆并取范数得到 (7.3)。□

该式不缺项，但依赖：同一 regularizer（包括线性项）、同一坐标、同一白化、同一 illumination weighting，以及 \(\mu>0\)。若改变 regularizer，必须增加 \(\|\Delta\Lambda\|\) 到 Hessian 扰动界，增加 \(\|\Delta g\|\) 到 gradient 扰动界。仅“都用 Tikhonov”而不同 \(\lambda\)，不满足同一正则器。

## 定理 7.2：可部署的后验 step certificate
令 \(s\) 为近似模型给出的任意候选，normal residual
\(\rho_{\rm solve}=\|\widetilde Hs+\widetilde h\|\)。
若有真实 Hessian 的认证下界 \(\mu_L>0\) 和
\(j_U\ge\|J\|\)，则
\[
 \boxed{\ 
 \|s-z^\star\|
 \le \frac{\rho_{\rm solve}
 +\epsilon\|\widetilde Js+\widetilde r\|
 +j_U(\epsilon\|s\|+\delta)}{\mu_L}.
 \ }                                                \tag{7.4}
\]
**证明。** 令 \(E=J-\widetilde J,d_r=r-\widetilde r\)。则
\[
 Hs+h=(\widetilde Hs+\widetilde h)
 +E^*(\widetilde Js+\widetilde r)+J^*(Es+d_r).
\]
再用 \(H(s-z^\star)=Hs+h\)。□

可取 \(j_U=\|\widetilde J\|_2+\epsilon\)。若 \(\Lambda\succeq\lambda I,\lambda>0\)，直接用 \(\mu_L=\lambda\)。否则可用认证的
\[
 \widetilde\mu_L\le\lambda_{\min}(\widetilde H),\quad
 \mu_L=\widetilde\mu_L-
        \big(2\|\widetilde J\|\epsilon+\epsilon^2\big)>0.
\]
单个最小 Ritz value 通常不是所需的全局特征值下界；应使用可验证的 enclosure、足够完整的小 reduced Hessian 分解、严格下界或 regularizer 下界。缺失时输出 `stability_uncertified`。

(7.4) 不依赖未知 \(\|z^\star\|\)，适合作为代码主验收；(7.3) 保留作理论比较。相对证书可用：若右端为 \(b<\|s\|\)，则
\(\|s-z^\star\|/\|z^\star\|\le b/(\|s\|-b)\)。否则只给绝对界。

## 推论 7.3：正交 spectral projection 有更锐利的界
只针对精确投影模型 \(\widetilde J=\Pi J\)，其中 \(\Pi\) 为正交测量投影，目标可继续使用 \(r\)（或用 \(\Pi r\)，两者最优解相同）。令
\(E=(I-\Pi)J,\epsilon=\|E\|\)。则
\[
 H-\widetilde H=E^*E,\quad
 h-\widetilde h=E^*(I-\Pi)r.                         \tag{7.5}
\]
若 \(\epsilon^2<\mu\)，
\[
 \|\widetilde z^\star-z^\star\|
 \le\frac{\epsilon\|(I-\Pi)r\|+\epsilon^2\|z^\star\|}
           {\mu-\epsilon^2}.                        \tag{7.6}
\]
**证明。** 正交性使保留行与丢弃行交叉项为零，然后重复定理 7.1 的 normal equation 推导。□
这揭示投影情况下 Hessian loss 是二阶 \(\epsilon^2\)，但真实 reduced dynamics 一般失去该正交结构，不能套用此更强式。

## 7.4 Inexact Krylov、约束、TR/line search 与非线性

**Inexact solve。** 若解的近似 normal residual 为 \(\rho_{\rm solve}\)，(7.3) 额外加 \(\rho_{\rm solve}/(\mu-\beta)\)；更推荐直接使用 (7.4)。必须检查相同 metric 的 normal equation residual。固定 GMRES/CG 迭代次数通常并不让输出成为 RHS 的线性映射，因为 Krylov 系数依赖 RHS；因此不能自动把“固定次数迭代”当作固定矩阵 \(\widetilde J\)。应使用一致的线性 reduced operator/adjoint、固定系数线性迭代，或另外分析 inexact oracle 误差。禁止仅通过一次 finite-difference 测试就认定任意 JVP/VJP 对一致。

**固定凸约束。** 在相同闭凸可行集（包括同一 TR ball/冻结线性约束）上，用 variational inequality 代替 normal equation。若近似 quadratic Hessian 下界为 \(\mu-\beta\)，
\[
 \|\widetilde z-z\|\le
 \frac{\|\Delta h\|+\beta\|z\|}{\mu-\beta}.
\]
证明：分别把另一解代入两套最优性不等式并相加，以强单调性控制解差。相同凸非光滑正则项也可用其 subgradient 的单调性抵消，但必须另外确保强凸性；不能把非凸或变化 active-set 的局部二阶导数未经检查代入。

**TR 的模型预测。** 对相同 \(s\)，二次模型 predicted reduction 的差满足
\[
 |{\rm pred}(s)-\widetilde{\rm pred}(s)|
 \le \eta_h\|s\|+\tfrac12\beta\|s\|^2,\quad
 \eta_h=\epsilon\|r\|+(\|J\|+\epsilon)\delta.          \tag{7.7}
\]
若用近似前向值评价接受率，还需分别在基点/候选点认证目标值误差。若前向残差向量误差 \(\le\eta_f\)，对应数据目标值误差
\[
 |f-\widetilde f|\le \eta_f\|\widetilde r\|+\tfrac12\eta_f^2
\]
（同一 regularizer 精确评价）。用区间算出 actual decrease 下界、predicted decrease 上界，只有实际下降下界大于指定接受阈值乘预测下降上界，且预测下降下界为正时，才能 `accept_certified`。其他情况交由高保真 objective 检查或拒绝，并非数学上已经接受。

**Line search。** 需要真实目标下降，或带误差区间的 Armijo 检查。近似模型的下降方向不自动是真实下降方向；gradient error 对方向导数的影响至多 \(\eta_h\|s\|\)。

**非线性重新线性化。** 若邻域中真实 \(DF\) 的 Lipschitz 常数为 \(L_J\)，则
\[
 \|F(x+Q_Ts)-F(x)-J_Ts\|\le\tfrac12 L_J\|s\|^2
\]
（固定线性材料坐标；非线性 retraction 还需其二阶项）。这给局部模型误差，不是全局重建保证。进入下一次 GN/TR 后，应重新计算材料局部场、切空间/metric、必要的 probes 和证书。模型选择 basis 在一次线性化内冻结；若要对 \(P(x)\) 自适应变化的整个算法求导，必须加入 \(DP\) 项并处理排序切换不可微。局部 step bound 不得包装成 global nonlinear reconstruction convergence。

---

# 8. 随机 tangent probes：无偏、置信区间、排序和停止

## 定义 8.1：探针观测
在固定线性化上取独立 \(z_k\in\mathbb R^d\)，\(\mathbb E[z_kz_k^*]=I_d\)。默认可部署证书使用 Gaussian \(N(0,I_d)\)；Rademacher 也给出合法但通常更保守的备选。
\[
 Y_{ik}=\sum_p|v_i^*K_pz_k|^2=z_k^*A_i z_k,\quad
 A_i=\sum_pK_p^*v_iv_i^*K_p\succeq0,\quad
 \widehat q_i=\frac1M\sum_{k=1}^M Y_{ik}.             \tag{8.1}
\]
所有 illumination 使用同一个 \(z_k\)，以对应共享物理材料 tangent。统计项先求各 illumination 能量之和，**不能先相干相加再平方**。

**复 current、实 material probe 的实现补充。** 若不先将 current/measurement 实化，而是用复 \(K_{c,p}\) 和单个复 SOM mode \(v_i\)，第8节方差及 Gaussian \(\chi^2\) 推导必须使用实对称 PSD
\[
 A_i^{\mathbb R}=\operatorname{Re}\sum_p
 K_{c,p}^Hv_iv_i^HK_{c,p}.
\]
于是 \(Y_{ik}=z_k^TA_i^{\mathbb R}z_k\)，\(q_i=\operatorname{tr}A_i^{\mathbb R}\)，方差中的 \(A_i\) 也应替换为 \(A_i^{\mathbb R}\)。等价实现是对复 mode 的二维实相位平面投影。不能把 real Gaussian probes 的方差直接写成未取实部的 Hermitian covariance 的 \(2\operatorname{tr}(A_i^2)/M\)。

## 定理 8.2：无偏与有限样本方差
\[
 \mathbb E\widehat q_i={\rm tr}A_i=q_i.              \tag{8.2}
\]
Gaussian 情况：
\[
 {\rm Var}(\widehat q_i)=\frac2M{\rm tr}(A_i^2)
 \le\frac{2q_i^2}{M}.
\]
实 Rademacher 情况：
\[
 {\rm Var}(\widehat q_i)=\frac2M\sum_{a\ne b}(A_i)_{ab}^2
 \le\frac{2q_i^2}{M}.                               \tag{8.3}
\]
若 \(u_k\) 在单位实球面均匀，必须用 \(z_k=\sqrt d\,u_k\)，或将
\(\sum_p|v_i^*K_pu_k|^2\) 乘 \(d\)。其单样本方差为
\[
 \frac{2}{d+2}\left(d\,{\rm tr}(A_i^2)-q_i^2\right). \tag{8.4}
\]
**证明。** 无偏由 trace 与 \(\mathbb Ezz^*=I\)。Gaussian 用四阶矩
\(\mathbb E z_a z_b z_c z_d=\delta_{ab}\delta_{cd}+\delta_{ac}\delta_{bd}+\delta_{ad}\delta_{bc}\)；
Rademacher 的对角平方恒为 1，只有非对角交叉项贡献方差。球面
\(\mathbb E u_au_b=\delta_{ab}/d\)，四阶矩为 Gaussian 配对项除 \(d(d+2)\)，代入并乘 \(d^2\)。□

## 定理 8.3：可直接求根的 Gaussian 置信区间
设 \(q_i>0\)，\(f(t)=t-1-\log t\)。对 \(0<a<1<b\)，
\[
 {\mathbb P}(\widehat q_i/q_i\le a)\le
 e^{-Mf(a)/2},\quad
 {\mathbb P}(\widehat q_i/q_i\ge b)\le e^{-Mf(b)/2}.    \tag{8.5}
\]
**证明。** 对 \(A_i\) 谱分解，令 \(w_j=\lambda_j/q_i\ge0,\sum w_j=1\)，
则 \(M\widehat q_i/q_i\) 是独立 \(\chi_M^2\) 的加权和。
对于 \(x>0\)，\(\prod_j(1+xw_j)\ge1+x\)；
对于 \(0<x<1\)，\(\prod_j(1-xw_j)\ge1-x\)。
所以其正/负指数矩均被 rank-one \(\chi_M^2\) 的相应指数矩控制。Markov 不等式及分别优化 \(x\) 给 (8.5)。□

给定总失败概率 \(\alpha\)、预定义阶段数 \(L\)、每阶段同时保护 \(N_{\rm test}\) 个 PSD quadratic traces，令
\[
 t=\log(2LN_{\rm test}/\alpha),\qquad
 f(a_M)=f(b_M)=2t/M,\quad 0<a_M<1<b_M.
\]
则所有阶段、所有 trace 同时满足的置信区间为
\[
 \boxed{\quad q_i\in
 [\,\widehat q_i/b_M,\ \widehat q_i/a_M\,].\quad}      \tag{8.6}
\]
无需完整 Jacobian、未知 \(q_i\) 或未知 effective rank。区间可能很宽，这是诚实的有限预算代价。\(q_i=0\) 时所有精确样本恒零；Gaussian 下若 \(q_i>0\)，全部精确样本为零的概率为零。但数值下溢/求解误差会改变判断，不能把浮点零作为证明。

**数值 JVP 误差。** 若对第 \(k\) 次 mode/illumination 聚合输出向量，已认证误差范数 \(\le e_{ik}\)，观测平方为 \(\widetilde Y_{ik}\)，则真实样本在
\[
 Y_{ik}\in
 [\,\max(0,\sqrt{\widetilde Y_{ik}}-e_{ik})^2,\ 
    (\sqrt{\widetilde Y_{ik}}+e_{ik})^2\,].
\]
先平均得到 \(\widehat q_{i,L},\widehat q_{i,U}\)，再用
\([\widehat q_{i,L}/b_M,\widehat q_{i,U}/a_M]\)。
\(e_{ik}\) 可由线性残差和认证 \(\gamma_L\) 得到；仅指定 iterative-solver tolerance、但没有 forward-error 或 resolvent 控制时，不能标为严格概率证书。

**Rademacher 备选。** 由 (8.3) 与 Chebyshev，
\(\mathbb P(|\widehat q_i-q_i|\ge u q_i)\le2/(Mu^2)\)。
令 \(u=\sqrt{2LN_{\rm test}/(M\alpha)}\)。若 \(u<1\)，可用
\([\widehat q/(1+u),\widehat q/(1-u)]\)；若 \(u\ge1\)，只能给上界 ∞ 等无信息区间。若已知 \(q_i\le b_i\)，因 Rademacher \(\|z\|^2=d\)，可改用有界变量 Hoeffding 的绝对半宽
\[
 d b_i\sqrt{\log(2LN_{\rm test}/\alpha)/(2M)}.
\]
两者常比 Gaussian 方案保守。Rademacher 得到零样本不排除非零 \(q_i\)：例如 rank-one 方向可与若干抽到的符号向量正交。

## 定理 8.4：具有 margin 的 top-\(r\) 恢复
设 \(\Delta=\tau_{(r)}-\tau_{(r+1)}>0\)，
\(0<u<\Delta/(\tau_{(r)}+\tau_{(r+1)})\)。若所有模态满足
\(|\widehat\tau_i-\tau_i|\le u\tau_i\)，则 top-\(r\) 集合正确。
Gaussian 对 \(0<u<1\) 有
\[
 {\mathbb P}(\text{任一模态相对误差}>u)
 \le 2s\exp(-Mu^2/8).                               \tag{8.7}
\]
故充分样本数为
\[
 M\ge \frac8{u^2}\log(2s/\alpha).
\]
例如取 \(u=\Delta/[2(\tau_{(r)}+\tau_{(r+1)})]\)，得到
\[
 M\ge32\left(\frac{\tau_{(r)}+\tau_{(r+1)}}{\Delta}\right)^2
          \log(2s/\alpha).                          \tag{8.8}
\]
**证明。** 最弱保留项的估计下界 \((1-u)\tau_{(r)}\) 大于最强丢弃项的估计上界 \((1+u)\tau_{(r+1)}\)。
由 (8.5)，\(f(1+u)\ge u^2/4,\ f(1-u)\ge u^2/2\)，得单项尾界并 union bound。□

(8.8) 是事前充分数，不是算法必须达到的固定数；在线算法用以下 interval separation，更容易在实际大 gap 上提前停止。

## 8.5 有原则的 8→16→32→64 停止
预先固定四个累计样本量，分配总 ranking 失败概率 \(\alpha_{\rm rank}\)，在 (8.6) 取 \(L=4\)。
每阶段：
1. 用累计样本求 \(\widehat\tau_i=\sigma_i^2\widehat q_i\) 及区间 \([L_i,U_i]\)；
2. 令 \(S\) 为估计 top-\(r\)，仅当
   \(\min_{i\in S}L_i>\max_{j\notin S}U_j+{\tt tol}_{\rm rank}\)
   时输出 `ranking_certified`；
3. 否则进入下一阶段。64 后仍不分离，输出 `unresolved/fallback`。

预定义各阶段 union bound 允许这种可选停止，尽管各阶段使用累计、彼此相关的样本。预算外新增阶段或无限检查需重新分配 alpha，不能沿用固定四次公式。

可对相等/近相等谱按组选择或允许近最优能量损失，但必须改写验收目标；不得在无 margin 时声称恢复唯一 top-\(r\)。Fallback 可以采用保守扩充、上一轮已验证模型、较高保真算子或直接拒绝本次降阶；**不得秘密要求完整 Jacobian**。

在联合区间事件上，任何自适应选择的固定 SOM subset 都有
\[
 E_{\rm trunc,F}^2\le\sum_{i\notin S}U_i.            \tag{8.9}
\]
这是固定基底有限集合的统一控制，不需要另做 selection holdout。

如果要控制电流尾而非测量尾，还需
\[
 q_{\rm null}=\sum_p\|(I-VV^*)K_p\|_F^2.
\]
它可作为额外一个 PSD quadratic trace 由 probes 估计，不需要枚举整个 null basis。
\(\sum_{i\notin S}q_i+q_{\rm null}\) 控制该 SOM projector 的总电流遗漏。

## 定理 8.6：实际 reduced solver 的独立 holdout 证书
固定由训练 probes 得到的 reduced solver/basis，随后使用与选择数据独立的新 Gaussian probes。令
\(E=J_T-\widetilde J_T\) 是这个固定线性算子的真实差。
\[
 \widehat e_F^2=M_v^{-1}\sum_k\|Ez_k\|^2.
\]
以 (8.6) 的下侧因子 \(a_{M_v}\) 得到
\[
 \|E\|_2\le\|E\|_F
 \le \sqrt{\widehat e_F^2/a_{M_v}}                  \tag{8.10}
\]
的指定概率上界（有数值误差时用均方样本上包络）。
**证明。** 对 PSD 矩阵 \(E^*E\) 用定理 8.3 的下尾界，再利用 operator≤Frobenius。□

同法可估计第 9 节的 tangent equation residual，而不做 full JVP。若使用 residual route，必须另乘经过认证的传播/测量放大常数。
选择/旋转/训练后的投影不能用同一批训练样本给非统一的 (8.10) 证书；这会引入 selection bias。重复修改模型并用同一验证集接受也须独立新验证或明确的同时置信分配。

少量探针对对角 relevance 的概率估计有理论支持；任意小 fixed budget 的确定性谱范数识别没有普遍保证。给定少于 \(d\) 个确定性 probes，可在其共同正交补构造不可见的非零误差算子。这不否定随机有置信/允许 unresolved 的方案。

---

# 9. 精确投影与真正 reduced solver：必须补齐的动力学误差

## 定义 9.1：冻结 basis 的实际 Galerkin solver
令 \(V_r\) 的列为固定正交 retained current basis，
\(P=V_rV_r^*,Q=I-P\)。假设
\(L_r=I_r-V_r^*AV_r\) 可逆，定义提升算子
\[
 \mathcal R_P=V_rL_r^{-1}V_r^*.
\]
基准前向系统 \(Lj_p=f_p,\ f_p=Xe_p^{\rm inc}\)。实际 reduced state
\[
 j_p^{(r)}=\mathcal R_P f_p .
\]
将同一个冻结 basis 的 reduced forward map 对真实材料求导，得到
\[
 K_p^{(r)}=\mathcal R_P\widetilde B_p,\qquad
 \widetilde B_pz=DX[Q_Tz](e_p^{\rm inc}+G_Dj_p^{(r)}). \tag{9.1}
\]
这里必须与第 1 节采用相同坐标变换。即使 \(B_p\) 是 full-state 的真实注入，也不能把 \(\mathcal R_P B_p\) 自动当作实际 reduced forward map 的导数；其真实导数使用 \(\widetilde B_p\)。

## 定理 9.2：投影/动力学/局部场误差的精确分解
令 \(\Delta B_p=B_p-\widetilde B_p\)。则
\[
 PK_p-K_p^{(r)}
 =\mathcal R_P\big(PAQ\,QK_p+P\Delta B_p\big),       \tag{9.2}
\]
\[
 K_p-K_p^{(r)}
 =QK_p+\mathcal R_P PAQ\,QK_p+\mathcal R_P P\Delta B_p,
\]
\[
 \boxed{\ 
 J_{T,p}-\widetilde J_{T,p}^{(r)}
 =G_SQK_p+G_S\mathcal R_P PAQ\,QK_p+
              G_S\mathcal R_PP\Delta B_p .
 \ }                                                \tag{9.3}
\]
**证明。** full \(P\) 块为
\((I_P-PAP)PK_p=PB_p+PAQ QK_p\)；
reduced 块为
\((I_P-PAP)K_p^{(r)}=P\widetilde B_p\)。
相减、乘提升逆得到 (9.2)，再分解 full error 并左乘 \(G_S\)。□

因此一般 \(PK_p\ne K_p^{(r)}\)。第二项是被删除的 current 经反馈重新进入保留状态后的测量贡献；第三项是 reduced base field 使 material injection 发生变化。

## 推论 9.3：总 Jacobian 误差上界
定义
\[
 E_{\rm trunc,F}
  =\left(\sum_p\|G_SQK_p\|_F^2\right)^{1/2},
\quad L_{Q,F}
  =\left(\sum_p\|QK_p\|_F^2\right)^{1/2},
\]
\[
 D_{B,F}
  =\left(\sum_p\|\Delta B_p\|_F^2\right)^{1/2},
\quad c_{\rm fb}=\|G_S\mathcal R_P PAQ\|_2,\quad
 c_{\rm inj}=\|G_S\mathcal R_PP\|_2 .
\]
则
\[
 \boxed{\ 
 \|J_T-\widetilde J_T^{(r)}\|_2
 \le \|J_T-\widetilde J_T^{(r)}\|_F
 \le E_{\rm trunc,F}+c_{\rm fb}L_{Q,F}+c_{\rm inj}D_{B,F}.
 \ }                                                \tag{9.4}
\]
也可以对每个 block 先使用 operator bounds、再将平方和开方，以避免不必要的 \(\sqrt d\) 因子；但不能把 operator bound 未经说明当成 Frobenius bound。
**证明。** 对 (9.3) 的纵向堆叠使用三角不等式、相容范数子乘性。□

对于 SOM 子集的 \(P\)，第 3/8 节给第一项和电流尾；对于一般旋转 TSOM/ROM basis，必须直接评估这些残余算子或用独立 probes，而不能继续使用固定 SOM 对角求和公式。

若 \(G_SQ=0\)，第一项为零，但第二项仍可能非零。若 \(PAQ=0\)、且 \(\Delta B=0\)，则精确投影和实际 reduced solver 才对材料 tangent 相同。TSOM 的 \(\|PAQ\|\le\|X\|\gamma_D\) 有助于控制反馈，但 \(\|\mathcal R_P\|\) 大时仍会放大。

## 9.4 不形成 Jacobian 的误差常数

设基准 reduced state residual
\(e_p^{\rm state}=f_p-Lj_p^{(r)}\)。则
\[
 j_p-j_p^{(r)}=R e_p^{\rm state},\quad
 \|j_p-j_p^{(r)}\|\le\gamma_L\|e_p^{\rm state}\|.       \tag{9.5}
\]
因材料注入是双线性，若在所用 metric 中有认证常数
\[
 \|DX[Q_Tz]w\|\le c_X\|z\|\|w\|,
\]
则
\[
 \|\Delta B_p\|_2
 \le c_X\|G_D\|\,\gamma_L\|e_p^{\rm state}\|,\quad
 \|\Delta B_p\|_F\le\sqrt d\,c_X\|G_D\|\gamma_L
                          \|e_p^{\rm state}\|.       \tag{9.6}
\]
可用真实 \(\|G_D(j-j_r)\|\) 或局部 block metric 得到更紧界。DDA 要把 \(Da\)、体积和场/电流 metric 纳入 \(c_X\)，不能只用 \(\|\chi\|\)。

近似 residual 与真实 residual 的误差满足
\[
 \delta\le
 \left[\sum_p(\|G_S\|\gamma_L\|e_p^{\rm state}\|)^2\right]^{1/2}.
                                                               \tag{9.7}
\]
这控制离散 full model 与 reduced model 的差；数据生成物理与离散 full model 的 discrepancy 另列，不由 (9.7) 覆盖。

### 定理 9.5：tangent equation residual 证书
对任意固定线性近似 tangent \(\widetilde K_p\)，包括非 Galerkin、N-Port、hyper-reduction，定义
\[
 T_p=B_p-L\widetilde K_p.
\]
则
\[
 K_p-\widetilde K_p=RT_p,\qquad
 \|J_T-\widetilde J_T\|_2
 \le\|G_SR\|_2\left(\sum_p\|T_p\|_F^2\right)^{1/2}.   \tag{9.8}
\]
若只有近似状态注入 \(\widetilde B_p\)，实际可计算
\(\widetilde T_p=\widetilde B_p-L\widetilde K_p\)，需用
\[
 \left(\sum_p\|T_p\|_F^2\right)^{1/2}
 \le \left(\sum_p\|\widetilde T_p\|_F^2\right)^{1/2}
       +D_{B,F}.                                   \tag{9.9}
\]
**证明。** 左乘 full equation residual 的逆；(9.9) 是 stacked triangle inequality。□

用独立 Gaussian probes 估计 \(\sum_p\|\widetilde T_p\|_F^2\)，每 probe 只需 reduced tangent solve 和 full \(L\) matvec，不要求 full Jacobian，也不要求 full tangent solve。代价是必须掌握认证的 \(\|G_SR\|\) 或 \(\|G_S\|\gamma_L\)；若没有，可改用独立 full JVP 验证 (8.10)。这两条路线任选其一，不得把未认证的 resolvent 估计当 rigorous upper bound。

\(\gamma_L\) 可来自 \(\|A\|<1\)、validated inf-sup/coercivity 下界、完整小系统奇异值或其他认证方法。Forward solver residual 本身不是 forward error，亦不能自行认证 \(\gamma_L\)。完整系统稳定也不保证 \(L_r\) 稳定；必须另检 reduced inf-sup/smallest singular value。

## 反例 9.6：spectral projection error 恰为零，实际 solver 仍错
\[
 G_S=(1,0),\quad P=\begin{pmatrix}1&0\\0&0\end{pmatrix},\quad
 A=\begin{pmatrix}0&1/2\\1/2&0\end{pmatrix},\quad B=e_1.
\]
\[
 K=(I-A)^{-1}B=(4/3,2/3)^T,\quad PK=(4/3,0)^T,
 \quad K^{(r)}=(1,0)^T
\]
（此处规定同一注入 \(B\)，隔离动力学机制）。因此
\[
 G_SK=G_SPK=4/3,\quad G_SK^{(r)}=1.
\]
投影误差为零，actual Jacobian error 为 \(1/3\)。丢掉的 \(e_2\) 虽不直接辐射到接收器，却通过 \(A_{12}\) 反馈改变可测的 \(e_1\)。
该矩阵是机制 unit test，不宣称其 \(G_S,A\) 已由特定天线几何实现。

### 反例 9.6b：连 reduced 基准场导数也一致的材料系统
令 \(X=\operatorname{diag}(\chi_1,\chi_2)\)，
\(G_D=\begin{pmatrix}0&1/2\\1/2&0\end{pmatrix}\)，
\(e^{\rm inc}=e_1\)，\(G_S=(1,0)\)，\(P=e_1e_1^*\)。
基点 \(\chi_1=\chi_2=1/2\)，材料 tangent 为 \(\delta\chi=e_1z\)。
真实测量和实际 reduced 测量分别为
\[
 F(\chi_1)=\frac{\chi_1}{1-\chi_1\chi_2/4},\qquad
 F_r(\chi_1)=\chi_1 .
\]
因此
\[
 J=256/225,\quad J_r=1,\quad G_S(I-P)K=0.
\]
在这个基点，\(B=(16/15,0)^T,\widetilde B=(1,0)^T\)，
\(K=(256/225,64/225)^T\)。反馈项贡献 \(16/225\)，
基准局部场差项贡献 \(15/225\)，总差 \(31/225\)。
这验证 (9.3) 的三项分解确实适用于物理材料形式的实际 reduced-map derivative，而不仅是人为指定的固定 RHS。该例仍为两单元代数系统，连续 Maxwell 实现另行验证。

## 9.7 完整理论—算法闭环
按第 8 节选择/检测，按 (9.4) 或 (9.8)–(9.9) 或 (8.10) 得到总 \(\epsilon_U\)，按 (9.7)/独立前向检查得到 \(\delta_U\)，按第 7 节得到 \(\mu_L\)。将这些**总误差上界**代入 (7.3)/(7.4)，才有 actual reduced inverse-update certificate。

若 `ranking_certified` 但 `dynamics_uncertified`，只能说模式排序已知；不能说 inverse update 已保真。若总误差证书通过但 top-\(r\) 排序未唯一，仍可以接受该 reduced model 的更新：**排序证书是解释/选择工具，不是实际更新保真的逻辑必要条件。** 这个分离避免把无关的微小排序 gap 变成模型使用的强制瓶颈。

---

# 9A. 必须纳入对照的已知 primal–dual 精确性条件

本节不是 A9 的首创贡献。参数化 interpolatory ROM 已经有 response gradient/Hessian matching 理论，尤其是 Baur–Beattie–Benner–Gugercin (2011) 的 Theorems 4.1–4.3 [R07]。这里把所需的一阶充分条件写成当前 Maxwell 接口，以免误导为“必须保留全部 tangent current 才能保留材料 Jacobian”。

## 定理 9A.1：基准 primal state + 测量 dual state 足以匹配方向导数
对 \(L(x)j_p=f_p(x)\)、固定 \(G_S\)，取冻结 trial/test bases \(V_r,W_r\)，
\(L_r=W_r^*LV_r\) 可逆，\(j_{r,p}=V_rL_r^{-1}W_r^*f_p\)。
给定测量方向 \(c\)，令 \(\psi_c=L^{-*}G_S^*c\)。如果在当前基点
\[
 j_p\in{\rm ran}V_r,\qquad \psi_c\in{\rm ran}W_r,
\]
则 \(j_{r,p}=j_p\)，并对每个 admissible material tangent \(z\)，
\[
 c^*J_{T,p}z=c^*\widetilde J_{T,p}^{(r)}z.             \tag{9A.1}
\]
**证明。** full primal state 已在 trial space且 reduced equation唯一，故其系数满足 reduced equation，得到 primal equality。同理测试空间中的真实 dual 满足 reduced adjoint equation，因此 lifted reduced dual等于 \(\psi_c\)。真实与降阶方向导数分别为
\(\psi_c^*(Df_p[z]-DL[z]j_p)\) 和
\(\psi_{r,c}^*(Df_p[z]-DL[z]j_{r,p})\)，两者相等。□

若所有 illuminations 的 primal states及所有测量方向的dual states均被覆盖，则完整材料 \(J_T\) 在当前点匹配，而无需将全部 \(K_p\) 的列放入trial space。若只覆盖某一个 residual方向，只能保证相应gradient作用，不能自动保证完整Hessian。Galerkin可取trial=test为primal和dual的并空间，但仍需检查reduced可逆性。

**成本和贡献边界。** 这提供不靠完整材料Jacobian的已知正面基线：P个primal states和指定数量的dual solves。但basis大小、测量有效秩、误差状态与真正工程成本并无免费保证。A9 的low-probe current ranking应和此类sensitivity-preserving ROM比较，而不是仅击败field-POD。

## 反例 9A.2：大 current leakage仍可有完全准确的实际ROM
\[
 L=\begin{pmatrix}1&-1/2\\-1/2&1\end{pmatrix},\
 G_S=(1,-1/2),\quad P=e_1e_1^*,\quad B=I_2.
\]
这里 \(K=L^{-1}\)，\(K_r=P\)，
\[
 G_SK=(1,0)=G_SK_r,\qquad
 G_SPK=(4/3,2/3)\ne G_SK.
\]
所以实际ROM的Jacobian可以比“投影真实K”更准；\(\|QK\|_F=\sqrt{20}/3\)并不小。
这不是固定SOM奇异子集的投影（其P不是G_S的singular projector），因而不违背第3节；它强调不同模型类别必须分清。第6/9节的范数加法界会忽略相消，是充分安全证书，不是必要性或最优性定理。

---

# 10. 最终数学裁决与实验边界

**核心 mismatch：成立。** 固定物理 metric 后，\(\sigma_i\) 描述任意 current 的接收可观测性，\(q_i\) 描述允许材料切向经局部场和全波传播能够实际激励多少该模态，\(\tau_i\) 是其测量能量；第 3–4 节给精确恒等式、等价条件和严格反例。它不是“完整 Jacobian 永远更好”的泛泛表述。

**inverse-scattering-specific 内容：存在。** 主要在 \(DX[\delta\chi]E_p\) 的矢量/材料支持限制，DDA 的材料链式导数，第 6 节 TSOM right-input 与 material-mixed output 的单向性，以及第 9 节隐藏电流对可测更新的反馈。这些比通用 trace weighting 更接近可投稿主贡献。

**strong TAP 理论核心：已形成可证明、可实现的候选核心，但不是期刊录用或已证首创的声明。** 单独 (3.1)、Davis–Kahan、Hutchinson 型估计、normal-equation perturbation 都不够。应把“Maxwell 材料注入—TSOM 单向保护缺口—实际 reduced dynamics—更新证书”作为主线，把这些通用工具作为支撑。最近邻 S-DBIM、参数导数保持的 primal–dual/interpolatory ROM、DINO、随机 ROM 误差估计、operator-valued inverse reduction 均需正面比较。第4A节还给出一个连续被动 Maxwell 的精确反转族，补上任意矩阵反例之外的物理存在性。

**仍需 3D Maxwell 验证的内容：**
- 在真实而非任意矩阵 \(G_D,G_S\)、物理介电常数 metric 与实际 \(Da\) 下，反转/transfer 是否具有非微小绝对测量能量和实际材料意义；
- 非共振与强多重散射两组中，TSOM bounds 是否够紧，weighted mixing 是否能预报真实 invariant-subspace rotation；
- 使用**真正** reduced forward 与匹配 adjoint，(9.4)/(9.8) 是否能在有用成本下通过 (7.4)，以及和 field-only、full-Jacobian、tangent-aware alternatives 的公平比较；
- 8–64 probes 的认证率、unresolved 率、误选率、总 PDE/adjoint/matvec 成本和节省是否值得；不得只报成功子集；
- 独立 discretization 或 refinement 下，物理排序/材料更新结论是否保持。有限 DDA 证明与 mesh convergence 要分开。

这些实验决定工程强度和论文说服力，不改变已经证明的有限维恒等式；未完成实验不能被写成已实现的加速/精度提升。

## 引用使用声明
(1.2) 属于经典 contrast-source / Newton–Kantorovich / DBIM sensitivity 结构；(6.1) 来自真实 TSOM retained-span 构造，不是 A9 首创；第 5、7、8 节分别运用标准矩阵扰动、强凸二次优化和随机二次型工具。具体来源、研究日期、已读/未读层级以及逐贡献重合度见 `A9_THEORY_CLAIM_LEDGER.md`，该文件也是本文的参考文献目录。
