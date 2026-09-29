# A15 后验参考验收（参考不足）

**状态：未通过。锁定测试未开启，Flow gate 关闭。**

两个全波表均完整生成：33³=35,937、41³=68,921 个合法材料节点，每节点三个频率，合计 314,574 次正演。没有删除失败节点重新归一化。源代码、协议及每个表的原始文件 SHA256 保存在各表 status.json；本地保留完整表。

独立连续 truth 的 validation 为三 seed ×128 queries。主噪声、基础 acquisition 的严格逐 query 比较只有 94/384 通过全部容差；最大 Q 均值差 0.11836、最大区间端点差 0.15384。新增 acquisition 同样未收敛。详见 `results/reference_validation_gate.json`。有效积分节点数不是 ESS 收敛证明。

已追加两个预定 query 的独立 tempered SMC（各两个独立 256 粒子运行，共 32,768 次全波求解），采用正确 conditional ESS、自适应温度、pCN 和祖先诊断。温度到达 1 不等于参考收敛；重复运行的方差差异仍不可忽略。

正权自适应 cubature 保留全部 prior 支持，对两个相同 query 分别使用约 9.7k 全波求解。每个叶 cell 使用 3³/5³ Gauss–Legendre 规则，所有叶节点及权重保留。没有只截取 MAP 周围区域。预算停止时嵌套 raw-moment 相对误差估计仍约 4.3% 和 8.6%；区间与方差也不能仅用这些估计认证。完整轨迹与节点位于 `adaptive_reference.json` / `adaptive_q*.npz`。

这是积分分辨率和有限预算的失败，**不是物理不存在的证据**。表上 population risk、prevalence、采集排序只作 validation 诊断，不能使用它们宣告 Gate A/B PASS。禁止开启 locked test，禁止训练 Flow。其余预算用于规定的 coupling 消融和数值/机制边界检查，不修改 prior 或阈值。
