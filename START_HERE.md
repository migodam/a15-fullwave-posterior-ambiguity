# A15 — 网页版 GPT 阅读入口

这是公开、无需登录的研究档案。请先读下面的结论边界，再按顺序读取报告和证据。理论和prompt是待验证输入，不是实验结论。

**本轮 verdict：NO-GO，针对当前 canonical family 和证据强度；后验参考不足不等于证明物理机制不存在。**

- Gate A：UNRESOLVED_REFERENCE_INSUFFICIENT；基础观测下仅94/384验证query通过全部33³/41³积分差异容差。
- Gate B：FAIL（未满足通过条件，含参考未决）；当前两模态解释未达到噪声精度。
- Locked test 未开启；FLOW NOT RUN — GATE CLOSED。不能写FLOW NOT NEEDED。
- 表上risk、prevalence、采集排序和相关图均为PROVISIONAL。44项测试通过只是实现证据。

## 推荐读取顺序（都是直接可读的公开文件）

- [完整报告，含Q1–Q10](https://raw.githubusercontent.com/migodam/a15-fullwave-posterior-ambiguity/main/A15/reports/A15_FULLWAVE_POSTERIOR_REPORT.md)
- [单文件报告合集，适合一次读取](https://raw.githubusercontent.com/migodam/a15-fullwave-posterior-ambiguity/main/ANALYSIS_BUNDLE.md)
- [冻结协议](https://raw.githubusercontent.com/migodam/a15-fullwave-posterior-ambiguity/main/A15/EXPERIMENT_PROTOCOL.md)
- [精确配置](https://raw.githubusercontent.com/migodam/a15-fullwave-posterior-ambiguity/main/A15/configs/protocol.json)
- [参考不足与积分验收](https://raw.githubusercontent.com/migodam/a15-fullwave-posterior-ambiguity/main/A15/reports/POSTERIOR_REFERENCE.md)
- [Gate A](https://raw.githubusercontent.com/migodam/a15-fullwave-posterior-ambiguity/main/A15/reports/GATE_A.md)
- [Gate B](https://raw.githubusercontent.com/migodam/a15-fullwave-posterior-ambiguity/main/A15/reports/GATE_B.md)
- [物理求解器验收](https://raw.githubusercontent.com/migodam/a15-fullwave-posterior-ambiguity/main/A15/reports/PHYSICS_VALIDATION.md)
- [等预算采集比较](https://raw.githubusercontent.com/migodam/a15-fullwave-posterior-ambiguity/main/A15/reports/ACQUISITION_COMPARISON.md)
- [11类图和视觉审查](https://raw.githubusercontent.com/migodam/a15-fullwave-posterior-ambiguity/main/FIGURES.md)
- [代码、数据、文献限制完整索引](https://raw.githubusercontent.com/migodam/a15-fullwave-posterior-ambiguity/main/FILE_INDEX.md)
- [原始实验要求](https://raw.githubusercontent.com/migodam/a15-fullwave-posterior-ambiguity/main/A15/A15_prompt.md)
- [待验证理论](https://raw.githubusercontent.com/migodam/a15-fullwave-posterior-ambiguity/main/A15/A15.md)

## 深入核查

- [全波圆柱物理后端](https://raw.githubusercontent.com/migodam/a15-fullwave-posterior-ambiguity/main/A15/src/cylinders.py)
- [复用VIE与材料切线](https://raw.githubusercontent.com/migodam/a15-fullwave-posterior-ambiguity/main/A15/src/physics.py)
- [后验工具](https://raw.githubusercontent.com/migodam/a15-fullwave-posterior-ambiguity/main/A15/src/posterior.py)
- [参考gate逐query原始记录](https://raw.githubusercontent.com/migodam/a15-fullwave-posterior-ambiguity/main/A15/results/reference_validation_gate.json)
- [population诊断](https://raw.githubusercontent.com/migodam/a15-fullwave-posterior-ambiguity/main/A15/results/population_diagnostics.json)
- [机器可读最终gate](https://raw.githubusercontent.com/migodam/a15-fullwave-posterior-ambiguity/main/A15/results/gates.json)
- [公开档案清单与SHA256](https://raw.githubusercontent.com/migodam/a15-fullwave-posterior-ambiguity/main/PUBLIC_MANIFEST.json)
- [发布说明与脱敏差异](https://raw.githubusercontent.com/migodam/a15-fullwave-posterior-ambiguity/main/PUBLICATION_NOTES.md)

分析时请区分：求解器准确性、后验积分准确性、物理机制解释、总体统计、采集设计、Flow必要性。不要把未认证的bootstrap区间当作包含积分误差，不要把对称性暗方向自动解释成task-relevant uncertainty。大型NPY/NPZ可下载复算；网页版阅读优先使用报告和JSON。
