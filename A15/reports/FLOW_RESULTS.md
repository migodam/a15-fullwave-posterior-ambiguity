# Flow 状态

**FLOW NOT RUN — GATE CLOSED**

Gate B 未通过，reference 不足；没有训练 FM、normalizing flow，没有生成50k/200k数据，没有使用任意材料latent或独立物理heads。

Laplace/1–4 Gaussian mixture 的正式精度与成本比较尚未执行，故不能写 FLOW NOT NEEDED。SMC 只承担独立参考诊断，不是已认证 posterior truth。Flow必要性、分支概率、amortization break-even均未验证。所有已生成场与电流来自同一材料realization的真实solver。
