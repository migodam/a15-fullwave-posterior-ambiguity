# 与 A9/A12 的接口

固定空间 VIE 沿用 M=I−DχGD；对允许参数方向 h，材料注入为 D_E Jχ h，电流变化为 M⁻¹D_EJχh，外部过滤为 GS。本轮 Fig11 把注入、放大、泄漏、任务导数四项显式分开，这是 material-current transfer 的自然数学接口。

圆柱质量归一化坐标随参数变化，不能把其 assembled RHS derivative 直接冒称固定空间材料注入。为避免该混淆，物理注入诊断在固定VIE网格上计算并标明其数值局限。

当前只足以 discussion/methodological connection；没有确认 population mechanism 或独立TAP论文贡献。未调用N-Port后端，不强行统一成N-Port成果，也不把标量2D结论升级为3D Maxwell。
