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
