# Gate A：未建立（参考不足）

**未通过执行 gate；物理存在性结论未决。** locked test 未开启。

主噪声 validation 41³ 的未认证诊断：R_Q=0.01292496，90%覆盖率=86.1979%，posterior-mean MSE=0.01345738；ambiguous query 为0/384，三个 seed 均为0。33³ 为1/384。它们既未接近10%门槛，也未通过 reference 收敛；不能把0/384当作真实总体零概率。Wilson 和 bootstrap 区间仅包括有限查询抽样误差，不包括积分偏差。

1/2/5%噪声和1/2/4/8通道的全部 validation 诊断均保存。5%基础 acquisition 的未认证比例上升，说明门槛结果对噪声敏感；不得据此把主噪声改成5%。

参考失败明细见 POSTERIOR_REFERENCE.md。没有只保留通过积分的94个 query 重算总体，也没有把它们当新 population。后续 positive claim 均阻断。
