# Fig1 修订与验证

本修订只针对Fig1。原来的第一列是带玻璃/器件掩模的零相位RIS；现已替换为 **No RIS (MC)**。
原有的其他图、相位优化与固定接收机比较继续保留各自协议。

- 重画入口：在RIS_phase目录运行`run_fig1_revision`。
- 成品：`plots/stage2/Fig1_irradiance.png`、同名FIG。
- 原始结果：`Fig1_data.mat`（完整MC透射光子、12幅主图和3幅全MC尾部图）、`Fig1_summary.csv`。
- 物理与显示检查：本目录`verification.json`、`mc_verification.csv`、`display_verification.csv`。
- 测试：`tests.csv`，5项新物理测试及19项原有测试全部通过；代码分析见`code_analyzer.csv`。

## 1. 束腰椭圆

从主图光强的功率归一化协方差矩阵求主轴：

$$
C=V\operatorname{diag}(\lambda_1,\lambda_2)V^T,\qquad
w_1=2\sqrt{\lambda_1},\quad w_2=2\sqrt{\lambda_2},\qquad
\boldsymbol r(\theta)=\boldsymbol c+V
\begin{pmatrix}w_1\cos\theta\\w_2\sin\theta\end{pmatrix}.
$$

虚线椭圆满足 $(\boldsymbol r-\boldsymbol c)^TC^{-1}(\boldsymbol r-\boldsymbol c)=4$。
不再用 $\sqrt{2\operatorname{tr}C}$ 的圆代替椭圆。主轴使用细点虚线。
这里是二阶束宽约定；对于非高斯场，不将其误称为逐点精确的1/e²等照度轮廓。

## 2. 无RIS的MC传播

直接读取`study.mat`中各距离保存的下表面MC光子：位置、方向、功率与弹道/散射标记。
三个距离的输入均为300000次发射、132881个到达共用1 m×1 m下表面评分窗的光子包；
评分窗是既有输入数据的范围。新的裸海面传播不再施加器件透明区、填充因子或玻璃板遮挡，
也不人为补齐保存输入范围外未记录的水下光子。

每个光子按Snell定律从水折射入空气，权重只乘水—空气Fresnel功率透射率：

$$
\boldsymbol q=n_w\boldsymbol d_{w,\parallel},\qquad
\boldsymbol d_a=R_R
\begin{pmatrix}\boldsymbol q/n_a\\\sqrt{1-|\boldsymbol q|^2/n_a^2}\end{pmatrix},
\quad P_{a,k}=T_{wa}(\boldsymbol q_k)P_{w,k}.
$$

全反射分量透射率为0。空气传播到与其他三组完全相同的世界水平接收面：

$$
t_k=\frac{z_R-z_{s,k}}{d_{a,k,z}},\qquad
\boldsymbol x_{R,k}=\boldsymbol x_{s,k}+t_k\boldsymbol d_{a,k}.
$$

没有LC分片相干叠加，因此不会产生器件阵列的周期格子。图像由MC落点的加权密度估计，
Scott带宽由弹道样本协方差与有效样本数给出；保存原始像素功率，未用高斯拟合替换MC结果。
1 W发射、水参数、入射方向、接收高度均沿用输入配置。图1显示湍流均值1时的平均光强，
与其他三组原光强图的含义一致。

## 3. 完整显示与尾部

其他三组主图显示全部原始波动场像素边界并留3%余量，椭圆和主轴端点也纳入坐标范围。
无RIS主图用四倍弹道坐标标准差加六倍密度带宽确定图域，覆盖本次输入中100%的弹道光子。
主图包括约88.3%–88.7%的总透射功率，其余是远处散射尾部。
每个直接通信面板的“Full MC”小图用可逆的asinh坐标压缩远场范围，包含**100%透射MC包功率**；
强度由每格真实物理面积计算，没有丢弃远处光子或把图外尾部当成零。
全部面板和小图统一使用同一个W/m²色标，未逐图归一化。

## 4. 数值结果

|接收高度(m)|No RIS峰值(W/m²)|Steering-only峰值|Flat top峰值|Robust ABER峰值|
|---:|---:|---:|---:|---:|
|5|1.1412|3.6923|4.0919|3.9641|
|7.5|0.64885|2.5684|4.9400|3.3171|
|10|0.41792|1.8914|5.2568|4.3140|

直接通信峰值均低于三组RIS；该排序在带宽减半、独立光子奇偶半样本中均保持。
独立半样本峰值差小于0.90%，二阶等效宽度差小于1.63%；带宽减半时峰值变化8.6%–9.2%，
这是有限MC密度估计的带宽敏感性，原始结果随数据保存。密度平滑边界功率损失小于0.084%。
水—空气能量账本残差约1.33e-15 W，全部透射光子的接收总功率为0.24558 W。
**较低的峰值光强不等于较低的整面总功率**：无RIS没有器件填充与玻璃损耗，功率分布更宽，
整面功率高于带器件的三组是模型的自然结果，未为强制排序而缩放。

旧的等效半径实线圆、粗黑轮廓/白色十字以及Fig1零相位器件板分支已从绘图实现中删除。
`rp_zero_phase`等仍被其他图使用的算法不是本次冗余代码，继续保留。

初版全阶段`audit/stage2/audit_artifacts.py`锁定的是当时README哈希；README已于本轮任务开始前的19:05修改，因此其历史哈希断言不用于此次Fig1验收。该断言未被放宽或替换；本次独立来源、结果与输入比对记录于本目录manifest。
