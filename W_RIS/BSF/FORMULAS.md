# RIS 下表面：弹道光与散射光的分解

## 1. 本次定义与物理参数

使用已保存的 `W_RIS/results/wris_results.mat` 中的 **18 个原始快照**，逐幅沿用其配置、Tx/RIS 三维位置与姿态、1 cm 网格、首次到达边界和有限板面。以下 (I) 均是板面辐照度，单位 W/m²，不是复振幅。

原始配置：(P_0=1\,\mathrm W\)，真空波长 (\lambda_0=532\,\mathrm{nm}\)，水折射率 (n_w=1.33\)，束散参数 (\vartheta=0.01\,\mathrm{rad}\)，(w_0=\lambda_0/(\pi n_w\vartheta)=1.273239\times10^{-5}\,\mathrm m\)。水质 c0.151：(a=0.114\,\mathrm{m^{-1}}\)，(b=0.037\,\mathrm{m^{-1}}\)，(c=a+b=0.151\,\mathrm{m^{-1}}\)，HG 参数 (g=0.924\)。RIS 面积 1 m × 1 m。

定义

\[
\boxed{I_{\rm BSF}(x_R,y_R;t)=I_b(x_R,y_R;t)+I_s(x_R,y_R;t).}
\]

(I_b\)：到达板面前 **零次散射**；(I_s\)：到达板面前 **至少一次散射**。两者是互斥路径集合，散射光中不再重复计入弹道光。

本次比较使用原结果的 `noTurbulence`。若比较某一冻结 EGG 增益 (H(t)\) 下的 `irradiance`，必须两边同时乘以同一个 (H(t)\)，即 (I^H_b=HI_b\)、(I^H_s=HI_s\)。不能把未经增益处理的公式与含增益的 MC 比较。

## 2. 将论文式（3-18）直接拆开

参考来源：`3yimiao/3yimiao.pdf` 第 7 个 PDF 页面（论文印刷页 33）的式（3-18），以及 `3yimiao/fig3-10/M1_bsf_240729_7_2.m`。本次只计算 Tx → RIS 下表面这一段，故不照搬 fig3-10 中的空气段、折射段和接收端变量。

RIS 上一点的世界坐标和 Tx 光学坐标为

\[
\mathbf q=\mathbf p_R+R_R[x_R,y_R,0]^T,\qquad
(\xi,\eta,z)^T=R_T^T(\mathbf q-\mathbf p_T),\qquad \rho=\sqrt{\xi^2+\eta^2}.
\]

这里 (R_T,R_R\) 均是局部坐标到世界坐标的主动旋转；沿用项目定义，无需重新定义欧拉角。令

\[
\mu(\mathbf q)=\frac{\mathbf n_R\cdot(\mathbf q-\mathbf p_T)}{|\mathbf q-\mathbf p_T|},\quad
W(z)=\sqrt{w_0^2+\vartheta^2z^2},\quad
H(\kappa,z)=bk_1\frac{e^{k_2z\kappa}-1}{k_2\kappa},\quad H(0,z)=bk_1z.
\]

在论文的特殊坐标下，(\xi=x_R,\eta=-y_R\cos\theta,z=L_R-y_R\sin\theta\)，与原式一致。原代码中的 `ek2=1.05` 对应这里的 (k_1\)，`k1=-0.07973` 对应这里的 (k_2\)，切勿按变量名字对调。

对式（3-18）的散射指数使用恒等式 (e^H=1+(e^H-1)\)，并做径向傅里叶逆变换，得到

\[
\boxed{I_b^{(3-18)}=\frac{2P_0\mu}{\pi W(z)^2}
\exp\!\left[-cz-\frac{2\rho^2}{W(z)^2}\right]}
\]

\[
\boxed{I_s^{(3-18)}=\frac{P_0\mu e^{-cz}}{2\pi}
\int_0^\infty e^{-W(z)^2\kappa^2/8}
\left[e^{H(\kappa,z)}-1\right]J_0(\kappa\rho)\,\kappa\,d\kappa.}
\]

两者相加就是式（3-18）的二维 S-BSF，有限束腰仅用于匹配原 MC 的源。`bsf_paper.m` 和 `bsf_paper_map.m` 实现这些公式，没有调整 (k_1,k_2\) 去拟合热力图。`expm1` 和零频极限避免原积分核的 (0/0\)，复合 Simpson 积分代替逐点无穷域二维复积分；不使用 `abs` 掩盖负值或积分误差。

原式的小角近似、拟合散射核和单一几何投影因子并不等同于完整 HG 输运。尤其 (k_1=1.05\) 使其零频散射核不是严格归一化的相函数；在固定轴向截面上，该公式预测总功率 (P_0e^{(-c+bk_1)z}\)，不能将其当作严格能量守恒的精确模型。本项目的实际 RIS 是斜面并可能截光，因此还须对真实板面积分。

## 3. 与原 MC 同模型的弹道项

原 MC 的源不是只在一个点发光。它独立抽样束腰位置 \(\mathbf u=(u,v)\) 和横向斜率 \(\boldsymbol\alpha\)：

\[
f_u(\mathbf u)=\frac{2}{\pi w_0^2}e^{-2|\mathbf u|^2/w_0^2},\qquad
f_\alpha(\boldsymbol\alpha)=\frac{2}{\pi\vartheta^2}e^{-2|\boldsymbol\alpha|^2/\vartheta^2},\qquad
f_{u,\alpha}=f_uf_\alpha.
\]

对给定发射位置 \(\mathbf p_0=\mathbf p_T+R_T[u,v,0]^T\)，记

\[
\mathbf r_0=R_T^T(\mathbf q-\mathbf p_0)=(\xi-u,\eta-v,z)^T,\quad
\ell_0=|\mathbf r_0|,\quad \mathbf n_T=R_T^T\mathbf n_R.
\]

斜率到 RIS 真面积的 Jacobian 是

\[
\left|\frac{\partial(\alpha_x,\alpha_y)}{\partial A_R}\right|
=\frac{\mathbf n_T\cdot\mathbf r_0}{z^3},\qquad
\boldsymbol\alpha=\frac{(\xi-u,\eta-v)}z.
\]

所以本次用于高精度结果的弹道项为

\[
\boxed{I_b(\mathbf q)=P_0\!\int_{\mathbb R^2}\! f_u(\mathbf u)
f_\alpha\!\left(\frac{\xi-u}{z},\frac{\eta-v}{z}\right)
e^{-c\ell_0}\frac{[\mathbf n_T\cdot\mathbf r_0]_+}{z^3}\,
\mathbf 1_{z>0,\ell_0\le L_{\max}}\,d^2\mathbf u.}
\]

`bsf_ballistic.m` 用 3 × 3 高斯求积计算极小有限束腰积分，板面单元平均使用 6 × 6 Gauss-Legendre 求积，并以 4 × 4 结果检查收敛。源束腰仅 12.7 μm，远小于约 0.1 m 的接收束宽。

若取点源极限 \(w_0\to0\)，上式化为可直接使用的闭式：

\[
\boxed{I_b^{\rm point}(\mathbf q)=\frac{2P_0}{\pi\vartheta^2}
\exp\!\left[-\frac{2\rho^2}{\vartheta^2z^2}-c\sqrt{z^2+\rho^2}\right]
\frac{[\mathbf n_R\cdot(\mathbf q-\mathbf p_T)]_+}{z^3}.}
\]

它保留实际光程衰减与精确落点 Jacobian。与论文弹道项相比，不能再在计数结果外额外乘一次入射余弦。

## 4. 完整 HG 散射项：精确输运公式

归一化 HG 相函数为

\[
p_g(\mathbf s\cdot\mathbf s')=
\frac{1-g^2}{4\pi[1+g^2-2g(\mathbf s\cdot\mathbf s')]^{3/2}},\qquad
\int_{S^2}p_g\,d\Omega=1.
\]

定义 \(\Gamma_m\) 为恰好经过 \(m\) 次散射、之前未穿越 RIS 首次计数平面或源平面、最后首次到达 RIS 下表面的路径集合。第 \(j\) 段散射前光程为 \(s_j\)，末次散射后到平面的剩余距离为 \(\tau_m\)，最终板面落点为 \(\mathbf H_m\)。令 \(\delta_A\) 表示按板面真面积定义的二维 Dirac 密度。与原 MC 同物理模型的散射 BSF 是

\[
\boxed{
I_s(\mathbf q)=P_0\sum_{m=1}^{\infty}
\int_{\Gamma_m}f_u(\mathbf u)f_\alpha(\boldsymbol\alpha)
\left[\prod_{j=1}^{m} b\,e^{-cs_j}
p_g(\mathbf s_{j-1}\cdot\mathbf s_j)\,ds_j\,d\Omega_j\right]
e^{-c\tau_m}\,\delta_A(\mathbf q-\mathbf H_m)
\,d^2\mathbf u\,d^2\boldsymbol\alpha .}
\]

这里 \(\mathbf s_0=R_T[\alpha_x,\alpha_y,1]^T/\sqrt{1+|\boldsymbol\alpha|^2}\)。有限 RIS 的裁剪直接作用于 \(\mathbf H_m\)，世界坐标、6 自由度姿态与偏轴位置全部进入路径边界，不用单一入射角替代整个几何。

这也是辐射输运方程的碰撞次数展开：\(L_s=\sum_{m\ge1}(\mathcal T\mathcal S)^mL_b\)，\(I_s=\int_{\mathbf s\cdot\mathbf n_R>0}L_s(\mathbf q,\mathbf s)(\mathbf s\cdot\mathbf n_R)d\Omega\)。它是积分形式的精确 BSF，**不是新的闭式近似**。计算结果中的 \(I_s\) 仍由独立随机数值积分估计；不能把它宣传为“无需 MC 的纯解析散射公式”。

## 5. 高精度求值与原 MC 的区别

`bsf_scatter_transport.m` 独立求解上面的 \(I_s\)。先按原源分布抽样 \((\mathbf u,\boldsymbol\alpha)\)，设无散射初始光线到首个 RIS 平面的距离为 \(L\)（并限制总光程）：

\[
Q=1-e^{-bL},\quad
s_1=-\frac1b\ln(1-UQ),\quad
w_1=Qe^{-as_1}.
\]

这是强制首次散射，保证样本用于散射项。随后散射间距取 Exp(\(b\))，每段连续吸收乘 \(e^{-as}\)，方向按完整 HG 更新。因

\[
(b e^{-bs})e^{-as}=b e^{-cs},
\]

它与原 MC 的 Exp(\(a+b\)) 自由程加隐式吸收权重具有相同路径积分期望，权重和随机实现不同。边界前的段按真实段长处理；不发生虚构的末段碰撞。对第 \(j\) 个板面单元：

\[
\widehat I_{s,j}=\frac{P_0}{N\Delta A}\sum_{k=1}^N
w_k\,\mathbf1_{\mathbf H_k\in A_j},\qquad
\widehat I_{{\rm BSF},j}=\overline I_{b,j}+\widehat I_{s,j}.
\]

没有读取原 MC 的强度来拟合 \(I_s\)，没有重用原图的随机种子，没有以 \(I_{\rm MC}-I_b\) 定义散射项。原结果始终保留，用于独立验证。

与原配置一致，计算有 200 m 总光程上限和 1000 次碰撞上限。代码报告被截断的散射权重；它量化该计算边界下的损失，不把有限次模拟称为已计算无限多阶。原 MC 的俄罗斯轮盘赌是无偏的；本实现用连续吸收权重，不需要轮盘赌。

对当前 \(a=0.114\,\mathrm{m^{-1}}\)，总光程截断导致的未计功率还具有直接上界 \(P_0e^{-aL_{\max}}=1.2534\times10^{-10}\,\mathrm W\)。在 200 m 内，纯散射泊松过程的均值仅为 \(bL_{\max}=7.4\)，1000 次散射上限的尾部远小于该光程截断误差。这两项均远小于本次数值目标。

## 6. 0.1% 的确切含义与验收

你确认的目标是 **公式数值求值的二维 RMS 标准误差 ≤0.1%，同时保留原 MC 噪声**：

\[
\epsilon_{\rm num}=\frac{\|\mathrm{SE}(\widehat I_s)\|_F}
{\|\overline I_b+\widehat I_s\|_F}\le10^{-3}.
\]

这是整体二维图的估计 1σ 标准误差，不是逐个暗像素相对误差的上界，也不是 95% 置信误差 ≤0.1%。原始 \(10^6\) 光子 MC 的像素噪声不会因此消失。原始像素比较不平滑、不重标功率、不对齐平移、不缩放峰值，直接报告：

- 总功率相对误差、真实单元峰值相对误差。
- 全板面相对 L2/L1 误差、质心、二维协方差、固定半径内功率。
- 原 MC SE 与模型求值 SE；仅对原 MC 有效计数 ≥25 的单元作近似 z 分数诊断。
- 论文 BSF 相对原 MC 的误差，以及论文 BSF 相对高精度输运模型的误差，后者用于观察小角模型的系统偏差。

随机数值 SE 按每个发射样本的一、二阶矩计算，跨独立批次合并时保留批次间均值项；不能把随机采样误差除以批次数代替标准误差公式。自适应停止仍使用估计 SE，因此额外独立种子复验用于检查该精度是否可复现。

## 7. 后续 RIS 调相接口

采用图中的互不相干散射集合假设时：

\[
U_b(\mathbf q)=\sqrt{I_b(\mathbf q)}e^{i\phi_b(\mathbf q)},\quad
U_b^+(\mathbf q)=U_b(\mathbf q)e^{i\phi_{\rm RIS}(\mathbf q)},\qquad
\langle I_{RX}\rangle=|\mathcal P[U_b^+]|^2+I_{s,RX}^{\rm MC}.
\]

`bsf_phase_input.m` 仅提供可替换的标量基模高斯相位：

接口同时施加有限 RIS 孔径；在扩展计算网格上，板外的 RIS 复场为零。入射面密度 `bsf_ballistic` 仍可用于检查板外光斑，实际 RIS 总功率只积分板面区域。

\[
\phi_b=k_wz+\frac{k_w\rho^2}{2R(z)}-\arctan(z/z_R),\quad
k_w=2\pi n_w/\lambda_0,\quad z_R=\pi n_ww_0^2/\lambda_0,\quad
R(z)^{-1}=z/(z^2+z_R^2).
\]

**0.1% 辐照度求值精度并不验证相位精度。** 强度热力图不能唯一确定复场；上述相位是明确添加的相干高斯模型假设。若以后采用精确角谱传播或其他发射波前，应替换这部分相位并单独检验传播采样与相位收敛。

特别地，原 1 cm 面元是统计强度的网格，不能直接作为光学 FFT 的复场采样网格。倾斜载波沿 RIS 两轴的梯度是

\[
\boldsymbol\gamma=k_w\mathbf d_T^T R_R(:,1\!:\!2),\qquad
\Delta x_j<\pi/|\gamma_j|.
\]

仅 30° 入射时，对有载波梯度的轴已需要约 0.4 μm 间距；大角度时约 0.2 μm。这还没有计入高斯曲率和 RIS 调相梯度。`bsf_phase_input` 返回这项诊断及 `gridPhaseSamplingCertified=false`；复场数组是给定坐标处的点值，绝不意味着原 MC 网格满足传播采样定理。可以采用已知连续相位再采样、解析高斯传播或经过验证的倾斜载波包络传播。

散射接口 `bsf_scatter_transport(...,SavePackets=true)` 返回世界坐标 `pos`、单位入射方向 `dir`、每包实际功率 `powerW` 和 RIS 局部落点 `localXY`。这些权重已经含 \(P_0/N\)，后续不可再次除以光子数。光程只作输运记录，不能将随机路径长度直接解释为已知的相干散射相位。不能把 \(\sqrt{I_b+I_s}\) 全部赋同一个高斯相位。

本次只建立并验证 RIS 下表面的输入；尚未验证 RIS 相位梯度改变散射出射方向的具体光学模型及 RX 波动/光线耦合传播。
