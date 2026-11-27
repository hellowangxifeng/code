# 自由漂浮 RIS 六自由度模型与复现说明

2026-11-05：实体更新为 1×1×0.01 m 玻璃与四周浮框。当前质量预算、安装高度、干湿接触边界和水动力范围见 [RIS浮框配平说明.md](RIS浮框配平说明.md)。下文的六自由度积分和辐射记忆仍使用，BEM 参考域按静态浸水状态选择；当前使用填充等效水线域。

本次实现是**完整六自由度刚体运动学、刚体惯性耦合，以及有明确边界的近似水动力模型**。它不是已验证的六自由度水动力数据库。偏航来自侧壁入射压力力矩和刚体姿态动力学；没有人为指定 `yaw(t)`，没有锚泊、水平恢复弹簧或偏航恢复弹簧。

当前运行和出图方法见 [README.md](README.md)。旧输出与审查报告只作为历史资料。

## 1. 审查结论与修改范围

| 模块           | 审查发现                                                                                   | 本次处理                                                                                                     |
| -------------- | ------------------------------------------------------------------------------------------ | ------------------------------------------------------------------------------------------------------------ |
| `RIS_float.m`  | 原`q=[heave,pitch,roll,surge,sway]`；三个竖向模态与两个水平解析解相互独立；没有 yaw 方程   | 替换拼接运动解为六自由度 RK4 积分，保留原薄板 BEM 核、奇偶对称化、辐射／绕射求解及几何缓存                   |
| 坐标与刚体惯性 | 原 pitch 为 X 高度斜率，右手 pitch 应取其负值；没有 Euler／刚体耦合                        | 统一右手 ZYX，使用单位四元数；当前惯量按玻璃/框体/附件/配重组合及平行轴定理计算；保留刚体陀螺、平移坐标耦合 |
| 波浪载荷       | 侧壁压力只用于水平合力，未积分偏航力矩及水平力作用高度                                     | 同一入射压力解析积分产生 Fx、Fy、侧壁 Mx、My、Mz；Mz 经独立四壁高阶数值积分校验                              |
| 辐射           | 原模型依靠每频率稳态 RAO，不适于直接处理任意初值与耦合运动                                 | 从原 BEM 的 A、B 拟合被动辐射记忆状态空间，输出拟合误差；原频域 RAO 仍保留为诊断量                           |
| 波浪生成       | 两处重复谱叠加，端点权重不同；原生成器改写全局 RNG，曲面逐点四重循环                       | 共享`RIS_wave_components` 与分块向量化波场；统一梯形端点半权重，删除旧端点加倍分支；使用局部随机流并核对连续谱积分 |
| 平面拟合       | 使用正规方程`X'*X`，条件数较差                                                             | 改用`X\z`，保留原接口、法向及单平面输出                                                                      |
| 代码清理       | 不可达绘图函数、只用于深水公式但不起作用的水深变量、未使用 sinc 辅助函数、调色函数名不匹配 | 删除不可达／重复逻辑及未调用的调色、四相位部署代码                              |
| 参数与绘图     | 物理参数集中程度不足；图表只有五自由度；输出文档与源码不同步                               | 统一默认参数，按四风速、三风向与 30° 质心轨迹生成表格图、CSV、FIG 和 MAT                               |

历史海浪数据、PDF 和既有图件保留；当前批量输出在 `FigOutput_RIS_WindDirection/`。

## 2. 坐标、符号与接口

惯性系 W：静水面 Z=0，Z 向上，X、Y 为水平右手系。本体系 B：原点在质心 G，X 沿 Lx，Y 沿 Ly，Z 指向板上表面。初始水平航向可由 `initialEulerAngles(3)` 指定；这是初始条件，后续角度由动力学积分。

标准输出：

\[
\eta=[x,y,z,\phi,\theta,\psi]^T
=[surge,sway,heave,roll,pitch,yaw]^T,\qquad
\nu_B=[u,v,w,p,q,r]^T.
\]

z 是相对于静水浮力平衡位置的质心升沉，正方向向上。三个角度均服从右手法则：正 roll 抬高 +Y 边；正 pitch 降低 +X 边；正 yaw 在俯视图使 +X 转向 +Y。风向参数保持原项目定义，传播主方向相对 +X 为 `90°-windAngleDeg`。

位置、方向的转换：

\[
R=R_z(\psi)R_y(\theta)R_x(\phi),\qquad
\boldsymbol x_W=\boldsymbol x_G+R\boldsymbol r_B,
\qquad\boldsymbol x_G=[center_x+x,center_y+y,z_G+z]^T.
\]

\[
R=\begin{bmatrix}
c_\psi c_\theta&c_\psi s_\theta s_\phi-s_\psi c_\phi&c_\psi s_\theta c_\phi+s_\psi s_\phi\\
s_\psi c_\theta&s_\psi s_\theta s_\phi+c_\psi c_\phi&s_\psi s_\theta c_\phi-c_\psi s_\phi\\
-s_\theta&c_\theta s_\phi&c_\theta c_\phi
\end{bmatrix}.
\]

向量由 W 到 B 使用 Rᵀ。旋转矩阵、质心、板上表面中心、四角和法向均取同一姿态；`normal=R*[0;0;1]`。状态积分使用 wxyz 单位四元数 Q：

\[
\dot Q=\tfrac12 Q\otimes[0,\boldsymbol\omega_B],\qquad
\dot{\boldsymbol x}_G=R\boldsymbol v_B.
\]

每个 RK4 完整步后归一化 Q。输出角度使用 ZYX 分解，yaw 作连续展开。Euler 角在 pitch=±90° 有表示奇异性，内部四元数没有该奇异性；大角运动本来也超出水动力近似范围。

兼容说明：`surge/sway/heave/roll/pitch` 字段继续存在；`q` 前五行位置保留，新增第六行 yaw。但 **pitch 已改为右手正方向**，旧程序若把 `pitch` 当作 X 高度斜率，应改用 `-pitch`。推荐新下游直接使用 `eta`、`rotationMatrix`、`bodyVelocity`。

`velocity` 是 `q` 的时间导数；`bodyVelocity(4:6,:)` 是本体角速度，并不等于 Euler 角速度。`massMatrix/stiffnessMatrix` 按 `qOrder` 排列；`rigidBodyMassMatrix/bodyMassMatrix` 按标准 `dofOrder` 排列，后者随姿态变化，为 6×6×Nt。原 `responseAmplitudeOperator` 三行保留为 `[heave,right-hand pitch,roll]` 的薄板频域诊断，`responseAmplitudeOperator6` 是完整参考线性载荷的近似六自由度 RAO；非线性时程不由这个 RAO 直接指定。

## 3. 质量、浮力与恢复力

\[
A_g=L_xL_y,\quad A_f=(L_x+2b)(L_y+2b)-A_g,\quad
m=\rho_g A_g h+\rho_f A_f H+m_{payload}+m_{ballast}.
\]

浮框吃水 d 由 `m/ρ=Af*d+Ag*clip(d−a,0,h)` 解出；a 是玻璃底面相对框底的固定安装高度，要求 d<H。玻璃底面 z=a−d，顶面 z=a−d+h。默认 a=0.055 m，m=29.871 kg，d=0.058432249322 m。质心和惯量用各部件质量、位置和几何计算，不再等同于均质玻璃板的厚度中心与惯量。

\[
K_z=\rho g A_{wp},\quad
K_\phi=\rho g I_{wp,x}+mg(z_B-z_G),\quad
K_\theta=\rho g I_{wp,y}+mg(z_B-z_G).
\]

玻璃恰好贴水时使用干侧环形水线面，`Awp=Af`，二次矩为外矩形减内矩形；玻璃部分静态浸没时，水线面含玻璃。接触点存在左右导数差异，动态模型不切换干湿状态。薄板小姿态静水恢复关系延伸为连续势能：

\[
U=\tfrac12K_z z^2+\tfrac12K_\phi n_{By}^2+\tfrac12K_\theta n_{Bx}^2,
\quad n_B=R^Te_z.
\]

\[
F_{hs,W}=-K_z z e_z,\qquad
\tau_{hs,B}=n_B\times[K_\theta n_{Bx},K_\phi n_{By},0]^T.
\]

这在小角度恢复为 `[-Kφ φ,-Kθ θ,0]`，且 `e_zᵀ R τhs,B=0`，不会引入惯性竖轴上的静水偏航恢复。surge、sway、yaw 没有弹簧约束。这个势能是小姿态近似的能量一致延伸，并非完整瞬时湿表面浮力积分。

## 4. 六自由度动力学及耦合

为避免把线性薄板辐射模型错误地变成未经求解的随体大转角辐射模型，附加惯性、激励及记忆均定义在初始**水平参考系 H**，其姿态 `R_H=Rz(ψ0)`。刚体惯性和姿态保留完整非线性。

设 `A∞,z`、`A∞,φ`、`A∞,θ` 来自 BEM 阻抗拟合，定义：

\[
D_t=mI+A_{\infty,z}e_ze_z^T,\quad
I_W=RI_BR^T,\quad
A_r=R_H\operatorname{diag}(A_{\infty,\phi},A_{\infty,\theta},0)R_H^T.
\]

惯性系动力学：

\[
D_t\dot v_W=F_{wave,W}+F_{hs,W}+F_{rad,W},
\]
\[
(I_W+A_r)\dot\omega_W+
\omega_W\times(I_W\omega_W)
=\tau_{wave,W}+\tau_{hs,W}+\tau_{rad,W}.
\]

对应代码中的本体系方程：

\[
M_B(R)\dot\nu_B+c_B(R,\nu_B)
=f_{wave,B}+f_{hs,B}+f_{rad,B},
\]
\[
M_B=\operatorname{blkdiag}(R^TD_tR,R^T(I_W+A_r)R),\qquad
c_B=\begin{bmatrix}
(R^TD_tR)(\omega_B\times v_B)\\
\omega_B\times(I_B\omega_B)
\end{bmatrix}.
\]

其中包含平移坐标耦合、完整 Euler 陀螺项及随姿态变化的有效惯量。附加质量固定于参考系，不能套用“常数随体附加质量”的 Coriolis 矩阵；那是另一种假设。初始参考水动力的小运动范围应单独检查。

偏航不是指定曲线：三个角速度由上式共同求解，再通过四元数得到 yaw。一般姿态下 `yawRate != r`。即使惯性竖轴力矩接近零，耦合横纵摇也可能产生小的 Euler yaw 几何变化；轴向规则波的镜像对称算例仍应保持零 yaw。

## 5. 波浪激励与偏航力矩的来源

保留原深水 Neumann 方向谱、cos² 方向扩散和 MT19937 随机相位：

\[
\zeta=\Re\sum_j c_j e^{i(\omega_jt-k_{xj}X-k_{yj}Y)},\quad
k_j=\omega_j^2/g,
\]
\[
S(\omega)=3.05\frac\pi4\omega^{-6}\exp[-2g^2/(U^2\omega^2)],\quad
D(\vartheta)=\frac2\pi\cos^2\vartheta,\quad |\vartheta|\le\pi/2.
\]

当前默认波段为 **0.05–4 Hz**，频率与方向均使用梯形端点半权重。海面生成器和 RIS 求解共用 `RIS_wave_components`；已删除频率端点加倍的旧分支。`wave` 输入直接复用给定分量，风速标签不会对其重复缩放。默认 50×36 离散和 400 点加密相对连续谱积分的误差由 `audit_sea_spectrum` 检查。

竖向合力、roll／pitch 力矩：复用原薄板辐射／绕射离散方程，广义形状函数为 `[1,x,y]`。玻璃不浸水时只保留环框区域的面元，内孔的自由液面不计为实体承重面；静态浸没玻璃使用填充的等效参考域。环框参考域网格必须对齐内外边缘；填充域没有内孔对齐要求。当前默认 panelCount=80。映射到右手坐标后为 `[Fz,Mx,My]=[F1,F3,-F2]`。环形参考域保留双轴镜像对称，三个奇偶扇区仍解耦。这是薄水线近似，不是实际阶梯湿表面的三维求解。

侧壁载荷：下列公式计算单个矩形四壁压力。当前总载荷取外矩形四壁减去内孔四壁；如玻璃底面低于水线，则内孔压力积分仅减去框底至玻璃底面之间的部分。内外壁都以同一组合重心计算力臂。法向外指：

\[
\widehat p=\rho g a e^{kz}e^{-i(k_xx+k_yy)},\quad
\widehat F=-\int_{S_{side}}\widehat p\,n\,dS,
\quad\widehat M_G=-\int_{S_{side}}r_G\times(\widehat p n)\,dS.
\]

记 `sinc(x)=sin(x)/x`，定义：

\[
J=L_xL_y\,sinc(k_xL_x/2)\,sinc(k_yL_y/2),\quad
D_z=\frac{1-e^{-kd}}k.
\]

每米波幅的水平力、偏航力矩为：

\[
\widehat F_x=i\rho gD_z k_xJ,\qquad
\widehat F_y=i\rho gD_z k_yJ,
\]
\[
\boxed{\widehat M_z=\rho gD_z\left(k_x\frac{\partial J}{\partial k_y}
-k_y\frac{\partial J}{\partial k_x}\right)}.
\]

这个式子由压力的力臂积分导出，保留有限板长与波相位差。不能把合力始终作用在板中心而丢掉 Mz。沿对称轴入射时 Mz=0；方板对角线入射也为零；一般斜向有限波长不必为零。长波极限的矩形主项正比 `kx*ky*(Lx²-Ly²)/12`；方板该项抵消，只有更高阶有限波长项，偏航更弱。

压力加权作用高度：

\[
\ell_z=-\frac1k+\frac d{e^{kd}-1}-z_G,
\qquad M_{x,side}=-\ell_z F_y,\quad M_{y,side}=\ell_z F_x.
\]

代码在 k≈0 使用级数，避免相消。该侧壁水平压力力矩补到薄板竖向载荷上；未再次添加底部入射竖向力，避免与 BEM 重复。

全部激励在初始参考点和水平航向计算，随后作为惯性系时变力／力矩输入刚体方程。**没有通过移动位置重新计算激励、遭遇频率或瞬时航向散射。** 它们属于需一致补齐的高阶／大航向模型，当前用 `maxHeadingChangeDeg` 等诊断暴露适用性。

## 6. 辐射记忆与初值

原 BEM 给出三模态 `A(ω),B(ω)`；用非负最小二乘拟合：

\[
B(\omega)+i\omega A(\omega)
\simeq i\omega A_\infty+
\sum_l c_l\frac{i\omega+a_l}{(i\omega+a_l)^2+b_l^2},
\quad A_\infty,c_l\ge0,\ a_l>0.
\]

每个模态的记忆状态满足：

\[
\dot x_l=-a_lx_l-b_ly_l+u,\quad
\dot y_l=b_lx_l-a_ly_l,\quad f_{mem}=-\sum_lc_lx_l.
\]

输入为参考系的 `[vW,z,ωH,x,ωH,y]`；水平及 yaw 辐射没有求解。对应记忆能量 `Σcl(xl²+yl²)/2` 的耗散为 `Σal cl(xl²+yl²)≥0`。完整系统无波浪时应满足：

\[
\frac d{dt}(T_{RB}+T_{added}+U+E_{memory})
=-\sum_l a_lc_l(x_l^2+y_l^2)\le0.
\]

这使拟合被动性和能量测试有明确依据。极点网格、模态数是数值近似设置，不是物理阻尼标定。默认拟合并非精确 BEM：`radiationFitRelativeRmsError` 必须随结果检查；更细的时间步不能消除此误差。`radiationModes` 可用于单独加密，不能把拟合残差隐藏为零。

默认 `initialization='waveSteady'` 从拟合线性稳态初始化 heave、roll、pitch 和记忆，水平初速度保持原波致轨道速度的默认约定；yaw 与 yaw 初始角速度默认零。完整初始状态可通过 `initialPosition/initialEulerAngles/initialBodyVelocity` 指定。`initialization='rest'` 将状态及辐射记忆初始化为零，再施加已有波浪。显式覆盖运动初值时默认波场记忆仍保留；要做由静水启动的释放试验请选 `rest`。

## 7. 运行、数据与图

```matlab
suite = run_RIS_wind_comparison();
report = audit_sea_spectrum();
out = RIS_float(struct('windSpeed',9,'windAngleDeg',30,'duration',10,'plot',true));
```

批量入口按 3、5、9、12 m/s 四行与 30°、60°、90° 三列输出海面和 RIS 位置；第四列为 30° 质心三维轨迹。位置图统一为绝对时间 27 s，轨迹为 25–35 s。所有三维图 z 轴固定 ±0.5 m，使用固定图框比例；范围外高度不改写，完整保存在 MAT、CSV 和 FIG 数组。详细文件索引与显示说明见 [README.md](README.md)。

旧的独立四相位部署、重复运动和极值对比绘图已从批量入口移除；`RIS_wave_states.m` 与未调用的调色函数已删除。纯海面入口保留 `main3_windseasurface_0326.m`，可供原论文工况重新生成修正后的海面。

本次验证包括连续谱积分、默认及加密频率网格方差、海面与运动输入的一致性、12 工况有限数值与姿态正交性、16 个图框的轴范围及图件视觉核对。水动力模型本身的实验适用性不由这些检查证明。

## 8. 模型局限与结果解释

1. **水平及偏航水动力不完整**：只有有限吃水侧壁 Froude–Krylov 压力，没有相应绕射、辐射附加质量和黏性阻力。缺失项以未求解标志明确输出；没有填入经验系数冒充已验证模型。
2. **参考线性水动力**：波浪激励及辐射关于初始位置／水平航向线性化；未推导瞬时随体六模态附加质量的输运及非线性流体Coriolis力矩，只保留参考系线性加速度和记忆项。完整的刚体陀螺项不代表这些流体非线性项已经求解。明显水平位移、大航向改变需要移动湿表面六模态 BEM／CFD 或试验系数。
3. **竖向薄板假设**：复用文献薄板模型，具有有限厚度的侧壁只作入射压力补充。不是有限吃水长方体的全耦合六模态势流解。
4. **辐射拟合误差**：正实有理近似有非零误差，需查看阻抗误差和响应误差；时间积分的四阶收敛不等于水动力模型收敛。
5. **静水／湿表面**：未求解越浪、拍击、出水、淹没、非线性水线变化。`maxIncidentWaveBoardMismatchM > draft` 是需要关注的失配指标，尤其薄板自由干舷很小。它使用入射波，没有包含局部散射波面，不能直接等同于真实越浪水深。
6. **高海况范围**：U=3–25 m/s 的扫描是数值稳定性检查，不能把所有这些海况认定为物理适用。小姿态范围、大航向改变、波板失配等应同时报告。
7. **外力与高阶机制**：没有直接风压、海流、系泊、推进、完整二阶漂移／慢漂、结构弹性及主动控制。长期无约束漂移和无偏航阻尼的偏航率不会自动回到零。
8. **无实测验证**：文献验证的是其自身模型和试验工况；本项目组合模型的数值测试、积分一致性、对称性及回归通过，不构成该 RIS 的试验验证。

## 9. 一手依据

- [Fossen: Marine Craft Model](https://www.fossen.biz/html/marineCraftModel.html)：六自由度刚体运动、惯性与坐标耦合、ZYX 与四元数。该网页使用 NED；本项目采用 Z 向上，本文已自行重列符号。
- [Michele 等：Wave energy extraction from rigid rectangular compound floating plates](https://doi.org/10.1016/j.jfluidstructs.2024.104193)：原三个竖向薄板辐射／绕射算法来源；[作者机构记录](https://art.torvergata.it/handle/2108/394478)。
- [WEC-Sim 理论文档](https://wec-sim.github.io/WEC-Sim/main/theory/theory.html)：Cummins 辐射卷积与状态空间近似的定义。项目的非负阻抗拟合由本次自行实现，不声称复刻 WEC-Sim 拟合器。
- 本文侧壁 Mz 解析式由入射压力表面积分直接推导，并用独立四壁高阶求积校验。

