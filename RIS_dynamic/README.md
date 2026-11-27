# RIS_dynamic：动态海况下浮动 RIS 与粗糙海面直达链路的 MC 对比（v3）

## 0. 概述

`RIS_phase` 中"无 RIS 直接传播"（No RIS direct）用的是**静止海面**：Fig1 的直接通信列用裸水—气平面折射 MC 包（`rp_mc_water_air`），湍流/入射角扫描用连续弹道模型（`rp_direct_wave`）。RIS 板本身也是静止的水平板（`wris_reference`，$R_{\mathrm{RIS}}=I_3$）。

本目录在**同一海况**下重新比较两条链路：

| 链路 | RIS_phase（原） | RIS_dynamic（本目录） |
|---|---|---|
| 无 RIS 直接传播 | 静止平面 $z=0$，弹道或已存 MC 包 | **逐帧 MC**：水下 HG 多次散射 → 冻结粗糙海面 $\eta(x,y,t_f)$（局部法向、矢量 Snell、Fresnel、全反射）→ 波峰再遮挡 → 空气 → 接收机 |
| RIS 链路 | 静止水平板 | 板随 `sea_wave/RIS_float` 六自由度运动，**LC 相位冻结**为静水设计值，逐帧 MC |
| 海况 | 无 | 每个风速一组方向波谱分量，**同时**驱动板的运动和直达链路穿过的海面 |
| 接收机 | 各方案对准自身光斑质心 | 同样对准静水质心，但**一次对准、不随波浪跟踪** |
| FOV | 主光线判据（$\psi\le\Theta$） | **逐光线**判据（粗糙海面使到达方向发散，主光线判据不再成立） |
| 统计 | 抖动 + 角误差 + EGG×GG | 上述全部 + **时间平均**（海况作为又一个随机变量） |

除此之外，光源、水体、玻璃/LC 阵列、器件角度量化、相位族（Steering-only / Flat top / Robust）、接收机抖动/角误差/孔径、EGG×GG 湍流、OOK 热噪声和 SNR 门限全部沿用 `RIS_phase/rp_stage2_config`，并直接调用 `RIS_phase` 的函数（`rp_stage2_scene`、`rp_phase`、`rp_finite_phase`、`rp_trace`、`rp_link_metrics` 等）。

**没有修改 `RIS_phase`、`W_RIS`、`sea_wave`、`yimiao` 下的任何文件**：本目录只调用它们，所有缓存、checkpoint 和结果都写在 `RIS_dynamic/results/`。

**交付状态：代码写完后没有运行过**（编写环境没有 MATLAB）。已做静态检查：代码块与括号配平、所有被调用的项目函数都存在、调用签名与 `RIS_phase`/`W_RIS`/`sea_wave` 的源码一致。请按第 4 节先跑测试和 smoke，再跑 production。上一版（v2，基于 W_RIS 的功率—时间对比）已整体替换。

---

## 1. 比较的七条链路

结果数组中链路维 `k` 的顺序固定为 `rd_schemes()`：

| k | key | 图例 | 说明 |
|---|---|---|---|
| 1 | `direct` | No RIS direct | 无 RIS，直接穿过（静止或粗糙）海面 |
| 2 | `steering` | Steering-only | `rp_phase`：每个单元尽量指向世界 $+z$ |
| 3 | `flat_top` | Flat top | `rp_finite_phase('flat_top')`，目标宽度最接近 10 cm |
| 4 | `robust_aber` | Robust ABER | robust 族中静水 ABER 最小 |
| 5 | `robust_op` | Robust OP | robust 族中静水 OP 最小 |
| 6 | `max_power` | Max mean power | 可行设计中静水平均接收功率最大 |
| 7 | `sea_aware` | Sea-aware ABER | **新增**：在独立海面实现上按时间平均 ABER 选宽度（第 2.6 节） |

静水选择规则与 `rp_stage2_distance` + `rp_stage2_select` 相同：每个（族, 宽度）在方向余量内取 ABER 最小的中心余量；可行性要求短轴半径和等效半径 $\ge 4$ cm、轴比 $\le 1.35$、$R_{50}\ge 2$ cm、倾角在最小可行倾角 $+1.5^\circ$ 以内。

---

## 2. 模型与公式

### 2.1 同一海况驱动海面与 RIS

每个风速 $U$ 用 `RIS_wave_components` 生成 Neumann 谱、$\cos^2$ 方向分布的分量 $\{a_j,\omega_j,\mathbf k_j\}$（随机相位由 `c.ris.seed` 固定）：

$$
\eta(\mathbf x,t)=\operatorname{Re}\sum_j a_j\,e^{i(\omega_j (t_0+t)-\mathbf k_j\cdot\mathbf x)},\qquad
H_s=4\sqrt{\tfrac12\sum_j |a_j|^2}.
$$

同一组分量传给 `RIS_float` 求板的六自由度运动（位姿 $R(t)$、下表面中心 $\mathbf p_b(t)$），也传给 `rd_surface` 生成直达链路在第 $f$ 帧穿过的冻结海面 $\eta(\mathbf x,t_f)$。海面在光束穿越点附近 $\pm 3$ m 内用 12.5 mm 网格三次插值（最短波长约 0.17 m），网格外用精确 Fourier 求和。

### 2.2 直达链路 Monte Carlo（`rd_direct_mc`）

光源与水下输运和 `W_RIS/wris_mc` 完全一致：束腰位置 $\mathbf u_0\sim\mathcal N(0,w_0^2/4)$、斜率 $\boldsymbol\tau\sim\mathcal N(0,\theta_d^2/4)$，自由程 $\ell=-\ln U/c_w$，碰撞时吸收 $w\,a_w/c_w$、剩余 $w\,b_w/c_w$，HG 散射（$g=0.924$ 时直接调用 `yimiao/Hen_Green`），俄罗斯轮盘，源平面背向逃逸，路径/散射次数上限。

光子在其直线段上**首次**向上穿过 $z=\eta(x,y,t_f)$ 的位置（上界步进 + 二分，容差 2 μm）发生折射。局部法向和切向波矢为

$$
\hat{\mathbf n}=\frac{(-\partial_x\eta,-\partial_y\eta,1)}{\sqrt{1+|\nabla\eta|^2}},\qquad
\mathbf q=n_w\big[\hat{\mathbf d}-(\hat{\mathbf d}\cdot\hat{\mathbf n})\hat{\mathbf n}\big],\qquad q=|\mathbf q|.
$$

$q\ge n_a$ 时全反射；否则空气中方向和透射率为

$$
\hat{\mathbf d}_a=\frac{\mathbf q+\sqrt{n_a^2-q^2}\,\hat{\mathbf n}}{n_a},\qquad
T=1-\tfrac12\left(|r_s|^2+|r_p|^2\right)\quad(\texttt{rp\_fresnel}).
$$

反射/全反射功率终止并记账（没有再追踪回水中的多次界面作用）。透射光线若在离开波浪层 $|z|\le\max|\eta|$ 之前再次进入水中（波峰遮挡），也记为损失（`c.surface.crestBlocking`）。能量账本

$$
P_{\mathrm{trans}}+P_{\mathrm{refl}}+P_{\mathrm{TIR}}+P_{\mathrm{abs}}+P_{\mathrm{back}}+P_{\mathrm{trunc}}=N+\Delta_{\mathrm{roulette}}
$$

逐帧检查，残差 $>10^{-9}$ 会给出警告。

静水基线用同一函数、平面 $z=0$（`rd_flat_surface`），因此"静止 → 粗糙"的变化只来自海面，不来自算法差异。

### 2.3 浮动 RIS 链路（`rd_ris_packets` + `rd_trace`）

1. 第 $f$ 帧的板位姿由 `rp_geometry` 从 `RIS_float` 读取；TX **固定不动**（两条链路共用）。
2. `wris_mc` 把光子输运到**运动中**的 1 m × 1 m 下表面（同一种子同时用于直达链路，公共随机数）。
3. 包转到板坐标系，水→玻璃 Snell + Fresnel，到 LC 屏高度 $h_i$，只保留落在 14.4 mm 透光区内的包（与 `rp_illuminate` 同一套取整和掩模）。
4. 每个单元施加**静水设计时的**量化仿射相位（`rd_slim_phase` 只保留 `rp_trace` 需要的字段），由 `rp_trace` 完成玻璃→空气折射并直线传播到世界水平面 $Z_R$。

包功率始终是原始 MC 权重 $P_0w/N$，不做确定性重加权，因此 RIS 和直达两条链路使用同一种无偏估计。

单元口径 14.4 mm、传播 5–10 m 时 $\lambda z/D\approx0.3$ mm，远小于 5 cm 接收孔径，几何光线足以描述动态下的光斑漂移和截获。`RIS_phase` 波动模型给出的像素内干涉细节在这里不重现。

### 2.4 固定接收机与逐光线 FOV（`rd_receiver` + `rd_capture`）

每条链路的接收机在静水下对准自身弹道光斑质心 $\mathbf C$（`RIS_phase` 约定），之后**不随波浪跟踪**。名义视轴 $\hat{\mathbf r}$ 从板上表面中心（直达链路用 TX 光轴与 $z=0$ 的交点）指向 $\mathbf C$。抖动与角误差样本与 `rp_aperture_channel` 完全相同：

$$
\boldsymbol\delta_s=\sigma_p\mathbf z^{\mathrm{pos}}_s,\qquad
\mathbf e_s=\boldsymbol\mu_\theta+\sigma_\theta\mathbf z^{\mathrm{ang}}_s,\qquad
\psi_s=|\mathbf e_s|,\qquad
\hat{\mathbf n}_s=\cos\psi_s\,\hat{\mathbf r}+\sin\psi_s\,\frac{e_{s1}\hat{\mathbf e}_1+e_{s2}\hat{\mathbf e}_2}{\psi_s}.
$$

第 $j$ 个包在接收面上的落点为 $\mathbf X_j$、方向为 $\hat{\mathbf d}_j$、功率为 $w_j$。它与样本 $s$ 的圆盘平面（过 $\mathbf P_s=(\mathbf C+\boldsymbol\delta_s,Z_R)$、法向 $\hat{\mathbf n}_s$）交于 $\mathbf Y=\mathbf X_j+t\hat{\mathbf d}_j$，其中

$$
t=\frac{\hat{\mathbf n}_s\cdot(\mathbf P_s-\mathbf X_j)}{\hat{\mathbf n}_s\cdot\hat{\mathbf d}_j},\qquad
G_s=\frac1{P_0}\sum_j w_j\,\mathbf 1\{|\mathbf Y-\mathbf P_s|\le a_R\}\,\mathbf 1\{\hat{\mathbf d}_j\cdot\hat{\mathbf n}_s\ge\cos\Theta\}.
$$

这是 `RIS_phase` 倾斜圆盘积分的光子包形式，区别在于 FOV 对**每条光线**判断。粗糙海面使直达光的到达方向显著发散，主光线判据会高估直达链路。为保证公平，RIS 链路也用同一判据，所以静水数值与 `RIS_phase` 会略有不同（见第 8 节第 1 条）。

为提速，先把落点相差 < 2 mm、方向余弦相差 < 2×10⁻³ 的包合并（功率和质心精确守恒），再按样本块做包围盒剪枝。

### 2.5 通信指标：时间作为又一个随机变量（`rd_metrics`）

每帧的增益样本 $\{G_{s,f}\}$ 交给 `rp_link_metrics`：均值为 1 的 EGG（海洋）× Gamma-Gamma（大气）衰落 $H$，热噪声 OOK，$\mathrm{SNR}=(\mathcal R P_tG H/\sigma_n)^2$。再对帧等权平均：

$$
\mathrm{ABER}(P_t)=\frac1{N_f}\sum_{f=1}^{N_f}\frac1{N_s}\sum_{s=1}^{N_s}\mathbb E_H\!\left[\tfrac12\operatorname{erfc}\!\frac{\mathcal R P_t G_{s,f}H}{2\sqrt2\,\sigma_n}\right],\qquad
\mathrm{OP}(P_t)=\frac1{N_f}\sum_f\frac1{N_s}\sum_s\Pr\!\left(H<\frac{P_{\mathrm{th}}}{P_tG_{s,f}}\right).
$$

即对时间、接收机抖动/角误差和湍流的联合平均。$G=0$ 的样本（光被全反射、光斑漂出孔径或落在 FOV 外）给出 BER = 1/2、中断。对 RIS 链路，帧间波动来自板的倾斜/升沉造成的光斑漂移；对直达链路，帧间波动来自海面斜率造成的折射偏转和全反射。

### 2.6 海况感知宽度选择（新增，`rd_train`）

静水最优宽度在波浪下不一定最优：光斑太窄，板一倾斜就漂出孔径；太宽则功率密度不足。`rd_train` 把静水候选库（每个族、每个宽度一个设计）放到**另一组独立的海面实现**上（波相位种子 `c.ris.seed + 1000`），按时间平均 ABER（或 OP）重选，结果作为第 7 条链路，再在评估海面上和其他方案一起比较。训练与评估用不同的海面实现，可以避免"在测试集上调参"带来的乐观偏差。

### 2.7 入射角扫描

与 `RIS_phase` Fig10 规则相同：每个入射角重算场景和载波，沿用 40° 下选定的族/宽度/余量。超过静水临界角 $\arcsin(n_a/n_w)=48.75^\circ$ 时，静水直达链路没有弹道光，接收机改用粗糙海面的长期平均光斑质心对准（`rd_pilot`）。这时可以直接看到粗糙海面如何把被全反射截断的直达链路部分"放出来"。

---

## 3. 目录结构

```text
RIS_dynamic/
├─ main_ris_dynamic.m        一键入口（参数区 → 计算 → 画图）
├─ README.md
├─ algorithm/                计算，不画图
│  ├─ rd_run.m               主流程：设计 → 训练 → 各海况 → 入射角 → 统计，checkpoint
│  ├─ rd_config.m            全部参数（production / smoke）
│  ├─ rd_setup.m / rd_check.m  路径、版本、工具箱与参数检查
│  ├─ rd_design.m            静水设计与选择（调用 RIS_phase）
│  ├─ rd_train.m             海况感知宽度选择
│  ├─ rd_condition.m         单个海况的逐帧配对 MC
│  ├─ rd_motion.m            波谱分量 + RIS_float（缓存）
│  ├─ rd_surface.m / rd_flat_surface.m / rd_crossing.m   粗糙/静止海面与首次穿越
│  ├─ rd_direct_mc.m         直达链路 MC
│  ├─ rd_ris_packets.m / rd_trace.m / rd_slim_phase.m    RIS 链路
│  ├─ rd_to_plane.m / rd_compress.m / rd_spot.m / rd_map.m
│  ├─ rd_receiver.m / rd_capture.m / rd_metrics.m        接收机与通信指标
│  ├─ rd_links.m / rd_pilot.m / rd_evaluate.m
│  └─ rd_cached.m / rd_code_hash.m / rd_physical.m / ...  checkpoint 工具
├─ plots/                    画图，只读 results/ris_dynamic.mat
│  ├─ plot_ris_dynamic.m     画图入口
│  ├─ rd_plot_options.m      所有画图选项（可覆盖）
│  └─ rd_fig_*.m             各图与公共样式
├─ tests/RdPhysicsTest.m     9 项快速物理/记账测试
└─ results/                  运行后生成（*.mat / *.csv / *.fig 被 .gitignore 排除，PNG 会入库）
   ├─ ris_dynamic.mat        轻量结果结构 R（画图只读它）
   ├─ ris_dynamic_*.csv      汇总表
   ├─ run_log.txt
   ├─ figures/               fig1–fig8
   ├─ checkpoints/           各阶段原始结果（含逐帧增益样本）
   └─ cache/                 RIS_phase 设计场景 MC 与 RIS_float 运动缓存
```

---

## 4. 运行步骤

环境：MATLAB **R2022a 及以上**，**Statistics and Machine Learning Toolbox**（`sobolset`、`ncx2cdf`，与 `RIS_phase` 相同）。`RIS_dynamic` 必须和 `RIS_phase`、`W_RIS`、`sea_wave`、`yimiao` 放在同一个 `code` 目录下。

在 MATLAB 中把当前目录切到 `...\code\RIS_dynamic`。

**步骤 0：测试（约 1 分钟）**

```matlab
addpath('algorithm','plots'); rd_setup();
r = runtests('tests/RdPhysicsTest.m'); table(r)
```

9 项应全部通过：平行光截获 $=I_h\pi a^2/d_z$、逐光线 FOV、抖动移位、平面 Snell + Fresnel、超临界全反射、`rd_surface` 与 `RIS_wave_field` 一致、粗糙海面能量账本闭合、包合并守恒、时间平均指标。

**步骤 1：smoke（几分钟）**

`main_ris_dynamic.m` 中设 `profile = 'smoke'`，按 F5。smoke 只有 1 个高度（7.5 m）、1 个风速（3 m/s）、3 帧、少量光子，入射角扫描只有 40° 和 50°（50° 用于走一遍"超临界 → pilot 对准"的分支）。图面稀疏是正常的，这一步只检查流程。

**步骤 2：production**

设 `profile = 'production'`，按 F5。中途中断后再运行即可，已完成的阶段会从 `results/checkpoints/` 复用（日志显示 `reusing checkpoint`）。修改物理参数或 `algorithm/` 中的物理代码后，相关 checkpoint 自动作废；只改 `rd_evaluate`、画图或 `c.metrics` 不会触发重算。

**步骤 3：只改图**

`doCompute = false`，修改 `plotOptions` 后运行；或直接：

```matlab
plot_ris_dynamic('results/ris_dynamic.mat',struct('curveWinds',[0 4 6],'formats',{{'png','pdf'}}));
```

**需要交回的文件**：`results/figures/*.png`、`results/ris_dynamic_*.csv`、`results/run_log.txt`（MAT 太大可以不传）。

---

## 5. 主要参数（`algorithm/rd_config.m`）

| 参数 | 默认值 | 含义 |
|---|---|---|
| `c.scene.distancesM` | `[7.5 10]` | 接收面高度（板上表面以上，m） |
| `c.scene.incidenceDeg` | 40 | 水中入射角，方位角沿用 RIS_phase |
| `c.sea.windsMps` | `[3 4 5 6]` | 风速，`RIS_float` 只接受 3–25 m/s；静水 $U=0$ 自动加在最前 |
| `c.sea.windDirectionDeg` | 50 | 风向（与 W_RIS/RIS_phase 一致） |
| `c.ris.duration / dt` | 30 s / 0.05 s | 每个风速的运动时长与积分步长 |
| `c.dynamic.frameStepS` | 0.5 s | 光学帧间隔，61 帧 |
| `c.dynamic.photons` | 5×10⁵ | 每条链路每帧的光子包数 |
| `c.dynamic.samples` | 2048 | 每帧接收机抖动/角误差样本（Sobol） |
| `c.dynamic.calmFrames` | 3 | 静水独立 MC 重复次数 |
| `c.study.fovHalfAngleRad` | 0.130 | FOV 半角（RIS_phase Fig2/3 取值） |
| `c.study.positionSigmaM / angleSigmaRad / angleBiasRad` | 0.08 / 0.025 / [0.025 0] | 继承自 RIS_phase |
| `c.study.turbulence / atmosphere` | weak / weak | 设计和训练使用的湍流 |
| `c.metrics.turbulence` | {weak,weak; strong,strong} | 后处理中计算的湍流组合（第一行为主图） |
| `c.surface.gridStep / gridHalfWidth` | 12.5 mm / 3 m | 海面插值网格 |
| `c.surface.crestBlocking` | true | 直达光被波峰再遮挡 |
| `c.training.*` | ABER，独立海面，1 s 帧，2×10⁵ 光子 | 海况感知选择 |
| `c.incidence.anglesDeg / windMps` | [30 … 54] / 5 | 入射角扫描 |
| `c.metrics.fixedPowerDBm` | 30 | 汇总表和 Fig6/7 的发射功率 |

---

## 6. 输出图

| 文件 | 内容 |
|---|---|
| `fig1_irradiance_maps` | z = 7.5 m 接收面光强：静水 / U = 5 m/s 某一帧 / U = 5 m/s 时间平均；列为直达、Steering-only、Flat top、Robust ABER、Sea-aware。实线圆为孔径，虚线圆为 RMS 抖动半径 |
| `fig2_aber_vs_power` | ABER–$P_t$，行为高度，列为静水 / 3 m/s / 5 m/s；点线为指向 FOV 地板 |
| `fig3_op_vs_power` | OP–$P_t$，同上 |
| `fig4_time_series` | U = 5 m/s：各链路逐帧增益、逐帧 OP、RIS 倾角 / 光束处海面坡度 / 直达弹道光全反射比例 |
| `fig5_power_cdf` | 瞬时接收功率 CDF（实线粗糙，虚线静水），显示深衰落尾部 |
| `fig6_wind_speed` | 固定 $P_t$ 下 ABER、OP 随风速变化 |
| `fig7_incidence_angle` | 固定 $P_t$ 下 ABER、OP 随入射角变化（实线 U = 5 m/s，虚线静水，点线为临界角） |
| `fig8_sea_aware_width` | 训练海面上时间平均 ABER/OP 随静水等效光斑半径的变化，实心点为选中的设计 |

样式按 RIS_phase 第二阶段：Times New Roman、10 pt、每个子图 2.67 × 2.37 in，图例放在坐标区上方，不遮挡曲线。字号、尺寸、格式（png/pdf/fig）、展示的风速/高度都在 `plots/rd_plot_options.m` 中修改，不需要重算。

---

## 7. 计算量（估计，未实测）

production 的耗时主要来自：4 个风速 × 61 帧 × 2 条 MC（各 5×10⁵ 包，直达链路还要做海面穿越搜索），以及每帧 14 次接收机截获。按经验估计为数小时量级，入射角扫描再加约 1–2 小时。建议先看 smoke 的 `run_log.txt` 中每帧耗时再外推。时间不够时可以优先减少 `c.sea.windsMps`、增大 `c.dynamic.frameStepS`，或者关掉 `c.incidence.enabled`。

---

## 8. 已知限制与建议（论文相关）

1. **FOV 判据变了。** RIS_phase 用主光线判据（$\psi_s\le\Theta$ 即全收），这里改为逐光线判据。对静止的 RIS 光斑，二者差别来自 1 m 口径造成的到达角展宽（约 ±50 mrad）；对粗糙海面直达光，差别是本质性的。建议论文统一采用逐光线判据，必要时 RIS_phase 静态图也按此重算一遍，避免静/动两部分口径不一致。
2. **公平性。** 两条链路共用 TX、光束、水体、每帧 MC 种子、接收机模型和湍流；接收机都对准各自静水质心，且都不跟踪。两点不对称：直达链路计入波峰遮挡，RIS 链路没有计入浮框侧壁遮挡、上浪和板附近波峰遮挡（RIS 出射光离竖直方向约 30°，直达光约 59°，前者掠过波峰的机会小得多）。另外未射到 1 m 板上的光直接丢弃，没有当作旁路直达光。论文中应说明这些假设。
3. **TX 固定。** W_RIS 中的 TX 航位推算运动没有引入，这样两条链路的差异只来自海面和板。如需加入，两条链路应共用同一 TX 轨迹。
4. **统计置信度。** 每个风速只有一组 30 s 的波浪实现。正式出图前建议增加时长或多组 `c.ris.seed`，并用批均值给出置信区间（结果中已保存逐帧指标 `R.main.frames{iw}.aberFrames/opFrames`，可以直接计算）。
5. **EGG/GG 与显式波浪不要重复计入。** EGG 表示气泡/温盐引起的水下闪烁，GG 用的是路径等效形状参数；波浪引起的衰落已经显式建模，正文中应说明 GG 参数不包含海面波动的贡献。
6. **RIS_float 的适用范围。** 模型是小角度近似六自由度模型，`R.main.maxTiltDeg` 和 `R.main.smallAngleValid` 记录了每个风速的最大倾角和是否小于 10°。风速较高时请先检查这两个量，再解读结果。
7. **海况感知选择可以作为一个贡献点。** Fig8 直接回答"波浪越大，最优光斑应该多宽"。由于训练和评估用的是不同的海面实现，Sea-aware 曲线不是在测试数据上调出来的。
8. **相位冻结**是保守设定：LC 不随波浪重新指向。如果后续要讨论闭环指向（例如用板的 IMU 姿态逐帧更新命令），可以在 `rd_condition` 中按 `g.RRis` 重算命令；当前接口已经支持。
