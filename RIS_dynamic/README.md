# RIS_dynamic：动态海况下浮动 RIS 与粗糙海面直达链路的 MC 对比（交付说明）

## 0. 概述

`W_RIS` / `RIS_phase` 中的“无 RIS 直接传播”基线用的是**静止海面**（平面 $z=0$）。本目录把它换成**同一海况下的粗糙海面直达链路**：

- 每个风况只生成一组方向波谱分量，同时驱动 RIS 六自由度运动（`sea_wave/RIS_float`）和直达链路穿过的水–气界面。
- 每一帧中，两条链路共用同一个 TX 位姿、光束、水体、MC 种子和物理接收机（RX）。
- 直达链路逐帧用 Monte Carlo 计算：水下 HG 多次散射，粗糙界面上做 Snell 折射、Fresnel 透射和全反射，再进入空气到达 RX。
- 旧的静止海面直达链路作为第三条曲线一起算出，方便说明“静止→粗糙”带来的变化。

**交付状态：代码写完后没有运行过。** 结果图和 CSV 需要运行者按第 4 节生成。`results/v1_previous/` 中的三张 PNG 来自本目录上一版代码，只供对照（见 §8 第 4 条）。

**没有修改的部分**：`W_RIS`、`RIS_phase`、`sea_wave`、`yimiao` 下的文件一律没动，本目录只调用它们。

---

## 1. 相对旧基线和上一版（v1）的改动

| 项目 | 旧基线（W_RIS/RIS_phase） | 本目录 v2 |
|---|---|---|
| 直达链路的海面 | 静止平面 $z=0$ | 同一海况、同一时刻的冻结粗糙海面（与 RIS 运动用同一组波分量） |
| 直达链路的方法 | 静止界面 | 逐帧 MC：HG 散射，按光线顺序找首次穿越，矢量 Snell、Fresnel、TIR，透射后被波峰再遮挡 |
| 静止海面基线 | — | 保留为 `flat` 链路，同 TX 运动、同种子逐帧计算（`c.comparison.flatBaseline`） |
| RIS 链路 | — | 沿用 W_RIS 的 MC 和 RIS_phase 的 LC 调相光线追迹，算法不变 |

v2 相对 v1（本目录的上一版）的工程改动：

1. 拆成 `algorithm/`（计算）、`plots/`（画图）和 `tests/`（回归测试），画图代码只读 MAT 文件。
2. 新增一键入口 `main_ris_dynamic.m`，参数集中在一个区块里。
3. 新增静止海面基线 `flat`。
4. 新增时间统计：均值、中位数、最值、闪烁指数、中断概率、零命中比例、RIS 相对两条直达链路的增益（dB），输出 `wris_dynamic_cases.csv`。
5. 新增两张图：功率 CDF（fig2）和随风速的统计图（fig5）。原来的三张图都加了 `flat` 曲线。
6. RIS 链路新增两项诊断：板上表面出射点低于入射海面的功率（`risOverwashedExitW`），以及可选的波峰遮挡（`c.comparison.risWaveBlocking`，默认关闭）。
7. 按风况保存 checkpoint，中断后重跑会接着算。config 或 `algorithm/*.m` 有任何改动，旧 checkpoint 自动作废。
8. 计算前先检查 MATLAB 版本和依赖函数路径。所有日志写入 `run_log.txt`。
9. 探测面光强图改存 `single`，MAT 体积减半。功率等标量仍为 `double`。

在相同参数下，`direct` 和 `ris` 两条链路的物理算法、随机种子与 v1 完全一致（见 §8 第 4 条）。

---

## 2. 目录结构

```text
RIS_dynamic/
├─ main_ris_dynamic.m          一键入口（参数区 → 计算 → 画图）
├─ README.md                   本文档
├─ algorithm/                  计算，不画图
│  ├─ run_wris_dynamic.m       主计算：逐风况/逐帧/逐角度，保存 MAT/CSV/checkpoint
│  ├─ wrd_config.m             全部默认参数（production / smoke）
│  ├─ wrd_setup.m              添加路径（W_RIS、RIS_phase、sea_wave、yimiao/fig3-7）
│  ├─ wrd_check_dependencies.m 版本与依赖检查
│  ├─ wrd_surface.m            当前帧冻结的 Fourier 海面：网格三次插值 + 网格外精确求和
│  ├─ wrd_flat_surface.m       静止海面 z=0（旧基线），与 wrd_surface 同一接口
│  ├─ wrd_crossing.m           光线与海面的首次穿越（上界步进 + 二分）
│  ├─ wrd_direct_mc.m          直达链路 MC（水下 HG + 海面折射/TIR + 空气段）
│  ├─ wrd_ris_mc.m             RIS 链路（W_RIS 的 wris_mc + RIS_phase 的 rp_trace）
│  ├─ wrd_receive.m            与公共 RX 平面求交，可选波峰遮挡
│  ├─ wrd_score.m              孔径 + FOV 计分、MC 标准误、探测面光强图
│  └─ wrd_summary.m            帧表与时间统计表
├─ plots/                      画图，只读 results/*.mat
│  ├─ plot_wris_dynamic.m      画图入口
│  ├─ wrd_fig_power.m          fig1 功率-时间
│  ├─ wrd_fig_cdf.m            fig2 功率 CDF
│  ├─ wrd_fig_tir.m            fig3 全反射占比与海面坡度
│  ├─ wrd_fig_maps.m           fig4 探测面光强
│  ├─ wrd_fig_summary.m        fig5 统计量-风速
│  └─ wrd_fig_*.m              公共样式、保存、图例、坐标范围
├─ tests/WrdPhysicsTest.m      9 项物理/记账回归测试（不调用 RIS_float，耗时很短）
└─ results/                    输出目录（运行后生成；*.mat/*.csv/*.fig 被 .gitignore 排除）
   └─ v1_previous/             上一版的三张 PNG，仅供对照
```

`RIS_dynamic` 必须与 `W_RIS`、`RIS_phase`、`sea_wave`、`yimiao`（旧副本中叫 `3yimiao` 也可以）放在同一个 `code` 目录下。

---

## 3. 运行环境

- MATLAB **R2022a 或更新版本**（原工程使用 R2024b），代码用到 `clim`、`tiledlayout`、`exportgraphics`、`arguments`。
- 只需要基础 MATLAB，**不需要任何工具箱**。
- 磁盘：production 结果 MAT 估计在百 MB 量级。另有每个风况一个 checkpoint。
- 运行时间：**没有实测过**。先跑 smoke，看 `run_log.txt` 里的耗时，再决定是否在空闲机器上跑 production。耗时主要来自 `RIS_float`、每帧的海面网格和 3 条链路 × 每帧 $10^5$ 个光子包。光子包数由 `c.mc.photons` 设置。

---

## 4. 运行步骤（交给运行者）

在 MATLAB 里把当前目录切到 `...\code\RIS_dynamic`。

### 步骤 0：回归测试（预计 1 分钟内）

```matlab
addpath('algorithm'); wrd_setup();
r = runtests('tests/WrdPhysicsTest.m'); table(r)
```

9 项应全部 `Passed`：Fourier 海面插值、多次穿越取第一次、静止平面、平面 Snell + Fresnel、超临界角全反射、斜面让超临界光线透射、孔径/FOV 计分、静止 RIS 链路能达到 RX、时间统计公式。

### 步骤 1：smoke 试跑（检查环境，预计数分钟）

打开 `main_ris_dynamic.m`，设置 `profile = 'smoke'`，按 F5 运行。也可以在命令行执行：

```matlab
addpath('algorithm','plots'); wrd_setup();
run_wris_dynamic('smoke');
plot_wris_dynamic(fullfile('results','smoke','wris_dynamic.mat'));
```

`results/smoke/` 下应生成 `wris_dynamic.mat`、两个 CSV、`run_log.txt` 和 `fig1`–`fig5` 的 PNG/FIG。smoke 只有 2000 个光子包、3 帧，图面稀疏是正常的，这一步只用来检查流程。

### 步骤 2：production 正式计算

在 `main_ris_dynamic.m` 中设置 `profile = 'production'`，按 F5 运行。输出在 `results/`。

如果中途中断，直接再运行一次即可，已完成的风况会从 `results/checkpoints/` 复用，日志里会显示 `reusing checkpoint`。

### 步骤 3：只改图、不重算

在 `main_ris_dynamic.m` 中设置 `doCompute = false`，修改 `plotOptions`（字体、字号、分辨率、光强图时刻），然后运行。也可以直接执行：

```matlab
plot_wris_dynamic('results/wris_dynamic.mat', struct('fontSize',12,'mapTimeS',15));
```

### 需要交回的文件

`results/` 下的 `fig1`–`fig5`（PNG + FIG）、`wris_dynamic_cases.csv`、`wris_dynamic_frames.csv`、`run_log.txt`、`wris_dynamic.mat`。

---

## 5. 主要参数（`algorithm/wrd_config.m`）

| 参数 | 默认值 | 含义 |
|---|---|---|
| `c.comparison.conditions` | `[3 50; 5 50]` | 风况 [风速 m/s, 风向 °]，取自 `wris_config` |
| `c.comparison.rotationXYDeg` | `[0 40; 0 49]` | TX 光轴：+z 先绕世界 x、再绕世界 y 旋转（°）；水中名义入射角 40°、49°（临界角约 48.75°） |
| `c.ris.duration / dt / startTime` | 20 s / 0.05 s / 25 s | RIS 运动积分区间（与 W_RIS 相同） |
| `c.comparison.frameStep` | 0.5 s | 光学帧间隔，必须是 `dt` 的整数倍，0–20 s 共 41 帧 |
| `c.comparison.flatBaseline` | `true` | 是否同时计算静止海面直达基线 |
| `c.comparison.risWaveBlocking` | `false` | RIS→RX 光线是否也检查被波峰遮挡（敏感性分析） |
| `c.mc.photons` | 100000 | 每条链路每帧的光子包数 |
| `c.receiver.center` | `[8.3 0 5]` m | 固定 RX，约为 40° 入射在静水中折射后、在 $z=5$ m 处的落点 |
| `c.receiver.normal` | 指向 RX 的单位向量 | RX 视轴，即沿名义水面原点到 RX 的到达方向 |
| `c.receiver.radius / fovDeg` | 0.05 m / 15° | 孔径半径和 FOV 半角 |
| `c.receiver.mapEdges` | ±6 m，240 格 | 探测面光强图的范围 |
| `c.surface.gridStep / gridHalfWidth` | 0.075 m / 40 m | 海面插值网格，网格外用精确 Fourier 求和 |
| `c.surface.marchStep / rootTolerance` | 0.04 m / 2 μm | 首次穿越的最小步长和二分容差 |
| `c.metrics.outageThresholdW` | 1e-6 W（−30 dBm） | 中断门限，**应按实际接收机灵敏度设置**。只影响统计和画图，不影响 MC |
| `c.run.resume` | `true` | 是否复用 checkpoint |

光源（1 W、532 nm、10 mrad 发散、高斯束腰）、水体（`c0.151`，HG $g=0.924$）、TX 运动（W_RIS 的 OU 航位推算）、RIS 器件（40×40 单元，节距 25 mm，透光 14.4 mm，±15.426°，量化 0.021°）全部继承自 `wris_config` / `rp_config`，保证与 W_RIS / RIS_phase 一致。

---

## 6. 计算模型

### 6.1 共同输入（每个风况一次）

1. `RIS_wave_components` 生成 Neumann 谱 + $\cos^2$ 方向分布的分量 $\{a_j,\omega_j,\mathbf k_j\}$，种子为 `c.ris.seed`。
2. 同一组分量传给 `RIS_float`，得到 RIS 六自由度运动。直达链路的海面也由同一组分量重建：
   $$\eta(\mathbf x,t)=\operatorname{Re}\sum_j a_j\,e^{\,i(\omega_j (t_0+t)-\mathbf k_j\cdot\mathbf x)}$$
   其中 $t_0$ 为 `startTime`。这个海面是**入射海面**，不含浮板的扰动。
3. TX 运动沿用 `wris_tx_motion`（W_RIS 的模型和种子），各角度分别计算。RX 固定不动。
4. 每一帧（时刻 $t_k$），三条链路使用同一个 `rp_geometry` 快照和同一个 MC 种子（公共随机数）。各帧也复用这个种子，以减小时间曲线上的 MC 抖动。

光传播远快于海面变化，所以每帧把海面和 RIS 位姿当作冻结。

### 6.2 粗糙海面直达链路（`wrd_direct_mc`）

- **发射**：高斯 Wigner 采样，位置和斜率相互独立，与 `wris_mc` 一致。每个光子包初始权重为 1，功率为 $P_0/N$。
- **水下输运**：自由程 $s=-\ln U/c$，隐式俘获（碰撞后权重乘以 $b/c$，吸收量记入能量账），HG 相函数（$g=0.924$，直接调用原工程的 `Hen_Green`），俄罗斯轮盘赌，路径上限，源平面背向逃逸。这些都与 W_RIS 相同。
- **首次穿越**（`wrd_crossing`）：在碰撞**之前**，先在本段自由程内搜索 $z-\eta(x,y)=0$ 的第一个根：
  1. 用振幅上界 $|\eta|\le A$ 把搜索范围限制在 $|z|\le A$ 的板层内。
  2. 如果 $d_z>|d_x|S_x+|d_y|S_y$（$S$ 为坡度上界），间隙函数单调，只检查端点。否则按“安全步长” $\max(\Delta_{\min},-\text{gap}/\text{bound})$ 逐步推进，保证不会跳过根，最小步长为 `marchStep`。
  3. 括住根后用二分法求到 `rootTolerance` 精度。
- **界面**：局部法向 $\mathbf n=(-\eta_x,-\eta_y,1)/\|\cdot\|$。切向分量 $\mathbf q=n_w(\mathbf d-(\mathbf d\cdot\mathbf n)\mathbf n)$。
  - 若 $|\mathbf q|\ge n_a$，发生全反射（TIR）。
  - 否则透射方向为 $\mathbf d_a=(\mathbf q+\sqrt{n_a^2-|\mathbf q|^2}\,\mathbf n)/n_a$，透射率为非偏振 Fresnel 系数 $T=1-(R_s+R_p)/2$（`rp_fresnel`）。
  - 首次界面上的反射和 TIR 能量记入能量账后终止。这部分光回到水中，对空中 RX 的贡献可以忽略。
- **空气段**：沿 $\mathbf d_a$ 直线传到 RX 平面。如果途中再次进入海面（被波峰挡住），这部分功率计入 `directWaveBlockedW`，不进入 RX。
- **能量账**：透射 + 反射 + TIR + 吸收 + 背逃 + 截断 − 轮盘赌修正 = $N$。每帧断言相对残差小于 $10^{-9}$。

### 6.3 静止海面基线（`flat`）

调用同一个 `wrd_direct_mc`，只是把海面换成 `wrd_flat_surface`（$\eta\equiv0$）。TX 运动、种子和 RX 都相同，这就是旧的“无 RIS 直接传播”在动态 TX 下的版本。

### 6.4 浮动 RIS 链路（`wrd_ris_mc`）

- 水下部分直接调用 `W_RIS/wris_mc`，在当前帧的动态板坐标系中求下表面落点、方向和权重。
- 水→玻璃 Fresnel 透射后，按到液晶层中面的位置分配到单元，单元透光区外的光被挡掉。
- 每个受照单元的目标方向由“TX→单元中心”和“单元中心→RX 中心”两段多层中心光线确定（`rp_layer_path`），再经过器件梯度限幅和 0.021° 量化（`rp_command`）。这相当于**每帧理想更新调相**。
- `rp_trace` 计算玻璃→空气传播和 Fresnel 透射，然后用与直达链路相同的 RX 判据计分。
- 板上方按“空气”处理。出射点低于入射海面（覆水）的功率单独记为 `risOverwashedExitW`，RIS 的露底/覆水样点比例也随结果保存。

### 6.5 RX 判据（三条链路完全相同）

光线与 RX 平面（法向为 `c.receiver.normal`）求交：

- 交点到 RX 中心的距离 ≤ `radius`，且到达方向满足 $\mathbf d\cdot\mathbf n_{rx}\ge\cos(\text{FOV})$ 的包，计入 RX 功率 `RxW`。
- 只满足孔径、不检查 FOV 的功率另存为 `ApertureNoFovW`。

单帧 MC 标准误为 $\mathrm{SE}=P_0\sqrt{(\sum w^2-(\sum w)^2/N)/(N(N-1))}$。探测面光强图是 FOV 门之前的辐照度。

### 6.6 时间统计（`wrd_summary`）

对每个（风况，角度，链路），用 $K$ 个帧等权统计：

- 均值 $\bar P$、标准差、中位数、最值
- 闪烁指数 $\sigma_I^2=\langle P^2\rangle/\bar P^2-1$
- 中断概率 $P_{out}=\Pr(P<P_{th})$
- 零命中比例（该帧没有任何光子包被接收）
- RIS 增益 $10\log_{10}(\bar P_{RIS}/\bar P_{direct})$，分别相对粗糙海面和静止海面直达链路

结果中**没有**叠加 EGG 湍流闪烁，功率随时间的起伏只来自海面、RIS 位姿和 TX 运动。

---

## 7. 输出文件

| 文件 | 内容 |
|---|---|
| `wris_dynamic.mat` | 变量 `results`：`config`；`cases{角度行, 风况列}`，每个元素含 `frames(k).direct/.flat/.ris/.geometry/...`；`frameTable`；`caseTable`；`codeHash`；`elapsedSeconds` |
| `wris_dynamic_frames.csv` | 每帧一行：海面高度/坡度，RIS 倾角和湿干比例；各链路的 `RxW/RxSEW/RxHits/RxBallisticW/RxScatterW/ApertureNoFovW`；直达链路的 `Tir/SeaIncident/Transmitted/WaveBlocked`；RIS 的 `BoardIncident/Active/Transmitted/Saturated/OverwashedExit`；能量残差 |
| `wris_dynamic_cases.csv` | 每个（风况，角度，链路）一行：`meanW, stdW, medianW, minW, maxW, meanDbm, scintillationIndex, outageProbability, zeroHitFraction, rmsFrameSEW, ballisticShare, risOverRoughDirectDb, risOverFlatDirectDb, significantWaveHeightM` |
| `fig1_rx_power_time` | 各链路 RX 功率随时间变化，误差棒为 95% MC 区间，× 表示该帧零命中，点线为中断门限 |
| `fig2_rx_power_cdf` | RX 功率的经验 CDF，曲线与门限线交点的纵坐标即中断概率 |
| `fig3_sea_tir` | 直达链路到达海面的功率中 TIR 所占比例（粗糙 vs 静止），以及板中心处的海面坡度 |
| `fig4_rx_maps` | 最大风况、$t$ = `mapTimeS` 时，三条链路的探测面光强（共用色标），白圈为 RX 孔径 |
| `fig5_summary_wind` | 时间平均功率（误差棒为最小–最大帧）和中断概率随风速的变化，每个入射角一列 |
| `run_log.txt` | 完整控制台日志，含进度和耗时 |
| `checkpoints/condition_*.mat` | 每个风况的中间结果，可以删除，删除后会重算 |

---

## 8. 验收检查清单

1. 步骤 0 的 9 项测试全部通过。
2. 运行过程中没有报 `wrd:Energy` 错误，即每帧能量账都已闭合。
3. `wris_dynamic_frames.csv` 中 `directEnergyResidual`、`risWaterEnergyResidual` 的绝对值都小于 1e-9。
4. **与 v1 对照**：`direct`、`ris` 两条链路的算法和种子与 v1 相同，所以 production 的 `roughDirect` / `ris` 时间均值应复现下表（v1 在原机器上的结果）。如果差异明显超过 MC 误差，说明运行环境或依赖代码有变化，需要排查。

   | 风速 (m/s) | 名义入射角 | 直达平均 RX 功率 (mW) | RIS 平均 RX 功率 (mW) | 直达平均 TIR 比例 |
   | ---: | ---: | ---: | ---: | ---: |
   | 3 | 40° | 2.334 | 67.042 | 4.1% |
   | 5 | 40° | 1.355 | 66.634 | 4.5% |
   | 3 | 49° | 0.006 | 66.841 | 53.0% |
   | 5 | 49° | 0.007 | 66.717 | 47.3% |

5. 物理上应满足的定性规律：
   - 49° 超过临界角，`flat` 链路应几乎全部 TIR，RX 功率接近 0 或零命中。粗糙海面的坡度会打开部分帧，所以 `roughDirect` 不为 0 但很低。
   - 40° 时 `flat` 是稳定的基准，`roughDirect` 会随坡度明显起伏，风速越大均值越低、闪烁指数越大。
   - RIS 链路逐帧重新对准 RX，均值高且起伏小。

---

## 9. 已知限制

- 几何光学粒子模型给出的是辐照度和功率，没有相干衍射或散斑。光斑亚结构需要另用波动传播计算。
- 直达链路在首次界面上的反射和 TIR 直接终止，不追踪“反射→水下再散射→再透射”这类高阶路径（对空中 RX 的贡献很小）。
- 首次穿越的最小步长为 4 cm，更细的多重穿越可能被跳过；网格海面在 ±40 m 内为三次插值。临界角附近功率很低的尾部，应结合 `RxHits` 和 MC 误差一起解读。
- RIS 链路假设板上方是空气、板下方是水（三层光路）。覆水时刻的结果仍按名义三层光路计算，应结合 `risTopWetSampleFraction` 和 `risOverwashedExitW` 解读。海面模型是**入射海面**，不含浮板扰动。
- 默认不考虑 RIS→RX 段的波峰遮挡（RIS 位于海面上），可以用 `risWaveBlocking=true` 做敏感性分析。直达链路始终考虑波峰遮挡。
- 每帧使用同一个 MC 种子：曲线更平滑，但各帧的 MC 误差相互关联。时间统计主要反映海况变化，`rmsFrameSEW` 给出单帧 MC 误差的量级。
- 默认中断门限 1e-6 W 只是占位值，应按实际接收机灵敏度修改。改后只需重画图，`caseTable` 需要重算统计（重新运行时会复用 checkpoint，很快）。

---

## 10. 常见问题

- **`wrd:Setup` / `wrd:Dependencies`**：目录相对位置不对，或者路径上有同名函数覆盖了原工程函数。按报错提示核对 `code/` 下的文件夹。
- **`RIS_float:Option` Unknown option**：`c.ris` 里多了 `get_default_RIS_params` 不认识的字段，不要往 `c.ris` 里添加自定义字段。
- **改了参数但结果没变**：checkpoint 按“物理参数 + 风况 + 算法代码哈希”匹配，参数变了会自动重算。如果仍有疑问，删除 `results/checkpoints/` 后重跑。
- **只想改图**：设置 `doCompute = false`，或直接调用 `plot_wris_dynamic`，不会触发任何 MC 计算。
- **打开 .fig 看不到窗口**：保存时已设置 `CreateFcn` 强制可见。如果仍不可见，用 `openfig('xxx.fig','visible')` 打开。
