# FOWT_RIS_eval — 浮式平台六自由度实测数据驱动的 RIS 链路评估与改进

本目录是本次工作的全部新增内容，位于工作区 `C:\Users\Spidey\Desktop\code\wyj`
内，不修改工作区外的任何文件。目标是：用波流水槽实测的浮式平台 6-DOF 数据
评估现有 RIS_model / RIS_phase / sea_wave / W_RIS 链路，并据此改进相位算法，
形成一篇可独立投稿的论文级工作。

## 1. 一句话结论

现有"逐单元各自尽量对准世界 +z"的相位算法在 60° 水入射角下会把 **1600 个
单元全部推到指令盒的同一个角点**，阵列因此散焦 3.0°；5 cm 接收孔径只收到
70 mW 中的 4.58 mW。改用"功率加权、可行域内的公共方向"相位合成后，99.94%
的功率被准直，5 cm 孔径收到 14.82 mW（+5.1 dB）。把该相位冻结在平台上仍会
因实测 0.14–1.10° 的横摇/纵摇产生深衰落（最差可用率 0.59）；用实测姿态做
逐帧点目标自适应控制后，全部 9 条实测记录的平均功率 14.90 mW、10% 分位
14.71 mW、质心误差 2.7 mm、可用率 100%，相对现有算法 **+5.2 dB 均值 /
+5.5 dB 10% 分位**。

## 2. 数据

原始数据（未改动）：
`../浮式平台六自由度实验数据/FOWT Platform Testing/`

- 9 条规则波：0.4 Hz 的 10/15/30 mm，0.8 Hz 的 15/30/45 mm，
  1.2 Hz 的 15/30/45 mm（含 181 s 的 10 mm/0.4 Hz 长记录）；
- 3 条衰减试验：heave、pitch、surge；
- Qualisys 6-DOF 183 Hz；6 个 EDL 波高仪 128 Hz。

处理结果（`data/processed/`）：

- `fowt_motion.npz`：统一到 100 Hz 的 t / pos / rpy / eta；
- `conditions.csv`：每条工况的稳态统计、RAO、相干性；
- `decay_summary.csv`：衰减周期与阻尼比；
- `nominal_design.mat`：RIS 几何、单元功率、全波锚点光斑（MATLAB 导出）。

关键实测值：静水纵摇配平 **+1.40°**（所有记录一致）、横摇 +0.45°；heave
固有周期 2.158 s、阻尼 1.08%；pitch 2.88 s、0.32%；0.4 Hz 的 heave RAO
1.6–2.1，0.8 Hz 0.27–0.29，1.2 Hz 0.06–0.07；roll RAO 14–70 deg/m。

## 3. 评估方法（多保真度回放）

```
实测 6-DOF ──> 姿态/平移 ──> 分层 Fermat 路径 + 离散指令 + 玻璃折射
                                    │
                                    ├─> 单元参考落点 / 功率加权质心
                                    │
全波 BSF/波动光学锚点 ──> 光斑图 ──> 按质心平移 ──> 5 cm 孔径功率
```

- **降阶模型**（`tools/fx_optics.py`）：1600 个单元、水/玻璃两层、48 次二分、
  0.021° 指令量化、±15.414° 限幅、`|q_out|<1` 修复；与 MATLAB 的
  `rp_trace`/`rp_carrier_landing` 在标称位姿上逐单元一致（质心差 <0.1 mm）。
- **全波锚点**（`src/fx_phase_design.m` + `tools/probe_improved.m`）：
  用现有 BSF/Fresnel/波动光学链验证新相位；生产相位与 `rp_phase` 的 +z 结果
  逐位一致（指令、梯度、相位活塞误差为 0）。
- **验证**：降阶模型 vs 全波，6 个缓存位姿的光斑中心差 7–51 mm（≤24% 光斑
  半径）；新相位的光斑中心差 2 mm；5 cm 孔径功率 4.57 vs 4.58 mW。

五个回放场景（`tools/replay.py`）：

| 场景 | 含义 |
|---|---|
| production static | 现有 +z 相位，冻结 |
| production point-adaptive | 现有器件逐帧瞄准自己的接收点（受饱和限制） |
| collimated static | 新相位，冻结 |
| collimated adaptive direction | 新相位，公共世界方向逐帧更新 |
| collimated point-adaptive | 新相位，逐单元瞄准固定接收点 |

## 4. 主要结果

### 4.1 生产相位的饱和

- 1600/1600 单元位于 (−15.414°, +15.414°)；
- 所需偏折 \|q_in\|=1.14，器件最多 0.376（对角），+z 不可达；
- 逐单元输出方向 5–95% 展宽 3.03°，加权标准差 0.31°；
- 全波光斑 90% 半径 0.211 m，峰值 2.01 W/m²；
- 5 cm 孔径 4.58 mW（6.5%），10 cm 16.70 mW（23.8%）；
- 天顶接收机（当前 `rx` 配置）偏离 6.2 m。

### 4.2 新相位

- 功率加权公共方向 `q*=[0.420, −0.690]`，53.9°，覆盖 99.94%；
- 加权方向标准差 0.0105°；
- 5 cm 孔径 14.82 mW（+5.1 dB），10 cm 40.70 mW（+3.9 dB）；
- 峰值 5.36 W/m²（曲率补偿版 7.57 W/m²）。

### 4.3 实测回放（50 Hz，5 cm）

| 场景 | 均值 mW | p10 mW | 平均损耗 dB | 最差可用率 | 质心 rms mm |
|---|---:|---:|---:|---:|---:|
| production static | 4.52 | 4.14 | −0.05 | 1.00 | 44.5 |
| production point | 4.70 | 4.48 | +0.12 | 1.00 | 24.9 |
| collimated static | 12.22 | 9.61 | −0.84 | 0.59 | 48.3 |
| collimated adaptive | 14.37 | 13.63 | −0.13 | 1.00 | 19.2 |
| **collimated point** | **14.90** | **14.71** | **+0.02** | **1.00** | **2.7** |

最差工况 30 mm/0.4 Hz：production 3.95/3.28 mW；collimated static 8.36/4.61，
可用率 0.59；point adaptive 14.93/14.58，可用率 1.00。

### 4.4 控制需求

- 一阶预测把 100 ms 延迟误差从 113 mm 降到 60 mm（−2.54 → −0.67 dB）；
- 20–50 ms 指令更新足够（100 ms 时 61 mm/−0.75 dB）；
- 0.05° 姿态噪声 → 14.8 mm/−0.04 dB；
- 二阶预测被实测 ~0.013° 高频抖动放大，不推荐；先用 5 Hz 低通再考虑。

### 4.5 sea_wave 桥接

- 实测波高仪 → 深水波分量 → `RIS_float` 的 `p.wave`；
- 30/45 mm、0.8 Hz 实测波下，1 m RIS 模型 pitch 2.41°/3.77°、最大 3.62°/5.67°，
  而水槽平台实测仅 0.75°/1.10°；
- 用 RIS 模型自身运动回放：30 mm/0.8 static 10.39 mW（avail 0.87）→ point
  14.85 mW（1.00）；45 mm/0.8 static 6.86 mW（avail 0.42）→ point 14.71 mW
  （1.00）。**实测水槽运动是保守的替代；真实 1 m 板的高海况收益更大。**

### 4.6 实测接口（可直接接回现有模块）

`src/fx_motion_from_experiment.m` 把 Qualisys 记录转换成 `RIS_float` 风格的
运动结构（`t`、`rotationMatrix`、`bottomCenter`、`normal`、角度与位移），
已用 `rp_geometry` 验证：45 mm/0.8 Hz 记录（6451 帧、100 Hz）在 k=1000 处
入射角 59.0°，旋转矩阵正交误差 2e-16。

```matlab
m = fx_motion_from_experiment('.../FOWTmotion_11_3_45mmp8Hz.mat');
g = rp_geometry(c,m,tx,rx,k);   % 现有 RIS_phase 入口
```

`src/fx_wave_from_gauge.m` 把波高仪记录转换成 `RIS_float` 的 `p.wave` 分量：

```matlab
w = fx_wave_from_gauge(t,etaGauge8,tStart=15,tEnd=35);
p.wave=w; o=RIS_float(p);
```

45 mm/0.8 Hz 记录给出 42 个分量、Hs=0.1204 m，合成误差 1.19 mm（波 std
30.1 mm 的 4%）。这两个函数就是"实测数据 → sea_wave / RIS_phase / W_RIS"
的接口，不需要改现有模块的默认路径。

## 5. 目录

```
FOWT_RIS_eval/
├─ CHECKPOINT.md              研究检查点日志
├─ README.md                  本文件
├─ src/
│  └─ fx_phase_design.m       可行域公共方向 / 点目标相位合成（MATLAB）
├─ tools/
│  ├─ fowt_data.py            原始 .mat/.txt 解析
│  ├─ build_dataset.py        清洗、统计、RAO、衰减
│  ├─ fx_optics.py            降阶分层光线 + 指令 + 孔径
│  ├─ replay.py               九记录五场景回放
│  ├─ control_study.py        延迟/预测/更新率/噪声
│  ├─ model_replay.py         实测波驱动 RIS 模型的回放
│  ├─ analyze_results.py      孔径权衡与聚合
│  ├─ make_figures_all.py     Fig1–4
│  ├─ make_figure5.py         Fig5
│  ├─ probe_improved.m        新相位全波锚点
│  ├─ probe_curv.m            曲率补偿全波锚点
│  ├─ export_nominal.m        导出 MAT 给 Python
│  ├─ probe_measured_wave.m   实测波 → RIS_float
│  └─ probe_*.m               诊断脚本
├─ data/processed/            清洗后的数据与标称设计
├─ results/                   所有 CSV / NPZ / MAT / RESULTS.md
├─ figures/                   Fig1–5 PNG（论文用）
├─ figure_atlas_matlab/       18 幅 MATLAB 中文科研图 + 图注 + 数据
└─ paper/
   ├─ manuscript.md           论文全文（Markdown + LaTeX 公式）
   └─ figures/                论文用图副本
```

## 6. 复现

在 `FOWT_RIS_eval/tools` 下：

```text
python build_dataset.py
python replay.py
python control_study.py
python model_replay.py
python analyze_results.py
python make_figures_all.py
python make_figure5.py
```

MATLAB 全波锚点（需要 `RIS_phase` 与 `sea_wave` 路径，脚本内已加）：

```matlab
cd FOWT_RIS_eval/tools
probe_improved      % 新相位全波
probe_curv          % 曲率补偿全波
export_nominal      % 导出 Python 用的标称设计
probe_measured_wave % 实测波驱动 RIS_float（约 75 s/工况）
```

## 7. 与四个模块的关系

- **RIS_model**：给出 40×40、25 mm 间距、1 μm 电极、±15.414° 的真实限制；
  评估发现 60° 入射下器件范围不足，提出"入射角-指令余量"设计规则。
- **RIS_phase**：新增 `fx_phase_design.m`，把逐单元 +z 目标改为功率加权公共
  方向/点目标，解决饱和散焦；保留原 `rp_phase` 作为 +z 特例的回归基准。
- **sea_wave**：新增实测波分量输入（`p.wave`）与实测 6-DOF 回放接口；给出
  实测衰减/RAO 与 1 m RIS 模型的周期、阻尼、倾角包络对比。
- **W_RIS**：把孔径功率、质心、可用率、衰落持续时间作为链路指标；用实测
  运动统计替代纯合成海况，给出 5 cm/2.5 cm/10 cm 孔径的增益曲线。

## 8. 重要边界

1. 水槽平台不是 1 m RIS 板：heave 固有周期 2.16 s vs 1.15 s，阻尼 1.08% vs
   28%（辐射阻尼）。因此**不能用实测数据标定 RIS 水动力系数**；本文用实测
   姿态作为低海况替代，用实测波驱动 RIS 模型给出高海况对照。
2. 没有硬件闭环：LC 响应时间、驱动电子、标定误差未实测；控制研究用指定
   延迟/噪声的仿真。
3. 全波锚点只有一个位姿；回放用质心平移近似光斑，已用 7–51 mm 的中心一致
   性和缓慢形状变化验证。
4. 链路指标是孔径功率，不含湍流/闪烁/背景光/BER；这些需要额外的接收机与
   信道模型。
5. RIS 模型无锚泊；水槽平台有锚泊。实测水平漂移作为平移输入，但 RIS 模型
   自身的漂移未建模。

## 9. 结论

现有 +z 相位在 60° 入射下是**不可行**的：器件指令范围不足，全部单元饱和，
阵列散焦，5 cm 孔径只收到 6.5% 功率，天顶接收机偏离 6.2 m。把目标改成
"功率加权可行公共方向"可把 99.94% 功率准直，5 cm 孔径 +5.1 dB；再把实测
姿态引入逐帧点目标控制，可在全部 9 条实测记录上保持 14.9 mW、100% 可用率，
相对现有算法 +5.2 dB。实测波驱动 1 m RIS 模型给出 2.4–3.8° 的更大倾角，
说明水槽回放是保守的，真实浮式 RIS 更需要这套协同设计。

## 10. MATLAB 中文科研图图集（18 幅）

`figure_atlas_matlab/` 按 `README.txt` 的最新要求，用 **MATLAB R2024b** 重绘了全部图，
不再用 Python 出图。连续曲线统一用 PCHIP 保形平滑，不使用柱状图；图例放在坐标区内的
空角（如 OP 图左下、递增曲线图右下、部分 ABER 图左下/右上），不压曲线、不在坐标区
正上方写字；坐标文字与 x/y 符号放大；热力图与对应曲线成组排版。

图目录：

| 文件 | 内容 |
|---|---|
| `fig01_flow` | 方法流程：实测数据 → 光学模型 → 相位设计 → 光斑/孔径 → 性能评估 |
| `fig02_slope_jpdf` | 海面斜率联合 PDF 等值线（U=3/9/12 m/s，全谱/足迹滤波） |
| `fig03_slope_pdf` | 顺/侧风向斜率 PDF 与 rms–风速曲线 |
| `fig05_link_conditions` | OP/ABER vs 发射功率、距离、入射角、孔径、风速、风向角（c=0.151/0.398） |
| `fig06_spots` | 线性/汇聚光斑、孔径功率、中心截面 |
| `fig07_motion_spots` | 110.8/19.7/7.4 mm 偏移光斑与功率时程 |
| `fig08_distance` | 5/10/15 m 光斑、功率–距离、捕获率与孔径对比 |
| `fig09_turbulence` | OP/ABER vs 湍流强度 |
| `fig10_dynamic_power` | 平均/10% 分位功率 vs 波幅（0.4/0.8/1.2 Hz） |
| `fig11_power_conditions` | 波幅、频率、入射角、风速、风向角、距离、孔径的功率对比 |
| `fig14_centroid_scatter` | 无跟踪/有跟踪质心散点 |
| `fig16_incidence` | 接收功率 vs 水入射角（曲线 + MC 菱形） |
| `fig17_aperture` | 圆斑/椭圆斑/宽斑的孔径功率差异（0.5–15 cm） |
| `fig18_op_conditions` | OP vs 波幅、频率、入射角、风速 |
| `fig19_aber_conditions` | ABER vs 波幅、频率、入射角、风速 |
| `fig20_distance_turb` | OP/ABER vs 距离与湍流 |
| `fig21_delay` | 质心、功率、归一化误差 vs 处理时延（c=0.151/0.398） |
| `fig22_update` | 更新周期 10–200 ms 与推荐 20–50 ms 区间 |

公共参数：532 nm、40×40 LC RIS（±15.414°）、10 m 水深、5–15 m 空气路径、5 cm 孔径、
OOK `R=0.5 A/W`、`σ_n=0.05 mA`、`P_th=1 mW`；条件 OP/ABER 取 `P_t=23 dBm`。
每幅图的详细图注、扫描范围与数据来源见 `figure_atlas_matlab/CAPTIONS.md`。

复现：

```matlab
cd FOWT_RIS_eval/figure_atlas_matlab
make_all
```

`make_all.m` 读取 `data/atlas_data.mat`（由 `export_atlas_data.py` 一次性整理，
只做数据准备，不做绘图），输出 300 dpi PNG 与可编辑 `.fig` 到 `figures/`。
旧的 `figure_atlas/` Python 图集与冗余图片已删除。
