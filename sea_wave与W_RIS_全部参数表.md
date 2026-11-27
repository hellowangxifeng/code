# sea_wave 与 W_RIS 全部参数表

> 2026-11-05 更新提示：本文下方为2026-10-30的历史参数快照，RIS实体、浮力、惯量及相关几何已更新，不能作为当前运行参数。当前完整输入见 `sea_wave/get_default_RIS_params.m`、`W_RIS/wris_config.m`；玻璃为1×1×0.01 m、浮框外尺寸1.2×1.2 m、实际总质量29.871 kg；静水下表面z=-3.432 mm、上表面z=+6.568 mm，详见 [RIS浮框配平说明](sea_wave/RIS浮框配平说明.md)。正式计算参数保存在本次MAT结果中。

核对日期：2026-10-30。依据本工作区当前源码；不以 README、旧 MAT/FIG/CSV 或历史图像中的值替代源码。

## 1. 范围、符号与读法

逐文件检查 `sea_wave/` 的 15 个、`W_RIS/` 的 17 个 MATLAB 文件（共 32 个），并补充实际调用的 `fig3-7/Table_abc.m`、`fig3-7/Hen_Green.m`。第 5 节按文件和函数作用域记录全部具名输入、输出、赋值变量、结构字段、循环索引、匿名函数输入及持久缓存；同名符号在不同作用域分别列出，共 1511 条符号记录。第 6 节补充没有独立变量名的调用选项和数值约束。

保留源码大小写和结构路径，例如 `sea_params.f_min`、`p.fMinHz`、`c.ris.fMinHz`。表中 `a`、`b`、`c`、`g`、`q`、`omega`、`R` 等必须连同文件/函数阅读；它们在不同函数有不同含义。标注“输入/输出”的符号，其具体值由调用方或上游计算确定；“计算式”保留源码表达式，不把派生量冒充独立可调参数。

世界坐标 z 向上。姿态采用右手主动旋转 `R=Rz(yaw)*Ry(pitch)*Rx(roll)`，`R` 将体坐标转为世界坐标。海浪传播风向从 +Y 向 +X 增大；该角与 RIS yaw、相对 +X 的光轴 azimuth 是不同变量。海浪相位为 `omega*absoluteTime-kx*x-ky*y`。没有注明角度单位的欧拉角内部量为 rad。

`eta=[surge;sway;heave;roll;pitch;yaw]`；兼容旧输出 `q=[heave;pitch;roll;surge;sway;yaw]`；`bodyVelocity=[u;v;w;p;q;r]` 中前三项 m/s、后三项 rad/s。体速度符号 `p/q/r` 是角速度分量。完整积分状态 `state=[position(3);quaternion(4);bodyVelocity(6);memoryX(3*Nr);memoryY(3*Nr)]`，四元数是 wxyz。`initialPosition(1:2)` 与 `horizontalInitialPosition` 是相对 `p.center` 的水平状态；输出世界质心再加 `p.center`，z 再加 `zG`。

RIS 的 `boardCenter`、`corners` 是上表面数据；光学接收中心 `g.pRis=centerOfMass-thickness/2*normal` 是下表面。第 2～4 节适合直接查配置；第 5～6 节用于逐源码核对。方括号内的数值严格保持源码顺序。

## 2. sea_wave 基础默认配置

### 2.1 海面默认参数：get_default_sea_params

| 源码符号 | 含义与单位 | 默认值/计算式 | 源码位置 |
| --- | --- | --- | --- |
| `sea_params.g` | 重力加速度（m/s²） | `9.81` | [sea_wave/get_default_sea_params.m:7](<C:/Users/Spidey/Desktop/code/sea_wave/get_default_sea_params.m:7>) |
| `sea_params.neumannConstant` | 源码 Neumann 频谱经验常数，默认 3.05；维度随源码谱表达式确定 | `3.05` | [sea_wave/get_default_sea_params.m:8](<C:/Users/Spidey/Desktop/code/sea_wave/get_default_sea_params.m:8>) |
| `sea_params.x_range` | 海面世界 x/y 轴显示与重建范围 [min,max]（m） | `[-5, 5]` | [sea_wave/get_default_sea_params.m:12](<C:/Users/Spidey/Desktop/code/sea_wave/get_default_sea_params.m:12>) |
| `sea_params.y_range` | 海面世界 x/y 轴显示与重建范围 [min,max]（m） | `[-5, 5]` | [sea_wave/get_default_sea_params.m:13](<C:/Users/Spidey/Desktop/code/sea_wave/get_default_sea_params.m:13>) |
| `sea_params.f_max` | 海浪最高频率（Hz） | `4.0` | [sea_wave/get_default_sea_params.m:16](<C:/Users/Spidey/Desktop/code/sea_wave/get_default_sea_params.m:16>) |
| `sea_params.grid_density` | 海面每个空间方向网格点数，由最高频率和最大跨度反推 | `ceil(L_max / dx_required) + 1` | [sea_wave/get_default_sea_params.m:29](<C:/Users/Spidey/Desktop/code/sea_wave/get_default_sea_params.m:29>) |
| `sea_params.f_min` | 海浪最低频率（Hz） | `0.05` | [sea_wave/get_default_sea_params.m:32](<C:/Users/Spidey/Desktop/code/sea_wave/get_default_sea_params.m:32>) |
| `sea_params.Fre_segment` | 对数频率离散点数（不是点数减一的区间数） | `50` | [sea_wave/get_default_sea_params.m:33](<C:/Users/Spidey/Desktop/code/sea_wave/get_default_sea_params.m:33>) |
| `sea_params.Dir_segment` | 方向离散点数，包含 ±π/2 端点 | `36` | [sea_wave/get_default_sea_params.m:34](<C:/Users/Spidey/Desktop/code/sea_wave/get_default_sea_params.m:34>) |
| `sea_params.U` | 风速（m/s） | `10` | [sea_wave/get_default_sea_params.m:35](<C:/Users/Spidey/Desktop/code/sea_wave/get_default_sea_params.m:35>) |
| `sea_params.t0` | 绝对海浪起算时刻（s） | `0` | [sea_wave/get_default_sea_params.m:36](<C:/Users/Spidey/Desktop/code/sea_wave/get_default_sea_params.m:36>) |
| `sea_params.seed` | 随机数种子；所属结构决定是海浪、TX、MC 或湍流随机流 | `42` | [sea_wave/get_default_sea_params.m:37](<C:/Users/Spidey/Desktop/code/sea_wave/get_default_sea_params.m:37>) |


由当前默认值计算，`dx_required≈0.04879094 m`，`L_x=L_y=L_max=10 m`，`sea_params.grid_density=206`。两端均包含，所以实际 206 点网格间距为 `10/205≈0.04878049 m`。`frequencyBins`/`Fre_segment` 与 `directionBins`/`Dir_segment` 实际用于 `logspace`/`linspace` 点数。

### 2.2 RIS/动力学默认值及 W_RIS 覆盖后的有效值

| 源码符号及继承字段 | 含义与单位 | sea_wave 默认值 | W_RIS 主流程有效值 | 默认定义位置 |
| --- | --- | --- | --- | --- |
| `p.Lx` / `c.ris.Lx` | RIS 体坐标 x 方向边长（m） | `1` | `1` | [sea_wave/get_default_RIS_params.m:4](<C:/Users/Spidey/Desktop/code/sea_wave/get_default_RIS_params.m:4>) |
| `p.Ly` / `c.ris.Ly` | RIS 体坐标 y 方向边长（m） | `1` | `1` | [sea_wave/get_default_RIS_params.m:4](<C:/Users/Spidey/Desktop/code/sea_wave/get_default_RIS_params.m:4>) |
| `p.mass` / `c.ris.mass` | RIS 质量（kg）；模型内 mass 为含附加质量的六分量广义质量 | `12` | `12` | [sea_wave/get_default_RIS_params.m:4](<C:/Users/Spidey/Desktop/code/sea_wave/get_default_RIS_params.m:4>) |
| `p.thickness` / `c.ris.thickness` | RIS 板厚（m） | `0.03` | `0.03` | [sea_wave/get_default_RIS_params.m:4](<C:/Users/Spidey/Desktop/code/sea_wave/get_default_RIS_params.m:4>) |
| `p.rho` / `c.ris.rho` | 海水密度（kg/m³） | `1025` | `1025` | [sea_wave/get_default_RIS_params.m:4](<C:/Users/Spidey/Desktop/code/sea_wave/get_default_RIS_params.m:4>) |
| `p.g` / `c.ris.g` | 重力加速度（m/s²） | `9.81` | `9.81` | [sea_wave/get_default_RIS_params.m:4](<C:/Users/Spidey/Desktop/code/sea_wave/get_default_RIS_params.m:4>) |
| `p.windSpeed` / `c.ris.windSpeed` | 风速（m/s） | `10` | `主图各列覆盖为 [3,3,5,9] m/s` | [sea_wave/get_default_RIS_params.m:5](<C:/Users/Spidey/Desktop/code/sea_wave/get_default_RIS_params.m:5>) |
| `p.windAngleDeg` / `c.ris.windAngleDeg` | 波浪传播方向角（°），从 +Y 朝 +X 增大 | `30` | `主图各列覆盖为 [30,50,50,30]°` | [sea_wave/get_default_RIS_params.m:5](<C:/Users/Spidey/Desktop/code/sea_wave/get_default_RIS_params.m:5>) |
| `p.duration` / `c.ris.duration` | 保存的运动持续时间（s） | `10` | `20` | [sea_wave/get_default_RIS_params.m:5](<C:/Users/Spidey/Desktop/code/sea_wave/get_default_RIS_params.m:5>) |
| `p.dt` / `c.ris.dt` | 运动输出采样间隔（s），不等于实际 RK4 内部步长 | `0.025` | `0.05` | [sea_wave/get_default_RIS_params.m:5](<C:/Users/Spidey/Desktop/code/sea_wave/get_default_RIS_params.m:5>) |
| `p.seed` / `c.ris.seed` | 随机数种子；所属结构决定是海浪、TX、MC 或湍流随机流 | `42` | `600` | [sea_wave/get_default_RIS_params.m:5](<C:/Users/Spidey/Desktop/code/sea_wave/get_default_RIS_params.m:5>) |
| `p.startTime` / `c.ris.startTime` | 绝对海浪起算时刻（s） | `0` | `25` | [sea_wave/get_default_RIS_params.m:6](<C:/Users/Spidey/Desktop/code/sea_wave/get_default_RIS_params.m:6>) |
| `p.center` / `c.ris.center` | RIS 的水平参考中心 [x;y]（m），或光斑加权中心 [x,y]（m）；以父结构为准 | `[0;0]` | `静水 [0;0]；各运动面板改为所选 site.center` | [sea_wave/get_default_RIS_params.m:6](<C:/Users/Spidey/Desktop/code/sea_wave/get_default_RIS_params.m:6>) |
| `p.fMinHz` / `c.ris.fMinHz` | 海浪最低频率（Hz） | `0.05` | `0.02` | [sea_wave/get_default_RIS_params.m:6](<C:/Users/Spidey/Desktop/code/sea_wave/get_default_RIS_params.m:6>) |
| `p.fMaxHz` / `c.ris.fMaxHz` | 海浪最高频率（Hz） | `4` | `3` | [sea_wave/get_default_RIS_params.m:6](<C:/Users/Spidey/Desktop/code/sea_wave/get_default_RIS_params.m:6>) |
| `p.frequencyBins` / `c.ris.frequencyBins` | 对数频率离散点数（不是点数减一的区间数） | `50` | `200` | [sea_wave/get_default_RIS_params.m:7](<C:/Users/Spidey/Desktop/code/sea_wave/get_default_RIS_params.m:7>) |
| `p.directionBins` / `c.ris.directionBins` | 方向离散点数，包含 ±π/2 端点 | `36` | `36` | [sea_wave/get_default_RIS_params.m:7](<C:/Users/Spidey/Desktop/code/sea_wave/get_default_RIS_params.m:7>) |
| `p.panelCount` / `c.ris.panelCount` | 薄板 BEM 每边面元数，须为 4..80 的偶数；完整网格 panelCount² | `24` | `24` | [sea_wave/get_default_RIS_params.m:7](<C:/Users/Spidey/Desktop/code/sea_wave/get_default_RIS_params.m:7>) |
| `p.horizontalInitialPosition` / `c.ris.horizontalInitialPosition` | 相对 center 的初始水平位移 [x;y]（m） | `[0;0]` | `[0;0]` | [sea_wave/get_default_RIS_params.m:8](<C:/Users/Spidey/Desktop/code/sea_wave/get_default_RIS_params.m:8>) |
| `p.horizontalInitialVelocity` / `c.ris.horizontalInitialVelocity` | 可选世界系初始水平速度 [vx;vy]（m/s），[] 保留初始化结果 | `[]` | `[]` | [sea_wave/get_default_RIS_params.m:8](<C:/Users/Spidey/Desktop/code/sea_wave/get_default_RIS_params.m:8>) |
| `p.initialPosition` / `c.ris.initialPosition` | 可选初始 [surge;sway;heave] 状态（m）；水平分量仍相对 center，[] 不覆盖 | `[]` | `[]` | [sea_wave/get_default_RIS_params.m:9](<C:/Users/Spidey/Desktop/code/sea_wave/get_default_RIS_params.m:9>) |
| `p.initialEulerAngles` / `c.ris.initialEulerAngles` | 可选初始 [roll;pitch;yaw]（rad），[] 不覆盖；默认初始 yaw=0 | `[]` | `[]` | [sea_wave/get_default_RIS_params.m:9](<C:/Users/Spidey/Desktop/code/sea_wave/get_default_RIS_params.m:9>) |
| `p.initialBodyVelocity` / `c.ris.initialBodyVelocity` | 可选体坐标初始 [u;v;w;p;q;r]（前 3 项 m/s，后 3 项 rad/s） | `[]` | `[]` | [sea_wave/get_default_RIS_params.m:9](<C:/Users/Spidey/Desktop/code/sea_wave/get_default_RIS_params.m:9>) |
| `p.initialization` / `c.ris.initialization` | 初始化方式：waveSteady（频域稳态）或 rest（静止状态） | `'waveSteady'` | `'waveSteady'` | [sea_wave/get_default_RIS_params.m:10](<C:/Users/Spidey/Desktop/code/sea_wave/get_default_RIS_params.m:10>) |
| `p.radiationModes` / `c.ris.radiationModes` | 被动辐射记忆有理拟合的振子模态数 | `64` | `64` | [sea_wave/get_default_RIS_params.m:10](<C:/Users/Spidey/Desktop/code/sea_wave/get_default_RIS_params.m:10>) |
| `p.radiationFrequencyBins` / `c.ris.radiationFrequencyBins` | 辐射拟合额外对数频率采样点数 | `80` | `80` | [sea_wave/get_default_RIS_params.m:10](<C:/Users/Spidey/Desktop/code/sea_wave/get_default_RIS_params.m:10>) |
| `p.integrationStep` / `c.ris.integrationStep` | RK4 期望最大内步长（s），还按稳定性自动细分 | `0.01` | `0.01` | [sea_wave/get_default_RIS_params.m:11](<C:/Users/Spidey/Desktop/code/sea_wave/get_default_RIS_params.m:11>) |
| `p.wave` / `c.ris.wave` | 外部共享的海浪分量结构；[] 时按参数生成 | `[]` | `静水不求运动；运动面板传入各列共享 wave` | [sea_wave/get_default_RIS_params.m:11](<C:/Users/Spidey/Desktop/code/sea_wave/get_default_RIS_params.m:11>) |
| `p.plot` / `c.ris.plot` | 是否自动绘制 RIS 单次运动；RIS_float 无输出且未覆盖时设为 true | `false` | `false` | [sea_wave/get_default_RIS_params.m:11](<C:/Users/Spidey/Desktop/code/sea_wave/get_default_RIS_params.m:11>) |


默认板面积 `area=1 m²`，`draft=12/(1025*1)≈0.0117073171 m`，`freeboard≈0.0182926829 m`，`zG≈0.0032926829 m`，`inertia=[1.0009;1.0009;2] kg·m²`。`RIS_float` 首先继承默认结构，再以调用方 `cfg` 覆盖；`initialBodyVelocity` 最后覆盖速度，`initialPosition` 和 `initialEulerAngles` 覆盖对应状态。`horizontalInitialVelocity` 只覆盖世界水平速度。`waveSteady` 初始化会生成垂向稳态运动与辐射记忆；`rest` 不叠加这些稳态初值。

### 2.3 独立入口/批量脚本的设置

| 入口/函数 | 源码符号 | 取值 | 含义 |
| --- | --- | --- | --- |
| `run_RIS_wind_comparison` | `winds / directions` | `[3,5,9,12] m/s / [30,60,90]°` | 4×3 动力学工况；第 4 图列为风向 30° 轨迹 |
| `run_RIS_wind_comparison` | `startTime / duration / snapshotTime / seed` | `25 s / 10 s / 2 s / 600` | 绝对运动时间 25..35 s，快照为绝对 27 s |
| `run_RIS_wind_comparison` | `zLimits` | `[-0.5,0.5] m` | 仅显示裁切；不改变高度原始数组 |
| `main3_windseasurface_0326` | `params.t0 / params.seed / exp_cases` | `27 s / 600 / [9,30;9,45;12,30]` | 纯海面，工况每行 [风速 m/s,传播风向 °] |
| `audit_sea_spectrum` | `winds / fine.frequencyBins` | `[3,5,9,12] / 400` | 连续谱与默认/加密频率积分比较 |
| `audit_sea_spectrum` | `p.windSpeed / p.fMinHz` | `12 / 0.1` | 仅旧端点错误重现算例覆盖 |
| `summarize_RIS_motion_existing` | `已有轨迹要求` | `12 文件、401 样本、dt=0.025 s、绝对 25..35 s` | 总体协方差/标准差按 N 归一化；均值移除、趋势保留；另外输出去线性趋势 RMS |
| `oceanSurfaceModel_phiwind_260315` | `seed / output_mode / plane_index` | `1 / all / [1,1]` | 省略第 11/12/13 输入时的接口默认值，其他输入必填 |
| `RIS_wave_components` | `quadrature` | `trapezoid` | 唯一允许规则，频率和方向两端均半权重 |
| `RIS_draw_state` | `heightLimits` | `[-0.5,0.5] m` | 可选显示高度范围 |
| `wris_geometry` | `mode` | `joint` | 可选 static/ris-only/tx-only/joint |
| `wris_mc` | `gain` | `1` | 省略第 4 输入时无公共强度增益 |
| `wris_egg` | `kind / n / seed` | `weak / 1 / []` | 种子为空时沿用调用方 RNG |
| `wris_config / wris_run` | `profile` | `production` | smoke 的光子数/重复数见下表 |


## 3. W_RIS 总配置：wris_config

`c.ris` 的全部继承字段与覆盖值已在第 2.2 节逐项列出。下表保留配置入口的实际赋值顺序；同一字段多次赋值时用换行列出。`c.comparison.incidenceDeg` 在主流程逐行覆盖 `c.geometry.incidenceDeg`。

| 源码符号 | 含义与单位 | 默认值/计算式及分支值 | 源码位置 |
| --- | --- | --- | --- |
| `c.version` | 配置版本标记字符串 | `'wris-3-rotation-table'` | [W_RIS/wris_config.m:4](<C:/Users/Spidey/Desktop/code/W_RIS/wris_config.m:4>) |
| `c.root` | 源码 W_RIS 目录 | `fileparts(mfilename('fullpath'))` | [W_RIS/wris_config.m:5](<C:/Users/Spidey/Desktop/code/W_RIS/wris_config.m:5>) |
| `c.geometry.distance` | TX 至 RIS 光学中心的空间直线距离（m） | `10` | [W_RIS/wris_config.m:6](<C:/Users/Spidey/Desktop/code/W_RIS/wris_config.m:6>) |
| `c.geometry.incidenceDeg` | 静水初始光轴相对 RIS 法线的入射角（°） | `60` | [W_RIS/wris_config.m:7](<C:/Users/Spidey/Desktop/code/W_RIS/wris_config.m:7>) |
| `c.geometry.azimuthDeg` | TX 光轴水平投影相对世界 +X 的方位角（°） | `90` | [W_RIS/wris_config.m:8](<C:/Users/Spidey/Desktop/code/W_RIS/wris_config.m:8>) |
| `c.ris` | 继承 get_default_RIS_params 后由 wris_config 覆盖的 RIS/海浪/动力学配置 | `get_default_RIS_params()` | [W_RIS/wris_config.m:9](<C:/Users/Spidey/Desktop/code/W_RIS/wris_config.m:9>) |
| `c.comparison.incidenceDeg` | 静水初始光轴相对 RIS 法线的入射角（°） | `[60 45 30]` | [W_RIS/wris_config.m:15](<C:/Users/Spidey/Desktop/code/W_RIS/wris_config.m:15>) |
| `c.comparison.conditions` | 四个运动列的 [风速(m/s),传播风向(°)] 工况矩阵 | `[3 30;3 50;5 50;9 30]` | [W_RIS/wris_config.m:16](<C:/Users/Spidey/Desktop/code/W_RIS/wris_config.m:16>) |
| `c.comparison.phaseIndex` | 3×4 站位选择索引：1峰、2谷、3空间上升坡、4空间下降坡 | `[1 1 1 1;2 2 2 2;3 4 3 4]` | [W_RIS/wris_config.m:17](<C:/Users/Spidey/Desktop/code/W_RIS/wris_config.m:17>) |
| `c.sampling.comparisonTime` | 所有面板统一的相对采样时刻（s） | `20` | [W_RIS/wris_config.m:18](<C:/Users/Spidey/Desktop/code/W_RIS/wris_config.m:18>) |
| `c.sea.quadrature` | 频率/方向积分规则，当前仅支持 trapezoid | `'trapezoid'` | [W_RIS/wris_config.m:19](<C:/Users/Spidey/Desktop/code/W_RIS/wris_config.m:19>) |
| `c.sea.siteSamples` | 沿主传播方向的空间剖面采样点数 | `801` | [W_RIS/wris_config.m:19](<C:/Users/Spidey/Desktop/code/W_RIS/wris_config.m:19>) |
| `c.sea.neumannConstant` | 源码 Neumann 频谱经验常数，默认 3.05；维度随源码谱表达式确定 | `seaDefaults.neumannConstant` | [W_RIS/wris_config.m:21](<C:/Users/Spidey/Desktop/code/W_RIS/wris_config.m:21>) |
| `c.tx.seed` | 随机数种子；所属结构决定是海浪、TX、MC 或湍流随机流 | `111` | [W_RIS/wris_config.m:22](<C:/Users/Spidey/Desktop/code/W_RIS/wris_config.m:22>) |
| `c.tx.initialYaw` | TX 平台初始航向角（rad） | `pi/2` | [W_RIS/wris_config.m:22](<C:/Users/Spidey/Desktop/code/W_RIS/wris_config.m:22>) |
| `c.tx.meanUVR` | TX 体坐标 [u;v;r] 的 OU 平稳均值（m/s,m/s,rad/s） | `[0.008;0;deg2rad(0.03)]` | [W_RIS/wris_config.m:23](<C:/Users/Spidey/Desktop/code/W_RIS/wris_config.m:23>) |
| `c.tx.stdUVR` | TX 体坐标 [u;v;r] 的 OU 平稳标准差（m/s,m/s,rad/s） | `[0.012;0.008;deg2rad(0.12)]` | [W_RIS/wris_config.m:24](<C:/Users/Spidey/Desktop/code/W_RIS/wris_config.m:24>) |
| `c.tx.correlationTime` | TX 三个 OU 分量的相关时间（s） | `[2;2;2]` | [W_RIS/wris_config.m:25](<C:/Users/Spidey/Desktop/code/W_RIS/wris_config.m:25>) |
| `c.tx.pointing` | 指向方式，当前仅支持 body-fixed（光轴随平台偏航） | `'body-fixed'` | [W_RIS/wris_config.m:26](<C:/Users/Spidey/Desktop/code/W_RIS/wris_config.m:26>) |
| `c.source.power` | 发射光功率（W） | `1` | [W_RIS/wris_config.m:27](<C:/Users/Spidey/Desktop/code/W_RIS/wris_config.m:27>) |
| `c.source.wavelength` | 真空/模型使用的光波长（m），默认 532e-9 | `532e-9` | [W_RIS/wris_config.m:27](<C:/Users/Spidey/Desktop/code/W_RIS/wris_config.m:27>) |
| `c.source.divergence` | 源码高斯径向发射律的角尺度（rad） | `0.01` | [W_RIS/wris_config.m:28](<C:/Users/Spidey/Desktop/code/W_RIS/wris_config.m:28>) |
| `c.source.waist` | 光源发射面高斯束腰尺度 wavelength/(pi*divergence)（m） | `c.source.wavelength/(pi*c.source.divergence)` | [W_RIS/wris_config.m:29](<C:/Users/Spidey/Desktop/code/W_RIS/wris_config.m:29>) |
| `c.water.type` | Table_abc 水体档位字符串 | `'c0.151'` | [W_RIS/wris_config.m:30](<C:/Users/Spidey/Desktop/code/W_RIS/wris_config.m:30>) |
| `c.water.absorption` | 水体吸收系数（m⁻¹），Table_abc 返回 a=c-b | `Table_abc(c.water.type)` | [W_RIS/wris_config.m:31](<C:/Users/Spidey/Desktop/code/W_RIS/wris_config.m:31>) |
| `c.water.scattering` | 水体散射系数（m⁻¹），Table_abc 返回 b | `Table_abc(c.water.type)` | [W_RIS/wris_config.m:31](<C:/Users/Spidey/Desktop/code/W_RIS/wris_config.m:31>) |
| `c.water.attenuation` | 水体总衰减系数（m⁻¹），Table_abc 返回 c=a+b；MC 内重新用 a+b | `Table_abc(c.water.type)` | [W_RIS/wris_config.m:31](<C:/Users/Spidey/Desktop/code/W_RIS/wris_config.m:31>) |
| `c.water.anisotropy` | HG 散射相函数各向异性因子，-1<值<1 | `0.924` | [W_RIS/wris_config.m:32](<C:/Users/Spidey/Desktop/code/W_RIS/wris_config.m:32>) |
| `c.turbulence.kind` | EGG 档位 weak/strong/none；缓存函数中为 motion/optics 类型 | `'weak'` | [W_RIS/wris_config.m:33](<C:/Users/Spidey/Desktop/code/W_RIS/wris_config.m:33>) |
| `c.turbulence.seed` | 随机数种子；所属结构决定是海浪、TX、MC 或湍流随机流 | `260617` | [W_RIS/wris_config.m:33](<C:/Users/Spidey/Desktop/code/W_RIS/wris_config.m:33>) |
| `c.mc.photons` | 单帧发射光子包总数，或验证采用的光子数档位 | `1e6`<br>`10000` | [W_RIS/wris_config.m:34](<C:/Users/Spidey/Desktop/code/W_RIS/wris_config.m:34>) |
| `c.mc.chunkSize` | MC 每批最大光子包数 | `100000` | [W_RIS/wris_config.m:34](<C:/Users/Spidey/Desktop/code/W_RIS/wris_config.m:34>) |
| `c.mc.maxScatter` | 每个光子包允许的最大散射次数 | `1000` | [W_RIS/wris_config.m:35](<C:/Users/Spidey/Desktop/code/W_RIS/wris_config.m:35>) |
| `c.mc.maxPath` | 每个光子包允许的最大累计路径长度（m） | `200` | [W_RIS/wris_config.m:35](<C:/Users/Spidey/Desktop/code/W_RIS/wris_config.m:35>) |
| `c.mc.rouletteThreshold` | 俄罗斯轮盘赌启动权重门限（无量纲） | `1e-4` | [W_RIS/wris_config.m:36](<C:/Users/Spidey/Desktop/code/W_RIS/wris_config.m:36>) |
| `c.mc.rouletteSurvival` | 低权重包轮盘赌存活概率（无量纲） | `0.1` | [W_RIS/wris_config.m:36](<C:/Users/Spidey/Desktop/code/W_RIS/wris_config.m:36>) |
| `c.mc.seed` | 随机数种子；所属结构决定是海浪、TX、MC 或湍流随机流 | `20261014` | [W_RIS/wris_config.m:37](<C:/Users/Spidey/Desktop/code/W_RIS/wris_config.m:37>) |
| `c.grid.bins` | 板面 x/y 方向分箱数 [nx,ny] | `[80 80]` | [W_RIS/wris_config.m:38](<C:/Users/Spidey/Desktop/code/W_RIS/wris_config.m:38>) |
| `c.grid.planeScale` | 延展评分窗口相对板面的边长倍数；保持单元面积不变 | `6` | [W_RIS/wris_config.m:38](<C:/Users/Spidey/Desktop/code/W_RIS/wris_config.m:38>) |
| `c.validation.photons` | 单帧发射光子包总数，或验证采用的光子数档位 | `[10000 100000 1000000]`<br>`[2000 10000 30000]` | [W_RIS/wris_config.m:39](<C:/Users/Spidey/Desktop/code/W_RIS/wris_config.m:39>) |
| `c.validation.repeats` | 每个光子数档位独立重复次数 | `4`<br>`3` | [W_RIS/wris_config.m:40](<C:/Users/Spidey/Desktop/code/W_RIS/wris_config.m:40>) |
| `c.validation.seed` | 随机数种子；所属结构决定是海浪、TX、MC 或湍流随机流 | `9500` | [W_RIS/wris_config.m:40](<C:/Users/Spidey/Desktop/code/W_RIS/wris_config.m:40>) |
| `c.output` | 输出目录 | `fullfile(c.root,'results',profile)` | [W_RIS/wris_config.m:41](<C:/Users/Spidey/Desktop/code/W_RIS/wris_config.m:41>) |
| `c.cache.folder` | 缓存目录 | `fullfile(c.root,'cache')` | [W_RIS/wris_config.m:42](<C:/Users/Spidey/Desktop/code/W_RIS/wris_config.m:42>) |
| `c.cache.enabled` | 是否复用物理计算缓存 | `true` | [W_RIS/wris_config.m:42](<C:/Users/Spidey/Desktop/code/W_RIS/wris_config.m:42>) |
| `c.plot.resolution` | 导出 PNG 分辨率（dpi） | `300` | [W_RIS/wris_config.m:43](<C:/Users/Spidey/Desktop/code/W_RIS/wris_config.m:43>) |
| `c.plot.dynamicRange` | 共享对数色标的最大/最小辐照度比 | `1e4` | [W_RIS/wris_config.m:43](<C:/Users/Spidey/Desktop/code/W_RIS/wris_config.m:43>) |
| `c.plot.padding` | 板面边界外绘图区留白相对板边长的比例 | `0.15` | [W_RIS/wris_config.m:43](<C:/Users/Spidey/Desktop/code/W_RIS/wris_config.m:43>) |


`c.water.type='c0.151'` 实际得到 `absorption=0.114 m⁻¹`、`scattering=0.037 m⁻¹`、`attenuation=0.151 m⁻¹`；`source.waist≈1.69340859e-5 m`。`grid.bins=[80,80]`、板长 1 m 时单格边长 0.0125 m、`cellArea=1.5625e-4 m²`；`planeScale=6` 给出 480×480 格、6×6 m 的延展显示/评分窗，实际接受孔径仍为 1×1 m。

| 配置项 | production | smoke |
| --- | --- | --- |
| `c.mc.photons` | 1000000 | 10000 |
| `c.validation.photons` | [10000,100000,1000000] | [2000,10000,30000] |
| `c.validation.repeats` | 4 | 3 |
| `c.output` | `fullfile(c.root,'results','production')` | `fullfile(c.root,'results','smoke')` |


主表为 3×5：行入射角 `[60,45,30]°`；第 1 列静水、TX/RIS 固定；后四列 `[U,windAngleDeg]=[3,30;3,50;5,50;9,30]`。站位矩阵 `phaseIndex=[1,1,1,1;2,2,2,2;3,4,3,4]` 中 1=空间峰、2=空间谷、3=空间上升过零、4=空间下降过零；统一相对时刻 20 s、绝对海浪时刻 45 s。`c.sea.neumannConstant` 可单独改变共享波浪幅度，实际缩放为 `sqrt(c.sea.neumannConstant/defaults.neumannConstant)`。

## 4. EGG 参数、外部水体参数与内置常量

### 4.1 EGG 全部分支

`p=[w,lambda,a,b,c]`；`info.parameters` 的对应符号是 `omega,lambda_e,a_g,b_g,c_g`。EGG 的 `omega` 是混合概率，区别于海浪 `wave.omega`（rad/s）。原始变量服从 `w*Exp(lambda)+(1-w)*GG(a,b,c)`，输出 `h=x/mu`；每个快照的所有光子共享同一个 h。

| 源码符号 | weak | strong | none | 含义 |
| --- | --- | --- | --- | --- |
| `w / info.parameters.omega` | 0.2130 | 0.4539 | 0 | 指数分量混合概率 |
| `lambda / info.parameters.lambda_e` | 0.3291 | 0.2744 | NaN | 指数分量尺度 |
| `a / info.parameters.a_g` | 1.4299 | 0.3008 | NaN | 广义 Gamma 形状 |
| `b / info.parameters.b_g` | 1.1817 | 1.7053 | NaN | 广义 Gamma 尺度 |
| `c / info.parameters.c_g` | 17.1984 | 54.1422 | NaN | 广义 Gamma 幂指数 |
| `tabSI` | 0.1484 | 0.4769 | 0 | 原表闪烁指数；区别于按分布矩计算的 theoreticalScintillationIndex |
| `flow` | 2.4 | 4.7 | 0 | 气泡流量（L/min），元数据 |
| `gradient` | 0.05 | 0.1 | 0 | 温度梯度（°C/cm），元数据 |
| `k / info.momentOrders` | 0:4 | 0:4 | 0:4 | 解析矩阶数 |
| `h` | x/mu | x/mu | ones(n,1) | 均值为 1 的快照增益 |


### 4.2 Table_abc 的全部可选水体与 Hen_Green 分支

| Table_abc 源码符号 | water_type | b（m⁻¹） | c（m⁻¹） | a=c-b（m⁻¹） |
| --- | --- | --- | --- | --- |
| `c3_2` | `c3.2` | 2.6 | 3.2 | 0.6 |
| `c2_675` | `c2.675` | 2.1734 | 2.675 | 0.5016 |
| `c1_92` | `c1.92` | 1.56 | 1.92 | 0.36 |
| `c1_475` | `c1.475` | 1.1984 | 1.475 | 0.2766 |
| `c1_0` | `c1.0` | 0.8125 | 1 | 0.1875 |
| `c0_5` | `c0.5` | 0.406 | 0.5 | 0.094 |
| `c0_275` | `c0.275` | 0.2234 | 0.275 | 0.0516 |
| `c0_15` | `c0.15` | 0.08 | 0.15 | 0.07 |
| `c0_305` | `c0.305` | 0.216 | 0.305 | 0.089 |
| `c0_151` | `c0.151` | 0.037 | 0.151 | 0.114 |
| `c0_398` | `c0.398` | 0.219 | 0.398 | 0.179 |
| `c2_190` | `c2.190` | 1.824 | 2.19 | 0.366 |
| `c2_17` | `c2.17` | 1.875 | 2.17 | 0.295 |


`Hen_Green(length,judge)`：`judge==1` 取 `g=0.86`；其他取 `g=0.924`，W_RIS 默认调用 `Hen_Green(ns,2)`。当 `c.water.anisotropy` 不等于 0.924 时，`wris_mc` 自己按 HG 逆变换抽样；`abs(gg)<1e-12` 时用各向同性 `mu=2*u-1`。

### 4.3 MATLAB 内置量与关键无独立变量名常数

| 源码常量/表达式 | 含义 | 实际用途 |
| --- | --- | --- |
| `pi` | 圆周率 | `3.141592653589793…，用于频率/角频率与角度转换` |
| `1i` | 虚数单位 | 复谐波、频域阻抗及 Green 核的虚部 |
| `realmin` | 最小正规正 double，约 2.2250738585072014e-308 | 防止 log(0) 或分母为 0；由 MATLAB 精度决定 |
| `eps` | double 在 1 附近的机器间隔，约 2.220446049250313e-16 | 共享色限的非零下限 |
| `Inf / NaN` | 无穷/未定义 | 无前向交点用 Inf；无有效主轴或空统计用 NaN |
| `true / false` | 布尔常量 | 功能开关、掩码和适用性标记 |
| `nargin / nargout` | MATLAB 输入/输出个数 | 决定接口默认值和自动绘图开关 |
| `4 / 64` | solvePotentials 高斯积分阶数 | 面元每方向 4 点（二维 16 点）；H0 用 64 点 |
| `0.02` | RIS_radiation 加权尺度下限比例 | scale=max(abs(target),0.02*max(abs(target))) |
| `0.6 / 0.3` | RK4 稳定性步长系数 | safeStep=min(requestedStep,0.6/poleRate,0.3/naturalRate) |
| `1e-3` | 侧壁深度/无 pi sinc 小量切换门限 | kd>1e-3 用解析力臂；abs(x)<1e-3 用 Taylor 展开 |
| `1e-9` | 深水色散相对一致性容差 | abs(hypot(kx,ky)-omega²/g) <= 1e-9*max(1,omega²/g) |
| `pi/max(w.omega)` | 运动输出采样时间上限（s） | p.dt 不能超过入射波 Nyquist 限制 |
| `1e-6` | BEM 非被动负阻尼容差 | 任一 B<-1e-6 时报错 |
| `10` | 小角度诊断范围（°） | 法向倾角超过 10° 则 smallAngleUnder10Deg=false |
| `4` | plateResponse 最大缓存条目数 | 新增第 5 套前移除最早的一套 |
| `64 / 2e6` | RIS_wave_field 内存分块限制 | block=max(1,min(64,floor(2e6/Npoints))) |
| `1e-10` | 多个快照采样匹配容差（s） | RIS_draw_state、批量快照及 W_RIS 比较时间使用 |
| `1e-12` | 射线平行/朝向及 HG 各向同性门限 | 具体比较方向由对应源码行给出 |
| `1e-10` | MC 源平面逃逸最短距离（m） | back<1e-10 的源面自交排除 |
| `0.9` | MC 极点辅助轴切换门限 | abs(oldDir(:,3))>0.9 时 helper 从 z 轴改为 x 轴 |
| `4` | D4sigma 光斑直径系数 | diam=4*sqrt(协方差特征值) |
| `0.05` | 透明水体理论轴比绝对误差门限 | abs(MC ratio-sec(alpha))<0.05 |
| `1e-10` | W_RIS 能量残差/功率积分误差门限 | 分别是无量纲残差与 W 单位绝对误差 |
| `0.25 / 4 / 0.02` | MC 主轴统计容差参数 | 0.25° 基础项、4 倍抽样标准误、轴比差分母下限 0.02 |
| `5` | 最高两档功率均值收敛容差倍数 | 小于 5 倍合成均值标准误 |
| `100*j+k / 901+row` | 验证种子偏移 | 重复种子 validation.seed+100*j+k；透明水体种子 901+row |
| `0.002 / 0.00003 / 1e-12` | 海谱积分审计门限 | 默认方差相对误差、400 点方差相对误差、频率权重和绝对误差 |
| `1e-13 / 1e-11` | bandVariance 数值积分容差 | AbsTol/RelTol |
| `2^32-1` | 随机种子最大允许值 | RIS 与 EGG 输入验证 |
| `[5,10]` | 运动窗口汇总终点（s） | 仅输出不超过当前 duration 的窗口 |
| `[-90,90) / 180 / 90` | 无向椭圆主轴角折叠约定（°） | mod(angle+90,180)-90；法向入射几何主轴为 NaN |
| `1 / 0 / 0.5` | 积分、几何和高斯定义中的系数 | 单位/零阵初始化、半板尺寸、半角及 RK4 标准系数等保留在第 5 节计算式，不单独当可调物理参数 |


### 4.4 表格/CSV 中以字符串定义的列名

这些列名直接出现在 `VariableNames`、`names` 或 CSV 读取表达式中，虽然不是独立 MATLAB 变量，也按源码原拼写补充说明。`snapshots.csv` 与运动统计 CSV 的结构字段已在第 5 节 `records.*`、`risRows.*` 中逐项列明。

| 输出表 | 源码列名 | 含义与单位 |
| --- | --- | --- |
| `RIS_trajectory_*.csv` | `elapsedTimeS` | 相对运动时间（s） |
| `RIS_trajectory_*.csv` | `seaTimeS` | 绝对海浪时刻 startTime+t（s） |
| `RIS_trajectory_*.csv` | `centerXM / centerYM / centerZM` | 质心世界 x/y/z 坐标（m） |
| `RIS_trajectory_*.csv` | `rollRad / pitchRad / yawRad` | 右手 ZYX 欧拉角（rad） |
| `RIS_trajectory_*.csv` | `uMps / vMps / wMps` | 体坐标 x/y/z 平移速度（m/s） |
| `RIS_trajectory_*.csv` | `pRadps / qRadps / rRadps` | 体坐标 x/y/z 角速度（rad/s），不是欧拉角导数 |
| `RIS_case_summary.csv` | `windMps / directionDeg / HsM` | 风速（m/s）/传播风向（°）/显著波高（m） |
| `RIS_case_summary.csv` | `snapshotCenterXM / snapshotCenterYM / snapshotCenterZM` | 快照时刻质心世界坐标（m） |
| `RIS_case_summary.csv` | `snapshotRollDeg / snapshotPitchDeg / snapshotYawDeg` | 快照时刻欧拉角（°） |
| `RIS_case_summary.csv` | `seaMinM / seaMaxM` | 显示海面空间网格中的最小/最大原始高程（m） |
| `RIS_case_summary.csv` | `seaOutsideZFraction` | 快照海面网格高度超出 ±0.5 m 显示窗的格点比例 |
| `RIS_case_summary.csv` | `pathMinZM / pathMaxZM` | 整条质心轨迹最小/最大世界高度（m） |
| `RIS_case_summary.csv` | `maxTiltDeg` | 全过程最大板法向倾角（°） |
| `RIS_case_summary.csv` | `radiationFitRelativeRmsError` | 三个辐射模态拟合相对 RMS 误差的最大值 |
| `corrected_spectrum_audit.csv` | `windMps / HsM` | 风速（m/s）/默认离散显著波高（m） |
| `corrected_spectrum_audit.csv` | `retainedVariance` | 当前频带连续方差占完整连续谱方差的比例 |
| `corrected_spectrum_audit.csv` | `defaultVarianceRelativeError` | 默认离散方差/独立连续谱频带方差-1 |
| `corrected_spectrum_audit.csv` | `refinedVarianceRelativeError` | 400 点离散方差/独立连续谱频带方差-1 |
| `corrected_spectrum_audit.csv` | `sumFrequencyWeightsHz` | 全部梯形频率权重之和（Hz），应为 fMaxHz-fMinHz |
| `convergence.csv` | `photons` | 当前验证光子包数量 |
| `convergence.csv` | `meanPowerW` | 独立重复的板内无 EGG 功率均值（W） |
| `convergence.csv` | `powerStdW` | 独立重复接收功率的样本标准差，std 默认 N-1 归一化（W） |
| `convergence.csv` | `reportedPowerSEW` | 各次 MC 自报功率标准误的均值（W） |
| `motion.csv` | `row / column` | 3×5 面板位置，运动记录仅包含后四列 |
| `motion.csv` | `windMps / windDirectionDeg` | 风速（m/s）/传播风向（°） |
| `motion.csv` | `siteXM / siteYM` | 初始参考海浪站位世界坐标（m） |
| `motion.csv` | `timeS` | 相对运动时间（s） |
| `motion.csv` | `surgeM / swayM / heaveM` | eta 三个位移状态（m）；水平相对 c.ris.center，heave 相对平衡质心高度 |
| `motion.csv` | `rollRad / pitchRad / yawRad` | RIS 欧拉角（rad） |
| `motion.csv` | `localSeaM` | 运动质心水平位置的真实入射海面高程（m） |
| `motion.csv` | `txXM / txYM / txZM` | TX 世界坐标（m） |
| `motion.csv` | `txYawRad` | TX 平台绝对 yaw（rad） |
| `tests.csv` | `table(tests)` | 运行 MATLAB 测试结果对象转表的字段；当前缺少所引用测试文件，不推定额外参数 |


## 5. 逐源码、逐函数的完整具名符号表

每行符号保留源码名称；结构数组/元胞元素的索引仅在“取值/计算式”列体现，例如 `records{count}.windMps` 的字段路径简写为 `records.windMps`，不改字段名。相同符号在同一函数中多次赋值合并为一行并保留不同表达式；动态字段 `p.(fields{j})` 的合法字段全集见第 2.2 节，其他动态访问由对应 names/required/pairs 列表给出。为避免重复长结构遮挡参数，结构父行写 `struct(...)` 时其每个字段均在同节单列。

### 5.1 sea_wave/audit_sea_spectrum.m

**作用域：`audit_sea_spectrum`**

| 源码符号 | 含义与单位/维度 | 源码取值/计算式（各次赋值） | 角色 | 源码位置 |
| --- | --- | --- | --- | --- |
| `outputFolder` | 计算结果与图表输出目录 | `fullfile(folder,'audit')` | 输入/赋值 | [sea_wave/audit_sea_spectrum.m:1](<C:/Users/Spidey/Desktop/code/sea_wave/audit_sea_spectrum.m:1>) |
| `report` | 独立连续谱积分、离散积分和旧端点错误对比的审计结果 | `struct('cases',cases,'fMinHz',sea.f_min, 'formerU12VarianceRelativeError',formerError,'correctedU12VarianceRelativeError',fixedError)` | 输出/赋值 | [sea_wave/audit_sea_spectrum.m:1](<C:/Users/Spidey/Desktop/code/sea_wave/audit_sea_spectrum.m:1>) |
| `folder` | 当前源码所在目录 | `fileparts(mfilename('fullpath'))` | 赋值 | [sea_wave/audit_sea_spectrum.m:4](<C:/Users/Spidey/Desktop/code/sea_wave/audit_sea_spectrum.m:4>) |
| `sea` | 由 get_default_sea_params 读取的基础海面配置 | `get_default_sea_params()` | 赋值 | [sea_wave/audit_sea_spectrum.m:7](<C:/Users/Spidey/Desktop/code/sea_wave/audit_sea_spectrum.m:7>) |
| `p` | 当前作用域的参数结构 | `get_default_RIS_params()` | 赋值 | [sea_wave/audit_sea_spectrum.m:8](<C:/Users/Spidey/Desktop/code/sea_wave/audit_sea_spectrum.m:8>) |
| `winds` | 待验证风速 [3,5,9,12]（m/s） | `[3,5,9,12]` | 赋值 | [sea_wave/audit_sea_spectrum.m:9](<C:/Users/Spidey/Desktop/code/sea_wave/audit_sea_spectrum.m:9>) |
| `values` | 汇总表数值矩阵 | `zeros(4,6)`<br>`[winds(j),4*sqrt(variance),spectrum.retainedFullSpectrumVariance, variance/exact-1,denseVariance/exact-1,sum(spectrum.frequencyWeightsHz)]` | 赋值 | [sea_wave/audit_sea_spectrum.m:9](<C:/Users/Spidey/Desktop/code/sea_wave/audit_sea_spectrum.m:9>) |
| `j` | 当前循环、分量或样本索引（无量纲，MATLAB 从 1 开始） | `1:4` | 循环索引 | [sea_wave/audit_sea_spectrum.m:10](<C:/Users/Spidey/Desktop/code/sea_wave/audit_sea_spectrum.m:10>) |
| `p.windSpeed` | 风速（m/s） | `winds(j)`<br>`12` | 赋值 | [sea_wave/audit_sea_spectrum.m:11](<C:/Users/Spidey/Desktop/code/sea_wave/audit_sea_spectrum.m:11>) |
| `wave` | 海浪分量结构，含 omega、kx、ky、coefficient；具体重名按下方作用域解释 | `RIS_wave_components(p)` | 赋值 | [sea_wave/audit_sea_spectrum.m:12](<C:/Users/Spidey/Desktop/code/sea_wave/audit_sea_spectrum.m:12>) |
| `spectrum` | 海浪频谱及离散积分诊断结构 | `RIS_wave_components(p)` | 赋值 | [sea_wave/audit_sea_spectrum.m:12](<C:/Users/Spidey/Desktop/code/sea_wave/audit_sea_spectrum.m:12>) |
| `variance` | 当前默认离散海面方差（m²） | `sum(abs(wave.coefficient).^2)/2` | 赋值 | [sea_wave/audit_sea_spectrum.m:13](<C:/Users/Spidey/Desktop/code/sea_wave/audit_sea_spectrum.m:13>) |
| `exact` | 独立连续谱数值积分得到的频带方差（m²） | `bandVariance(p,sea.neumannConstant)` | 赋值 | [sea_wave/audit_sea_spectrum.m:14](<C:/Users/Spidey/Desktop/code/sea_wave/audit_sea_spectrum.m:14>) |
| `fine` | 将 frequencyBins 改为 400 的加密海浪参数 | `p` | 赋值 | [sea_wave/audit_sea_spectrum.m:15](<C:/Users/Spidey/Desktop/code/sea_wave/audit_sea_spectrum.m:15>) |
| `fine.frequencyBins` | 对数频率离散点数（不是点数减一的区间数） | `400` | 赋值 | [sea_wave/audit_sea_spectrum.m:15](<C:/Users/Spidey/Desktop/code/sea_wave/audit_sea_spectrum.m:15>) |
| `denseWave` | 400 个频率点生成的海浪分量结构 | `RIS_wave_components(fine)` | 赋值 | [sea_wave/audit_sea_spectrum.m:16](<C:/Users/Spidey/Desktop/code/sea_wave/audit_sea_spectrum.m:16>) |
| `denseVariance` | 加密离散海面方差（m²） | `sum(abs(denseWave.coefficient).^2)/2` | 赋值 | [sea_wave/audit_sea_spectrum.m:17](<C:/Users/Spidey/Desktop/code/sea_wave/audit_sea_spectrum.m:17>) |
| `p.fMinHz` | 海浪最低频率（Hz） | `0.1` | 赋值 | [sea_wave/audit_sea_spectrum.m:24](<C:/Users/Spidey/Desktop/code/sea_wave/audit_sea_spectrum.m:24>) |
| `fixedError` | U=12、fMinHz=0.1 工况修正后方差相对误差 | `sum(abs(wave.coefficient).^2)/(2*exact)-1` | 赋值 | [sea_wave/audit_sea_spectrum.m:27](<C:/Users/Spidey/Desktop/code/sea_wave/audit_sea_spectrum.m:27>) |
| `f` | 去重的海浪离散频率（Hz） | `unique(wave.omega)/(2*pi)` | 赋值 | [sea_wave/audit_sea_spectrum.m:29](<C:/Users/Spidey/Desktop/code/sea_wave/audit_sea_spectrum.m:29>) |
| `S` | U=12 对比算例中的 Neumann 角频率谱值 | `sea.neumannConstant*pi/4*(2*pi*f).^(-6).*exp(-2*p.g^2./(p.windSpeed*2*pi*f).^2)` | 赋值 | [sea_wave/audit_sea_spectrum.m:30](<C:/Users/Spidey/Desktop/code/sea_wave/audit_sea_spectrum.m:30>) |
| `extraVariance` | 旧实现错误加倍频率两端点引入的额外方差（m²） | `2*pi*sum(S([1,end]).*spectrum.frequencyWeightsHz([1,end]))` | 赋值 | [sea_wave/audit_sea_spectrum.m:31](<C:/Users/Spidey/Desktop/code/sea_wave/audit_sea_spectrum.m:31>) |
| `formerError` | 重建旧错误后的方差相对误差 | `fixedError+extraVariance/exact` | 赋值 | [sea_wave/audit_sea_spectrum.m:32](<C:/Users/Spidey/Desktop/code/sea_wave/audit_sea_spectrum.m:32>) |
| `cases` | 四风速离散谱与连续谱的审计结果表 | `array2table(values,'VariableNames',{'windMps','HsM','retainedVariance', 'defaultVarianceRelativeError','refinedVarianceRelativeError','sumFrequencyWeightsHz'})` | 赋值 | [sea_wave/audit_sea_spectrum.m:34](<C:/Users/Spidey/Desktop/code/sea_wave/audit_sea_spectrum.m:34>) |
| `report.cases` | 连续谱积分与离散谱方差核验表 | `cases` | 结构字段 | [sea_wave/audit_sea_spectrum.m:36](<C:/Users/Spidey/Desktop/code/sea_wave/audit_sea_spectrum.m:36>) |
| `report.fMinHz` | 海浪最低频率（Hz） | `sea.f_min` | 结构字段 | [sea_wave/audit_sea_spectrum.m:36](<C:/Users/Spidey/Desktop/code/sea_wave/audit_sea_spectrum.m:36>) |
| `report.formerU12VarianceRelativeError` | 审计中重现旧频率端点加倍错误的 U=12 方差相对误差 | `formerError` | 结构字段 | [sea_wave/audit_sea_spectrum.m:37](<C:/Users/Spidey/Desktop/code/sea_wave/audit_sea_spectrum.m:37>) |
| `report.correctedU12VarianceRelativeError` | 修复端点半权重后的 U=12 方差相对误差 | `fixedError` | 结构字段 | [sea_wave/audit_sea_spectrum.m:37](<C:/Users/Spidey/Desktop/code/sea_wave/audit_sea_spectrum.m:37>) |


**作用域：`bandVariance`**

| 源码符号 | 含义与单位/维度 | 源码取值/计算式（各次赋值） | 角色 | 源码位置 |
| --- | --- | --- | --- | --- |
| `p` | 当前作用域的参数结构 | 函数输入/输出；值随调用或上游结果确定 | 输入 | [sea_wave/audit_sea_spectrum.m:44](<C:/Users/Spidey/Desktop/code/sea_wave/audit_sea_spectrum.m:44>) |
| `C` | Neumann 频谱经验常数，调用时为 sea.neumannConstant | 函数输入/输出；值随调用或上游结果确定 | 输入 | [sea_wave/audit_sea_spectrum.m:44](<C:/Users/Spidey/Desktop/code/sea_wave/audit_sea_spectrum.m:44>) |
| `variance` | 指定频带连续角频率谱积分得到的海面方差（m²） | `integral(@(w) C*pi/4*w.^(-6).*exp(-2*p.g^2./(p.windSpeed*w).^2), 2*pi*p.fMinHz,2*pi*p.fMaxHz,'AbsTol',1e-13,'RelTol',1e-11)` | 输出/赋值 | [sea_wave/audit_sea_spectrum.m:44](<C:/Users/Spidey/Desktop/code/sea_wave/audit_sea_spectrum.m:44>) |
| `w` | 连续角频率积分自变量（rad/s）；此处不是波浪结构 | `variance = integral(@(w) C*pi/4*w.^(-6).*exp(-2*p.g^2./(p.windSpeed*w).^2), 2*pi*p.fMinHz,2*pi*p.fMaxHz,'AbsTol',1e-13,'RelTol',1e-11)` | 匿名函数输入 | [sea_wave/audit_sea_spectrum.m:46](<C:/Users/Spidey/Desktop/code/sea_wave/audit_sea_spectrum.m:46>) |


### 5.2 sea_wave/get_default_RIS_params.m

**作用域：`get_default_RIS_params`**

| 源码符号 | 含义与单位/维度 | 源码取值/计算式（各次赋值） | 角色 | 源码位置 |
| --- | --- | --- | --- | --- |
| `p` | 当前作用域的参数结构 | `struct(...)（所有字段在本节逐行列明）` | 输出/赋值 | [sea_wave/get_default_RIS_params.m:1](<C:/Users/Spidey/Desktop/code/sea_wave/get_default_RIS_params.m:1>) |
| `sea` | 由 get_default_sea_params 读取的基础海面配置 | `get_default_sea_params()` | 赋值 | [sea_wave/get_default_RIS_params.m:3](<C:/Users/Spidey/Desktop/code/sea_wave/get_default_RIS_params.m:3>) |
| `p.Lx` | RIS 体坐标 x 方向边长（m） | `1` | 结构字段 | [sea_wave/get_default_RIS_params.m:4](<C:/Users/Spidey/Desktop/code/sea_wave/get_default_RIS_params.m:4>) |
| `p.Ly` | RIS 体坐标 y 方向边长（m） | `1` | 结构字段 | [sea_wave/get_default_RIS_params.m:4](<C:/Users/Spidey/Desktop/code/sea_wave/get_default_RIS_params.m:4>) |
| `p.mass` | RIS 质量（kg）；模型内 mass 为含附加质量的六分量广义质量 | `12` | 结构字段 | [sea_wave/get_default_RIS_params.m:4](<C:/Users/Spidey/Desktop/code/sea_wave/get_default_RIS_params.m:4>) |
| `p.thickness` | RIS 板厚（m） | `0.03` | 结构字段 | [sea_wave/get_default_RIS_params.m:4](<C:/Users/Spidey/Desktop/code/sea_wave/get_default_RIS_params.m:4>) |
| `p.rho` | 海水密度（kg/m³） | `1025` | 结构字段 | [sea_wave/get_default_RIS_params.m:4](<C:/Users/Spidey/Desktop/code/sea_wave/get_default_RIS_params.m:4>) |
| `p.g` | 重力加速度（m/s²） | `sea.g` | 结构字段 | [sea_wave/get_default_RIS_params.m:4](<C:/Users/Spidey/Desktop/code/sea_wave/get_default_RIS_params.m:4>) |
| `p.windSpeed` | 风速（m/s） | `sea.U` | 结构字段 | [sea_wave/get_default_RIS_params.m:5](<C:/Users/Spidey/Desktop/code/sea_wave/get_default_RIS_params.m:5>) |
| `p.windAngleDeg` | 波浪传播方向角（°），从 +Y 朝 +X 增大 | `30` | 结构字段 | [sea_wave/get_default_RIS_params.m:5](<C:/Users/Spidey/Desktop/code/sea_wave/get_default_RIS_params.m:5>) |
| `p.duration` | 保存的运动持续时间（s） | `10` | 结构字段 | [sea_wave/get_default_RIS_params.m:5](<C:/Users/Spidey/Desktop/code/sea_wave/get_default_RIS_params.m:5>) |
| `p.dt` | 运动输出采样间隔（s），不等于实际 RK4 内部步长 | `0.025` | 结构字段 | [sea_wave/get_default_RIS_params.m:5](<C:/Users/Spidey/Desktop/code/sea_wave/get_default_RIS_params.m:5>) |
| `p.seed` | 随机数种子；所属结构决定是海浪、TX、MC 或湍流随机流 | `sea.seed` | 结构字段 | [sea_wave/get_default_RIS_params.m:5](<C:/Users/Spidey/Desktop/code/sea_wave/get_default_RIS_params.m:5>) |
| `p.startTime` | 绝对海浪起算时刻（s） | `sea.t0` | 结构字段 | [sea_wave/get_default_RIS_params.m:6](<C:/Users/Spidey/Desktop/code/sea_wave/get_default_RIS_params.m:6>) |
| `p.center` | RIS 的水平参考中心 [x;y]（m），或光斑加权中心 [x,y]（m）；以父结构为准 | `[0;0]` | 结构字段 | [sea_wave/get_default_RIS_params.m:6](<C:/Users/Spidey/Desktop/code/sea_wave/get_default_RIS_params.m:6>) |
| `p.fMinHz` | 海浪最低频率（Hz） | `sea.f_min` | 结构字段 | [sea_wave/get_default_RIS_params.m:6](<C:/Users/Spidey/Desktop/code/sea_wave/get_default_RIS_params.m:6>) |
| `p.fMaxHz` | 海浪最高频率（Hz） | `sea.f_max` | 结构字段 | [sea_wave/get_default_RIS_params.m:6](<C:/Users/Spidey/Desktop/code/sea_wave/get_default_RIS_params.m:6>) |
| `p.frequencyBins` | 对数频率离散点数（不是点数减一的区间数） | `sea.Fre_segment` | 结构字段 | [sea_wave/get_default_RIS_params.m:7](<C:/Users/Spidey/Desktop/code/sea_wave/get_default_RIS_params.m:7>) |
| `p.directionBins` | 方向离散点数，包含 ±π/2 端点 | `sea.Dir_segment` | 结构字段 | [sea_wave/get_default_RIS_params.m:7](<C:/Users/Spidey/Desktop/code/sea_wave/get_default_RIS_params.m:7>) |
| `p.panelCount` | 薄板 BEM 每边面元数，须为 4..80 的偶数；完整网格 panelCount² | `24` | 结构字段 | [sea_wave/get_default_RIS_params.m:7](<C:/Users/Spidey/Desktop/code/sea_wave/get_default_RIS_params.m:7>) |
| `p.horizontalInitialPosition` | 相对 center 的初始水平位移 [x;y]（m） | `[0;0]` | 结构字段 | [sea_wave/get_default_RIS_params.m:8](<C:/Users/Spidey/Desktop/code/sea_wave/get_default_RIS_params.m:8>) |
| `p.horizontalInitialVelocity` | 可选世界系初始水平速度 [vx;vy]（m/s），[] 保留初始化结果 | `[]` | 结构字段 | [sea_wave/get_default_RIS_params.m:8](<C:/Users/Spidey/Desktop/code/sea_wave/get_default_RIS_params.m:8>) |
| `p.initialPosition` | 可选初始 [surge;sway;heave] 状态（m）；水平分量仍相对 center，[] 不覆盖 | `[]` | 结构字段 | [sea_wave/get_default_RIS_params.m:9](<C:/Users/Spidey/Desktop/code/sea_wave/get_default_RIS_params.m:9>) |
| `p.initialEulerAngles` | 可选初始 [roll;pitch;yaw]（rad），[] 不覆盖；默认初始 yaw=0 | `[]` | 结构字段 | [sea_wave/get_default_RIS_params.m:9](<C:/Users/Spidey/Desktop/code/sea_wave/get_default_RIS_params.m:9>) |
| `p.initialBodyVelocity` | 可选体坐标初始 [u;v;w;p;q;r]（前 3 项 m/s，后 3 项 rad/s） | `[]` | 结构字段 | [sea_wave/get_default_RIS_params.m:9](<C:/Users/Spidey/Desktop/code/sea_wave/get_default_RIS_params.m:9>) |
| `p.initialization` | 初始化方式：waveSteady（频域稳态）或 rest（静止状态） | `'waveSteady'` | 结构字段 | [sea_wave/get_default_RIS_params.m:10](<C:/Users/Spidey/Desktop/code/sea_wave/get_default_RIS_params.m:10>) |
| `p.radiationModes` | 被动辐射记忆有理拟合的振子模态数 | `64` | 结构字段 | [sea_wave/get_default_RIS_params.m:10](<C:/Users/Spidey/Desktop/code/sea_wave/get_default_RIS_params.m:10>) |
| `p.radiationFrequencyBins` | 辐射拟合额外对数频率采样点数 | `80` | 结构字段 | [sea_wave/get_default_RIS_params.m:10](<C:/Users/Spidey/Desktop/code/sea_wave/get_default_RIS_params.m:10>) |
| `p.integrationStep` | RK4 期望最大内步长（s），还按稳定性自动细分 | `0.01` | 结构字段 | [sea_wave/get_default_RIS_params.m:11](<C:/Users/Spidey/Desktop/code/sea_wave/get_default_RIS_params.m:11>) |
| `p.wave` | 外部共享的海浪分量结构；[] 时按参数生成 | `[]` | 结构字段 | [sea_wave/get_default_RIS_params.m:11](<C:/Users/Spidey/Desktop/code/sea_wave/get_default_RIS_params.m:11>) |
| `p.plot` | 是否自动绘制 RIS 单次运动；RIS_float 无输出且未覆盖时设为 true | `false` | 结构字段 | [sea_wave/get_default_RIS_params.m:11](<C:/Users/Spidey/Desktop/code/sea_wave/get_default_RIS_params.m:11>) |


### 5.3 sea_wave/get_default_sea_params.m

**作用域：`get_default_sea_params`**

| 源码符号 | 含义与单位/维度 | 源码取值/计算式（各次赋值） | 角色 | 源码位置 |
| --- | --- | --- | --- | --- |
| `sea_params` | 海面配置结构 | 函数输入/输出；值随调用或上游结果确定 | 输出 | [sea_wave/get_default_sea_params.m:1](<C:/Users/Spidey/Desktop/code/sea_wave/get_default_sea_params.m:1>) |
| `sea_params.g` | 重力加速度（m/s²） | `9.81` | 赋值 | [sea_wave/get_default_sea_params.m:7](<C:/Users/Spidey/Desktop/code/sea_wave/get_default_sea_params.m:7>) |
| `sea_params.neumannConstant` | 源码 Neumann 频谱经验常数，默认 3.05；维度随源码谱表达式确定 | `3.05` | 赋值 | [sea_wave/get_default_sea_params.m:8](<C:/Users/Spidey/Desktop/code/sea_wave/get_default_sea_params.m:8>) |
| `g` | 局部重力加速度副本（m/s²） | `sea_params.g` | 赋值 | [sea_wave/get_default_sea_params.m:9](<C:/Users/Spidey/Desktop/code/sea_wave/get_default_sea_params.m:9>) |
| `sea_params.x_range` | 海面世界 x/y 轴显示与重建范围 [min,max]（m） | `[-5, 5]` | 赋值 | [sea_wave/get_default_sea_params.m:12](<C:/Users/Spidey/Desktop/code/sea_wave/get_default_sea_params.m:12>) |
| `sea_params.y_range` | 海面世界 x/y 轴显示与重建范围 [min,max]（m） | `[-5, 5]` | 赋值 | [sea_wave/get_default_sea_params.m:13](<C:/Users/Spidey/Desktop/code/sea_wave/get_default_sea_params.m:13>) |
| `sea_params.f_max` | 海浪最高频率（Hz） | `4.0` | 赋值 | [sea_wave/get_default_sea_params.m:16](<C:/Users/Spidey/Desktop/code/sea_wave/get_default_sea_params.m:16>) |
| `dx_required` | 最高频率对应的空间 Nyquist 最大网格间距 g/(4*pi*f_max²)（m） | `g / (4 * pi * sea_params.f_max^2)` | 赋值 | [sea_wave/get_default_sea_params.m:21](<C:/Users/Spidey/Desktop/code/sea_wave/get_default_sea_params.m:21>) |
| `L_x` | 海面 x/y 范围跨度（m） | `sea_params.x_range(2) - sea_params.x_range(1)` | 赋值 | [sea_wave/get_default_sea_params.m:24](<C:/Users/Spidey/Desktop/code/sea_wave/get_default_sea_params.m:24>) |
| `L_y` | 海面 x/y 范围跨度（m） | `sea_params.y_range(2) - sea_params.y_range(1)` | 赋值 | [sea_wave/get_default_sea_params.m:25](<C:/Users/Spidey/Desktop/code/sea_wave/get_default_sea_params.m:25>) |
| `L_max` | 两方向跨度的较大值（m），用于确定统一网格点数 | `max(L_x, L_y)` | 赋值 | [sea_wave/get_default_sea_params.m:26](<C:/Users/Spidey/Desktop/code/sea_wave/get_default_sea_params.m:26>) |
| `sea_params.grid_density` | 海面每个空间方向网格点数，由最高频率和最大跨度反推 | `ceil(L_max / dx_required) + 1` | 赋值 | [sea_wave/get_default_sea_params.m:29](<C:/Users/Spidey/Desktop/code/sea_wave/get_default_sea_params.m:29>) |
| `sea_params.f_min` | 海浪最低频率（Hz） | `0.05` | 赋值 | [sea_wave/get_default_sea_params.m:32](<C:/Users/Spidey/Desktop/code/sea_wave/get_default_sea_params.m:32>) |
| `sea_params.Fre_segment` | 对数频率离散点数（不是点数减一的区间数） | `50` | 赋值 | [sea_wave/get_default_sea_params.m:33](<C:/Users/Spidey/Desktop/code/sea_wave/get_default_sea_params.m:33>) |
| `sea_params.Dir_segment` | 方向离散点数，包含 ±π/2 端点 | `36` | 赋值 | [sea_wave/get_default_sea_params.m:34](<C:/Users/Spidey/Desktop/code/sea_wave/get_default_sea_params.m:34>) |
| `sea_params.U` | 风速（m/s） | `10` | 赋值 | [sea_wave/get_default_sea_params.m:35](<C:/Users/Spidey/Desktop/code/sea_wave/get_default_sea_params.m:35>) |
| `sea_params.t0` | 绝对海浪起算时刻（s） | `0` | 赋值 | [sea_wave/get_default_sea_params.m:36](<C:/Users/Spidey/Desktop/code/sea_wave/get_default_sea_params.m:36>) |
| `sea_params.seed` | 随机数种子；所属结构决定是海浪、TX、MC 或湍流随机流 | `42` | 赋值 | [sea_wave/get_default_sea_params.m:37](<C:/Users/Spidey/Desktop/code/sea_wave/get_default_sea_params.m:37>) |


### 5.4 sea_wave/main3_windseasurface_0326.m

**作用域：`script`**

| 源码符号 | 含义与单位/维度 | 源码取值/计算式（各次赋值） | 角色 | 源码位置 |
| --- | --- | --- | --- | --- |
| `params` | 海面配置结构 | `get_default_sea_params()` | 赋值 | [sea_wave/main3_windseasurface_0326.m:3](<C:/Users/Spidey/Desktop/code/sea_wave/main3_windseasurface_0326.m:3>) |
| `params.t0` | 绝对海浪起算时刻（s） | `27` | 赋值 | [sea_wave/main3_windseasurface_0326.m:4](<C:/Users/Spidey/Desktop/code/sea_wave/main3_windseasurface_0326.m:4>) |
| `params.seed` | 随机数种子；所属结构决定是海浪、TX、MC 或湍流随机流 | `600` | 赋值 | [sea_wave/main3_windseasurface_0326.m:5](<C:/Users/Spidey/Desktop/code/sea_wave/main3_windseasurface_0326.m:5>) |
| `exp_cases` | 纯海面入口工况 [9,30;9,45;12,30]，每行为 [风速 m/s,风向 °] | `[9,30;9,45;12,30]` | 赋值 | [sea_wave/main3_windseasurface_0326.m:6](<C:/Users/Spidey/Desktop/code/sea_wave/main3_windseasurface_0326.m:6>) |
| `output_dir` | 计算结果与图表输出目录 | `fullfile(fileparts(mfilename('fullpath')),'FigOutput_SeaSurface')` | 赋值 | [sea_wave/main3_windseasurface_0326.m:7](<C:/Users/Spidey/Desktop/code/sea_wave/main3_windseasurface_0326.m:7>) |
| `sim_data` | 各纯海面工况的 X/Y/Z/wave 数据元胞数组 | `cell(1,size(exp_cases,1))`<br>`struct('X',X,'Y',Y,'Z',Z,'wave',wave)` | 赋值 | [sea_wave/main3_windseasurface_0326.m:9](<C:/Users/Spidey/Desktop/code/sea_wave/main3_windseasurface_0326.m:9>) |
| `i` | 当前循环、分量或样本索引（无量纲，MATLAB 从 1 开始） | `1:size(exp_cases,1)` | 循环索引 | [sea_wave/main3_windseasurface_0326.m:10](<C:/Users/Spidey/Desktop/code/sea_wave/main3_windseasurface_0326.m:10>) |
| `p` | 当前作用域的参数结构 | `struct(...)（所有字段在本节逐行列明）` | 赋值 | [sea_wave/main3_windseasurface_0326.m:11](<C:/Users/Spidey/Desktop/code/sea_wave/main3_windseasurface_0326.m:11>) |
| `p.g` | 重力加速度（m/s²） | `params.g` | 结构字段 | [sea_wave/main3_windseasurface_0326.m:11](<C:/Users/Spidey/Desktop/code/sea_wave/main3_windseasurface_0326.m:11>) |
| `p.fMinHz` | 海浪最低频率（Hz） | `params.f_min` | 结构字段 | [sea_wave/main3_windseasurface_0326.m:11](<C:/Users/Spidey/Desktop/code/sea_wave/main3_windseasurface_0326.m:11>) |
| `p.fMaxHz` | 海浪最高频率（Hz） | `params.f_max` | 结构字段 | [sea_wave/main3_windseasurface_0326.m:11](<C:/Users/Spidey/Desktop/code/sea_wave/main3_windseasurface_0326.m:11>) |
| `p.frequencyBins` | 对数频率离散点数（不是点数减一的区间数） | `params.Fre_segment` | 结构字段 | [sea_wave/main3_windseasurface_0326.m:12](<C:/Users/Spidey/Desktop/code/sea_wave/main3_windseasurface_0326.m:12>) |
| `p.directionBins` | 方向离散点数，包含 ±π/2 端点 | `params.Dir_segment` | 结构字段 | [sea_wave/main3_windseasurface_0326.m:12](<C:/Users/Spidey/Desktop/code/sea_wave/main3_windseasurface_0326.m:12>) |
| `p.windSpeed` | 风速（m/s） | `exp_cases(i,1)` | 结构字段 | [sea_wave/main3_windseasurface_0326.m:13](<C:/Users/Spidey/Desktop/code/sea_wave/main3_windseasurface_0326.m:13>) |
| `p.windAngleDeg` | 波浪传播方向角（°），从 +Y 朝 +X 增大 | `exp_cases(i,2)` | 结构字段 | [sea_wave/main3_windseasurface_0326.m:13](<C:/Users/Spidey/Desktop/code/sea_wave/main3_windseasurface_0326.m:13>) |
| `p.seed` | 随机数种子；所属结构决定是海浪、TX、MC 或湍流随机流 | `params.seed` | 结构字段 | [sea_wave/main3_windseasurface_0326.m:13](<C:/Users/Spidey/Desktop/code/sea_wave/main3_windseasurface_0326.m:13>) |
| `wave` | 海浪分量结构，含 omega、kx、ky、coefficient；具体重名按下方作用域解释 | `RIS_wave_components(p)` | 赋值 | [sea_wave/main3_windseasurface_0326.m:14](<C:/Users/Spidey/Desktop/code/sea_wave/main3_windseasurface_0326.m:14>) |
| `field` | 提供给海面重建函数的 wave/parameters 封装，或当前块复波场 | `struct('wave',wave,'parameters',struct('startTime',params.t0))` | 赋值 | [sea_wave/main3_windseasurface_0326.m:15](<C:/Users/Spidey/Desktop/code/sea_wave/main3_windseasurface_0326.m:15>) |
| `field.wave` | 外部共享的海浪分量结构；[] 时按参数生成 | `wave` | 结构字段 | [sea_wave/main3_windseasurface_0326.m:15](<C:/Users/Spidey/Desktop/code/sea_wave/main3_windseasurface_0326.m:15>) |
| `field.parameters` | 当前计算实际使用的参数结构；EGG 中为混合分布参数 | `struct('startTime',params.t0)` | 结构字段 | [sea_wave/main3_windseasurface_0326.m:15](<C:/Users/Spidey/Desktop/code/sea_wave/main3_windseasurface_0326.m:15>) |
| `field.parameters.startTime` | 绝对海浪起算时刻（s） | `params.t0` | 结构字段 | [sea_wave/main3_windseasurface_0326.m:15](<C:/Users/Spidey/Desktop/code/sea_wave/main3_windseasurface_0326.m:15>) |
| `X` | 海面网格的世界 x 坐标矩阵（m） | `RIS_sea_surface(field, linspace(params.x_range(1),params.x_range(2),params.grid_density), linspace(params.y_range(1),params.y_range(2),params.grid_density),0)` | 赋值 | [sea_wave/main3_windseasurface_0326.m:16](<C:/Users/Spidey/Desktop/code/sea_wave/main3_windseasurface_0326.m:16>) |
| `Y` | 海面网格的世界 y 坐标矩阵（m） | `RIS_sea_surface(field, linspace(params.x_range(1),params.x_range(2),params.grid_density), linspace(params.y_range(1),params.y_range(2),params.grid_density),0)` | 赋值 | [sea_wave/main3_windseasurface_0326.m:16](<C:/Users/Spidey/Desktop/code/sea_wave/main3_windseasurface_0326.m:16>) |
| `Z` | 海面高度矩阵（m） | `RIS_sea_surface(field, linspace(params.x_range(1),params.x_range(2),params.grid_density), linspace(params.y_range(1),params.y_range(2),params.grid_density),0)` | 赋值 | [sea_wave/main3_windseasurface_0326.m:16](<C:/Users/Spidey/Desktop/code/sea_wave/main3_windseasurface_0326.m:16>) |
| `sim_data.X` | 海面世界 x/y 坐标与高程矩阵（m），按字段区分 | `X` | 结构字段 | [sea_wave/main3_windseasurface_0326.m:19](<C:/Users/Spidey/Desktop/code/sea_wave/main3_windseasurface_0326.m:19>) |
| `sim_data.Y` | 海面世界 x/y 坐标与高程矩阵（m），按字段区分 | `Y` | 结构字段 | [sea_wave/main3_windseasurface_0326.m:19](<C:/Users/Spidey/Desktop/code/sea_wave/main3_windseasurface_0326.m:19>) |
| `sim_data.Z` | 海面世界 x/y 坐标与高程矩阵（m），按字段区分 | `Z` | 结构字段 | [sea_wave/main3_windseasurface_0326.m:19](<C:/Users/Spidey/Desktop/code/sea_wave/main3_windseasurface_0326.m:19>) |
| `sim_data.wave` | 外部共享的海浪分量结构；[] 时按参数生成 | `wave` | 结构字段 | [sea_wave/main3_windseasurface_0326.m:19](<C:/Users/Spidey/Desktop/code/sea_wave/main3_windseasurface_0326.m:19>) |
| `fig` | 图窗句柄；频率等重名按所在函数解释 | `figure('Color','w','Visible','off','Position',[50,50,880,720])` | 赋值 | [sea_wave/main3_windseasurface_0326.m:20](<C:/Users/Spidey/Desktop/code/sea_wave/main3_windseasurface_0326.m:20>) |
| `ax` | 绘图坐标轴句柄 | `nexttile(tiledlayout(fig,1,1,'Padding','compact'))` | 赋值 | [sea_wave/main3_windseasurface_0326.m:21](<C:/Users/Spidey/Desktop/code/sea_wave/main3_windseasurface_0326.m:21>) |
| `cb` | 颜色条句柄 | `colorbar(ax)` | 赋值 | [sea_wave/main3_windseasurface_0326.m:26](<C:/Users/Spidey/Desktop/code/sea_wave/main3_windseasurface_0326.m:26>) |
| `cb.Label.String` | 颜色条标题文本 | `'Elevation (m)'` | 赋值 | [sea_wave/main3_windseasurface_0326.m:26](<C:/Users/Spidey/Desktop/code/sea_wave/main3_windseasurface_0326.m:26>) |
| `name` | 当前输出文件的基础名称 | `sprintf('Surf_U%02d_Phi%02d',exp_cases(i,:))` | 赋值 | [sea_wave/main3_windseasurface_0326.m:29](<C:/Users/Spidey/Desktop/code/sea_wave/main3_windseasurface_0326.m:29>) |


### 5.5 sea_wave/oceanSurfaceModel_phiwind_260315.m

**作用域：`oceanSurfaceModel_phiwind_260315`**

| 源码符号 | 含义与单位/维度 | 源码取值/计算式（各次赋值） | 角色 | 源码位置 |
| --- | --- | --- | --- | --- |
| `f_min` | 海浪最低/最高频率（Hz），必填输入 | 函数输入/输出；值随调用或上游结果确定 | 输入 | [sea_wave/oceanSurfaceModel_phiwind_260315.m:1](<C:/Users/Spidey/Desktop/code/sea_wave/oceanSurfaceModel_phiwind_260315.m:1>) |
| `f_max` | 海浪最低/最高频率（Hz），必填输入 | 函数输入/输出；值随调用或上游结果确定 | 输入 | [sea_wave/oceanSurfaceModel_phiwind_260315.m:1](<C:/Users/Spidey/Desktop/code/sea_wave/oceanSurfaceModel_phiwind_260315.m:1>) |
| `Fre_segment` | 频率/方向离散点数，必填输入 | 函数输入/输出；值随调用或上游结果确定 | 输入 | [sea_wave/oceanSurfaceModel_phiwind_260315.m:1](<C:/Users/Spidey/Desktop/code/sea_wave/oceanSurfaceModel_phiwind_260315.m:1>) |
| `Dir_segment` | 频率/方向离散点数，必填输入 | 函数输入/输出；值随调用或上游结果确定 | 输入 | [sea_wave/oceanSurfaceModel_phiwind_260315.m:1](<C:/Users/Spidey/Desktop/code/sea_wave/oceanSurfaceModel_phiwind_260315.m:1>) |
| `U` | 风速（m/s），必填输入 | 函数输入/输出；值随调用或上游结果确定 | 输入 | [sea_wave/oceanSurfaceModel_phiwind_260315.m:1](<C:/Users/Spidey/Desktop/code/sea_wave/oceanSurfaceModel_phiwind_260315.m:1>) |
| `phi_wind` | 从 +Y 朝 +X 的传播风向角（rad），接口中转成 windAngleDeg | 函数输入/输出；值随调用或上游结果确定 | 输入 | [sea_wave/oceanSurfaceModel_phiwind_260315.m:1](<C:/Users/Spidey/Desktop/code/sea_wave/oceanSurfaceModel_phiwind_260315.m:1>) |
| `x_range` | 世界 x/y 空间范围 [min,max]（m），必填输入 | 函数输入/输出；值随调用或上游结果确定 | 输入 | [sea_wave/oceanSurfaceModel_phiwind_260315.m:1](<C:/Users/Spidey/Desktop/code/sea_wave/oceanSurfaceModel_phiwind_260315.m:1>) |
| `y_range` | 世界 x/y 空间范围 [min,max]（m），必填输入 | 函数输入/输出；值随调用或上游结果确定 | 输入 | [sea_wave/oceanSurfaceModel_phiwind_260315.m:1](<C:/Users/Spidey/Desktop/code/sea_wave/oceanSurfaceModel_phiwind_260315.m:1>) |
| `grid_density` | 每个空间方向的网格点数，必填输入 | 函数输入/输出；值随调用或上游结果确定 | 输入 | [sea_wave/oceanSurfaceModel_phiwind_260315.m:1](<C:/Users/Spidey/Desktop/code/sea_wave/oceanSurfaceModel_phiwind_260315.m:1>) |
| `t0` | 本次静态海面的绝对时间（s），必填输入 | 函数输入/输出；值随调用或上游结果确定 | 输入 | [sea_wave/oceanSurfaceModel_phiwind_260315.m:1](<C:/Users/Spidey/Desktop/code/sea_wave/oceanSurfaceModel_phiwind_260315.m:1>) |
| `seed` | 随机数种子（非负整数） | `1` | 输入/赋值 | [sea_wave/oceanSurfaceModel_phiwind_260315.m:1](<C:/Users/Spidey/Desktop/code/sea_wave/oceanSurfaceModel_phiwind_260315.m:1>) |
| `output_mode` | 输出模式 'all'（默认）或 'single' | `'all'` | 输入/赋值 | [sea_wave/oceanSurfaceModel_phiwind_260315.m:1](<C:/Users/Spidey/Desktop/code/sea_wave/oceanSurfaceModel_phiwind_260315.m:1>) |
| `plane_index` | 单平面模式选择的 [i,j] 或线性索引，默认 [1,1] | `[1, 1]` | 输入/赋值 | [sea_wave/oceanSurfaceModel_phiwind_260315.m:1](<C:/Users/Spidey/Desktop/code/sea_wave/oceanSurfaceModel_phiwind_260315.m:1>) |
| `XX` | 海面网格的世界 x 坐标矩阵（m） | `meshgrid(linspace(x_range(1),x_range(2),grid_density), linspace(y_range(1),y_range(2),grid_density))` | 输出/赋值 | [sea_wave/oceanSurfaceModel_phiwind_260315.m:1](<C:/Users/Spidey/Desktop/code/sea_wave/oceanSurfaceModel_phiwind_260315.m:1>) |
| `YY` | 海面网格的世界 y 坐标矩阵（m） | `meshgrid(linspace(x_range(1),x_range(2),grid_density), linspace(y_range(1),y_range(2),grid_density))` | 输出/赋值 | [sea_wave/oceanSurfaceModel_phiwind_260315.m:1](<C:/Users/Spidey/Desktop/code/sea_wave/oceanSurfaceModel_phiwind_260315.m:1>) |
| `Z_surf` | 海面高度矩阵（m） | `reshape(RIS_wave_field(field,[XX(:),YY(:)],0),size(XX))` | 输出/赋值 | [sea_wave/oceanSurfaceModel_phiwind_260315.m:1](<C:/Users/Spidey/Desktop/code/sea_wave/oceanSurfaceModel_phiwind_260315.m:1>) |
| `plane_data` | 所有网格单元拟合小平面的方程/法向/中心/误差/顶点结构 | `calculatePlaneParameters(XX,YY,Z_surf)` | 输出/赋值 | [sea_wave/oceanSurfaceModel_phiwind_260315.m:1](<C:/Users/Spidey/Desktop/code/sea_wave/oceanSurfaceModel_phiwind_260315.m:1>) |
| `selected_plane_data` | single 模式所选小平面结构；all 模式为 [] | `[]`<br>`getSelectedPlaneData(plane_data,plane_index)` | 输出/赋值 | [sea_wave/oceanSurfaceModel_phiwind_260315.m:1](<C:/Users/Spidey/Desktop/code/sea_wave/oceanSurfaceModel_phiwind_260315.m:1>) |
| `wave` | 海浪分量结构，含 omega、kx、ky、coefficient；具体重名按下方作用域解释 | `RIS_wave_components(p)` | 输出/赋值 | [sea_wave/oceanSurfaceModel_phiwind_260315.m:1](<C:/Users/Spidey/Desktop/code/sea_wave/oceanSurfaceModel_phiwind_260315.m:1>) |
| `sea` | 由 get_default_sea_params 读取的基础海面配置 | `get_default_sea_params()` | 赋值 | [sea_wave/oceanSurfaceModel_phiwind_260315.m:40](<C:/Users/Spidey/Desktop/code/sea_wave/oceanSurfaceModel_phiwind_260315.m:40>) |
| `p` | 当前作用域的参数结构 | `struct(...)（所有字段在本节逐行列明）` | 赋值 | [sea_wave/oceanSurfaceModel_phiwind_260315.m:41](<C:/Users/Spidey/Desktop/code/sea_wave/oceanSurfaceModel_phiwind_260315.m:41>) |
| `p.g` | 重力加速度（m/s²） | `sea.g` | 结构字段 | [sea_wave/oceanSurfaceModel_phiwind_260315.m:41](<C:/Users/Spidey/Desktop/code/sea_wave/oceanSurfaceModel_phiwind_260315.m:41>) |
| `p.fMinHz` | 海浪最低频率（Hz） | `f_min` | 结构字段 | [sea_wave/oceanSurfaceModel_phiwind_260315.m:41](<C:/Users/Spidey/Desktop/code/sea_wave/oceanSurfaceModel_phiwind_260315.m:41>) |
| `p.fMaxHz` | 海浪最高频率（Hz） | `f_max` | 结构字段 | [sea_wave/oceanSurfaceModel_phiwind_260315.m:41](<C:/Users/Spidey/Desktop/code/sea_wave/oceanSurfaceModel_phiwind_260315.m:41>) |
| `p.frequencyBins` | 对数频率离散点数（不是点数减一的区间数） | `Fre_segment` | 结构字段 | [sea_wave/oceanSurfaceModel_phiwind_260315.m:41](<C:/Users/Spidey/Desktop/code/sea_wave/oceanSurfaceModel_phiwind_260315.m:41>) |
| `p.directionBins` | 方向离散点数，包含 ±π/2 端点 | `Dir_segment` | 结构字段 | [sea_wave/oceanSurfaceModel_phiwind_260315.m:42](<C:/Users/Spidey/Desktop/code/sea_wave/oceanSurfaceModel_phiwind_260315.m:42>) |
| `p.windSpeed` | 风速（m/s） | `U` | 结构字段 | [sea_wave/oceanSurfaceModel_phiwind_260315.m:42](<C:/Users/Spidey/Desktop/code/sea_wave/oceanSurfaceModel_phiwind_260315.m:42>) |
| `p.windAngleDeg` | 波浪传播方向角（°），从 +Y 朝 +X 增大 | `rad2deg(phi_wind)` | 结构字段 | [sea_wave/oceanSurfaceModel_phiwind_260315.m:42](<C:/Users/Spidey/Desktop/code/sea_wave/oceanSurfaceModel_phiwind_260315.m:42>) |
| `p.seed` | 随机数种子；所属结构决定是海浪、TX、MC 或湍流随机流 | `seed` | 结构字段 | [sea_wave/oceanSurfaceModel_phiwind_260315.m:42](<C:/Users/Spidey/Desktop/code/sea_wave/oceanSurfaceModel_phiwind_260315.m:42>) |
| `field` | 提供给海面重建函数的 wave/parameters 封装，或当前块复波场 | `struct('wave',wave,'parameters',struct('startTime',t0))` | 赋值 | [sea_wave/oceanSurfaceModel_phiwind_260315.m:46](<C:/Users/Spidey/Desktop/code/sea_wave/oceanSurfaceModel_phiwind_260315.m:46>) |
| `field.wave` | 外部共享的海浪分量结构；[] 时按参数生成 | `wave` | 结构字段 | [sea_wave/oceanSurfaceModel_phiwind_260315.m:46](<C:/Users/Spidey/Desktop/code/sea_wave/oceanSurfaceModel_phiwind_260315.m:46>) |
| `field.parameters` | 当前计算实际使用的参数结构；EGG 中为混合分布参数 | `struct('startTime',t0)` | 结构字段 | [sea_wave/oceanSurfaceModel_phiwind_260315.m:46](<C:/Users/Spidey/Desktop/code/sea_wave/oceanSurfaceModel_phiwind_260315.m:46>) |
| `field.parameters.startTime` | 绝对海浪起算时刻（s） | `t0` | 结构字段 | [sea_wave/oceanSurfaceModel_phiwind_260315.m:46](<C:/Users/Spidey/Desktop/code/sea_wave/oceanSurfaceModel_phiwind_260315.m:46>) |


**作用域：`calculatePlaneParameters`**

| 源码符号 | 含义与单位/维度 | 源码取值/计算式（各次赋值） | 角色 | 源码位置 |
| --- | --- | --- | --- | --- |
| `XX` | 海面网格的世界 x 坐标矩阵（m） | 函数输入/输出；值随调用或上游结果确定 | 输入 | [sea_wave/oceanSurfaceModel_phiwind_260315.m:53](<C:/Users/Spidey/Desktop/code/sea_wave/oceanSurfaceModel_phiwind_260315.m:53>) |
| `YY` | 海面网格的世界 y 坐标矩阵（m） | 函数输入/输出；值随调用或上游结果确定 | 输入 | [sea_wave/oceanSurfaceModel_phiwind_260315.m:53](<C:/Users/Spidey/Desktop/code/sea_wave/oceanSurfaceModel_phiwind_260315.m:53>) |
| `Z_surf` | 海面高度矩阵（m） | 函数输入/输出；值随调用或上游结果确定 | 输入 | [sea_wave/oceanSurfaceModel_phiwind_260315.m:53](<C:/Users/Spidey/Desktop/code/sea_wave/oceanSurfaceModel_phiwind_260315.m:53>) |
| `plane_data` | 所有网格单元的小平面参数结构 | `struct()` | 输出/赋值 | [sea_wave/oceanSurfaceModel_phiwind_260315.m:53](<C:/Users/Spidey/Desktop/code/sea_wave/oceanSurfaceModel_phiwind_260315.m:53>) |
| `rows` | 海面网格行/列数；单元数组尺寸为 (rows-1)×(cols-1) | `size(XX)` | 赋值 | [sea_wave/oceanSurfaceModel_phiwind_260315.m:56](<C:/Users/Spidey/Desktop/code/sea_wave/oceanSurfaceModel_phiwind_260315.m:56>) |
| `cols` | 海面网格行/列数；单元数组尺寸为 (rows-1)×(cols-1) | `size(XX)` | 赋值 | [sea_wave/oceanSurfaceModel_phiwind_260315.m:56](<C:/Users/Spidey/Desktop/code/sea_wave/oceanSurfaceModel_phiwind_260315.m:56>) |
| `plane_data.A` | 小平面 z=A*x+B*y+C 的系数；A/B 无量纲，C 为 m；BEM 重名见作用域 | `zeros(rows-1, cols-1)`<br>`A` | 赋值 | [sea_wave/oceanSurfaceModel_phiwind_260315.m:60](<C:/Users/Spidey/Desktop/code/sea_wave/oceanSurfaceModel_phiwind_260315.m:60>) |
| `plane_data.B` | 小平面 z=A*x+B*y+C 的系数；A/B 无量纲，C 为 m；BEM 重名见作用域 | `zeros(rows-1, cols-1)`<br>`B` | 赋值 | [sea_wave/oceanSurfaceModel_phiwind_260315.m:61](<C:/Users/Spidey/Desktop/code/sea_wave/oceanSurfaceModel_phiwind_260315.m:61>) |
| `plane_data.C` | 小平面 z=A*x+B*y+C 的系数；A/B 无量纲，C 为 m；BEM 重名见作用域 | `zeros(rows-1, cols-1)`<br>`C` | 赋值 | [sea_wave/oceanSurfaceModel_phiwind_260315.m:62](<C:/Users/Spidey/Desktop/code/sea_wave/oceanSurfaceModel_phiwind_260315.m:62>) |
| `plane_data.normal_x` | 拟合单元平面的向上单位法向量各分量 | `zeros(rows-1, cols-1)`<br>`normal(1)` | 赋值 | [sea_wave/oceanSurfaceModel_phiwind_260315.m:63](<C:/Users/Spidey/Desktop/code/sea_wave/oceanSurfaceModel_phiwind_260315.m:63>) |
| `plane_data.normal_y` | 拟合单元平面的向上单位法向量各分量 | `zeros(rows-1, cols-1)`<br>`normal(2)` | 赋值 | [sea_wave/oceanSurfaceModel_phiwind_260315.m:64](<C:/Users/Spidey/Desktop/code/sea_wave/oceanSurfaceModel_phiwind_260315.m:64>) |
| `plane_data.normal_z` | 拟合单元平面的向上单位法向量各分量 | `zeros(rows-1, cols-1)`<br>`normal(3)` | 赋值 | [sea_wave/oceanSurfaceModel_phiwind_260315.m:65](<C:/Users/Spidey/Desktop/code/sea_wave/oceanSurfaceModel_phiwind_260315.m:65>) |
| `plane_data.centroid_x` | 拟合单元四顶点均值中心的世界坐标（m） | `zeros(rows-1, cols-1)`<br>`mean(vertices(:,1))` | 赋值 | [sea_wave/oceanSurfaceModel_phiwind_260315.m:66](<C:/Users/Spidey/Desktop/code/sea_wave/oceanSurfaceModel_phiwind_260315.m:66>) |
| `plane_data.centroid_y` | 拟合单元四顶点均值中心的世界坐标（m） | `zeros(rows-1, cols-1)`<br>`mean(vertices(:,2))` | 赋值 | [sea_wave/oceanSurfaceModel_phiwind_260315.m:67](<C:/Users/Spidey/Desktop/code/sea_wave/oceanSurfaceModel_phiwind_260315.m:67>) |
| `plane_data.centroid_z` | 拟合单元四顶点均值中心的世界坐标（m） | `zeros(rows-1, cols-1)`<br>`mean(vertices(:,3))` | 赋值 | [sea_wave/oceanSurfaceModel_phiwind_260315.m:68](<C:/Users/Spidey/Desktop/code/sea_wave/oceanSurfaceModel_phiwind_260315.m:68>) |
| `plane_data.fit_error` | 小平面四顶点最小二乘高度残差 RMS（m） | `zeros(rows-1, cols-1)`<br>`fitError` | 赋值 | [sea_wave/oceanSurfaceModel_phiwind_260315.m:69](<C:/Users/Spidey/Desktop/code/sea_wave/oceanSurfaceModel_phiwind_260315.m:69>) |
| `plane_data.vertices` | 小平面四顶点世界坐标（m）或各单元顶点元胞数组 | `cell(rows-1, cols-1)`<br>`vertices` | 赋值 | [sea_wave/oceanSurfaceModel_phiwind_260315.m:70](<C:/Users/Spidey/Desktop/code/sea_wave/oceanSurfaceModel_phiwind_260315.m:70>) |
| `idx_i` | 当前小平面单元的行/列索引 | `1:rows-1` | 循环索引 | [sea_wave/oceanSurfaceModel_phiwind_260315.m:73](<C:/Users/Spidey/Desktop/code/sea_wave/oceanSurfaceModel_phiwind_260315.m:73>) |
| `idx_j` | 当前小平面单元的行/列索引 | `1:cols-1` | 循环索引 | [sea_wave/oceanSurfaceModel_phiwind_260315.m:74](<C:/Users/Spidey/Desktop/code/sea_wave/oceanSurfaceModel_phiwind_260315.m:74>) |
| `vertices` | 当前单元四个世界坐标顶点（4×3，m） | `[ XX(idx_i,idx_j), YY(idx_i,idx_j), Z_surf(idx_i,idx_j); XX(idx_i,idx_j+1), YY(idx_i,idx_j+1), Z_surf(idx_i,idx_j+1); XX(idx_i+1,idx_j), YY(idx_i+1,idx_j), Z_surf(idx_i+1,idx_j); XX(idx_i+1,idx_j+1), YY(idx_i+1,idx_j+1), Z_surf(idx_i+1,idx_j+1) ]` | 赋值 | [sea_wave/oceanSurfaceModel_phiwind_260315.m:75](<C:/Users/Spidey/Desktop/code/sea_wave/oceanSurfaceModel_phiwind_260315.m:75>) |
| `A` | 拟合 z=A*x+B*y+C 的系数（A/B 无量纲，C 为 m） | `fitPlaneLS(vertices)` | 赋值 | [sea_wave/oceanSurfaceModel_phiwind_260315.m:82](<C:/Users/Spidey/Desktop/code/sea_wave/oceanSurfaceModel_phiwind_260315.m:82>) |
| `B` | 拟合 z=A*x+B*y+C 的系数（A/B 无量纲，C 为 m） | `fitPlaneLS(vertices)` | 赋值 | [sea_wave/oceanSurfaceModel_phiwind_260315.m:82](<C:/Users/Spidey/Desktop/code/sea_wave/oceanSurfaceModel_phiwind_260315.m:82>) |
| `C` | 拟合 z=A*x+B*y+C 的系数（A/B 无量纲，C 为 m） | `fitPlaneLS(vertices)` | 赋值 | [sea_wave/oceanSurfaceModel_phiwind_260315.m:82](<C:/Users/Spidey/Desktop/code/sea_wave/oceanSurfaceModel_phiwind_260315.m:82>) |
| `fitError` | 四顶点高度拟合残差的 RMS（m） | `fitPlaneLS(vertices)` | 赋值 | [sea_wave/oceanSurfaceModel_phiwind_260315.m:82](<C:/Users/Spidey/Desktop/code/sea_wave/oceanSurfaceModel_phiwind_260315.m:82>) |
| `normal` | 当前板面或平面单位法向量 | `fitPlaneLS(vertices)` | 赋值 | [sea_wave/oceanSurfaceModel_phiwind_260315.m:82](<C:/Users/Spidey/Desktop/code/sea_wave/oceanSurfaceModel_phiwind_260315.m:82>) |


**作用域：`fitPlaneLS`**

| 源码符号 | 含义与单位/维度 | 源码取值/计算式（各次赋值） | 角色 | 源码位置 |
| --- | --- | --- | --- | --- |
| `vertices` | 待拟合的四个三维顶点（m） | 函数输入/输出；值随调用或上游结果确定 | 输入 | [sea_wave/oceanSurfaceModel_phiwind_260315.m:97](<C:/Users/Spidey/Desktop/code/sea_wave/oceanSurfaceModel_phiwind_260315.m:97>) |
| `A` | 最小二乘平面 z=A*x+B*y+C 的系数（A/B 无量纲，C 为 m） | `coefficients(1)` | 输出/赋值 | [sea_wave/oceanSurfaceModel_phiwind_260315.m:97](<C:/Users/Spidey/Desktop/code/sea_wave/oceanSurfaceModel_phiwind_260315.m:97>) |
| `B` | 最小二乘平面 z=A*x+B*y+C 的系数（A/B 无量纲，C 为 m） | `coefficients(2)` | 输出/赋值 | [sea_wave/oceanSurfaceModel_phiwind_260315.m:97](<C:/Users/Spidey/Desktop/code/sea_wave/oceanSurfaceModel_phiwind_260315.m:97>) |
| `C` | 最小二乘平面 z=A*x+B*y+C 的系数（A/B 无量纲，C 为 m） | `coefficients(3)` | 输出/赋值 | [sea_wave/oceanSurfaceModel_phiwind_260315.m:97](<C:/Users/Spidey/Desktop/code/sea_wave/oceanSurfaceModel_phiwind_260315.m:97>) |
| `fitError` | 高度预测残差均方根（m） | `sqrt(mean((z - z_pred).^2))` | 输出/赋值 | [sea_wave/oceanSurfaceModel_phiwind_260315.m:97](<C:/Users/Spidey/Desktop/code/sea_wave/oceanSurfaceModel_phiwind_260315.m:97>) |
| `normal` | 当前板面或平面单位法向量 | `normal_vec / normal_length`<br>`[0, 0, -1]`<br>`-normal` | 输出/赋值 | [sea_wave/oceanSurfaceModel_phiwind_260315.m:97](<C:/Users/Spidey/Desktop/code/sea_wave/oceanSurfaceModel_phiwind_260315.m:97>) |
| `x` | 当前坐标分量（m）；重名参数按所在函数解释 | `vertices(:,1)` | 赋值 | [sea_wave/oceanSurfaceModel_phiwind_260315.m:99](<C:/Users/Spidey/Desktop/code/sea_wave/oceanSurfaceModel_phiwind_260315.m:99>) |
| `y` | 当前坐标分量（m）；重名参数按所在函数解释 | `vertices(:,2)` | 赋值 | [sea_wave/oceanSurfaceModel_phiwind_260315.m:99](<C:/Users/Spidey/Desktop/code/sea_wave/oceanSurfaceModel_phiwind_260315.m:99>) |
| `z` | 当前坐标分量（m）；重名参数按所在函数解释 | `vertices(:,3)` | 赋值 | [sea_wave/oceanSurfaceModel_phiwind_260315.m:99](<C:/Users/Spidey/Desktop/code/sea_wave/oceanSurfaceModel_phiwind_260315.m:99>) |
| `X_design` | 平面最小二乘设计矩阵 [x,y,1] | `[x, y, ones(size(x))]` | 赋值 | [sea_wave/oceanSurfaceModel_phiwind_260315.m:100](<C:/Users/Spidey/Desktop/code/sea_wave/oceanSurfaceModel_phiwind_260315.m:100>) |
| `coefficients` | 最小二乘平面系数列向量 [A;B;C] | `X_design \ z` | 赋值 | [sea_wave/oceanSurfaceModel_phiwind_260315.m:101](<C:/Users/Spidey/Desktop/code/sea_wave/oceanSurfaceModel_phiwind_260315.m:101>) |
| `z_pred` | 拟合平面在四顶点水平位置处预测的高度（m） | `X_design * coefficients` | 赋值 | [sea_wave/oceanSurfaceModel_phiwind_260315.m:105](<C:/Users/Spidey/Desktop/code/sea_wave/oceanSurfaceModel_phiwind_260315.m:105>) |
| `normal_vec` | 未归一化平面法向 [A,B,-1] | `[A, B, -1]` | 赋值 | [sea_wave/oceanSurfaceModel_phiwind_260315.m:108](<C:/Users/Spidey/Desktop/code/sea_wave/oceanSurfaceModel_phiwind_260315.m:108>) |
| `normal_length` | 未归一化法向的欧氏模长 | `norm(normal_vec)` | 赋值 | [sea_wave/oceanSurfaceModel_phiwind_260315.m:109](<C:/Users/Spidey/Desktop/code/sea_wave/oceanSurfaceModel_phiwind_260315.m:109>) |


**作用域：`getSelectedPlaneData`**

| 源码符号 | 含义与单位/维度 | 源码取值/计算式（各次赋值） | 角色 | 源码位置 |
| --- | --- | --- | --- | --- |
| `plane_data` | 全部小平面参数结构 | 函数输入/输出；值随调用或上游结果确定 | 输入 | [sea_wave/oceanSurfaceModel_phiwind_260315.m:120](<C:/Users/Spidey/Desktop/code/sea_wave/oceanSurfaceModel_phiwind_260315.m:120>) |
| `plane_index` | 待选单元二维或线性索引 | 函数输入/输出；值随调用或上游结果确定 | 输入 | [sea_wave/oceanSurfaceModel_phiwind_260315.m:120](<C:/Users/Spidey/Desktop/code/sea_wave/oceanSurfaceModel_phiwind_260315.m:120>) |
| `selected_data` | 选定小平面参数及可读方程结构 | `struct()` | 输出/赋值 | [sea_wave/oceanSurfaceModel_phiwind_260315.m:120](<C:/Users/Spidey/Desktop/code/sea_wave/oceanSurfaceModel_phiwind_260315.m:120>) |
| `rows` | 小平面单元数组的行/列数 | `size(plane_data.A)` | 赋值 | [sea_wave/oceanSurfaceModel_phiwind_260315.m:122](<C:/Users/Spidey/Desktop/code/sea_wave/oceanSurfaceModel_phiwind_260315.m:122>) |
| `cols` | 小平面单元数组的行/列数 | `size(plane_data.A)` | 赋值 | [sea_wave/oceanSurfaceModel_phiwind_260315.m:122](<C:/Users/Spidey/Desktop/code/sea_wave/oceanSurfaceModel_phiwind_260315.m:122>) |
| `i` | 当前循环、分量或样本索引（无量纲，MATLAB 从 1 开始） | `ind2sub([rows, cols], plane_index)`<br>`plane_index(1)` | 赋值 | [sea_wave/oceanSurfaceModel_phiwind_260315.m:124](<C:/Users/Spidey/Desktop/code/sea_wave/oceanSurfaceModel_phiwind_260315.m:124>) |
| `j` | 当前循环、分量或样本索引（无量纲，MATLAB 从 1 开始） | `ind2sub([rows, cols], plane_index)`<br>`plane_index(2)` | 赋值 | [sea_wave/oceanSurfaceModel_phiwind_260315.m:124](<C:/Users/Spidey/Desktop/code/sea_wave/oceanSurfaceModel_phiwind_260315.m:124>) |
| `selected_data.index` | 所选网格平面的二维 [i,j] 索引 | `[i, j]` | 赋值 | [sea_wave/oceanSurfaceModel_phiwind_260315.m:130](<C:/Users/Spidey/Desktop/code/sea_wave/oceanSurfaceModel_phiwind_260315.m:130>) |
| `selected_data.linear_index` | 所选网格平面的 MATLAB 线性索引 | `sub2ind([rows, cols], i, j)` | 赋值 | [sea_wave/oceanSurfaceModel_phiwind_260315.m:131](<C:/Users/Spidey/Desktop/code/sea_wave/oceanSurfaceModel_phiwind_260315.m:131>) |
| `selected_data.equation` | 小平面方程 A/B/C 系数结构 | `struct('A', plane_data.A(i,j), 'B', plane_data.B(i,j), 'C', plane_data.C(i,j))` | 赋值 | [sea_wave/oceanSurfaceModel_phiwind_260315.m:132](<C:/Users/Spidey/Desktop/code/sea_wave/oceanSurfaceModel_phiwind_260315.m:132>) |
| `selected_data.equation.A` | 小平面 z=A*x+B*y+C 的系数；A/B 无量纲，C 为 m；BEM 重名见作用域 | `plane_data.A(i,j)` | 结构字段 | [sea_wave/oceanSurfaceModel_phiwind_260315.m:132](<C:/Users/Spidey/Desktop/code/sea_wave/oceanSurfaceModel_phiwind_260315.m:132>) |
| `selected_data.equation.B` | 小平面 z=A*x+B*y+C 的系数；A/B 无量纲，C 为 m；BEM 重名见作用域 | `plane_data.B(i,j)` | 结构字段 | [sea_wave/oceanSurfaceModel_phiwind_260315.m:132](<C:/Users/Spidey/Desktop/code/sea_wave/oceanSurfaceModel_phiwind_260315.m:132>) |
| `selected_data.equation.C` | 小平面 z=A*x+B*y+C 的系数；A/B 无量纲，C 为 m；BEM 重名见作用域 | `plane_data.C(i,j)` | 结构字段 | [sea_wave/oceanSurfaceModel_phiwind_260315.m:132](<C:/Users/Spidey/Desktop/code/sea_wave/oceanSurfaceModel_phiwind_260315.m:132>) |
| `selected_data.normal_vector` | 所选小平面的向上单位法向量 | `[plane_data.normal_x(i,j), plane_data.normal_y(i,j), plane_data.normal_z(i,j)]` | 赋值 | [sea_wave/oceanSurfaceModel_phiwind_260315.m:133](<C:/Users/Spidey/Desktop/code/sea_wave/oceanSurfaceModel_phiwind_260315.m:133>) |
| `selected_data.centroid` | 所选单元四顶点均值中心 [x,y,z]（m） | `[plane_data.centroid_x(i,j), plane_data.centroid_y(i,j), plane_data.centroid_z(i,j)]` | 赋值 | [sea_wave/oceanSurfaceModel_phiwind_260315.m:134](<C:/Users/Spidey/Desktop/code/sea_wave/oceanSurfaceModel_phiwind_260315.m:134>) |
| `selected_data.fit_error` | 小平面四顶点最小二乘高度残差 RMS（m） | `plane_data.fit_error(i,j)` | 赋值 | [sea_wave/oceanSurfaceModel_phiwind_260315.m:135](<C:/Users/Spidey/Desktop/code/sea_wave/oceanSurfaceModel_phiwind_260315.m:135>) |
| `selected_data.vertices` | 小平面四顶点世界坐标（m）或各单元顶点元胞数组 | `plane_data.vertices{i,j}` | 赋值 | [sea_wave/oceanSurfaceModel_phiwind_260315.m:136](<C:/Users/Spidey/Desktop/code/sea_wave/oceanSurfaceModel_phiwind_260315.m:136>) |
| `selected_data.equation_string` | 平面方程的可读文本 | `sprintf('z = %.6f*x + %.6f*y + %.6f', selected_data.equation.A, selected_data.equation.B, selected_data.equation.C)` | 赋值 | [sea_wave/oceanSurfaceModel_phiwind_260315.m:137](<C:/Users/Spidey/Desktop/code/sea_wave/oceanSurfaceModel_phiwind_260315.m:137>) |


### 5.6 sea_wave/RIS_draw_state.m

**作用域：`RIS_draw_state`**

| 源码符号 | 含义与单位/维度 | 源码取值/计算式（各次赋值） | 角色 | 源码位置 |
| --- | --- | --- | --- | --- |
| `ax` | 绘图坐标轴句柄 | 函数输入/输出；值随调用或上游结果确定 | 输入 | [sea_wave/RIS_draw_state.m:1](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_draw_state.m:1>) |
| `motion` | 含上表面四顶点、旋转矩阵和板尺寸的 RIS 运动结果 | 函数输入/输出；值随调用或上游结果确定 | 输入 | [sea_wave/RIS_draw_state.m:1](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_draw_state.m:1>) |
| `time` | 相对 parameters.startTime 的时间（s）或其采样向量 | 函数输入/输出；值随调用或上游结果确定 | 输入 | [sea_wave/RIS_draw_state.m:1](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_draw_state.m:1>) |
| `scene` | 海面显示快照结构（X/Y/Z）或其集合 | 函数输入/输出；值随调用或上游结果确定 | 输入 | [sea_wave/RIS_draw_state.m:1](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_draw_state.m:1>) |
| `heightLimits` | 海面高度、轨迹与色轴共同显示范围，默认 [-0.5,0.5]（m） | `[-0.5,0.5]` | 输入/声明默认值 | [sea_wave/RIS_draw_state.m:1](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_draw_state.m:1>) |
| `index` | 被选中的分量、样本或数组位置索引 | `find(abs(motion.t-time)<1e-10,1)` | 赋值 | [sea_wave/RIS_draw_state.m:11](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_draw_state.m:11>) |
| `normal` | 当前板面或平面单位法向量 | `motion.normal(:,index)` | 赋值 | [sea_wave/RIS_draw_state.m:19](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_draw_state.m:19>) |
| `center` | 当前 RIS 上表面中心世界坐标（m） | `motion.boardCenter(:,index)` | 赋值 | [sea_wave/RIS_draw_state.m:19](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_draw_state.m:19>) |
| `R` | 体坐标到世界坐标的主动旋转矩阵（无量纲） | `motion.rotationMatrix(:,:,index)`<br>`[(board(2,:)-board(1,:)).'/motion.parameters.Lx, (board(4,:)-board(1,:)).'/motion.parameters.Ly,normal]` | 赋值 | [sea_wave/RIS_draw_state.m:23](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_draw_state.m:23>) |
| `board` | 旧格式运动结果中的四个上表面顶点（m） | `motion.corners(:,:,index)` | 赋值 | [sea_wave/RIS_draw_state.m:26](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_draw_state.m:26>) |
| `colors` | 板体 x/y 方向箭头的两种 RGB 颜色 | `[0.8,0.12,0.12;0.1,0.55,0.2]` | 赋值 | [sea_wave/RIS_draw_state.m:30](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_draw_state.m:30>) |
| `bodyAxis` | 要显示的板体轴索引 1 或 2 | `1:2` | 循环索引 | [sea_wave/RIS_draw_state.m:31](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_draw_state.m:31>) |
| `direction` | 板体方向箭头的世界位移向量，长度 0.4*min(Lx,Ly)（m） | `0.4*min(motion.parameters.Lx,motion.parameters.Ly)*R(:,bodyAxis)` | 赋值 | [sea_wave/RIS_draw_state.m:32](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_draw_state.m:32>) |
| `cb` | 颜色条句柄 | `colorbar(ax)` | 赋值 | [sea_wave/RIS_draw_state.m:42](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_draw_state.m:42>) |
| `cb.Label.String` | 颜色条标题文本 | `'Elevation (m)'` | 赋值 | [sea_wave/RIS_draw_state.m:42](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_draw_state.m:42>) |


### 5.7 sea_wave/RIS_float.m

**作用域：`RIS_float`**

| 源码符号 | 含义与单位/维度 | 源码取值/计算式（各次赋值） | 角色 | 源码位置 |
| --- | --- | --- | --- | --- |
| `cfg` | 调用方传入的 RIS 参数覆盖结构；仅接受默认参数结构中已有字段 | `struct()` | 输入/声明默认值 | [sea_wave/RIS_float.m:1](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:1>) |
| `out` | 当前函数输出结果结构 | `assembleOutput(p,model,t,states,wave,spectrum,referenceCenter,draft,zG,step)` | 输出/赋值 | [sea_wave/RIS_float.m:1](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:1>) |
| `p` | 当前作用域的参数结构 | `get_default_RIS_params()`<br>`cfg.(fields{j})`<br>`validateParameters(p)` | 赋值 | [sea_wave/RIS_float.m:18](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:18>) |
| `p.plot` | 是否自动绘制 RIS 单次运动；RIS_float 无输出且未覆盖时设为 true | `nargout==0` | 赋值 | [sea_wave/RIS_float.m:19](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:19>) |
| `fields` | 字段名、文件名或表格列名集合（用途见计算式） | `fieldnames(cfg)` | 赋值 | [sea_wave/RIS_float.m:20](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:20>) |
| `j` | 当前循环、分量或样本索引（无量纲，MATLAB 从 1 开始） | `1:numel(fields)` | 循环索引 | [sea_wave/RIS_float.m:21](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:21>) |
| `area` | 当前板面或网格单元面积（m²）；见所在函数计算式 | `p.Lx*p.Ly` | 赋值 | [sea_wave/RIS_float.m:28](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:28>) |
| `draft` | 静水平衡吃水 mass/(rho*Lx*Ly)（m） | `p.mass/(p.rho*area)` | 赋值 | [sea_wave/RIS_float.m:29](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:29>) |
| `freeboard` | 静水干舷 thickness-draft（m） | `p.thickness-draft` | 赋值 | [sea_wave/RIS_float.m:33](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:33>) |
| `zG` | 静水平衡质心高度 thickness/2-draft（m） | `p.thickness/2-draft` | 赋值 | [sea_wave/RIS_float.m:34](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:34>) |
| `inertia` | 实体矩形板三主转动惯量 [Ix;Iy;Iz]（kg·m²） | `p.mass/12*[p.Ly^2+p.thickness^2;p.Lx^2+p.thickness^2;p.Lx^2+p.Ly^2]` | 赋值 | [sea_wave/RIS_float.m:35](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:35>) |
| `rigidMass` | 六分量刚体广义质量 [m;m;m;Ix;Iy;Iz]（kg / kg·m²） | `[p.mass;p.mass;p.mass;inertia]` | 赋值 | [sea_wave/RIS_float.m:36](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:36>) |
| `restoring` | 垂向模态 [heave;roll;pitch] 的静水恢复刚度（N/m / N·m/rad） | `[p.rho*p.g*area; p.rho*p.g*p.Lx*p.Ly^3/12-p.mass*p.g*freeboard/2; p.rho*p.g*p.Ly*p.Lx^3/12-p.mass*p.g*freeboard/2]` | 赋值 | [sea_wave/RIS_float.m:37](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:37>) |
| `wave` | 海浪分量结构，含 omega、kx、ky、coefficient；具体重名按下方作用域解释 | `makeWave(p)` | 赋值 | [sea_wave/RIS_float.m:43](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:43>) |
| `spectrum` | 海浪频谱及离散积分诊断结构 | `makeWave(p)` | 赋值 | [sea_wave/RIS_float.m:43](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:43>) |
| `heading` | 水动力参考航向角，来自初始 yaw，默认 0（rad） | `0`<br>`p.initialEulerAngles(3)` | 赋值 | [sea_wave/RIS_float.m:44](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:44>) |
| `headingR` | 初始水平航向的参考系到世界系旋转矩阵 | `RIS_rotation([0;0;heading])` | 赋值 | [sea_wave/RIS_float.m:46](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:46>) |
| `bodyWave` | 将波数旋转到初始航向体轴后的同一组入射海浪分量 | `wave` | 赋值 | [sea_wave/RIS_float.m:47](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:47>) |
| `bodyK` | 初始体坐标波数水平分量矩阵（2×Nwave，m⁻¹） | `headingR(1:2,1:2).'*[wave.kx.';wave.ky.']` | 赋值 | [sea_wave/RIS_float.m:48](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:48>) |
| `bodyWave.kx` | 海浪波数的世界 x/y 分量（m⁻¹）；bodyWave 中为初始体坐标分量 | `bodyK(1,:).'` | 赋值 | [sea_wave/RIS_float.m:49](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:49>) |
| `bodyWave.ky` | 海浪波数的世界 x/y 分量（m⁻¹）；bodyWave 中为初始体坐标分量 | `bodyK(2,:).'` | 赋值 | [sea_wave/RIS_float.m:49](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:49>) |
| `fitFrequency` | 额外用于辐射记忆拟合的对数角频率（rad/s） | `2*pi*logspace(log10(p.fMinHz),log10(p.fMaxHz),p.radiationFrequencyBins).'` | 赋值 | [sea_wave/RIS_float.m:50](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:50>) |
| `frequency` | 海浪频率与拟合频率去重合并的角频率向量（rad/s） | `unique([wave.omega;fitFrequency])` | 赋值 | [sea_wave/RIS_float.m:51](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:51>) |
| `index` | 被选中的分量、样本或数组位置索引 | `ismember(wave.omega,frequency)` | 赋值 | [sea_wave/RIS_float.m:52](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:52>) |
| `modalMass` | 原 BEM 顺序 [heave;x-height slope;y-height slope] 的刚体质量/惯量 | `rigidMass([3,5,4])` | 赋值 | [sea_wave/RIS_float.m:54](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:54>) |
| `modalK` | 同一 BEM 顺序的静水恢复刚度 | `restoring([1,3,2])` | 赋值 | [sea_wave/RIS_float.m:55](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:55>) |
| `oldRao` | 原三模态频域响应算子；x-height slope 与右手 pitch 反号 | `plateResponse(p,frequency,bodyWave.kx,bodyWave.ky,index,modalMass,modalK)` | 赋值 | [sea_wave/RIS_float.m:56](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:56>) |
| `A` | 原 BEM 三模态附加广义质量频响（kg / kg·m²） | `plateResponse(p,frequency,bodyWave.kx,bodyWave.ky,index,modalMass,modalK)` | 赋值 | [sea_wave/RIS_float.m:56](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:56>) |
| `B` | 原 BEM 三模态辐射广义阻尼（kg/s / kg·m²/s） | `plateResponse(p,frequency,bodyWave.kx,bodyWave.ky,index,modalMass,modalK)` | 赋值 | [sea_wave/RIS_float.m:56](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:56>) |
| `F` | 原三模态单位波幅激励复幅（N/m / N·m/m） | `plateResponse(p,frequency,bodyWave.kx,bodyWave.ky,index,modalMass,modalK)` | 赋值 | [sea_wave/RIS_float.m:56](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:56>) |
| `radiation` | 转换为 [heave;roll;pitch] 顺序后的辐射记忆参数结构 | `RIS_radiation(frequency,A,B,p.radiationModes)` | 赋值 | [sea_wave/RIS_float.m:57](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:57>) |
| `modes` | 六维广义状态中受垂向辐射影响的索引 [3,4,5] | `[3,4,5]` | 赋值 | [sea_wave/RIS_float.m:59](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:59>) |
| `radiation.residues` | 各辐射振子对三个垂向模态的非负留数矩阵 | `radiation.residues(:,[1,3,2])` | 赋值 | [sea_wave/RIS_float.m:60](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:60>) |
| `radiation.infiniteAddedMass` | 无穷频率附加广义质量；升沉项 kg、转动项 kg·m² | `radiation.infiniteAddedMass([1,3,2])` | 赋值 | [sea_wave/RIS_float.m:61](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:61>) |
| `radiation.relativeImpedanceRmsError` | 按频响幅值加权的阻抗拟合相对均方根误差（无量纲，三个垂向模态） | `radiation.relativeImpedanceRmsError([1,3,2])` | 赋值 | [sea_wave/RIS_float.m:62](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:62>) |
| `radiation.fittedAddedMass` | 辐射有理模型拟合的附加广义质量频响 | `radiation.fittedAddedMass([1,3,2],:)` | 赋值 | [sea_wave/RIS_float.m:63](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:63>) |
| `radiation.fittedDamping` | 辐射有理模型拟合的广义阻尼频响 | `radiation.fittedDamping([1,3,2],:)` | 赋值 | [sea_wave/RIS_float.m:64](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:64>) |
| `totalMass` | 刚体广义质量叠加无穷频附加质量后的六分量质量/惯量 | `rigidMass`<br>`totalMass(modes)+radiation.infiniteAddedMass` | 赋值 | [sea_wave/RIS_float.m:65](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:65>) |
| `excitation` | 最终六维单位波幅波浪激励，合并底面 BEM 与侧壁压力 | `RIS_sidewall_load(p,bodyWave,draft,zG)`<br>`F(1,:)`<br>`excitation(4,:)+F(3,:)`<br>`excitation(5,:)-F(2,:)` | 赋值 | [sea_wave/RIS_float.m:67](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:67>) |
| `initialXY` | 相对 p.center 的初始水平状态（m） | `p.horizontalInitialPosition`<br>`p.initialPosition(1:2)` | 赋值 | [sea_wave/RIS_float.m:71](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:71>) |
| `referenceCenter` | 用于线性载荷相位的初始世界水平位置 p.center+initialXY（m） | `p.center+initialXY` | 赋值 | [sea_wave/RIS_float.m:73](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:73>) |
| `c` | 包含初始绝对时间与参考位置相位的海浪复幅（m） | `wave.coefficient.*exp(1i*(wave.omega*p.startTime- wave.kx*referenceCenter(1)-wave.ky*referenceCenter(2)))` | 赋值 | [sea_wave/RIS_float.m:74](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:74>) |
| `forceHat` | 六维各海浪分量的激励复幅（N / N·m） | `excitation.*c.'` | 赋值 | [sea_wave/RIS_float.m:76](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:76>) |
| `model` | 积分器使用的广义质量、恢复力、辐射、参考旋转与激励结构 | `struct(...)（所有字段在本节逐行列明）` | 赋值 | [sea_wave/RIS_float.m:77](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:77>) |
| `model.mass` | 含无穷频附加质量的六分量广义质量/惯量 | `totalMass` | 结构字段 | [sea_wave/RIS_float.m:77](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:77>) |
| `model.rigidMass` | 六分量刚体广义质量 [m;m;m;Ix;Iy;Iz]（kg / kg·m²） | `rigidMass` | 结构字段 | [sea_wave/RIS_float.m:77](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:77>) |
| `model.restoring` | 垂向模态 [heave;roll;pitch] 的静水恢复刚度（N/m / N·m/rad） | `restoring` | 结构字段 | [sea_wave/RIS_float.m:77](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:77>) |
| `model.radiation` | 辐射记忆广义力，或模型的辐射参数结构；由父结构与表达式确定 | `radiation` | 结构字段 | [sea_wave/RIS_float.m:78](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:78>) |
| `model.headingR` | 初始水平航向的参考系到世界系旋转矩阵 | `headingR` | 结构字段 | [sea_wave/RIS_float.m:78](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:78>) |
| `model.forceHat` | 包含绝对起时与参考位置相位的六维波浪激励复幅 | `forceHat` | 结构字段 | [sea_wave/RIS_float.m:78](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:78>) |
| `model.omega` | 海浪角频率（rad/s）；EGG parameters.omega 则是混合概率 | `wave.omega` | 结构字段 | [sea_wave/RIS_float.m:78](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:78>) |
| `state` | 初始状态：[3位置;4四元数;6体速度;3Nr记忆x;3Nr记忆y] | `initialState(p,model,initialXY,heading)` | 赋值 | [sea_wave/RIS_float.m:79](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:79>) |
| `admittance` | 六自由度动态阻抗逐元素倒数 | `initialState(p,model,initialXY,heading)` | 赋值 | [sea_wave/RIS_float.m:79](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:79>) |
| `rao` | 六维单位波幅响应算子 admittance.*excitation | `admittance.*excitation` | 赋值 | [sea_wave/RIS_float.m:80](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:80>) |
| `t` | 相对 parameters.startTime 的时间（s）或其采样向量 | `unique([0:p.dt:p.duration,p.duration])` | 赋值 | [sea_wave/RIS_float.m:81](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:81>) |
| `states` | 所有保存时间点的完整状态矩阵 | `integrate(model,state,t,p.integrationStep)` | 赋值 | [sea_wave/RIS_float.m:82](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:82>) |
| `step` | RK4 实际最大内部步长（s） | `integrate(model,state,t,p.integrationStep)` | 赋值 | [sea_wave/RIS_float.m:82](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:82>) |
| `out.freeboard` | 静水干舷 thickness-draft（m） | `freeboard` | 赋值 | [sea_wave/RIS_float.m:84](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:84>) |
| `out.addedMass` | 原 BEM 三模态附加广义质量，顺序 [heave;x-height slope;y-height slope] | `A` | 赋值 | [sea_wave/RIS_float.m:85](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:85>) |
| `out.radiationDamping` | 原 BEM 三模态辐射广义阻尼，顺序同 addedMass | `B` | 赋值 | [sea_wave/RIS_float.m:86](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:86>) |
| `out.hydroFrequencyHz` | BEM/辐射计算的合并频率向量（Hz） | `frequency/(2*pi)` | 赋值 | [sea_wave/RIS_float.m:87](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:87>) |
| `out.excitationPerMetre` | 每米入射波幅对应的六维激励复幅（N/m / N·m/m） | `excitation` | 赋值 | [sea_wave/RIS_float.m:88](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:88>) |
| `out.responseAmplitudeOperator` | 旧三模态响应算子，经符号处理后顺序 [heave;pitch;roll] | `[oldRao(1,:);-oldRao(2,:);oldRao(3,:)]` | 赋值 | [sea_wave/RIS_float.m:89](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:89>) |
| `out.responseAmplitudeOperator6` | 六自由度响应算子 admittance.*excitation | `rao` | 赋值 | [sea_wave/RIS_float.m:90](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:90>) |
| `out.responseComponents` | 各海浪分量对应的六维响应复幅 | `rao.*c.'` | 赋值 | [sea_wave/RIS_float.m:91](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:91>) |
| `out.windows` | 5 s、10 s 窗口汇总；仅保留不超过 duration 的窗口 | `windowSummary(out,[5,10])` | 赋值 | [sea_wave/RIS_float.m:92](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:92>) |


**作用域：`initialState`**

| 源码符号 | 含义与单位/维度 | 源码取值/计算式（各次赋值） | 角色 | 源码位置 |
| --- | --- | --- | --- | --- |
| `p` | 当前作用域的参数结构 | 函数输入/输出；值随调用或上游结果确定 | 输入 | [sea_wave/RIS_float.m:95](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:95>) |
| `m` | 六自由度动力学模型结构 | 函数输入/输出；值随调用或上游结果确定 | 输入 | [sea_wave/RIS_float.m:95](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:95>) |
| `initialXY` | 初始水平状态（m），相对 p.center | 函数输入/输出；值随调用或上游结果确定 | 输入 | [sea_wave/RIS_float.m:95](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:95>) |
| `heading` | 初始参考 yaw（rad） | 函数输入/输出；值随调用或上游结果确定 | 输入 | [sea_wave/RIS_float.m:95](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:95>) |
| `state` | 完整初始状态 [位置3;四元数4;体速度6;辐射记忆6Nr] | `zeros(13+6*numel(rad.width),1)`<br>`[initialXY;0]`<br>`displacement(3)`<br>`rotation.'*worldVelocity`<br>`angular`<br>`[memoryX(:);memoryY(:)]`<br>`p.initialPosition`<br>`quaternion`<br>`p.initialBodyVelocity` | 输出/赋值 | [sea_wave/RIS_float.m:95](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:95>) |
| `rao` | 六维动态柔度 1./Z；不同于乘入激励后的 responseAmplitudeOperator6 | `1./Z` | 输出/赋值 | [sea_wave/RIS_float.m:95](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:95>) |
| `s` | 虚频率 1i*omega 行向量（s⁻¹） | `1i*m.omega.'` | 赋值 | [sea_wave/RIS_float.m:97](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:97>) |
| `rad` | 辐射记忆模型结构的局部副本 | `m.radiation` | 赋值 | [sea_wave/RIS_float.m:98](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:98>) |
| `basis` | 辐射振子在虚频率上的响应基函数矩阵 | `(s+rad.width)./((s+rad.width).^2+rad.frequency.^2)` | 赋值 | [sea_wave/RIS_float.m:99](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:99>) |
| `impedance` | 三个垂向模态的辐射阻抗 | `rad.residues.'*basis` | 赋值 | [sea_wave/RIS_float.m:100](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:100>) |
| `stiffness` | 六维恢复刚度 [0;0;Kheave;Kroll;Kpitch;0] | `[0;0;m.restoring;0]` | 赋值 | [sea_wave/RIS_float.m:101](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:101>) |
| `Z` | 质量、恢复力及辐射项组成的六维动态阻抗 | `m.mass.*s.^2+stiffness`<br>`Z([3,4,5],:)+s.*impedance` | 赋值 | [sea_wave/RIS_float.m:102](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:102>) |
| `response` | 每个海浪分量的六维稳态响应复幅 | `m.forceHat./Z` | 赋值 | [sea_wave/RIS_float.m:104](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:104>) |
| `angles` | 右手 ZYX 欧拉角 [roll;pitch;yaw]（rad） | `[0;0;heading]`<br>`displacement(4:5)`<br>`p.initialEulerAngles` | 赋值 | [sea_wave/RIS_float.m:108](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:108>) |
| `displacement` | 初始时刻各谐波响应实部之和（m / rad） | `real(sum(response,2))` | 赋值 | [sea_wave/RIS_float.m:110](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:110>) |
| `velocity` | 初始时刻各谐波一阶导数实部之和（m/s / rad/s） | `real(sum(response.*s,2))` | 赋值 | [sea_wave/RIS_float.m:111](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:111>) |
| `worldVelocity` | 初始质心平移速度的世界系分量（m/s） | `m.headingR*velocity(1:3)`<br>`rotation*state(8:10)`<br>`p.horizontalInitialVelocity` | 赋值 | [sea_wave/RIS_float.m:114](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:114>) |
| `angular` | 初始体坐标角速度（rad/s）；稳态初始化时偏航角速度强制为 0 | `velocity(4:6)`<br>`0` | 赋值 | [sea_wave/RIS_float.m:115](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:115>) |
| `rotation` | 体坐标到世界坐标的主动旋转矩阵（无量纲） | `RIS_rotation(angles)` | 赋值 | [sea_wave/RIS_float.m:116](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:116>) |
| `modalVelocity` | 三个垂向模态的频域速度响应 | `response([3,4,5],:).*s` | 赋值 | [sea_wave/RIS_float.m:119](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:119>) |
| `memoryX` | 被动辐射振子的两组实初始记忆状态（各 Nr×3） | `real(basis*modalVelocity.')` | 赋值 | [sea_wave/RIS_float.m:120](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:120>) |
| `memoryY` | 被动辐射振子的两组实初始记忆状态（各 Nr×3） | `real((rad.frequency./((s+rad.width).^2+rad.frequency.^2))*modalVelocity.')` | 赋值 | [sea_wave/RIS_float.m:121](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:121>) |
| `quaternion` | 姿态单位四元数，顺序 [w;x;y;z] | `RIS_rotation(angles)` | 赋值 | [sea_wave/RIS_float.m:126](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:126>) |


**作用域：`integrate`**

| 源码符号 | 含义与单位/维度 | 源码取值/计算式（各次赋值） | 角色 | 源码位置 |
| --- | --- | --- | --- | --- |
| `m` | 六自由度动力学模型结构 | 函数输入/输出；值随调用或上游结果确定 | 输入 | [sea_wave/RIS_float.m:135](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:135>) |
| `state` | 当前完整积分状态向量 | `state+h/6*(k1+2*k2+2*k3+k4)`<br>`state(4:7)/norm(state(4:7))` | 输入/赋值 | [sea_wave/RIS_float.m:135](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:135>) |
| `t` | 相对 parameters.startTime 的时间（s）或其采样向量 | 函数输入/输出；值随调用或上游结果确定 | 输入 | [sea_wave/RIS_float.m:135](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:135>) |
| `requestedStep` | 用户要求的最大 RK4 内部步长（s） | 函数输入/输出；值随调用或上游结果确定 | 输入 | [sea_wave/RIS_float.m:135](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:135>) |
| `states` | 全部输出时刻保存的完整状态矩阵 | `zeros(numel(state),numel(t))`<br>`state` | 输出/赋值 | [sea_wave/RIS_float.m:135](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:135>) |
| `largestStep` | 实际采用过的最大 RK4 内部步长（s） | `0`<br>`max(largestStep,h)` | 输出/赋值 | [sea_wave/RIS_float.m:135](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:135>) |
| `poleRate` | 所有辐射振子的最大极点模长（s⁻¹） | `max(hypot(m.radiation.width,m.radiation.frequency))` | 赋值 | [sea_wave/RIS_float.m:138](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:138>) |
| `naturalRate` | 按垂向恢复刚度/有效质量计算的最大固有角频率（s⁻¹） | `sqrt(max(m.restoring./m.mass([3,4,5])))` | 赋值 | [sea_wave/RIS_float.m:139](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:139>) |
| `safeStep` | 稳定性限制后的最大内步长 min(requestedStep,0.6/poleRate,0.3/naturalRate)（s） | `min([requestedStep,0.6/poleRate,0.3/naturalRate])` | 赋值 | [sea_wave/RIS_float.m:140](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:140>) |
| `j` | 当前循环、分量或样本索引（无量纲，MATLAB 从 1 开始） | `2:numel(t)` | 循环索引 | [sea_wave/RIS_float.m:143](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:143>) |
| `steps` | 两个保存时刻之间所需的 RK4 子步数 | `ceil((t(j)-t(j-1))/safeStep)` | 赋值 | [sea_wave/RIS_float.m:144](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:144>) |
| `h` | 本输出区间实际内部积分步长（s） | `(t(j)-t(j-1))/steps` | 赋值 | [sea_wave/RIS_float.m:145](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:145>) |
| `stageTimes` | 当前输出区间 RK4 全部半步激励采样时刻（s） | `t(j-1)+(0:2*steps)*h/2` | 赋值 | [sea_wave/RIS_float.m:148](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:148>) |
| `forces` | RK4 子阶段预计算的六维实谐波激励（N / N·m） | `real(m.forceHat*exp(1i*(m.omega*stageTimes)))` | 赋值 | [sea_wave/RIS_float.m:149](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:149>) |
| `k` | 当前循环、分量或样本索引（无量纲，MATLAB 从 1 开始） | `1:steps` | 循环索引 | [sea_wave/RIS_float.m:150](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:150>) |
| `first` | 当前分块/积分子步的起始索引 | `2*k-1` | 赋值 | [sea_wave/RIS_float.m:151](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:151>) |
| `k1` | 经典四阶 Runge–Kutta 四个阶段的完整状态导数 | `derivative(state,m,forces(:,first))` | 赋值 | [sea_wave/RIS_float.m:152](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:152>) |
| `k2` | 经典四阶 Runge–Kutta 四个阶段的完整状态导数 | `derivative(state+h/2*k1,m,forces(:,first+1))` | 赋值 | [sea_wave/RIS_float.m:153](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:153>) |
| `k3` | 经典四阶 Runge–Kutta 四个阶段的完整状态导数 | `derivative(state+h/2*k2,m,forces(:,first+1))` | 赋值 | [sea_wave/RIS_float.m:154](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:154>) |
| `k4` | 经典四阶 Runge–Kutta 四个阶段的完整状态导数 | `derivative(state+h*k3,m,forces(:,first+2))` | 赋值 | [sea_wave/RIS_float.m:155](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:155>) |


**作用域：`derivative`**

| 源码符号 | 含义与单位/维度 | 源码取值/计算式（各次赋值） | 角色 | 源码位置 |
| --- | --- | --- | --- | --- |
| `state` | 完整当前积分状态 [位置;四元数;体速度;两组记忆] | 函数输入/输出；值随调用或上游结果确定 | 输入 | [sea_wave/RIS_float.m:165](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:165>) |
| `m` | 动力学模型结构 | 函数输入/输出；值随调用或上游结果确定 | 输入 | [sea_wave/RIS_float.m:165](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:165>) |
| `force` | 当前时刻初始参考轴六维波浪力/力矩（N / N·m） | 函数输入/输出；值随调用或上游结果确定 | 输入 | [sea_wave/RIS_float.m:165](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:165>) |
| `d` | 完整状态时间导数 | `[worldVelocity;dq;acceleration;dx(:);dy(:)]` | 输出/赋值 | [sea_wave/RIS_float.m:165](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:165>) |
| `load` | 当前波浪、静水、辐射、耦合广义力与加速度及残差结构 | `struct(...)（所有字段在本节逐行列明）` | 输出/赋值 | [sea_wave/RIS_float.m:165](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:165>) |
| `n` | 世界竖直方向在当前体坐标中的表示 R(3,:).' | `R(3,:).';`<br>`R(3,:).'` | 声明默认值/赋值 | [sea_wave/RIS_float.m:178](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:178>) |
| `R` | 体坐标到世界坐标的主动旋转矩阵（无量纲） | `RIS_rotation(state(4:7))` | 赋值 | [sea_wave/RIS_float.m:168](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:168>) |
| `q` | 当前姿态单位四元数 wxyz；此处不是六自由度 q | `RIS_rotation(state(4:7))` | 赋值 | [sea_wave/RIS_float.m:168](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:168>) |
| `nu` | 体坐标六速度 [u;v;w;p;q;r]（m/s / rad/s） | `state(8:13)` | 赋值 | [sea_wave/RIS_float.m:169](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:169>) |
| `v` | 体坐标三平移速度（m/s） | `nu(1:3)` | 赋值 | [sea_wave/RIS_float.m:169](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:169>) |
| `omega` | 体坐标三角速度（rad/s）；此处不是海浪角频率 | `nu(4:6)` | 赋值 | [sea_wave/RIS_float.m:169](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:169>) |
| `ref` | 初始航向参考轴到世界轴的旋转矩阵 | `m.headingR` | 赋值 | [sea_wave/RIS_float.m:170](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:170>) |
| `worldVelocity` | 质心世界平移速度（m/s） | `R*v` | 赋值 | [sea_wave/RIS_float.m:171](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:171>) |
| `worldOmega` | 世界坐标角速度（rad/s） | `R*omega` | 赋值 | [sea_wave/RIS_float.m:171](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:171>) |
| `worldMass` | 平移有效质量矩阵，垂向含附加质量（kg） | `diag([m.rigidMass(1),m.rigidMass(2),m.mass(3)])` | 赋值 | [sea_wave/RIS_float.m:172](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:172>) |
| `worldInertia` | 姿态变换后的刚体世界转动惯量矩阵（kg·m²） | `R*diag(m.rigidMass(4:6))*R.'` | 赋值 | [sea_wave/RIS_float.m:173](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:173>) |
| `addedInertia` | 初始参考轴固定的横摇/纵摇附加转动惯量矩阵（kg·m²） | `ref*diag([m.radiation.infiniteAddedMass(2:3);0])*ref.'` | 赋值 | [sea_wave/RIS_float.m:174](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:174>) |
| `effectiveInertia` | 世界系刚体惯量与附加转动惯量之和（kg·m²） | `worldInertia+addedInertia` | 赋值 | [sea_wave/RIS_float.m:175](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:175>) |
| `effectiveBody` | 当前体坐标有效 6×6 广义质量矩阵 | `blkdiag(R.'*worldMass*R,R.'*effectiveInertia*R)` | 赋值 | [sea_wave/RIS_float.m:176](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:176>) |
| `coupling` | 体轴旋转产生的平移耦合与陀螺力矩（N / N·m） | `[effectiveBody(1:3,1:3)*cross(omega,v);cross(omega,m.rigidMass(4:6).*omega)]` | 赋值 | [sea_wave/RIS_float.m:177](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:177>) |
| `grad` | 恢复势能对板体竖直方向分量的梯度，用于静水恢复力矩 | `[m.restoring(3)*n(1);m.restoring(2)*n(2);0]` | 赋值 | [sea_wave/RIS_float.m:179](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:179>) |
| `hydrostatic` | 体坐标静水恢复广义力（N / N·m） | `[R.'*[0;0;-m.restoring(1)*state(3)];cross(n,grad)]` | 赋值 | [sea_wave/RIS_float.m:180](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:180>) |
| `waveLoad` | 转换到当前体坐标的波浪广义力（N / N·m） | `[R.'*ref*force(1:3);R.'*ref*force(4:6)]` | 赋值 | [sea_wave/RIS_float.m:181](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:181>) |
| `nr` | 每个垂向模态的辐射振子数 | `numel(m.radiation.width)` | 赋值 | [sea_wave/RIS_float.m:182](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:182>) |
| `x` | Nr×3 辐射振子的两组实记忆状态；此处不是空间坐标 | `reshape(state(14:13+3*nr),nr,3)` | 赋值 | [sea_wave/RIS_float.m:183](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:183>) |
| `y` | Nr×3 辐射振子的两组实记忆状态；此处不是空间坐标 | `reshape(state(14+3*nr:end),nr,3)` | 赋值 | [sea_wave/RIS_float.m:184](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:184>) |
| `memory` | 三个垂向模态的记忆广义力 sum(residues.*x) | `sum(m.radiation.residues.*x,1).'` | 赋值 | [sea_wave/RIS_float.m:185](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:185>) |
| `radiationLoad` | 当前体坐标辐射记忆广义力（N / N·m） | `[R.'*ref*[0;0;-memory(1)];R.'*ref*[-memory(2:3);0]]` | 赋值 | [sea_wave/RIS_float.m:186](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:186>) |
| `acceleration` | 六维体速度时间导数（m/s² / rad/s²） | `effectiveBody\(waveLoad+hydrostatic+radiationLoad-coupling)` | 赋值 | [sea_wave/RIS_float.m:187](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:187>) |
| `refVelocity` | 辐射参考系的 [垂向速度;横摇角速度;纵摇角速度] | `[worldVelocity(3);ref(:,1:2).'*worldOmega]` | 赋值 | [sea_wave/RIS_float.m:188](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:188>) |
| `dx` | 两组辐射振子记忆状态的时间导数 | `-m.radiation.width.*x-m.radiation.frequency.*y+refVelocity.'` | 赋值 | [sea_wave/RIS_float.m:189](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:189>) |
| `dy` | 两组辐射振子记忆状态的时间导数 | `m.radiation.frequency.*x-m.radiation.width.*y` | 赋值 | [sea_wave/RIS_float.m:190](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:190>) |
| `dq` | 单位四元数 wxyz 的时间导数（s⁻¹） | `0.5*[-q(2:4).'*omega;q(1)*omega+cross(q(2:4),omega)]` | 赋值 | [sea_wave/RIS_float.m:191](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:191>) |
| `load.wave` | 当前时刻体坐标波浪激励广义力（N / N·m） | `waveLoad` | 结构字段 | [sea_wave/RIS_float.m:194](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:194>) |
| `load.hydrostatic` | 静水恢复广义力（前三项 N，后三项 N·m） | `hydrostatic` | 结构字段 | [sea_wave/RIS_float.m:194](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:194>) |
| `load.radiation` | 辐射记忆广义力，或模型的辐射参数结构；由父结构与表达式确定 | `radiationLoad` | 结构字段 | [sea_wave/RIS_float.m:194](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:194>) |
| `load.coriolis` | 体坐标陀螺及参考附加惯性耦合广义力 [N;N;N;N·m;N·m;N·m] | `coupling` | 结构字段 | [sea_wave/RIS_float.m:195](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:195>) |
| `load.acceleration` | 六维体速度的时间导数（m/s² / rad/s²） | `acceleration` | 结构字段 | [sea_wave/RIS_float.m:195](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:195>) |
| `load.effectiveMass` | 当前姿态下有效六维广义质量矩阵（体坐标） | `effectiveBody` | 结构字段 | [sea_wave/RIS_float.m:195](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:195>) |
| `load.residual` | 牛顿—欧拉方程广义力残差（N / N·m） | `effectiveBody*acceleration+coupling-waveLoad-hydrostatic-radiationLoad` | 结构字段 | [sea_wave/RIS_float.m:196](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:196>) |


**作用域：`assembleOutput`**

| 源码符号 | 含义与单位/维度 | 源码取值/计算式（各次赋值） | 角色 | 源码位置 |
| --- | --- | --- | --- | --- |
| `p` | 当前作用域的参数结构 | 函数输入/输出；值随调用或上游结果确定 | 输入 | [sea_wave/RIS_float.m:199](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:199>) |
| `m` | 六自由度动力学模型结构 | 函数输入/输出；值随调用或上游结果确定 | 输入 | [sea_wave/RIS_float.m:199](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:199>) |
| `t` | 相对 parameters.startTime 的时间（s）或其采样向量 | 函数输入/输出；值随调用或上游结果确定 | 输入 | [sea_wave/RIS_float.m:199](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:199>) |
| `states` | 输出时刻的完整积分状态矩阵 | 函数输入/输出；值随调用或上游结果确定 | 输入 | [sea_wave/RIS_float.m:199](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:199>) |
| `wave` | 海浪分量结构，含 omega、kx、ky、coefficient；具体重名按下方作用域解释 | 函数输入/输出；值随调用或上游结果确定 | 输入 | [sea_wave/RIS_float.m:199](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:199>) |
| `spectrum` | 海浪频谱及离散积分诊断结构 | 函数输入/输出；值随调用或上游结果确定 | 输入 | [sea_wave/RIS_float.m:199](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:199>) |
| `referenceCenter` | 生成载荷时采用的初始世界水平中心（m） | 函数输入/输出；值随调用或上游结果确定 | 输入 | [sea_wave/RIS_float.m:199](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:199>) |
| `draft` | 静水平衡吃水 mass/(rho*Lx*Ly)（m） | 函数输入/输出；值随调用或上游结果确定 | 输入 | [sea_wave/RIS_float.m:199](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:199>) |
| `zG` | 静水平衡质心高度 thickness/2-draft（m） | 函数输入/输出；值随调用或上游结果确定 | 输入 | [sea_wave/RIS_float.m:199](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:199>) |
| `step` | 实际采用过的最大内部积分步长（s） | 函数输入/输出；值随调用或上游结果确定 | 输入 | [sea_wave/RIS_float.m:199](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:199>) |
| `out` | 当前函数输出结果结构 | `struct(...)（所有字段在本节逐行列明）` | 输出/赋值 | [sea_wave/RIS_float.m:199](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:199>) |
| `count` | 输出时间样本总数 | `numel(t)` | 赋值 | [sea_wave/RIS_float.m:201](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:201>) |
| `rotation` | 所有样本体轴到世界轴旋转矩阵（3×3×Nt） | `zeros(3,3,count)`<br>`R` | 赋值 | [sea_wave/RIS_float.m:202](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:202>) |
| `angles` | 所有样本欧拉角 [roll;pitch;yaw]（rad） | `zeros(3,count)`<br>`RIS_rotation(states(4:7,j))`<br>`unwrap(angles(3,:))` | 赋值 | [sea_wave/RIS_float.m:202](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:202>) |
| `vel` | 所有样本世界平移速度（3×Nt，m/s） | `zeros(3,count)`<br>`R*states(8:10,j)` | 赋值 | [sea_wave/RIS_float.m:202](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:202>) |
| `effectiveMass` | 所有样本体坐标有效广义质量矩阵（6×6×Nt） | `zeros(6,6,count)`<br>`load.effectiveMass` | 赋值 | [sea_wave/RIS_float.m:203](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:203>) |
| `acceleration` | 所有样本体速度导数（6×Nt，m/s² / rad/s²） | `zeros(6,count)`<br>`load.acceleration` | 赋值 | [sea_wave/RIS_float.m:204](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:204>) |
| `loads` | 波浪/静水/辐射/耦合四类广义力数组（6×Nt×4） | `zeros(6,count,4)`<br>`reshape([load.wave,load.hydrostatic,load.radiation,load.coriolis],6,1,4)` | 赋值 | [sea_wave/RIS_float.m:204](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:204>) |
| `residual` | 所有时间点牛顿—欧拉方程残差（N / N·m） | `zeros(6,count)`<br>`load.residual` | 赋值 | [sea_wave/RIS_float.m:204](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:204>) |
| `corners` | 上表面四顶点世界坐标（4×3×Nt，m） | `zeros(4,3,count)`<br>`(R*local+center(:,j)).'` | 赋值 | [sea_wave/RIS_float.m:205](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:205>) |
| `center` | 质心世界坐标 states(1:3,:)+[p.center;zG]（m） | `states(1:3,:)+[p.center;zG]` | 赋值 | [sea_wave/RIS_float.m:205](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:205>) |
| `local` | 上表面四顶点相对质心的体坐标（3×4，m） | `[-p.Lx/2,-p.Ly/2,p.thickness/2;p.Lx/2,-p.Ly/2,p.thickness/2; p.Lx/2,p.Ly/2,p.thickness/2;-p.Lx/2,p.Ly/2,p.thickness/2].'` | 赋值 | [sea_wave/RIS_float.m:206](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:206>) |
| `normal` | 当前板面或平面单位法向量 | `zeros(3,count)`<br>`R(:,3)` | 赋值 | [sea_wave/RIS_float.m:208](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:208>) |
| `boardCenter` | 上表面中心世界坐标（m） | `normal`<br>`center(:,j)+R(:,3)*p.thickness/2` | 赋值 | [sea_wave/RIS_float.m:208](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:208>) |
| `eulerRate` | 由体角速度经 T 变换得到的 ZYX 欧拉角速度（rad/s） | `normal`<br>`T*states(11:13,j)` | 赋值 | [sea_wave/RIS_float.m:208](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:208>) |
| `forces` | 所有输出时刻的参考轴波浪激励（N / N·m） | `real(m.forceHat*exp(1i*(m.omega*t)))` | 赋值 | [sea_wave/RIS_float.m:209](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:209>) |
| `energy` | 刚体动能+恢复势能+辐射振子能量（J） | `zeros(1,count)`<br>`0.5*states(8:13,j).'*load.effectiveMass*states(8:13,j)+potential+memoryEnergy` | 赋值 | [sea_wave/RIS_float.m:210](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:210>) |
| `orthogonality` | 每时刻旋转矩阵正交误差 Frobenius 范数 | `zeros(1,count)`<br>`norm(R.'*R-eye(3),'fro')` | 赋值 | [sea_wave/RIS_float.m:210](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:210>) |
| `power` | 波浪输入功率减辐射记忆耗散（W） | `zeros(1,count)`<br>`states(8:13,j).'*load.wave- sum(m.radiation.width.*m.radiation.residues.*(x.^2+y.^2),'all')` | 赋值 | [sea_wave/RIS_float.m:210](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:210>) |
| `j` | 当前循环、分量或样本索引（无量纲，MATLAB 从 1 开始） | `1:count`<br>`1:4` | 循环索引 | [sea_wave/RIS_float.m:211](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:211>) |
| `R` | 体坐标到世界坐标的主动旋转矩阵（无量纲） | `RIS_rotation(states(4:7,j))` | 赋值 | [sea_wave/RIS_float.m:212](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:212>) |
| `a` | 当前 roll/pitch 欧拉角（rad） | `angles(1,j)` | 赋值 | [sea_wave/RIS_float.m:218](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:218>) |
| `b` | 当前 roll/pitch 欧拉角（rad） | `angles(2,j)` | 赋值 | [sea_wave/RIS_float.m:218](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:218>) |
| `T` | 将体角速度转换为 ZYX 欧拉角速度的 3×3 矩阵 | `[1,sin(a)*tan(b),cos(a)*tan(b);0,cos(a),-sin(a); 0,sin(a)/cos(b),cos(a)/cos(b)]` | 赋值 | [sea_wave/RIS_float.m:219](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:219>) |
| `load` | 当前时刻的动力学受力与加速度诊断结构 | `derivative(states(:,j),m,forces(:,j))` | 赋值 | [sea_wave/RIS_float.m:222](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:222>) |
| `nr` | 每个垂向模态的辐射振子数 | `numel(m.radiation.width)` | 赋值 | [sea_wave/RIS_float.m:227](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:227>) |
| `x` | 当前样本的两组 Nr×3 实辐射记忆状态 | `reshape(states(14:13+3*nr,j),nr,3)` | 赋值 | [sea_wave/RIS_float.m:228](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:228>) |
| `y` | 当前样本的两组 Nr×3 实辐射记忆状态 | `reshape(states(14+3*nr:end,j),nr,3)` | 赋值 | [sea_wave/RIS_float.m:228](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:228>) |
| `memoryEnergy` | 辐射振子储能（J） | `0.5*sum(m.radiation.residues.*(x.^2+y.^2),'all')` | 赋值 | [sea_wave/RIS_float.m:229](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:229>) |
| `potential` | 升沉与倾斜的静水恢复势能（J） | `0.5*(m.restoring(1)*states(3,j)^2+ m.restoring(2)*R(3,2)^2+m.restoring(3)*R(3,1)^2)` | 赋值 | [sea_wave/RIS_float.m:230](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:230>) |
| `eta` | 标准六自由度状态 [surge;sway;heave;roll;pitch;yaw] | `[states(1:3,:);angles]` | 赋值 | [sea_wave/RIS_float.m:238](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:238>) |
| `order` | 旧 q 行序对应 eta 的索引 [3,5,4,1,2,6] | `[3,5,4,1,2,6]` | 赋值 | [sea_wave/RIS_float.m:238](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:238>) |
| `q` | 按旧行序重排的六自由度状态 | `eta(order,:)` | 赋值 | [sea_wave/RIS_float.m:239](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:239>) |
| `qd` | 按旧 q 行序的平移/欧拉角时间导数 | `[vel;eulerRate]`<br>`qd(order,:)` | 赋值 | [sea_wave/RIS_float.m:239](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:239>) |
| `field` | 提供给海面重建函数的 wave/parameters 封装，或当前块复波场 | `struct('wave',wave,'parameters',p)` | 赋值 | [sea_wave/RIS_float.m:240](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:240>) |
| `field.wave` | 外部共享的海浪分量结构；[] 时按参数生成 | `wave` | 结构字段 | [sea_wave/RIS_float.m:240](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:240>) |
| `field.parameters` | 当前计算实际使用的参数结构；EGG 中为混合分布参数 | `p` | 结构字段 | [sea_wave/RIS_float.m:240](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:240>) |
| `etaCenter` | 当前质心水平位置的入射波面高程（m） | `RIS_wave_field(field,center(1:2,:).',t).'` | 赋值 | [sea_wave/RIS_float.m:241](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:241>) |
| `tilt` | 板法线与世界竖直的倾角（°） | `acosd(max(-1,min(1,normal(3,:))))` | 赋值 | [sea_wave/RIS_float.m:242](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:242>) |
| `k` | 深水海浪波数 omega²/g（m⁻¹） | `wave.omega.^2/p.g` | 赋值 | [sea_wave/RIS_float.m:243](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:243>) |
| `orbitalWeight` | 各波浪分量的水平轨道响应平方权重（m²） | `sum(abs(m.forceHat(1:2,:)./(m.rigidMass(1)*wave.omega.'.^2)).^2,1).'` | 赋值 | [sea_wave/RIS_float.m:244](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:244>) |
| `shortFraction` | 短波条件 k*max(Lx,Ly)>pi 对应的轨道权重比例 | `sum(orbitalWeight(k*max(p.Lx,p.Ly)>pi))/max(sum(orbitalWeight),realmin)` | 赋值 | [sea_wave/RIS_float.m:245](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:245>) |
| `phaseGradient` | 按海浪复幅加权的 kx/omega、ky/omega 矩阵 | `[wave.kx./wave.omega,wave.ky./wave.omega].*abs(wave.coefficient)` | 赋值 | [sea_wave/RIS_float.m:247](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:247>) |
| `encounterCovariance` | 遭遇频率诊断使用的水平相位梯度二阶矩阵 | `phaseGradient.'*phaseGradient` | 赋值 | [sea_wave/RIS_float.m:248](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:248>) |
| `dopplerVariance` | 各时刻速度对应的遭遇频率变化方差指标 | `sum(vel(1:2,:).*(encounterCovariance*vel(1:2,:)),1)` | 赋值 | [sea_wave/RIS_float.m:249](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:249>) |
| `encounterRatio` | 平均遭遇频率变化相对指标（无量纲） | `sqrt(max(0,mean(dopplerVariance))/max(sum(abs(wave.coefficient).^2),realmin))` | 赋值 | [sea_wave/RIS_float.m:250](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:250>) |
| `mismatchByTime` | 每时刻四顶点与参考入射海面高差的最大值（m） | `zeros(1,count)`<br>`max(mismatchByTime,abs(sea-(reshape(corners(j,3,:),1,[])-(p.thickness-draft))))` | 赋值 | [sea_wave/RIS_float.m:251](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:251>) |
| `points` | 世界坐标采样点数组；波面接口为 N×2 的 [x,y]（m） | `reshape(corners(j,1:2,:),2,[]).'` | 赋值 | [sea_wave/RIS_float.m:253](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:253>) |
| `sea` | 各顶点水平位置对应的海浪高度时序（m） | `RIS_wave_field(field,points,t).'` | 赋值 | [sea_wave/RIS_float.m:254](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:254>) |
| `mismatch` | 所有顶点与所有时刻最大海面/板面参考高差（m） | `max(mismatchByTime)` | 赋值 | [sea_wave/RIS_float.m:257](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:257>) |
| `validity` | 适用范围与数值质量诊断结构 | `struct(...)（所有字段在本节逐行列明）` | 赋值 | [sea_wave/RIS_float.m:258](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:258>) |
| `validity.maxTiltDeg` | 全过程最大法向倾角（°） | `max(tilt)` | 结构字段 | [sea_wave/RIS_float.m:258](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:258>) |
| `validity.smallAngleUnder10Deg` | 全过程倾角是否不超过 10° | `all(tilt<=10)` | 结构字段 | [sea_wave/RIS_float.m:258](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:258>) |
| `validity.maxHeadingChangeDeg` | 相对初始航向的最大偏航变化绝对值（°） | `max(abs(rad2deg(angles(3,:)-angles(3,1))))` | 结构字段 | [sea_wave/RIS_float.m:259](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:259>) |
| `validity.hydrodynamicsComputedWithBEM` | 是否在本次模型中使用薄板 BEM 水动力 | `true` | 结构字段 | [sea_wave/RIS_float.m:260](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:260>) |
| `validity.usesExternalHydrodynamicDatabase` | 是否使用外部水动力数据库（源码 false） | `false` | 结构字段 | [sea_wave/RIS_float.m:260](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:260>) |
| `validity.verticalRadiationDiffractionSolved` | 是否求解垂向辐射/绕射（true） | `true` | 结构字段 | [sea_wave/RIS_float.m:261](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:261>) |
| `validity.radiationMemoryUsesFit` | 辐射记忆是否来自被动有理拟合（true） | `true` | 结构字段 | [sea_wave/RIS_float.m:261](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:261>) |
| `validity.horizontalDiffractionSolved` | 对应水平绕射/偏航绕射/水平与偏航辐射/完整二阶漂移/直接风载/实测标定/本次网格收敛/水动力验证开关，源码均为 false；字段名区分项目 | `false` | 结构字段 | [sea_wave/RIS_float.m:262](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:262>) |
| `validity.yawDiffractionSolved` | 对应水平绕射/偏航绕射/水平与偏航辐射/完整二阶漂移/直接风载/实测标定/本次网格收敛/水动力验证开关，源码均为 false；字段名区分项目 | `false` | 结构字段 | [sea_wave/RIS_float.m:262](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:262>) |
| `validity.horizontalYawRadiationSolved` | 对应水平绕射/偏航绕射/水平与偏航辐射/完整二阶漂移/直接风载/实测标定/本次网格收敛/水动力验证开关，源码均为 false；字段名区分项目 | `false` | 结构字段 | [sea_wave/RIS_float.m:263](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:263>) |
| `validity.fullSecondOrderDriftSolved` | 对应水平绕射/偏航绕射/水平与偏航辐射/完整二阶漂移/直接风载/实测标定/本次网格收敛/水动力验证开关，源码均为 false；字段名区分项目 | `false` | 结构字段 | [sea_wave/RIS_float.m:263](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:263>) |
| `validity.directWindLoadIncluded` | 对应水平绕射/偏航绕射/水平与偏航辐射/完整二阶漂移/直接风载/实测标定/本次网格收敛/水动力验证开关，源码均为 false；字段名区分项目 | `false` | 结构字段 | [sea_wave/RIS_float.m:264](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:264>) |
| `validity.experimentallyCalibrated` | 对应水平绕射/偏航绕射/水平与偏航辐射/完整二阶漂移/直接风载/实测标定/本次网格收敛/水动力验证开关，源码均为 false；字段名区分项目 | `false` | 结构字段 | [sea_wave/RIS_float.m:264](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:264>) |
| `validity.horizontalShortWaveVarianceFraction` | 水平轨道响应权重中 k*max(Lx,Ly)>pi 的短波方差比例 | `shortFraction` | 结构字段 | [sea_wave/RIS_float.m:265](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:265>) |
| `validity.rmsEncounterFrequencyChange` | 按轨迹速度计算的遭遇频率变化均方根相对指标 | `encounterRatio` | 结构字段 | [sea_wave/RIS_float.m:266](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:266>) |
| `validity.maxIncidentWaveBoardMismatchM` | 全过程最大入射波面与板面吃水参考高差（m） | `mismatch` | 结构字段 | [sea_wave/RIS_float.m:267](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:267>) |
| `validity.incidentMismatchExceedsDraft` | 上述最大高差是否超过平衡吃水 | `mismatch>draft` | 结构字段 | [sea_wave/RIS_float.m:267](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:267>) |
| `validity.numericalSolutionFinite` | 所有积分状态是否均为有限值 | `all(isfinite(states(:)))` | 结构字段 | [sea_wave/RIS_float.m:268](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:268>) |
| `validity.maxRotationOrthogonalityError` | 旋转矩阵 RᵀR-I 的最大 Frobenius 范数 | `max(orthogonality)` | 结构字段 | [sea_wave/RIS_float.m:269](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:269>) |
| `validity.maxEquationResidual` | 全过程最大广义方程残差绝对值 | `max(abs(residual),[],'all')` | 结构字段 | [sea_wave/RIS_float.m:270](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:270>) |
| `validity.radiationFitRelativeRmsError` | 按频响幅值加权的阻抗拟合相对均方根误差（无量纲，三个垂向模态） | `m.radiation.relativeImpedanceRmsError` | 结构字段 | [sea_wave/RIS_float.m:271](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:271>) |
| `validity.referenceLinearHydrodynamics` | 是否采用初始参考轴的线性水动力（true） | `true` | 结构字段 | [sea_wave/RIS_float.m:272](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:272>) |
| `validity.meshConvergenceCheckedThisCall` | 对应水平绕射/偏航绕射/水平与偏航辐射/完整二阶漂移/直接风载/实测标定/本次网格收敛/水动力验证开关，源码均为 false；字段名区分项目 | `false` | 结构字段 | [sea_wave/RIS_float.m:272](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:272>) |
| `out.t` | 相对时间采样向量（s） | `t` | 结构字段 | [sea_wave/RIS_float.m:273](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:273>) |
| `out.absoluteTime` | 绝对海浪时间 startTime+t（s） | `p.startTime+t` | 结构字段 | [sea_wave/RIS_float.m:273](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:273>) |
| `out.eta` | 六自由度位移/姿态 [surge;sway;heave;roll;pitch;yaw]（m / rad） | `eta` | 结构字段 | [sea_wave/RIS_float.m:273](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:273>) |
| `out.q` | 兼容旧行序的六自由度状态 [heave;pitch;roll;surge;sway;yaw]（m / rad） | `q` | 结构字段 | [sea_wave/RIS_float.m:273](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:273>) |
| `out.dofOrder` | 对应 eta/q 数组的自由度名称和顺序 | `{{'surge','sway','heave','roll','pitch','yaw'}}` | 结构字段 | [sea_wave/RIS_float.m:274](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:274>) |
| `out.qOrder` | 对应 eta/q 数组的自由度名称和顺序 | `{{'heave','pitch','roll','surge','sway','yaw'}}` | 结构字段 | [sea_wave/RIS_float.m:275](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:275>) |
| `out.surge` | 世界坐标的纵荡、横荡、升沉状态（m），水平状态相对 center | `q(4,:)` | 结构字段 | [sea_wave/RIS_float.m:274](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:274>) |
| `out.sway` | 世界坐标的纵荡、横荡、升沉状态（m），水平状态相对 center | `q(5,:)` | 结构字段 | [sea_wave/RIS_float.m:274](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:274>) |
| `out.heave` | 世界坐标的纵荡、横荡、升沉状态（m），水平状态相对 center | `q(1,:)` | 结构字段 | [sea_wave/RIS_float.m:274](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:274>) |
| `out.roll` | 右手 ZYX 横摇、纵摇、偏航（rad）；yaw 时序已 unwrap | `q(3,:)` | 结构字段 | [sea_wave/RIS_float.m:274](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:274>) |
| `out.pitch` | 右手 ZYX 横摇、纵摇、偏航（rad）；yaw 时序已 unwrap | `q(2,:)` | 结构字段 | [sea_wave/RIS_float.m:274](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:274>) |
| `out.yaw` | 右手 ZYX 横摇、纵摇、偏航（rad）；yaw 时序已 unwrap | `q(6,:)` | 结构字段 | [sea_wave/RIS_float.m:277](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:277>) |
| `out.rollDeg` | 对应欧拉角的度数（°） | `rad2deg(q(3,:))` | 结构字段 | [sea_wave/RIS_float.m:278](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:278>) |
| `out.pitchDeg` | 对应欧拉角的度数（°） | `rad2deg(q(2,:))` | 结构字段 | [sea_wave/RIS_float.m:278](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:278>) |
| `out.yawDeg` | 对应欧拉角的度数（°） | `rad2deg(q(6,:))` | 结构字段 | [sea_wave/RIS_float.m:278](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:278>) |
| `out.velocity` | 与 q 相同行序的导数；平移为世界速度，转动为欧拉角速度（m/s / rad/s） | `qd` | 结构字段 | [sea_wave/RIS_float.m:279](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:279>) |
| `out.bodyVelocity` | 体坐标六速度 [u;v;w;p;q;r]（m/s / rad/s） | `states(8:13,:)` | 结构字段 | [sea_wave/RIS_float.m:279](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:279>) |
| `out.bodyAcceleration` | 体坐标六速度的时间导数（m/s² / rad/s²） | `acceleration` | 结构字段 | [sea_wave/RIS_float.m:279](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:279>) |
| `out.rotationMatrix` | 体坐标到世界坐标的 3×3×Nt 主动旋转矩阵 | `rotation` | 结构字段 | [sea_wave/RIS_float.m:280](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:280>) |
| `out.quaternion` | 姿态 wxyz 单位四元数时序（4×Nt） | `states(4:7,:)` | 结构字段 | [sea_wave/RIS_float.m:280](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:280>) |
| `out.normal` | 板体 +z 轴的世界单位法向量（3×Nt） | `normal` | 结构字段 | [sea_wave/RIS_float.m:280](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:280>) |
| `out.corners` | RIS 上表面四顶点的世界坐标（4×3×Nt，m） | `corners` | 结构字段 | [sea_wave/RIS_float.m:281](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:281>) |
| `out.boardCenter` | RIS 上表面中心的世界坐标（m）；W_RIS 接收面使用下表面 | `boardCenter` | 结构字段 | [sea_wave/RIS_float.m:281](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:281>) |
| `out.centerOfMass` | RIS 质心世界坐标（m），含静水平衡 zG 与 center 偏置 | `center` | 结构字段 | [sea_wave/RIS_float.m:281](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:281>) |
| `out.centerTrajectory` | 质心的世界水平轨迹（2×Nt，m） | `center(1:2,:)` | 结构字段 | [sea_wave/RIS_float.m:282](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:282>) |
| `out.horizontalVelocity` | 质心的世界水平速度（2×Nt，m/s） | `vel(1:2,:)` | 结构字段 | [sea_wave/RIS_float.m:282](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:282>) |
| `out.waveLoadBody` | 对应波浪/静水/辐射/耦合广义力的体坐标时序（N / N·m），名称决定分项 | `loads(:,:,1)` | 结构字段 | [sea_wave/RIS_float.m:283](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:283>) |
| `out.hydrostaticLoadBody` | 对应波浪/静水/辐射/耦合广义力的体坐标时序（N / N·m），名称决定分项 | `loads(:,:,2)` | 结构字段 | [sea_wave/RIS_float.m:283](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:283>) |
| `out.radiationLoadBody` | 对应波浪/静水/辐射/耦合广义力的体坐标时序（N / N·m），名称决定分项 | `loads(:,:,3)` | 结构字段 | [sea_wave/RIS_float.m:284](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:284>) |
| `out.coriolisLoadBody` | 对应波浪/静水/辐射/耦合广义力的体坐标时序（N / N·m），名称决定分项 | `loads(:,:,4)` | 结构字段 | [sea_wave/RIS_float.m:284](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:284>) |
| `out.waveLoadReference` | 初始航向参考轴中的六维波浪激励（N / N·m） | `forces` | 结构字段 | [sea_wave/RIS_float.m:285](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:285>) |
| `out.horizontalWaveForce` | 世界 x/y 分量的波浪水平力（N） | `m.headingR(1:2,1:2)*forces(1:2,:)` | 结构字段 | [sea_wave/RIS_float.m:285](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:285>) |
| `out.yawWaveMoment` | 参考系 z 轴波浪偏航力矩（N·m） | `forces(6,:)` | 结构字段 | [sea_wave/RIS_float.m:286](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:286>) |
| `out.equationResidual` | 动力学方程残差时序（N / N·m） | `residual` | 结构字段 | [sea_wave/RIS_float.m:286](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:286>) |
| `out.energy` | 动能、静水势能和辐射记忆能之和（J） | `energy` | 结构字段 | [sea_wave/RIS_float.m:287](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:287>) |
| `out.netPower` | 波浪输入功率减辐射记忆耗散功率（W） | `power` | 结构字段 | [sea_wave/RIS_float.m:287](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:287>) |
| `out.integrationStepUsed` | 本次积分实际使用过的最大内部步长（s） | `step` | 结构字段 | [sea_wave/RIS_float.m:287](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:287>) |
| `out.mass` | 旧三模态 [heave;pitch;roll] 刚体质量/转动惯量向量 | `m.rigidMass([3,5,4])` | 结构字段 | [sea_wave/RIS_float.m:288](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:288>) |
| `out.massMatrix` | 按旧 q 行序排列的刚体广义质量对角矩阵 | `diag(m.rigidMass(order))` | 结构字段 | [sea_wave/RIS_float.m:288](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:288>) |
| `out.bodyMassMatrix` | 随姿态变化的有效体坐标六维质量矩阵时序 | `effectiveMass` | 结构字段 | [sea_wave/RIS_float.m:289](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:289>) |
| `out.rigidBodyMassMatrix` | 按 [surge;sway;heave;roll;pitch;yaw] 排列的刚体质量矩阵 | `diag(m.rigidMass)` | 结构字段 | [sea_wave/RIS_float.m:289](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:289>) |
| `out.stiffness` | 按旧 q 行序的静水恢复刚度向量，水平/偏航项为 0 | `[m.restoring([1,3,2]);0;0;0]` | 结构字段 | [sea_wave/RIS_float.m:290](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:290>) |
| `out.stiffnessMatrix` | 由 stiffness 构造的对角矩阵 | `diag([m.restoring([1,3,2]);0;0;0])` | 结构字段 | [sea_wave/RIS_float.m:291](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:291>) |
| `out.radiationModel` | 重排为 [heave;roll;pitch] 的辐射记忆模型结构 | `m.radiation` | 结构字段 | [sea_wave/RIS_float.m:292](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:292>) |
| `out.etaCenter` | 运动质心水平位置处的真实入射波面高度（m） | `etaCenter` | 结构字段 | [sea_wave/RIS_float.m:292](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:292>) |
| `out.incidentWaveBoardMismatchM` | 各时刻四个上表面顶点扣除干舷后与入射波面的最大高差（m） | `mismatchByTime` | 结构字段 | [sea_wave/RIS_float.m:293](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:293>) |
| `out.etaReferenceCenter` | 固定参考中心处的入射波面高度时序（m） | `real((wave.coefficient.*exp(1i*(wave.omega*p.startTime- wave.kx*referenceCenter(1)-wave.ky*referenceCenter(2)))).'*exp(1i*(wave.omega*t)))` | 结构字段 | [sea_wave/RIS_float.m:294](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:294>) |
| `out.referenceCenter` | 生成载荷时采用的初始世界水平中心（m） | `referenceCenter` | 结构字段 | [sea_wave/RIS_float.m:296](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:296>) |
| `out.referenceHeading` | 线性水动力参考航向（rad） | `atan2(m.headingR(2,1),m.headingR(1,1))` | 结构字段 | [sea_wave/RIS_float.m:296](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:296>) |
| `out.wave` | 外部共享的海浪分量结构；[] 时按参数生成 | `wave` | 结构字段 | [sea_wave/RIS_float.m:297](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:297>) |
| `out.waveFrequencyHz` | 实际入射海浪分量频率 omega/(2*pi)（Hz） | `wave.omega/(2*pi)` | 结构字段 | [sea_wave/RIS_float.m:297](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:297>) |
| `out.spectrum` | 频谱积分、波高和斜率诊断结构 | `spectrum` | 结构字段 | [sea_wave/RIS_float.m:297](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:297>) |
| `out.parameters` | 当前计算实际使用的参数结构；EGG 中为混合分布参数 | `p` | 结构字段 | [sea_wave/RIS_float.m:297](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:297>) |
| `out.draft` | 静水平衡吃水 mass/(rho*Lx*Ly)（m） | `draft` | 结构字段 | [sea_wave/RIS_float.m:298](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:298>) |
| `out.equilibriumCenterOfMassZ` | 静水平衡质心世界高度（m） | `zG` | 结构字段 | [sea_wave/RIS_float.m:298](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:298>) |
| `out.validity` | 水动力适用范围与数值误差诊断结构 | `validity` | 结构字段 | [sea_wave/RIS_float.m:298](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:298>) |
| `out.hydroMode` | 水动力模型说明文本 | `'6DOF Newton-Euler + reference-axis vertical BEM memory + sidewall FK'` | 结构字段 | [sea_wave/RIS_float.m:299](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:299>) |


**作用域：`windowSummary`**

| 源码符号 | 含义与单位/维度 | 源码取值/计算式（各次赋值） | 角色 | 源码位置 |
| --- | --- | --- | --- | --- |
| `o` | 当前 RIS 运动结果或动力学模型；具体类型由所在函数输入确定 | 函数输入/输出；值随调用或上游结果确定 | 输入 | [sea_wave/RIS_float.m:301](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:301>) |
| `windows` | 待汇总的相对时间窗口终点，主函数传 [5,10]（s） | 函数输入/输出；值随调用或上游结果确定 | 输入 | [sea_wave/RIS_float.m:301](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:301>) |
| `summaries` | 各时间窗口的终点状态、倾角和水平位移结构数组 | `struct('duration',{},'finalQ',{},'maxTiltDeg',{},'horizontalDisplacementM',{})`<br>`struct(...)（所有字段在本节逐行列明）` | 输出/赋值 | [sea_wave/RIS_float.m:301](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:301>) |
| `summaries.duration` | 保存的运动持续时间（s） | `{}`<br>`window` | 结构字段 | [sea_wave/RIS_float.m:303](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:303>) |
| `summaries.finalQ` | 窗口结束时按旧 q 行序的六自由度状态（m / rad） | `{}`<br>`qend` | 结构字段 | [sea_wave/RIS_float.m:303](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:303>) |
| `summaries.maxTiltDeg` | 全过程最大法向倾角（°） | `{}`<br>`max(acosd(max(-1,min(1,o.normal(3,ids)))))` | 结构字段 | [sea_wave/RIS_float.m:303](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:303>) |
| `summaries.horizontalDisplacementM` | 窗口内起终点水平位移模长（m） | `{}`<br>`norm(qend(4:5)-o.q(4:5,1))` | 结构字段 | [sea_wave/RIS_float.m:303](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:303>) |
| `window` | 当前窗口终点（s） | `windows(windows<=o.t(end))` | 循环索引 | [sea_wave/RIS_float.m:304](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:304>) |
| `qend` | 在窗口终点以 pchip 插值得到的旧行序状态（m / rad） | `interp1(o.t.',o.q.',window,'pchip').'` | 赋值 | [sea_wave/RIS_float.m:305](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:305>) |
| `ids` | 当前窗口内的输出样本布尔掩码 | `o.t<=window` | 赋值 | [sea_wave/RIS_float.m:306](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:306>) |


**作用域：`plotMotion`**

| 源码符号 | 含义与单位/维度 | 源码取值/计算式（各次赋值） | 角色 | 源码位置 |
| --- | --- | --- | --- | --- |
| `o` | 当前 RIS 运动结果或动力学模型；具体类型由所在函数输入确定 | 函数输入/输出；值随调用或上游结果确定 | 输入 | [sea_wave/RIS_float.m:312](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:312>) |


**作用域：`validateParameters`**

| 源码符号 | 含义与单位/维度 | 源码取值/计算式（各次赋值） | 角色 | 源码位置 |
| --- | --- | --- | --- | --- |
| `p` | 当前作用域的参数结构 | `p.(names{j})(:)`<br>`p.(pairs{j})(:)` | 输入/输出/赋值 | [sea_wave/RIS_float.m:328](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:328>) |
| `positive` | 必须为正、有限、实标量的参数字段名列表 | `{'Lx','Ly','mass','thickness','rho','g','duration','dt','fMinHz','fMaxHz','integrationStep'}` | 赋值 | [sea_wave/RIS_float.m:330](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:330>) |
| `j` | 当前循环、分量或样本索引（无量纲，MATLAB 从 1 开始） | `1:numel(positive)`<br>`1:3`<br>`1:numel(pairs)` | 循环索引 | [sea_wave/RIS_float.m:331](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:331>) |
| `sizes` | 可选初始位置/欧拉角/体速度的长度 [3,3,6] | `[3,3,6]` | 赋值 | [sea_wave/RIS_float.m:344](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:344>) |
| `names` | 可选三类初始状态的字段名 | `{'initialPosition','initialEulerAngles','initialBodyVelocity'}` | 赋值 | [sea_wave/RIS_float.m:344](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:344>) |
| `pairs` | 水平中心、初始位移、初始速度三个二维字段名 | `{'center','horizontalInitialPosition','horizontalInitialVelocity'}` | 赋值 | [sea_wave/RIS_float.m:353](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:353>) |


**作用域：`makeWave`**

| 源码符号 | 含义与单位/维度 | 源码取值/计算式（各次赋值） | 角色 | 源码位置 |
| --- | --- | --- | --- | --- |
| `p` | 当前作用域的参数结构 | 函数输入/输出；值随调用或上游结果确定 | 输入 | [sea_wave/RIS_float.m:363](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:363>) |
| `w` | 外部传入或新生成的海浪分量结构 | `p.wave`<br>`w.(required{j})(:)`<br>`RIS_wave_components(p)` | 输出/赋值 | [sea_wave/RIS_float.m:363](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:363>) |
| `s` | 海浪频谱和方差/显著波高/斜率诊断结构 | `struct('name','supplied components','retainedFullSpectrumVariance',NaN)`<br>`RIS_wave_components(p)` | 输出/赋值 | [sea_wave/RIS_float.m:363](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:363>) |
| `required` | 外部波浪必须具备的字段 omega/kx/ky/coefficient | `{'omega','kx','ky','coefficient'}` | 赋值 | [sea_wave/RIS_float.m:367](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:367>) |
| `j` | 当前循环、分量或样本索引（无量纲，MATLAB 从 1 开始） | `1:4` | 循环索引 | [sea_wave/RIS_float.m:371](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:371>) |
| `s.name` | 模型名称或输出名称文本 | `'supplied components'` | 结构字段 | [sea_wave/RIS_float.m:382](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:382>) |
| `s.retainedFullSpectrumVariance` | 所取频带保留的完整连续 Neumann 谱方差比例（无量纲） | `NaN` | 结构字段 | [sea_wave/RIS_float.m:382](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:382>) |
| `s.waveVarianceM2` | 离散海面方差 sum(abs(coefficient).^2)/2（m²） | `sum(abs(w.coefficient).^2)/2` | 赋值 | [sea_wave/RIS_float.m:389](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:389>) |
| `s.significantWaveHeightM` | 显著波高 4*sqrt(waveVarianceM2)（m） | `4*sqrt(s.waveVarianceM2)` | 赋值 | [sea_wave/RIS_float.m:390](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:390>) |
| `s.rmsIncidentSlope` | 入射海面空间总斜率均方根（无量纲） | `sqrt(sum(abs(w.coefficient).^2.*(w.omega.^2/p.g).^2)/2)` | 赋值 | [sea_wave/RIS_float.m:391](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:391>) |


**作用域：`plateResponse`**

| 源码符号 | 含义与单位/维度 | 源码取值/计算式（各次赋值） | 角色 | 源码位置 |
| --- | --- | --- | --- | --- |
| `p` | 当前作用域的参数结构 | 函数输入/输出；值随调用或上游结果确定 | 输入 | [sea_wave/RIS_float.m:393](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:393>) |
| `omega` | 合并后的 BEM 计算角频率（rad/s） | 函数输入/输出；值随调用或上游结果确定 | 输入 | [sea_wave/RIS_float.m:393](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:393>) |
| `kx` | 初始航向体坐标海浪波数分量（m⁻¹） | 函数输入/输出；值随调用或上游结果确定 | 输入 | [sea_wave/RIS_float.m:393](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:393>) |
| `ky` | 初始航向体坐标海浪波数分量（m⁻¹） | 函数输入/输出；值随调用或上游结果确定 | 输入 | [sea_wave/RIS_float.m:393](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:393>) |
| `index` | 每个入射波分量在合并 omega 中的位置 | 函数输入/输出；值随调用或上游结果确定 | 输入 | [sea_wave/RIS_float.m:393](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:393>) |
| `M` | 旧三模态刚体广义质量/惯量 | 函数输入/输出；值随调用或上游结果确定 | 输入 | [sea_wave/RIS_float.m:393](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:393>) |
| `K` | 旧三模态静水恢复刚度 | 函数输入/输出；值随调用或上游结果确定 | 输入 | [sea_wave/RIS_float.m:393](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:393>) |
| `rao` | 旧三模态单位波幅频率响应算子 | `F./Z(:,index)` | 输出/赋值 | [sea_wave/RIS_float.m:393](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:393>) |
| `A` | 旧三模态 BEM 附加广义质量频响 | `entry.mass` | 输出/赋值 | [sea_wave/RIS_float.m:393](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:393>) |
| `B` | 旧三模态 BEM 辐射广义阻尼频响 | `entry.damping` | 输出/赋值 | [sea_wave/RIS_float.m:393](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:393>) |
| `F` | 各入射波对应的旧三模态单位波幅激励复幅 | `complex(zeros(3,numel(kx)))`<br>`p.rho*p.g*4*delta*[r(:,1).'* (cx.*cy); r(:,2).'* (sx.*cy); r(:,3).'* (cx.*sy)]` | 输出/赋值 | [sea_wave/RIS_float.m:393](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:393>) |
| `keys` | 最多保留 4 套几何/频率输入的持久缓存键集合 | `{}`<br>`[]`<br>`key` | 持久缓存/赋值 | [sea_wave/RIS_float.m:396](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:396>) |
| `entries` | 与 keys 对应的势函数、附加质量、阻尼和配点持久缓存 | `{}`<br>`[]`<br>`struct('response',response,'mass',mass,'damping',damping,'points',points)` | 持久缓存/赋值 | [sea_wave/RIS_float.m:396](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:396>) |
| `key` | 由 Lx/Ly/rho/g/panelCount/omega 组成的当前 BEM 缓存键 | `{p.Lx,p.Ly,p.rho,p.g,p.panelCount,omega}` | 赋值 | [sea_wave/RIS_float.m:398](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:398>) |
| `hit` | 命中缓存的索引；空表示需要重新求解 | `find(cellfun(@(item) isequal(item,key),keys),1)`<br>`numel(keys)` | 赋值 | [sea_wave/RIS_float.m:399](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:399>) |
| `item` | cellfun 匿名函数中的历史缓存键 | `hit = find(cellfun(@(item) isequal(item,key),keys),1)` | 匿名函数输入 | [sea_wave/RIS_float.m:399](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:399>) |
| `response` | 三个奇偶模态的面元势解频响 | `solvePotentials(p,omega)` | 赋值 | [sea_wave/RIS_float.m:401](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:401>) |
| `mass` | 旧三模态 BEM 附加广义质量频响 | `solvePotentials(p,omega)` | 赋值 | [sea_wave/RIS_float.m:401](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:401>) |
| `damping` | 旧三模态 BEM 辐射广义阻尼频响 | `solvePotentials(p,omega)` | 赋值 | [sea_wave/RIS_float.m:401](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:401>) |
| `points` | 一个象限的 BEM 面元中心坐标（m） | `solvePotentials(p,omega)` | 赋值 | [sea_wave/RIS_float.m:401](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:401>) |
| `entries.response` | 缓存的三个奇偶模态的面元势函数频响 | `response` | 结构字段 | [sea_wave/RIS_float.m:404](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:404>) |
| `entries.mass` | 缓存的旧 BEM 三模态附加广义质量频响（升沉 kg，转动 kg·m²） | `mass` | 结构字段 | [sea_wave/RIS_float.m:404](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:404>) |
| `entries.damping` | 缓存的三模态 BEM 辐射广义阻尼频响 | `damping` | 结构字段 | [sea_wave/RIS_float.m:404](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:404>) |
| `entries.points` | 缓存的一个象限 BEM 面元中心坐标（m） | `points` | 结构字段 | [sea_wave/RIS_float.m:404](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:404>) |
| `entry` | 当前选中的 BEM 缓存数据结构 | `entries{hit}` | 赋值 | [sea_wave/RIS_float.m:407](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:407>) |
| `xy` | 一个象限的 BEM 面元中心坐标（m） | `entry.points` | 赋值 | [sea_wave/RIS_float.m:410](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:410>) |
| `delta` | 完整单个面元面积 Lx*Ly/panelCount²（m²） | `p.Lx*p.Ly/p.panelCount^2` | 赋值 | [sea_wave/RIS_float.m:411](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:411>) |
| `j` | 当前循环、分量或样本索引（无量纲，MATLAB 从 1 开始） | `1:numel(omega)` | 循环索引 | [sea_wave/RIS_float.m:413](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:413>) |
| `ids` | 当前 BEM 频率对应的海浪分量索引 | `find(index==j)` | 赋值 | [sea_wave/RIS_float.m:414](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:414>) |
| `cx` | 面元上的 x/y 偶对称入射相位 cos(k*x/y) | `cos(xy(:,1)*kx(ids).')` | 赋值 | [sea_wave/RIS_float.m:416](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:416>) |
| `sx` | 面元上的 x/y 奇对称入射相位 -1i*sin(k*x/y) | `-1i*sin(xy(:,1)*kx(ids).')` | 赋值 | [sea_wave/RIS_float.m:416](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:416>) |
| `cy` | 面元上的 x/y 偶对称入射相位 cos(k*x/y) | `cos(xy(:,2)*ky(ids).')` | 赋值 | [sea_wave/RIS_float.m:417](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:417>) |
| `sy` | 面元上的 x/y 奇对称入射相位 -1i*sin(k*x/y) | `-1i*sin(xy(:,2)*ky(ids).')` | 赋值 | [sea_wave/RIS_float.m:417](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:417>) |
| `r` | 当前频率的三个面元势函数响应列 | `entry.response(:,:,j)` | 赋值 | [sea_wave/RIS_float.m:418](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:418>) |
| `Z` | 旧三模态动态阻抗 K-omega²*(M+A)+1i*omega*B | `K-omega.'.^2.*(M+A)+1i*omega.'.*B` | 赋值 | [sea_wave/RIS_float.m:422](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:422>) |


**作用域：`solvePotentials`**

| 源码符号 | 含义与单位/维度 | 源码取值/计算式（各次赋值） | 角色 | 源码位置 |
| --- | --- | --- | --- | --- |
| `p` | 当前作用域的参数结构 | 函数输入/输出；值随调用或上游结果确定 | 输入 | [sea_wave/RIS_float.m:425](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:425>) |
| `omega` | 要求解的角频率向量（rad/s） | 函数输入/输出；值随调用或上游结果确定 | 输入 | [sea_wave/RIS_float.m:425](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:425>) |
| `R` | 三个奇偶扇区的面元势响应数组；此处不是旋转矩阵 | `complex(zeros(size(P,1),3,numel(omega)))`<br>`rr` | 输出/赋值 | [sea_wave/RIS_float.m:425](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:425>) |
| `A` | 三个旧 BEM 模态的附加广义质量频响 | `zeros(3,numel(omega))`<br>`real(T)` | 输出/赋值 | [sea_wave/RIS_float.m:425](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:425>) |
| `B` | 三个旧 BEM 模态的辐射广义阻尼频响 | `A`<br>`-omega(f)*imag(T)` | 输出/赋值 | [sea_wave/RIS_float.m:425](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:425>) |
| `xy` | 一个象限中 (n/2)² 个面元中心的 [x,y] 坐标（m） | `[(ix(:)+0.5)*dx,(iy(:)+0.5)*dy]` | 输出/赋值 | [sea_wave/RIS_float.m:425](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:425>) |
| `quadrature` | 16 个二维高斯面积积分权重（m²） | `dx*dy*weights(i)*weights(j)/4;`<br>`zeros(1,16)`<br>`dx*dy*weights(i)*weights(j)/4` | 声明默认值/赋值 | [sea_wave/RIS_float.m:440](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:440>) |
| `n` | 完整矩形每边面元数 p.panelCount | `p.panelCount` | 赋值 | [sea_wave/RIS_float.m:427](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:427>) |
| `dx` | 单个面元沿 x/y 的边长（m） | `p.Lx/n` | 赋值 | [sea_wave/RIS_float.m:428](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:428>) |
| `dy` | 单个面元沿 x/y 的边长（m） | `p.Ly/n` | 赋值 | [sea_wave/RIS_float.m:428](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:428>) |
| `ix` | 一个象限内的面元整数索引网格 | `ndgrid(0:n/2-1,0:n/2-1)` | 赋值 | [sea_wave/RIS_float.m:429](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:429>) |
| `iy` | 一个象限内的面元整数索引网格 | `ndgrid(0:n/2-1,0:n/2-1)` | 赋值 | [sea_wave/RIS_float.m:429](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:429>) |
| `ox` | 完整网格的 x/y 偏移矩阵（m） | `ndgrid((0:n-1)*dx,(0:n-1)*dy)` | 赋值 | [sea_wave/RIS_float.m:431](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:431>) |
| `oy` | 完整网格的 x/y 偏移矩阵（m） | `ndgrid((0:n-1)*dx,(0:n-1)*dy)` | 赋值 | [sea_wave/RIS_float.m:431](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:431>) |
| `singular` | 自由面 Green 函数奇异核的解析面元积分表 | `zeros(n)`<br>`singular+sx*sy*(xx.*log(yy+rr)+yy.*log(xx+rr))/(2*pi)` | 赋值 | [sea_wave/RIS_float.m:432](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:432>) |
| `sx` | 镜像/面元角点积分使用的 ±1 符号；奇偶扇区也用该符号 | `[-1,1]`<br>`parity(mode,1)` | 循环索引/赋值 | [sea_wave/RIS_float.m:433](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:433>) |
| `sy` | 镜像/面元角点积分使用的 ±1 符号；奇偶扇区也用该符号 | `[-1,1]`<br>`parity(mode,2)` | 循环索引/赋值 | [sea_wave/RIS_float.m:434](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:434>) |
| `xx` | 构造奇异核时的坐标偏移（m），或镜像积分查表的整数索引 | `ox+sx*dx/2`<br>`round(abs(xy(:,1)-mirrors(j,1)*xy(:,1).')/dx)+1` | 赋值 | [sea_wave/RIS_float.m:435](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:435>) |
| `yy` | 构造奇异核时的坐标偏移（m），或镜像积分查表的整数索引 | `oy+sy*dy/2`<br>`round(abs(xy(:,2)-mirrors(j,2)*xy(:,2).')/dy)+1` | 赋值 | [sea_wave/RIS_float.m:435](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:435>) |
| `rr` | 积分点径向距离（m）；模态求解阶段重用为面元势解向量 | `hypot(xx,yy)`<br>`(eye(size(H))+k*H)\P(:,mode)` | 赋值 | [sea_wave/RIS_float.m:435](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:435>) |
| `nodes` | 4 点 Gauss–Legendre 一维节点/权重 | `gaussLegendre(4)` | 赋值 | [sea_wave/RIS_float.m:439](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:439>) |
| `weights` | 4 点 Gauss–Legendre 一维节点/权重 | `gaussLegendre(4)` | 赋值 | [sea_wave/RIS_float.m:439](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:439>) |
| `radius` | 各面元 4×4 高斯点的径向距离（n²×16，m） | `zeros(n*n,16)`<br>`reshape(hypot(ox+nodes(i)*dx/2,oy+nodes(j)*dy/2),[],1)` | 赋值 | [sea_wave/RIS_float.m:440](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:440>) |
| `count` | 当前二维高斯积分点序号 1..16 | `0`<br>`count+1` | 赋值 | [sea_wave/RIS_float.m:441](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:441>) |
| `i` | 当前循环、分量或样本索引（无量纲，MATLAB 从 1 开始） | `1:4` | 循环索引 | [sea_wave/RIS_float.m:442](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:442>) |
| `j` | 当前循环、分量或样本索引（无量纲，MATLAB 从 1 开始） | `1:4` | 循环索引 | [sea_wave/RIS_float.m:443](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:443>) |
| `mirrors` | 四象限镜像的 x/y 符号组合 | `[1,1;-1,1;1,-1;-1,-1]` | 赋值 | [sea_wave/RIS_float.m:449](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:449>) |
| `ids` | 四象限镜像点在积分表中的线性索引元胞 | `cell(1,4)`<br>`sub2ind([n,n],xx,yy)` | 赋值 | [sea_wave/RIS_float.m:450](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:450>) |
| `P` | 旧模态面元形函数 [1,x,y] | `[ones(size(xy,1),1),xy]` | 赋值 | [sea_wave/RIS_float.m:456](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:456>) |
| `parity` | 三个模态的 x/y 对称符号 [1,1;-1,1;1,-1] | `[1,1;-1,1;1,-1]` | 赋值 | [sea_wave/RIS_float.m:457](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:457>) |
| `gn` | 64 点 Gauss–Legendre 节点/权重，用于 Struve H0 积分 | `gaussLegendre(64)` | 赋值 | [sea_wave/RIS_float.m:461](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:461>) |
| `gw` | 64 点 Gauss–Legendre 节点/权重，用于 Struve H0 积分 | `gaussLegendre(64)` | 赋值 | [sea_wave/RIS_float.m:461](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:461>) |
| `angles` | H0 积分变量角度，映射到 [0,pi/2]（rad） | `(gn+1)*pi/4` | 赋值 | [sea_wave/RIS_float.m:462](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:462>) |
| `f` | 当前角频率的数组索引 | `1:numel(omega)` | 循环索引 | [sea_wave/RIS_float.m:463](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:463>) |
| `k` | 当前深水波数 omega(f)²/g（m⁻¹） | `omega(f)^2/p.g` | 赋值 | [sea_wave/RIS_float.m:464](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:464>) |
| `z` | 无量纲 Green 核函数自变量 k*radius | `k*radius` | 赋值 | [sea_wave/RIS_float.m:465](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:465>) |
| `H0` | 由 64 点高斯积分计算的 Struve H0 值 | `reshape(sin(z(:)*cos(angles).')*(gw/2),size(z))` | 赋值 | [sea_wave/RIS_float.m:466](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:466>) |
| `regular` | 自由面 Green 函数的正则部分复核值 | `-k/4*(bessely(0,z)+H0)-1i*k/2*besselj(0,z)` | 赋值 | [sea_wave/RIS_float.m:467](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:467>) |
| `table` | 奇异解析积分加正则高斯积分构成的完整面元 Green 表 | `singular+reshape(regular*quadrature.',n,n)` | 赋值 | [sea_wave/RIS_float.m:468](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:468>) |
| `mode` | 当前旧 BEM 模态索引 1..3 | `1:3` | 循环索引 | [sea_wave/RIS_float.m:469](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:469>) |
| `H` | 应用对应奇偶对称条件后形成的 BEM 核矩阵 | `table(ids{1})+sx*table(ids{2})+sy*table(ids{3})+sx*sy*table(ids{4})` | 赋值 | [sea_wave/RIS_float.m:471](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:471>) |
| `T` | 由形函数和势解计算的复水动力系数，其实部给 A，-omega*虚部给 B | `p.rho*4*dx*dy*(P(:,mode).'*(P(:,mode)-rr))/k` | 赋值 | [sea_wave/RIS_float.m:474](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:474>) |


**作用域：`gaussLegendre`**

| 源码符号 | 含义与单位/维度 | 源码取值/计算式（各次赋值） | 角色 | 源码位置 |
| --- | --- | --- | --- | --- |
| `n` | Gauss–Legendre 求积阶数 | 函数输入/输出；值随调用或上游结果确定 | 输入 | [sea_wave/RIS_float.m:483](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:483>) |
| `x` | 按升序排列的 [-1,1] Gauss–Legendre 节点 | `sort(diag(D))` | 输出/赋值 | [sea_wave/RIS_float.m:483](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:483>) |
| `w` | 对应高斯积分权重 | `2*V(1,order).'.^2` | 输出/赋值 | [sea_wave/RIS_float.m:483](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:483>) |
| `b` | Legendre Jacobi 对称三对角矩阵的副对角元素 | `(1:n-1)./sqrt(4*(1:n-1).^2-1)` | 赋值 | [sea_wave/RIS_float.m:485](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:485>) |
| `V` | Jacobi 矩阵特征向量/特征值矩阵 | `eig(diag(b,1)+diag(b,-1))` | 赋值 | [sea_wave/RIS_float.m:486](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:486>) |
| `D` | Jacobi 矩阵特征向量/特征值矩阵 | `eig(diag(b,1)+diag(b,-1))` | 赋值 | [sea_wave/RIS_float.m:486](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:486>) |
| `order` | 节点特征值升序排列索引 | `sort(diag(D))` | 赋值 | [sea_wave/RIS_float.m:487](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:487>) |


**作用域：`plotSea`**

| 源码符号 | 含义与单位/维度 | 源码取值/计算式（各次赋值） | 角色 | 源码位置 |
| --- | --- | --- | --- | --- |
| `o` | 当前 RIS 运动结果或动力学模型；具体类型由所在函数输入确定 | 函数输入/输出；值随调用或上游结果确定 | 输入 | [sea_wave/RIS_float.m:490](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:490>) |
| `high` | 同一轨迹升沉最大/最小值的样本索引 | `max(o.heave)` | 赋值 | [sea_wave/RIS_float.m:494](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:494>) |
| `low` | 同一轨迹升沉最大/最小值的样本索引 | `min(o.heave)` | 赋值 | [sea_wave/RIS_float.m:494](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:494>) |
| `rising` | 同一轨迹最大/最小升沉速度的样本索引 | `max(o.velocity(1,:))` | 赋值 | [sea_wave/RIS_float.m:495](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:495>) |
| `falling` | 同一轨迹最大/最小升沉速度的样本索引 | `min(o.velocity(1,:))` | 赋值 | [sea_wave/RIS_float.m:495](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:495>) |
| `indices` | 被选中的分量、样本或数组位置索引 | `[high,rising,low,falling]` | 赋值 | [sea_wave/RIS_float.m:496](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:496>) |
| `labels` | 工况、相位或图形标签文本 | `["Maximum heave","Maximum vertical velocity","Minimum heave","Minimum vertical velocity"]` | 赋值 | [sea_wave/RIS_float.m:497](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:497>) |
| `halfWidth` | 围绕板中心的海面显示半宽 2.5*max(Lx,Ly)（m） | `2.5*max(o.parameters.Lx,o.parameters.Ly)` | 赋值 | [sea_wave/RIS_float.m:498](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:498>) |
| `lambdaMin` | 最高海浪角频率对应的最短深水波长（m） | `2*pi*o.parameters.g/max(o.wave.omega)^2` | 赋值 | [sea_wave/RIS_float.m:499](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:499>) |
| `count` | 满足最短波长半波空间采样要求的显示网格点数 | `ceil(2*halfWidth/(lambdaMin/2))+1` | 赋值 | [sea_wave/RIS_float.m:500](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:500>) |
| `scenes` | 海面显示快照结构（X/Y/Z）或其集合 | `cell(1,4)`<br>`struct('X',X,'Y',Y,'Z',Z)` | 赋值 | [sea_wave/RIS_float.m:501](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:501>) |
| `j` | 当前循环、分量或样本索引（无量纲，MATLAB 从 1 开始） | `1:4` | 循环索引 | [sea_wave/RIS_float.m:502](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:502>) |
| `center` | 当前极值样本的世界水平质心坐标（m） | `o.centerTrajectory(:,indices(j))` | 赋值 | [sea_wave/RIS_float.m:503](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:503>) |
| `X` | 海面网格的世界 x 坐标矩阵（m） | `RIS_sea_surface(o,center(1)+linspace(-halfWidth,halfWidth,count), center(2)+linspace(-halfWidth,halfWidth,count),o.t(indices(j)))` | 赋值 | [sea_wave/RIS_float.m:504](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:504>) |
| `Y` | 海面网格的世界 y 坐标矩阵（m） | `RIS_sea_surface(o,center(1)+linspace(-halfWidth,halfWidth,count), center(2)+linspace(-halfWidth,halfWidth,count),o.t(indices(j)))` | 赋值 | [sea_wave/RIS_float.m:504](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:504>) |
| `Z` | 海面高度矩阵（m） | `RIS_sea_surface(o,center(1)+linspace(-halfWidth,halfWidth,count), center(2)+linspace(-halfWidth,halfWidth,count),o.t(indices(j)))` | 赋值 | [sea_wave/RIS_float.m:504](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:504>) |
| `scenes.X` | 海面世界 x/y 坐标与高程矩阵（m），按字段区分 | `X` | 结构字段 | [sea_wave/RIS_float.m:506](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:506>) |
| `scenes.Y` | 海面世界 x/y 坐标与高程矩阵（m），按字段区分 | `Y` | 结构字段 | [sea_wave/RIS_float.m:506](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:506>) |
| `scenes.Z` | 海面世界 x/y 坐标与高程矩阵（m），按字段区分 | `Z` | 结构字段 | [sea_wave/RIS_float.m:506](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:506>) |
| `layout` | 分面绘图布局句柄 | `tiledlayout(2,2,'TileSpacing','compact','Padding','compact')` | 赋值 | [sea_wave/RIS_float.m:509](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:509>) |
| `ax` | 绘图坐标轴句柄 | `nexttile(layout)` | 赋值 | [sea_wave/RIS_float.m:511](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:511>) |
| `index` | 被选中的分量、样本或数组位置索引 | `indices(j)` | 赋值 | [sea_wave/RIS_float.m:511](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:511>) |


### 5.8 sea_wave/RIS_radiation.m

**作用域：`RIS_radiation`**

| 源码符号 | 含义与单位/维度 | 源码取值/计算式（各次赋值） | 角色 | 源码位置 |
| --- | --- | --- | --- | --- |
| `omega` | 输入 BEM 角频率（rad/s） | 函数输入/输出；值随调用或上游结果确定 | 输入 | [sea_wave/RIS_radiation.m:1](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_radiation.m:1>) |
| `addedMass` | 三垂向模态的 BEM 附加广义质量频响（3×Nfreq） | 函数输入/输出；值随调用或上游结果确定 | 输入 | [sea_wave/RIS_radiation.m:1](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_radiation.m:1>) |
| `damping` | 三垂向模态的 BEM 广义阻尼频响（3×Nfreq） | 函数输入/输出；值随调用或上游结果确定 | 输入 | [sea_wave/RIS_radiation.m:1](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_radiation.m:1>) |
| `modeCount` | 被动拟合振子数，整数且至少 4 | 函数输入/输出；值随调用或上游结果确定 | 输入 | [sea_wave/RIS_radiation.m:1](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_radiation.m:1>) |
| `radiation` | 拟合完成的被动辐射记忆模型结构 | `struct(...)（所有字段在本节逐行列明）` | 输出/赋值 | [sea_wave/RIS_radiation.m:1](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_radiation.m:1>) |
| `frequency` | 从 min(omega)/2 到 2*max(omega) 的振子中心角频率（rad/s） | `logspace(log10(min(omega)/2),log10(max(omega)*2),modeCount).'` | 赋值 | [sea_wave/RIS_radiation.m:11](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_radiation.m:11>) |
| `width` | 按对数极点间隔确定的振子衰减率（s⁻¹），不是阻力经验系数 | `frequency*(exp(log(frequency(end)/frequency(1))/(modeCount-1))-1)` | 赋值 | [sea_wave/RIS_radiation.m:13](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_radiation.m:13>) |
| `s` | 拉普拉斯虚频率 1i*omega（s⁻¹） | `1i*omega` | 赋值 | [sea_wave/RIS_radiation.m:14](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_radiation.m:14>) |
| `basis` | 各振子的 (s+a)/((s+a)²+b²) 复频响基函数 | `(s+width.')./((s+width.').^2+frequency.'.^2)` | 赋值 | [sea_wave/RIS_radiation.m:15](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_radiation.m:15>) |
| `design` | 追加无穷频附加质量列 s 后的复设计矩阵 | `[basis,s]` | 赋值 | [sea_wave/RIS_radiation.m:16](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_radiation.m:16>) |
| `residues` | 各振子各模态的非负拟合留数 | `zeros(modeCount,3)`<br>`solution(1:end-1)` | 赋值 | [sea_wave/RIS_radiation.m:17](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_radiation.m:17>) |
| `infinity` | 三个模态非负拟合的无穷频附加广义质量 | `zeros(3,1)`<br>`solution(end)` | 赋值 | [sea_wave/RIS_radiation.m:17](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_radiation.m:17>) |
| `fittedMass` | 重建的附加广义质量/阻尼频响（按变量名分别对应） | `zeros(size(addedMass))`<br>`imag(fit)./omega` | 赋值 | [sea_wave/RIS_radiation.m:18](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_radiation.m:18>) |
| `fittedDamping` | 重建的附加广义质量/阻尼频响（按变量名分别对应） | `fittedMass`<br>`real(fit)` | 赋值 | [sea_wave/RIS_radiation.m:18](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_radiation.m:18>) |
| `errorRelative` | 三个模态的加权阻抗相对 RMS 拟合误差 | `zeros(3,1)`<br>`norm((fit-target)./scale)/sqrt(numel(omega))` | 赋值 | [sea_wave/RIS_radiation.m:19](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_radiation.m:19>) |
| `j` | 当前循环、分量或样本索引（无量纲，MATLAB 从 1 开始） | `1:3` | 循环索引 | [sea_wave/RIS_radiation.m:20](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_radiation.m:20>) |
| `target` | 某模态目标复阻抗 damping+1i*omega*addedMass | `damping(j,:).'+s.*addedMass(j,:).'` | 赋值 | [sea_wave/RIS_radiation.m:21](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_radiation.m:21>) |
| `scale` | 频响加权尺度，下限为该模态最大目标幅度的 0.02 倍 | `max(abs(target),0.02*max(abs(target)))` | 赋值 | [sea_wave/RIS_radiation.m:22](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_radiation.m:22>) |
| `weighted` | 逐频率加权后的复设计矩阵 | `design./scale` | 赋值 | [sea_wave/RIS_radiation.m:23](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_radiation.m:23>) |
| `realDesign` | 将实部与虚部纵向堆叠后的实非负最小二乘矩阵 | `[real(weighted);imag(weighted)]` | 赋值 | [sea_wave/RIS_radiation.m:24](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_radiation.m:24>) |
| `columnScale` | 实设计矩阵每列的欧氏范数，用于列归一化 | `vecnorm(realDesign)` | 赋值 | [sea_wave/RIS_radiation.m:25](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_radiation.m:25>) |
| `solution` | 非负最小二乘解，包含振子留数与无穷频附加质量 | `lsqnonneg(realDesign./columnScale,[real(target./scale);imag(target./scale)])`<br>`solution./columnScale.'` | 赋值 | [sea_wave/RIS_radiation.m:26](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_radiation.m:26>) |
| `fit` | 拟合得到的目标复阻抗 | `design*solution` | 赋值 | [sea_wave/RIS_radiation.m:30](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_radiation.m:30>) |
| `radiation.frequency` | 辐射振子的中心角频率（rad/s） | `frequency` | 结构字段 | [sea_wave/RIS_radiation.m:35](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_radiation.m:35>) |
| `radiation.width` | 辐射振子的衰减率/极点实部绝对值（s⁻¹） | `width` | 结构字段 | [sea_wave/RIS_radiation.m:35](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_radiation.m:35>) |
| `radiation.residues` | 各辐射振子对三个垂向模态的非负留数矩阵 | `residues` | 结构字段 | [sea_wave/RIS_radiation.m:35](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_radiation.m:35>) |
| `radiation.infiniteAddedMass` | 无穷频率附加广义质量；升沉项 kg、转动项 kg·m² | `infinity` | 结构字段 | [sea_wave/RIS_radiation.m:36](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_radiation.m:36>) |
| `radiation.fittedAddedMass` | 辐射有理模型拟合的附加广义质量频响 | `fittedMass` | 结构字段 | [sea_wave/RIS_radiation.m:36](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_radiation.m:36>) |
| `radiation.fittedDamping` | 辐射有理模型拟合的广义阻尼频响 | `fittedDamping` | 结构字段 | [sea_wave/RIS_radiation.m:37](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_radiation.m:37>) |
| `radiation.relativeImpedanceRmsError` | 按频响幅值加权的阻抗拟合相对均方根误差（无量纲，三个垂向模态） | `errorRelative` | 结构字段 | [sea_wave/RIS_radiation.m:37](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_radiation.m:37>) |
| `radiation.passive` | 非负留数和非负无穷频附加质量的被动性标志 | `all(residues(:)>=0) && all(infinity>=0)` | 结构字段 | [sea_wave/RIS_radiation.m:38](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_radiation.m:38>) |


### 5.9 sea_wave/RIS_rotation.m

**作用域：`RIS_rotation`**

| 源码符号 | 含义与单位/维度 | 源码取值/计算式（各次赋值） | 角色 | 源码位置 |
| --- | --- | --- | --- | --- |
| `value` | 输入姿态：3 元 [roll;pitch;yaw]（rad）或非零 4 元 wxyz 四元数 | 函数输入/输出；值随调用或上游结果确定 | 输入 | [sea_wave/RIS_rotation.m:1](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_rotation.m:1>) |
| `R` | 体坐标到世界坐标的主动旋转矩阵（无量纲） | `[1-2*(y*y+z*z),2*(x*y-w*z),2*(x*z+w*y); 2*(x*y+w*z),1-2*(x*x+z*z),2*(y*z-w*x); 2*(x*z-w*y),2*(y*z+w*x),1-2*(x*x+y*y)]` | 输出/赋值 | [sea_wave/RIS_rotation.m:1](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_rotation.m:1>) |
| `quaternion` | 姿态单位四元数，顺序 [w;x;y;z] | `[c(1)*c(2)*c(3)+s(1)*s(2)*s(3); s(1)*c(2)*c(3)-c(1)*s(2)*s(3); c(1)*s(2)*c(3)+s(1)*c(2)*s(3); c(1)*c(2)*s(3)-s(1)*s(2)*c(3)]`<br>`value/norm(value)` | 输出/赋值 | [sea_wave/RIS_rotation.m:1](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_rotation.m:1>) |
| `euler` | 右手 ZYX 欧拉角 [roll;pitch;yaw]（rad） | `[atan2(R(3,2),R(3,3));asin(max(-1,min(1,-R(3,1))));atan2(R(2,1),R(1,1))]` | 输出/赋值 | [sea_wave/RIS_rotation.m:1](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_rotation.m:1>) |
| `c` | 各欧拉角半角的余弦/正弦向量 | `cos(value/2)` | 赋值 | [sea_wave/RIS_rotation.m:7](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_rotation.m:7>) |
| `s` | 各欧拉角半角的余弦/正弦向量 | `sin(value/2)` | 赋值 | [sea_wave/RIS_rotation.m:7](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_rotation.m:7>) |
| `w` | 单位四元数的实部 w 和虚部 x/y/z 分量；此处不是空间坐标 | `quaternion(1)` | 赋值 | [sea_wave/RIS_rotation.m:17](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_rotation.m:17>) |
| `x` | 单位四元数的实部 w 和虚部 x/y/z 分量；此处不是空间坐标 | `quaternion(2)` | 赋值 | [sea_wave/RIS_rotation.m:17](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_rotation.m:17>) |
| `y` | 单位四元数的实部 w 和虚部 x/y/z 分量；此处不是空间坐标 | `quaternion(3)` | 赋值 | [sea_wave/RIS_rotation.m:17](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_rotation.m:17>) |
| `z` | 单位四元数的实部 w 和虚部 x/y/z 分量；此处不是空间坐标 | `quaternion(4)` | 赋值 | [sea_wave/RIS_rotation.m:17](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_rotation.m:17>) |


### 5.10 sea_wave/RIS_sea_surface.m

**作用域：`RIS_sea_surface`**

| 源码符号 | 含义与单位/维度 | 源码取值/计算式（各次赋值） | 角色 | 源码位置 |
| --- | --- | --- | --- | --- |
| `motion` | 当前 RIS 运动结果或动力学模型；具体类型由所在函数输入确定 | 函数输入/输出；值随调用或上游结果确定 | 输入 | [sea_wave/RIS_sea_surface.m:1](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_sea_surface.m:1>) |
| `x` | 当前坐标分量（m）；重名参数按所在函数解释 | 函数输入/输出；值随调用或上游结果确定 | 输入 | [sea_wave/RIS_sea_surface.m:1](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_sea_surface.m:1>) |
| `y` | 当前坐标分量（m）；重名参数按所在函数解释 | 函数输入/输出；值随调用或上游结果确定 | 输入 | [sea_wave/RIS_sea_surface.m:1](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_sea_surface.m:1>) |
| `time` | 相对 parameters.startTime 的时间（s）或其采样向量 | 函数输入/输出；值随调用或上游结果确定 | 输入 | [sea_wave/RIS_sea_surface.m:1](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_sea_surface.m:1>) |
| `X` | 海面网格的世界 x 坐标矩阵（m） | `meshgrid(x,y)` | 输出/赋值 | [sea_wave/RIS_sea_surface.m:1](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_sea_surface.m:1>) |
| `Y` | 海面网格的世界 y 坐标矩阵（m） | `meshgrid(x,y)` | 输出/赋值 | [sea_wave/RIS_sea_surface.m:1](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_sea_surface.m:1>) |
| `Z` | 海面高度矩阵（m） | `reshape(RIS_wave_field(motion,[X(:),Y(:)],time),size(X))` | 输出/赋值 | [sea_wave/RIS_sea_surface.m:1](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_sea_surface.m:1>) |


### 5.11 sea_wave/RIS_sidewall_load.m

**作用域：`RIS_sidewall_load`**

| 源码符号 | 含义与单位/维度 | 源码取值/计算式（各次赋值） | 角色 | 源码位置 |
| --- | --- | --- | --- | --- |
| `p` | 当前作用域的参数结构 | 函数输入/输出；值随调用或上游结果确定 | 输入 | [sea_wave/RIS_sidewall_load.m:1](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_sidewall_load.m:1>) |
| `wave` | 海浪分量结构，含 omega、kx、ky、coefficient；具体重名按下方作用域解释 | 函数输入/输出；值随调用或上游结果确定 | 输入 | [sea_wave/RIS_sidewall_load.m:1](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_sidewall_load.m:1>) |
| `draft` | 静水平衡吃水 mass/(rho*Lx*Ly)（m） | 函数输入/输出；值随调用或上游结果确定 | 输入 | [sea_wave/RIS_sidewall_load.m:1](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_sidewall_load.m:1>) |
| `zG` | 静水平衡质心高度 thickness/2-draft（m） | 函数输入/输出；值随调用或上游结果确定 | 输入 | [sea_wave/RIS_sidewall_load.m:1](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_sidewall_load.m:1>) |
| `load` | 侧壁每米波幅的 [Fx;Fy;Fz;Mx;My;Mz] 复激励（N/m / N·m/m） | `[fx;fy;zeros(size(fx));-lever.*fy;lever.*fx;mz]` | 输出/赋值 | [sea_wave/RIS_sidewall_load.m:1](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_sidewall_load.m:1>) |
| `kx` | 初始体坐标入射波数分量行向量（m⁻¹） | `wave.kx(:).'` | 赋值 | [sea_wave/RIS_sidewall_load.m:11](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_sidewall_load.m:11>) |
| `ky` | 初始体坐标入射波数分量行向量（m⁻¹） | `wave.ky(:).'` | 赋值 | [sea_wave/RIS_sidewall_load.m:11](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_sidewall_load.m:11>) |
| `k` | 深水波数 omega²/g（m⁻¹） | `wave.omega(:).'.^2/p.g` | 赋值 | [sea_wave/RIS_sidewall_load.m:11](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_sidewall_load.m:11>) |
| `sx` | 归一化矩形平面波积分 sinc(kx*Lx/2)、sinc(ky*Ly/2)，sinc=sin(x)/x | `sincDerivative(kx*p.Lx/2)` | 赋值 | [sea_wave/RIS_sidewall_load.m:12](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_sidewall_load.m:12>) |
| `dsx` | 上述无 pi 定义 sinc 对其自变量的一阶导数 | `sincDerivative(kx*p.Lx/2)` | 赋值 | [sea_wave/RIS_sidewall_load.m:12](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_sidewall_load.m:12>) |
| `sy` | 归一化矩形平面波积分 sinc(kx*Lx/2)、sinc(ky*Ly/2)，sinc=sin(x)/x | `sincDerivative(ky*p.Ly/2)` | 赋值 | [sea_wave/RIS_sidewall_load.m:13](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_sidewall_load.m:13>) |
| `dsy` | 上述无 pi 定义 sinc 对其自变量的一阶导数 | `sincDerivative(ky*p.Ly/2)` | 赋值 | [sea_wave/RIS_sidewall_load.m:13](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_sidewall_load.m:13>) |
| `area` | 当前板面或网格单元面积（m²）；见所在函数计算式 | `p.Lx*p.Ly` | 赋值 | [sea_wave/RIS_sidewall_load.m:14](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_sidewall_load.m:14>) |
| `footprint` | 矩形水平面积的复平面波积分 area*sx*sy（m²） | `area*sx.*sy` | 赋值 | [sea_wave/RIS_sidewall_load.m:15](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_sidewall_load.m:15>) |
| `dx` | footprint 对 kx/ky 的偏导数（m³）；此处不是网格间距 | `area*p.Lx/2*dsx.*sy` | 赋值 | [sea_wave/RIS_sidewall_load.m:16](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_sidewall_load.m:16>) |
| `dy` | footprint 对 kx/ky 的偏导数（m³）；此处不是网格间距 | `area*p.Ly/2*sx.*dsy` | 赋值 | [sea_wave/RIS_sidewall_load.m:17](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_sidewall_load.m:17>) |
| `depth` | 从 -draft 到 0 的 exp(k*z) 深度积分（m） | `-expm1(-k*draft)./k` | 赋值 | [sea_wave/RIS_sidewall_load.m:18](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_sidewall_load.m:18>) |
| `kd` | 波数乘吃水 k*draft（无量纲） | `k*draft` | 赋值 | [sea_wave/RIS_sidewall_load.m:20](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_sidewall_load.m:20>) |
| `lever` | 相对质心的压力加权垂向力臂（m） | `-draft/2+zeros(size(k))`<br>`1./k(large)-draft./expm1(kd(large))`<br>`-lever(large)`<br>`-draft/2+k(~large)*draft^2/12-k(~large).^3*draft^4/720`<br>`lever-zG` | 赋值 | [sea_wave/RIS_sidewall_load.m:21](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_sidewall_load.m:21>) |
| `large` | kd>1e-3 的布尔掩码，用于切换解析式和小量级数 | `kd>1e-3` | 赋值 | [sea_wave/RIS_sidewall_load.m:22](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_sidewall_load.m:22>) |
| `factor` | rho*g*depth，侧壁波浪压力深度积分系数 | `p.rho*p.g*depth` | 赋值 | [sea_wave/RIS_sidewall_load.m:27](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_sidewall_load.m:27>) |
| `fx` | 侧壁单位入射波幅的水平复力（N/m） | `1i*factor.*kx.*footprint` | 赋值 | [sea_wave/RIS_sidewall_load.m:28](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_sidewall_load.m:28>) |
| `fy` | 侧壁单位入射波幅的水平复力（N/m） | `1i*factor.*ky.*footprint` | 赋值 | [sea_wave/RIS_sidewall_load.m:29](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_sidewall_load.m:29>) |
| `mz` | 侧壁单位入射波幅的偏航复力矩（N·m/m） | `factor.*(kx.*dy-ky.*dx)` | 赋值 | [sea_wave/RIS_sidewall_load.m:30](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_sidewall_load.m:30>) |


**作用域：`sincDerivative`**

| 源码符号 | 含义与单位/维度 | 源码取值/计算式（各次赋值） | 角色 | 源码位置 |
| --- | --- | --- | --- | --- |
| `x` | 无 pi 归一化 sinc 的无量纲自变量 | 函数输入/输出；值随调用或上游结果确定 | 输入 | [sea_wave/RIS_sidewall_load.m:33](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_sidewall_load.m:33>) |
| `s` | sin(x)/x 的数值稳定计算结果 | `ones(size(x))`<br>`1-x(small).^2/6+x(small).^4/120`<br>`sin(x(~small))./x(~small)` | 输出/赋值 | [sea_wave/RIS_sidewall_load.m:33](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_sidewall_load.m:33>) |
| `ds` | sinc 对 x 的一阶导数 | `zeros(size(x))`<br>`-x(small)/3+x(small).^3/30-x(small).^5/840`<br>`(x(~small).*cos(x(~small))-sin(x(~small)))./x(~small).^2` | 输出/赋值 | [sea_wave/RIS_sidewall_load.m:33](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_sidewall_load.m:33>) |
| `small` | abs(x)<1e-3 的小量级数展开掩码 | `abs(x)<1e-3` | 赋值 | [sea_wave/RIS_sidewall_load.m:36](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_sidewall_load.m:36>) |


### 5.12 sea_wave/RIS_wave_components.m

**作用域：`RIS_wave_components`**

| 源码符号 | 含义与单位/维度 | 源码取值/计算式（各次赋值） | 角色 | 源码位置 |
| --- | --- | --- | --- | --- |
| `p` | 当前作用域的参数结构 | 函数输入/输出；值随调用或上游结果确定 | 输入 | [sea_wave/RIS_wave_components.m:1](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_wave_components.m:1>) |
| `quadrature` | 可选积分规则参数，默认及唯一合法值 'trapezoid' | `'trapezoid' end` | 输入/声明默认值 | [sea_wave/RIS_wave_components.m:1](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_wave_components.m:1>) |
| `wave` | 海浪分量结构，含 omega、kx、ky、coefficient；具体重名按下方作用域解释 | `struct('omega',om(:),'kx',reshape(k.*cos(beta),[],1), 'ky',reshape(k.*sin(beta),[],1),'coefficient',amplitude(:).*exp(1i*phase(:)))` | 输出/赋值 | [sea_wave/RIS_wave_components.m:1](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_wave_components.m:1>) |
| `spectrum` | 海浪频谱及离散积分诊断结构 | `struct(...)（所有字段在本节逐行列明）` | 输出/赋值 | [sea_wave/RIS_wave_components.m:1](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_wave_components.m:1>) |
| `f` | 对数离散波频率列向量（Hz） | `logspace(log10(p.fMinHz),log10(p.fMaxHz),p.frequencyBins).'` | 赋值 | [sea_wave/RIS_wave_components.m:9](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_wave_components.m:9>) |
| `th` | 相对主传播方向的离散方向角（rad，[-pi/2,pi/2]） | `linspace(-pi/2,pi/2,p.directionBins)` | 赋值 | [sea_wave/RIS_wave_components.m:10](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_wave_components.m:10>) |
| `df` | 频率非均匀梯形积分权重（Hz） | `[(f(2)-f(1))/2;(f(3:end)-f(1:end-2))/2;(f(end)-f(end-1))/2]` | 赋值 | [sea_wave/RIS_wave_components.m:11](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_wave_components.m:11>) |
| `dw` | 方向梯形积分权重（rad） | `repmat(th(2)-th(1),size(th))`<br>`dw([1,end])/2` | 赋值 | [sea_wave/RIS_wave_components.m:12](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_wave_components.m:12>) |
| `om` | 二维角频率网格（rad/s） | `ndgrid(2*pi*f,th)` | 赋值 | [sea_wave/RIS_wave_components.m:14](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_wave_components.m:14>) |
| `theta` | 二维相对传播方向角网格（rad） | `ndgrid(2*pi*f,th)` | 赋值 | [sea_wave/RIS_wave_components.m:14](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_wave_components.m:14>) |
| `sea` | 由 get_default_sea_params 读取的基础海面配置 | `get_default_sea_params()` | 赋值 | [sea_wave/RIS_wave_components.m:15](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_wave_components.m:15>) |
| `S` | Neumann 单边角频率海面高程谱，S=3.05*pi/4*om^(-6)*exp(-2*g²/(U²*om²)) | `sea.neumannConstant*pi/4*om.^(-6).*exp(-2*p.g^2./(p.windSpeed^2*om.^2))` | 赋值 | [sea_wave/RIS_wave_components.m:16](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_wave_components.m:16>) |
| `D` | 方向展宽密度 2/pi*cos(theta)²（按弧度积分） | `2/pi*cos(theta).^2` | 赋值 | [sea_wave/RIS_wave_components.m:17](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_wave_components.m:17>) |
| `amplitude` | 每个方向/频率分量的实波幅 sqrt(2*S*D*2*pi*df*dw)（m） | `sqrt(2*S.*D.*(2*pi*df*dw))` | 赋值 | [sea_wave/RIS_wave_components.m:18](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_wave_components.m:18>) |
| `stream` | 局部 Mersenne Twister 随机数流 | `RandStream('mt19937ar','Seed',p.seed)` | 赋值 | [sea_wave/RIS_wave_components.m:19](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_wave_components.m:19>) |
| `phase` | 独立均匀随机初相位 [0,2*pi)（rad） | `2*pi*rand(stream,size(amplitude))` | 赋值 | [sea_wave/RIS_wave_components.m:20](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_wave_components.m:20>) |
| `beta` | 相对世界 +X 的传播角 theta+pi/2-deg2rad(windAngleDeg)（rad） | `theta+pi/2-deg2rad(p.windAngleDeg)` | 赋值 | [sea_wave/RIS_wave_components.m:21](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_wave_components.m:21>) |
| `k` | 深水波数 om²/g（m⁻¹） | `om.^2/p.g` | 赋值 | [sea_wave/RIS_wave_components.m:22](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_wave_components.m:22>) |
| `wave.omega` | 海浪角频率（rad/s）；EGG parameters.omega 则是混合概率 | `om(:)` | 结构字段 | [sea_wave/RIS_wave_components.m:23](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_wave_components.m:23>) |
| `wave.kx` | 海浪波数的世界 x/y 分量（m⁻¹）；bodyWave 中为初始体坐标分量 | `reshape(k.*cos(beta),[],1)` | 结构字段 | [sea_wave/RIS_wave_components.m:23](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_wave_components.m:23>) |
| `wave.ky` | 海浪波数的世界 x/y 分量（m⁻¹）；bodyWave 中为初始体坐标分量 | `reshape(k.*sin(beta),[],1)` | 结构字段 | [sea_wave/RIS_wave_components.m:24](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_wave_components.m:24>) |
| `wave.coefficient` | 海浪复振幅 amplitude*exp(1i*phase)（m），参考世界原点与绝对 t=0 | `amplitude(:).*exp(1i*phase(:))` | 结构字段 | [sea_wave/RIS_wave_components.m:24](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_wave_components.m:24>) |
| `a` | 计算谱保留比例的参数 2*g²/windSpeed²（s⁻²） | `2*p.g^2/p.windSpeed^2` | 赋值 | [sea_wave/RIS_wave_components.m:25](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_wave_components.m:25>) |
| `retained` | 当前频带占完整连续谱方差的解析比例 | `gammainc(a/(2*pi*p.fMinHz)^2,2.5)-gammainc(a/(2*pi*p.fMaxHz)^2,2.5)` | 赋值 | [sea_wave/RIS_wave_components.m:26](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_wave_components.m:26>) |
| `spectrum.name` | 模型名称或输出名称文本 | `'project Neumann S(omega), cos^2 spreading'` | 结构字段 | [sea_wave/RIS_wave_components.m:27](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_wave_components.m:27>) |
| `spectrum.quadrature` | 频率/方向积分规则，当前仅支持 trapezoid | `quadrature` | 结构字段 | [sea_wave/RIS_wave_components.m:28](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_wave_components.m:28>) |
| `spectrum.retainedFullSpectrumVariance` | 所取频带保留的完整连续 Neumann 谱方差比例（无量纲） | `retained` | 结构字段 | [sea_wave/RIS_wave_components.m:28](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_wave_components.m:28>) |
| `spectrum.frequencyWeightsHz` | 非均匀频率梯形积分权重（Hz），两端半权重 | `df` | 结构字段 | [sea_wave/RIS_wave_components.m:29](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_wave_components.m:29>) |
| `spectrum.directionWeightsRad` | 方向梯形积分权重（rad），两端半权重 | `dw` | 结构字段 | [sea_wave/RIS_wave_components.m:29](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_wave_components.m:29>) |
| `spectrum.frequencyBandHz` | 海浪频率范围 [fMinHz,fMaxHz]（Hz） | `[p.fMinHz,p.fMaxHz]` | 结构字段 | [sea_wave/RIS_wave_components.m:29](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_wave_components.m:29>) |


### 5.13 sea_wave/RIS_wave_field.m

**作用域：`RIS_wave_field`**

| 源码符号 | 含义与单位/维度 | 源码取值/计算式（各次赋值） | 角色 | 源码位置 |
| --- | --- | --- | --- | --- |
| `motion` | 当前 RIS 运动结果或动力学模型；具体类型由所在函数输入确定 | 函数输入/输出；值随调用或上游结果确定 | 输入 | [sea_wave/RIS_wave_field.m:1](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_wave_field.m:1>) |
| `points` | 世界坐标采样点数组；波面接口为 N×2 的 [x,y]（m） | 函数输入/输出；值随调用或上游结果确定 | 输入 | [sea_wave/RIS_wave_field.m:1](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_wave_field.m:1>) |
| `time` | 相对 parameters.startTime 的时间（s）或其采样向量 | 函数输入/输出；值随调用或上游结果确定 | 输入 | [sea_wave/RIS_wave_field.m:1](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_wave_field.m:1>) |
| `elevation` | 海面高程（m） | `zeros(size(points,1),1)`<br>`elevation+real(sum(field,1)).'` | 输出/赋值 | [sea_wave/RIS_wave_field.m:1](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_wave_field.m:1>) |
| `gradient` | 海面 x/y 方向解析斜率（无量纲） | `zeros(size(points))`<br>`gradient+real([-1i*wave.kx(ids).';-1i*wave.ky(ids).']*field).'` | 输出/赋值 | [sea_wave/RIS_wave_field.m:1](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_wave_field.m:1>) |
| `wave` | 海浪分量结构，含 omega、kx、ky、coefficient；具体重名按下方作用域解释 | `motion.wave` | 赋值 | [sea_wave/RIS_wave_field.m:13](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_wave_field.m:13>) |
| `absoluteTime` | 绝对海浪时刻（s） | `motion.parameters.startTime+time(:).'` | 赋值 | [sea_wave/RIS_wave_field.m:14](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_wave_field.m:14>) |
| `block` | 一次计算的海浪分量数，上限 64，并控制中间矩阵约不超过 2e6 元素 | `max(1,min(64,floor(2e6/size(points,1))))` | 赋值 | [sea_wave/RIS_wave_field.m:17](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_wave_field.m:17>) |
| `first` | 当前分块/积分子步的起始索引 | `1:block:numel(wave.omega)` | 循环索引 | [sea_wave/RIS_wave_field.m:18](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_wave_field.m:18>) |
| `ids` | 被选中的分量、样本或数组位置索引 | `first:min(first+block-1,numel(wave.omega))` | 赋值 | [sea_wave/RIS_wave_field.m:19](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_wave_field.m:19>) |
| `phase` | omega*绝对时间-kx*x-ky*y 的时空相位（rad） | `wave.omega(ids)*absoluteTime-wave.kx(ids)*points(:,1).'- wave.ky(ids)*points(:,2).'` | 赋值 | [sea_wave/RIS_wave_field.m:20](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_wave_field.m:20>) |
| `field` | 当前海浪分块在所有点/时刻的复高程贡献（m） | `wave.coefficient(ids).*exp(1i*phase)` | 赋值 | [sea_wave/RIS_wave_field.m:22](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_wave_field.m:22>) |


### 5.14 sea_wave/run_RIS_wind_comparison.m

**作用域：`run_RIS_wind_comparison`**

| 源码符号 | 含义与单位/维度 | 源码取值/计算式（各次赋值） | 角色 | 源码位置 |
| --- | --- | --- | --- | --- |
| `outputFolder` | 计算结果与图表输出目录 | `fullfile(folder,'FigOutput_RIS_WindDirection')` | 输入/赋值 | [sea_wave/run_RIS_wind_comparison.m:1](<C:/Users/Spidey/Desktop/code/sea_wave/run_RIS_wind_comparison.m:1>) |
| `suite` | 12 个动力学工况及其显示、导出结果的总结构 | `struct('motions',{motions},'scenes',{scenes},'summary',summary, 'meta',meta,'panelFiles',panelFiles,'overviewFile',overviewFile)` | 输出/赋值 | [sea_wave/run_RIS_wind_comparison.m:1](<C:/Users/Spidey/Desktop/code/sea_wave/run_RIS_wind_comparison.m:1>) |
| `folder` | 当前源码所在目录 | `fileparts(mfilename('fullpath'))` | 赋值 | [sea_wave/run_RIS_wind_comparison.m:5](<C:/Users/Spidey/Desktop/code/sea_wave/run_RIS_wind_comparison.m:5>) |
| `sea` | 由 get_default_sea_params 读取的基础海面配置 | `get_default_sea_params()` | 赋值 | [sea_wave/run_RIS_wind_comparison.m:10](<C:/Users/Spidey/Desktop/code/sea_wave/run_RIS_wind_comparison.m:10>) |
| `winds` | 四行风速 [3,5,9,12]（m/s） | `[3,5,9,12]` | 赋值 | [sea_wave/run_RIS_wind_comparison.m:11](<C:/Users/Spidey/Desktop/code/sea_wave/run_RIS_wind_comparison.m:11>) |
| `directions` | 三列传播风向角 [30,60,90]（°） | `[30,60,90]` | 赋值 | [sea_wave/run_RIS_wind_comparison.m:11](<C:/Users/Spidey/Desktop/code/sea_wave/run_RIS_wind_comparison.m:11>) |
| `snapshotTime` | 统一海面快照的相对时刻，默认 2 s（绝对 27 s） | `2` | 赋值 | [sea_wave/run_RIS_wind_comparison.m:12](<C:/Users/Spidey/Desktop/code/sea_wave/run_RIS_wind_comparison.m:12>) |
| `startTime` | 动力学绝对海浪起时，默认 25 s | `25` | 赋值 | [sea_wave/run_RIS_wind_comparison.m:12](<C:/Users/Spidey/Desktop/code/sea_wave/run_RIS_wind_comparison.m:12>) |
| `duration` | 动力学记录持续时间，默认 10 s | `10` | 赋值 | [sea_wave/run_RIS_wind_comparison.m:12](<C:/Users/Spidey/Desktop/code/sea_wave/run_RIS_wind_comparison.m:12>) |
| `seed` | 随机数种子（非负整数） | `600` | 赋值 | [sea_wave/run_RIS_wind_comparison.m:12](<C:/Users/Spidey/Desktop/code/sea_wave/run_RIS_wind_comparison.m:12>) |
| `zLimits` | 海面、轨迹、色轴显示裁切范围 [-0.5,0.5]（m） | `[-0.5,0.5]` | 赋值 | [sea_wave/run_RIS_wind_comparison.m:13](<C:/Users/Spidey/Desktop/code/sea_wave/run_RIS_wind_comparison.m:13>) |
| `x` | 当前坐标分量（m）；重名参数按所在函数解释 | `linspace(sea.x_range(1),sea.x_range(2),sea.grid_density)` | 赋值 | [sea_wave/run_RIS_wind_comparison.m:14](<C:/Users/Spidey/Desktop/code/sea_wave/run_RIS_wind_comparison.m:14>) |
| `y` | 当前坐标分量（m）；重名参数按所在函数解释 | `linspace(sea.y_range(1),sea.y_range(2),sea.grid_density)` | 赋值 | [sea_wave/run_RIS_wind_comparison.m:15](<C:/Users/Spidey/Desktop/code/sea_wave/run_RIS_wind_comparison.m:15>) |
| `motions` | 4×3 各风速/风向 RIS 完整运动结果元胞数组 | `cell(4,3)`<br>`o` | 赋值 | [sea_wave/run_RIS_wind_comparison.m:16](<C:/Users/Spidey/Desktop/code/sea_wave/run_RIS_wind_comparison.m:16>) |
| `scenes` | 4×3 对应快照 X/Y/Z 海面结构元胞数组 | `cell(4,3)`<br>`struct('X',X,'Y',Y,'Z',Z)` | 赋值 | [sea_wave/run_RIS_wind_comparison.m:16](<C:/Users/Spidey/Desktop/code/sea_wave/run_RIS_wind_comparison.m:16>) |
| `values` | 12 行、16 列工况统计数值矩阵，列名见 summary 构造式 | `zeros(12,16)`<br>`[winds(j),directions(k),o.spectrum.significantWaveHeightM, center.',o.rollDeg(sample),o.pitchDeg(sample),o.yawDeg(sample), min(Z,[],'all'),max(Z,[],'all'),mean(abs(Z(:))>zLimits(2)), min(o.centerOfMass(3,:)),max(o.centerOfMass(3,:)), o.validity.maxTiltDeg,max(o.validity.radiationFitRelativeRmsError)]` | 赋值 | [sea_wave/run_RIS_wind_comparison.m:17](<C:/Users/Spidey/Desktop/code/sea_wave/run_RIS_wind_comparison.m:17>) |
| `row` | 工况统计矩阵当前行号 1..12 | `0`<br>`row+1` | 赋值 | [sea_wave/run_RIS_wind_comparison.m:17](<C:/Users/Spidey/Desktop/code/sea_wave/run_RIS_wind_comparison.m:17>) |
| `panelFiles` | 4×4 单独图形 FIG 文件路径数组 | `strings(4,4)`<br>`saveFigure(f,outputFolder,name)` | 赋值 | [sea_wave/run_RIS_wind_comparison.m:18](<C:/Users/Spidey/Desktop/code/sea_wave/run_RIS_wind_comparison.m:18>) |
| `j` | 当前循环、分量或样本索引（无量纲，MATLAB 从 1 开始） | `1:4` | 循环索引 | [sea_wave/run_RIS_wind_comparison.m:19](<C:/Users/Spidey/Desktop/code/sea_wave/run_RIS_wind_comparison.m:19>) |
| `k` | 当前循环、分量或样本索引（无量纲，MATLAB 从 1 开始） | `1:3`<br>`1:4` | 循环索引 | [sea_wave/run_RIS_wind_comparison.m:20](<C:/Users/Spidey/Desktop/code/sea_wave/run_RIS_wind_comparison.m:20>) |
| `o` | 当前 RIS 运动结果或动力学模型；具体类型由所在函数输入确定 | `RIS_float(struct('windSpeed',winds(j),'windAngleDeg',directions(k), 'startTime',startTime,'duration',duration,'seed',seed,'plot',false))` | 赋值 | [sea_wave/run_RIS_wind_comparison.m:21](<C:/Users/Spidey/Desktop/code/sea_wave/run_RIS_wind_comparison.m:21>) |
| `X` | 海面网格的世界 x 坐标矩阵（m） | `RIS_sea_surface(o,x,y,snapshotTime)` | 赋值 | [sea_wave/run_RIS_wind_comparison.m:23](<C:/Users/Spidey/Desktop/code/sea_wave/run_RIS_wind_comparison.m:23>) |
| `Y` | 海面网格的世界 y 坐标矩阵（m） | `RIS_sea_surface(o,x,y,snapshotTime)` | 赋值 | [sea_wave/run_RIS_wind_comparison.m:23](<C:/Users/Spidey/Desktop/code/sea_wave/run_RIS_wind_comparison.m:23>) |
| `Z` | 海面高度矩阵（m） | `RIS_sea_surface(o,x,y,snapshotTime)` | 赋值 | [sea_wave/run_RIS_wind_comparison.m:23](<C:/Users/Spidey/Desktop/code/sea_wave/run_RIS_wind_comparison.m:23>) |
| `scenes.X` | 海面世界 x/y 坐标与高程矩阵（m），按字段区分 | `X` | 结构字段 | [sea_wave/run_RIS_wind_comparison.m:25](<C:/Users/Spidey/Desktop/code/sea_wave/run_RIS_wind_comparison.m:25>) |
| `scenes.Y` | 海面世界 x/y 坐标与高程矩阵（m），按字段区分 | `Y` | 结构字段 | [sea_wave/run_RIS_wind_comparison.m:25](<C:/Users/Spidey/Desktop/code/sea_wave/run_RIS_wind_comparison.m:25>) |
| `scenes.Z` | 海面世界 x/y 坐标与高程矩阵（m），按字段区分 | `Z` | 结构字段 | [sea_wave/run_RIS_wind_comparison.m:25](<C:/Users/Spidey/Desktop/code/sea_wave/run_RIS_wind_comparison.m:25>) |
| `sample` | 运动结果中对应 snapshotTime 的已保存时间点索引 | `find(abs(o.t-snapshotTime)<1e-10,1)` | 赋值 | [sea_wave/run_RIS_wind_comparison.m:26](<C:/Users/Spidey/Desktop/code/sea_wave/run_RIS_wind_comparison.m:26>) |
| `center` | 快照时刻 RIS 质心世界坐标（m） | `o.centerOfMass(:,sample)` | 赋值 | [sea_wave/run_RIS_wind_comparison.m:27](<C:/Users/Spidey/Desktop/code/sea_wave/run_RIS_wind_comparison.m:27>) |
| `summary` | 工况或运动统计汇总表 | `array2table(values,'VariableNames',{'windMps','directionDeg','HsM', 'snapshotCenterXM','snapshotCenterYM','snapshotCenterZM','snapshotRollDeg', 'snapshotPitchDeg','snapshotYawDeg','seaMinM','seaMaxM','seaOutsideZFraction', 'pathMinZM','pathMaxZM','maxTiltDeg','radiationFitRelativeRmsError'})` | 赋值 | [sea_wave/run_RIS_wind_comparison.m:39](<C:/Users/Spidey/Desktop/code/sea_wave/run_RIS_wind_comparison.m:39>) |
| `f` | 图窗句柄；频率等重名按所在函数解释 | `newFigure(880,720)`<br>`newFigure(2000,1600)` | 赋值 | [sea_wave/run_RIS_wind_comparison.m:48](<C:/Users/Spidey/Desktop/code/sea_wave/run_RIS_wind_comparison.m:48>) |
| `ax` | 绘图坐标轴句柄 | `nexttile(tiledlayout(f,1,1,'Padding','compact'))`<br>`nexttile(layout)` | 赋值 | [sea_wave/run_RIS_wind_comparison.m:49](<C:/Users/Spidey/Desktop/code/sea_wave/run_RIS_wind_comparison.m:49>) |
| `name` | 当前输出文件的基础名称 | `sprintf('RIS_sea_U%02d_Phi%02d_t27',winds(j),directions(k))`<br>`sprintf('RIS_center_path_U%02d_Phi30',winds(j))` | 赋值 | [sea_wave/run_RIS_wind_comparison.m:52](<C:/Users/Spidey/Desktop/code/sea_wave/run_RIS_wind_comparison.m:52>) |
| `layout` | 分面绘图布局句柄 | `tiledlayout(f,4,4,'Padding','compact','TileSpacing','compact')` | 赋值 | [sea_wave/run_RIS_wind_comparison.m:62](<C:/Users/Spidey/Desktop/code/sea_wave/run_RIS_wind_comparison.m:62>) |
| `overviewFile` | 4×4 总览 FIG 文件路径 | `saveFigure(f,outputFolder,'RIS_wind_direction_table')` | 赋值 | [sea_wave/run_RIS_wind_comparison.m:74](<C:/Users/Spidey/Desktop/code/sea_wave/run_RIS_wind_comparison.m:74>) |
| `meta` | 记录实际时间、波谱、显示与工况设置的元数据结构 | `struct(...)（所有字段在本节逐行列明）` | 赋值 | [sea_wave/run_RIS_wind_comparison.m:76](<C:/Users/Spidey/Desktop/code/sea_wave/run_RIS_wind_comparison.m:76>) |
| `meta.winds` | 批量仿真的风速列表（m/s） | `winds` | 结构字段 | [sea_wave/run_RIS_wind_comparison.m:76](<C:/Users/Spidey/Desktop/code/sea_wave/run_RIS_wind_comparison.m:76>) |
| `meta.directionsDeg` | 波浪传播方向角（°），从 +Y 朝 +X 增大 | `directions` | 结构字段 | [sea_wave/run_RIS_wind_comparison.m:76](<C:/Users/Spidey/Desktop/code/sea_wave/run_RIS_wind_comparison.m:76>) |
| `meta.seed` | 随机数种子；所属结构决定是海浪、TX、MC 或湍流随机流 | `seed` | 结构字段 | [sea_wave/run_RIS_wind_comparison.m:76](<C:/Users/Spidey/Desktop/code/sea_wave/run_RIS_wind_comparison.m:76>) |
| `meta.startTimeS` | 批量运动的绝对海浪起始时刻（s） | `startTime` | 结构字段 | [sea_wave/run_RIS_wind_comparison.m:77](<C:/Users/Spidey/Desktop/code/sea_wave/run_RIS_wind_comparison.m:77>) |
| `meta.durationS` | 保存的运动持续时间（s） | `duration` | 结构字段 | [sea_wave/run_RIS_wind_comparison.m:77](<C:/Users/Spidey/Desktop/code/sea_wave/run_RIS_wind_comparison.m:77>) |
| `meta.snapshotElapsedTimeS` | 海面/RIS 快照的相对时间（s） | `snapshotTime` | 结构字段 | [sea_wave/run_RIS_wind_comparison.m:77](<C:/Users/Spidey/Desktop/code/sea_wave/run_RIS_wind_comparison.m:77>) |
| `meta.snapshotSeaTimeS` | 快照的绝对海浪时间（s） | `startTime+snapshotTime` | 结构字段 | [sea_wave/run_RIS_wind_comparison.m:78](<C:/Users/Spidey/Desktop/code/sea_wave/run_RIS_wind_comparison.m:78>) |
| `meta.zLimitsM` | 海面高程及轨迹图显示裁切范围（m） | `zLimits` | 结构字段 | [sea_wave/run_RIS_wind_comparison.m:78](<C:/Users/Spidey/Desktop/code/sea_wave/run_RIS_wind_comparison.m:78>) |
| `meta.frequencyBandHz` | 海浪频率范围 [fMinHz,fMaxHz]（Hz） | `[sea.f_min,sea.f_max]` | 结构字段 | [sea_wave/run_RIS_wind_comparison.m:79](<C:/Users/Spidey/Desktop/code/sea_wave/run_RIS_wind_comparison.m:79>) |
| `meta.frequencyBins` | 对数频率离散点数（不是点数减一的区间数） | `sea.Fre_segment` | 结构字段 | [sea_wave/run_RIS_wind_comparison.m:79](<C:/Users/Spidey/Desktop/code/sea_wave/run_RIS_wind_comparison.m:79>) |
| `meta.directionBins` | 方向离散点数，包含 ±π/2 端点 | `sea.Dir_segment` | 结构字段 | [sea_wave/run_RIS_wind_comparison.m:80](<C:/Users/Spidey/Desktop/code/sea_wave/run_RIS_wind_comparison.m:80>) |
| `meta.quadrature` | 频率/方向积分规则，当前仅支持 trapezoid | `'trapezoid'` | 结构字段 | [sea_wave/run_RIS_wind_comparison.m:80](<C:/Users/Spidey/Desktop/code/sea_wave/run_RIS_wind_comparison.m:80>) |
| `meta.amplitudeScale` | 海浪高度显示/生成的附加倍率，固定 1 | `1` | 结构字段 | [sea_wave/run_RIS_wind_comparison.m:80](<C:/Users/Spidey/Desktop/code/sea_wave/run_RIS_wind_comparison.m:80>) |
| `meta.plotBoxAspectRatio` | 三维坐标框显示比例 [1,1,0.7]，未缩放原始高度 | `[1,1,0.7]` | 结构字段 | [sea_wave/run_RIS_wind_comparison.m:81](<C:/Users/Spidey/Desktop/code/sea_wave/run_RIS_wind_comparison.m:81>) |
| `meta.trajectoryPoint` | 轨迹所代表的点，固定为质心 | `'centre of mass'` | 结构字段 | [sea_wave/run_RIS_wind_comparison.m:81](<C:/Users/Spidey/Desktop/code/sea_wave/run_RIS_wind_comparison.m:81>) |
| `meta.displayCropOnly` | 超出 z 显示范围只裁图、不修改数值的标志 | `true` | 结构字段 | [sea_wave/run_RIS_wind_comparison.m:82](<C:/Users/Spidey/Desktop/code/sea_wave/run_RIS_wind_comparison.m:82>) |
| `meta.hydrodynamicallyValidated` | 对应水平绕射/偏航绕射/水平与偏航辐射/完整二阶漂移/直接风载/实测标定/本次网格收敛/水动力验证开关，源码均为 false；字段名区分项目 | `false` | 结构字段 | [sea_wave/run_RIS_wind_comparison.m:82](<C:/Users/Spidey/Desktop/code/sea_wave/run_RIS_wind_comparison.m:82>) |
| `suite.motions` | 全部 RIS 运动结果/海面快照元胞数组 | `{motions}` | 结构字段 | [sea_wave/run_RIS_wind_comparison.m:85](<C:/Users/Spidey/Desktop/code/sea_wave/run_RIS_wind_comparison.m:85>) |
| `suite.scenes` | 全部 RIS 运动结果/海面快照元胞数组 | `{scenes}` | 结构字段 | [sea_wave/run_RIS_wind_comparison.m:85](<C:/Users/Spidey/Desktop/code/sea_wave/run_RIS_wind_comparison.m:85>) |
| `suite.summary` | 工况统计表 | `summary` | 结构字段 | [sea_wave/run_RIS_wind_comparison.m:85](<C:/Users/Spidey/Desktop/code/sea_wave/run_RIS_wind_comparison.m:85>) |
| `suite.meta` | 批量仿真的参数与显示元数据 | `meta` | 结构字段 | [sea_wave/run_RIS_wind_comparison.m:86](<C:/Users/Spidey/Desktop/code/sea_wave/run_RIS_wind_comparison.m:86>) |
| `suite.panelFiles` | 单元 FIG 文件路径数组/总览 FIG 文件路径 | `panelFiles` | 结构字段 | [sea_wave/run_RIS_wind_comparison.m:86](<C:/Users/Spidey/Desktop/code/sea_wave/run_RIS_wind_comparison.m:86>) |
| `suite.overviewFile` | 单元 FIG 文件路径数组/总览 FIG 文件路径 | `overviewFile` | 结构字段 | [sea_wave/run_RIS_wind_comparison.m:86](<C:/Users/Spidey/Desktop/code/sea_wave/run_RIS_wind_comparison.m:86>) |


**作用域：`drawCell`**

| 源码符号 | 含义与单位/维度 | 源码取值/计算式（各次赋值） | 角色 | 源码位置 |
| --- | --- | --- | --- | --- |
| `ax` | 绘图坐标轴句柄 | 函数输入/输出；值随调用或上游结果确定 | 输入 | [sea_wave/run_RIS_wind_comparison.m:88](<C:/Users/Spidey/Desktop/code/sea_wave/run_RIS_wind_comparison.m:88>) |
| `row` | 4×4 总图中的风速行/面板列；第 4 列是风向 30° 的轨迹 | 函数输入/输出；值随调用或上游结果确定 | 输入 | [sea_wave/run_RIS_wind_comparison.m:88](<C:/Users/Spidey/Desktop/code/sea_wave/run_RIS_wind_comparison.m:88>) |
| `column` | 4×4 总图中的风速行/面板列；第 4 列是风向 30° 的轨迹 | 函数输入/输出；值随调用或上游结果确定 | 输入 | [sea_wave/run_RIS_wind_comparison.m:88](<C:/Users/Spidey/Desktop/code/sea_wave/run_RIS_wind_comparison.m:88>) |
| `motions` | 全部 RIS 运动结果/海面快照元胞数组 | 函数输入/输出；值随调用或上游结果确定 | 输入 | [sea_wave/run_RIS_wind_comparison.m:88](<C:/Users/Spidey/Desktop/code/sea_wave/run_RIS_wind_comparison.m:88>) |
| `scenes` | 海面显示快照结构（X/Y/Z）或其集合 | 函数输入/输出；值随调用或上游结果确定 | 输入 | [sea_wave/run_RIS_wind_comparison.m:88](<C:/Users/Spidey/Desktop/code/sea_wave/run_RIS_wind_comparison.m:88>) |
| `time` | 相对 parameters.startTime 的时间（s）或其采样向量 | 函数输入/输出；值随调用或上游结果确定 | 输入 | [sea_wave/run_RIS_wind_comparison.m:88](<C:/Users/Spidey/Desktop/code/sea_wave/run_RIS_wind_comparison.m:88>) |
| `zLimits` | 海面与质心轨迹的显示裁切高度范围 [-0.5,0.5]（m） | 函数输入/输出；值随调用或上游结果确定 | 输入 | [sea_wave/run_RIS_wind_comparison.m:88](<C:/Users/Spidey/Desktop/code/sea_wave/run_RIS_wind_comparison.m:88>) |
| `o` | 当前 RIS 运动结果或动力学模型；具体类型由所在函数输入确定 | `motions{row,column}`<br>`motions{row,1}` | 赋值 | [sea_wave/run_RIS_wind_comparison.m:91](<C:/Users/Spidey/Desktop/code/sea_wave/run_RIS_wind_comparison.m:91>) |
| `p` | 当前轨迹的质心世界位置数组（3×Nt，m）；此处不是配置 | `o.centerOfMass` | 赋值 | [sea_wave/run_RIS_wind_comparison.m:97](<C:/Users/Spidey/Desktop/code/sea_wave/run_RIS_wind_comparison.m:97>) |
| `inRange` | 质心高度处于 zLimits 内的样本索引 | `find(p(3,:)>=zLimits(1) & p(3,:)<=zLimits(2))` | 赋值 | [sea_wave/run_RIS_wind_comparison.m:103](<C:/Users/Spidey/Desktop/code/sea_wave/run_RIS_wind_comparison.m:103>) |
| `outside` | 质心高度超出显示范围的时间样本比例 | `mean(abs(p(3,:))>zLimits(2))` | 赋值 | [sea_wave/run_RIS_wind_comparison.m:113](<C:/Users/Spidey/Desktop/code/sea_wave/run_RIS_wind_comparison.m:113>) |


**作用域：`newFigure`**

| 源码符号 | 含义与单位/维度 | 源码取值/计算式（各次赋值） | 角色 | 源码位置 |
| --- | --- | --- | --- | --- |
| `width` | 输出图像宽、高（pixel） | 函数输入/输出；值随调用或上游结果确定 | 输入 | [sea_wave/run_RIS_wind_comparison.m:118](<C:/Users/Spidey/Desktop/code/sea_wave/run_RIS_wind_comparison.m:118>) |
| `height` | 输出图像宽、高（pixel） | 函数输入/输出；值随调用或上游结果确定 | 输入 | [sea_wave/run_RIS_wind_comparison.m:118](<C:/Users/Spidey/Desktop/code/sea_wave/run_RIS_wind_comparison.m:118>) |
| `f` | 图窗句柄；频率等重名按所在函数解释 | `figure('Color','w','Visible','off','Units','pixels','Position',[40,40,width,height])` | 输出/赋值 | [sea_wave/run_RIS_wind_comparison.m:118](<C:/Users/Spidey/Desktop/code/sea_wave/run_RIS_wind_comparison.m:118>) |


**作用域：`saveFigure`**

| 源码符号 | 含义与单位/维度 | 源码取值/计算式（各次赋值） | 角色 | 源码位置 |
| --- | --- | --- | --- | --- |
| `f` | 图窗句柄；频率等重名按所在函数解释 | 函数输入/输出；值随调用或上游结果确定 | 输入 | [sea_wave/run_RIS_wind_comparison.m:122](<C:/Users/Spidey/Desktop/code/sea_wave/run_RIS_wind_comparison.m:122>) |
| `folder` | 当前源码所在目录 | 函数输入/输出；值随调用或上游结果确定 | 输入 | [sea_wave/run_RIS_wind_comparison.m:122](<C:/Users/Spidey/Desktop/code/sea_wave/run_RIS_wind_comparison.m:122>) |
| `name` | 当前输出文件的基础名称 | 函数输入/输出；值随调用或上游结果确定 | 输入 | [sea_wave/run_RIS_wind_comparison.m:122](<C:/Users/Spidey/Desktop/code/sea_wave/run_RIS_wind_comparison.m:122>) |
| `path` | 生成或读取的文件完整路径 | `string(fullfile(folder,[name,'.fig']))` | 输出/赋值 | [sea_wave/run_RIS_wind_comparison.m:122](<C:/Users/Spidey/Desktop/code/sea_wave/run_RIS_wind_comparison.m:122>) |


**作用域：`exportTrajectory`**

| 源码符号 | 含义与单位/维度 | 源码取值/计算式（各次赋值） | 角色 | 源码位置 |
| --- | --- | --- | --- | --- |
| `o` | 当前 RIS 运动结果或动力学模型；具体类型由所在函数输入确定 | 函数输入/输出；值随调用或上游结果确定 | 输入 | [sea_wave/run_RIS_wind_comparison.m:128](<C:/Users/Spidey/Desktop/code/sea_wave/run_RIS_wind_comparison.m:128>) |
| `folder` | 当前源码所在目录 | 函数输入/输出；值随调用或上游结果确定 | 输入 | [sea_wave/run_RIS_wind_comparison.m:128](<C:/Users/Spidey/Desktop/code/sea_wave/run_RIS_wind_comparison.m:128>) |
| `wind` | 当前轨迹风速（m/s），用于输出文件名 | 函数输入/输出；值随调用或上游结果确定 | 输入 | [sea_wave/run_RIS_wind_comparison.m:128](<C:/Users/Spidey/Desktop/code/sea_wave/run_RIS_wind_comparison.m:128>) |
| `direction` | 当前轨迹风向角（°），用于输出文件名 | 函数输入/输出；值随调用或上游结果确定 | 输入 | [sea_wave/run_RIS_wind_comparison.m:128](<C:/Users/Spidey/Desktop/code/sea_wave/run_RIS_wind_comparison.m:128>) |
| `values` | 14 列完整轨迹数据：相对/绝对时刻、质心、欧拉角、体速度 | `[o.t.',o.absoluteTime.',o.centerOfMass.',o.eta(4:6,:).',o.bodyVelocity.']` | 赋值 | [sea_wave/run_RIS_wind_comparison.m:130](<C:/Users/Spidey/Desktop/code/sea_wave/run_RIS_wind_comparison.m:130>) |
| `names` | 上述 14 列 CSV 表头；带 S/M/Rad/Mps/Radps 的后缀直接标示单位 | `{'elapsedTimeS','seaTimeS','centerXM','centerYM','centerZM', 'rollRad','pitchRad','yawRad','uMps','vMps','wMps','pRadps','qRadps','rRadps'}` | 赋值 | [sea_wave/run_RIS_wind_comparison.m:131](<C:/Users/Spidey/Desktop/code/sea_wave/run_RIS_wind_comparison.m:131>) |


### 5.15 sea_wave/summarize_RIS_motion_existing.m

**作用域：`script`**

| 源码符号 | 含义与单位/维度 | 源码取值/计算式（各次赋值） | 角色 | 源码位置 |
| --- | --- | --- | --- | --- |
| `risStatsFolder` | 当前源码所在目录 | `fileparts(mfilename('fullpath'))` | 赋值 | [sea_wave/summarize_RIS_motion_existing.m:4](<C:/Users/Spidey/Desktop/code/sea_wave/summarize_RIS_motion_existing.m:4>) |
| `risInputFolder` | 已有 12 条海面/RIS 轨迹 CSV 的输入目录 | `fullfile(risStatsFolder,'FigOutput_RIS_WindDirection')` | 赋值 | [sea_wave/summarize_RIS_motion_existing.m:5](<C:/Users/Spidey/Desktop/code/sea_wave/summarize_RIS_motion_existing.m:5>) |
| `risOutputFolder` | 计算结果与图表输出目录 | `fullfile(risStatsFolder,'motion_parameterization')` | 赋值 | [sea_wave/summarize_RIS_motion_existing.m:6](<C:/Users/Spidey/Desktop/code/sea_wave/summarize_RIS_motion_existing.m:6>) |
| `risFiles` | 匹配 RIS_trajectory_U*_Phi*.csv 的文件信息数组 | `dir(fullfile(risInputFolder,'RIS_trajectory_U*_Phi*.csv'))` | 赋值 | [sea_wave/summarize_RIS_motion_existing.m:8](<C:/Users/Spidey/Desktop/code/sea_wave/summarize_RIS_motion_existing.m:8>) |
| `risRows` | 逐轨迹计算的运动统计结构数组 | `struct([])` | 赋值 | [sea_wave/summarize_RIS_motion_existing.m:10](<C:/Users/Spidey/Desktop/code/sea_wave/summarize_RIS_motion_existing.m:10>) |
| `risIndex` | 当前输入轨迹 CSV 的序号 | `1:numel(risFiles)` | 循环索引 | [sea_wave/summarize_RIS_motion_existing.m:11](<C:/Users/Spidey/Desktop/code/sea_wave/summarize_RIS_motion_existing.m:11>) |
| `risToken` | 由文件名解析出的风速与风向字符串 | `regexp(risFiles(risIndex).name,'U(\d+)_Phi(\d+)','tokens','once')` | 赋值 | [sea_wave/summarize_RIS_motion_existing.m:12](<C:/Users/Spidey/Desktop/code/sea_wave/summarize_RIS_motion_existing.m:12>) |
| `risT` | 读入的一条原始轨迹数据表 | `readtable(fullfile(risFiles(risIndex).folder,risFiles(risIndex).name),TextType="string")` | 赋值 | [sea_wave/summarize_RIS_motion_existing.m:13](<C:/Users/Spidey/Desktop/code/sea_wave/summarize_RIS_motion_existing.m:13>) |
| `TextType` | readtable 的文本类型选项，取 string；不是物理参数 | `"string")` | 赋值 | [sea_wave/summarize_RIS_motion_existing.m:13](<C:/Users/Spidey/Desktop/code/sea_wave/summarize_RIS_motion_existing.m:13>) |
| `risTime` | 当前轨迹相对时间向量（s） | `risT.elapsedTimeS` | 赋值 | [sea_wave/summarize_RIS_motion_existing.m:14](<C:/Users/Spidey/Desktop/code/sea_wave/summarize_RIS_motion_existing.m:14>) |
| `risP` | 质心世界位置样本矩阵（Nt×3，m） | `risT{:,{'centerXM','centerYM','centerZM'}}` | 赋值 | [sea_wave/summarize_RIS_motion_existing.m:15](<C:/Users/Spidey/Desktop/code/sea_wave/summarize_RIS_motion_existing.m:15>) |
| `risAngles` | 欧拉角样本矩阵 [roll,pitch,yaw]（rad） | `risT{:,{'rollRad','pitchRad','yawRad'}}` | 赋值 | [sea_wave/summarize_RIS_motion_existing.m:16](<C:/Users/Spidey/Desktop/code/sea_wave/summarize_RIS_motion_existing.m:16>) |
| `risC` | 按 N 归一化的质心位置总体协方差矩阵（m²） | `cov(risP,1)` | 赋值 | [sea_wave/summarize_RIS_motion_existing.m:20](<C:/Users/Spidey/Desktop/code/sea_wave/summarize_RIS_motion_existing.m:20>) |
| `risCentered` | 质心坐标减各列时间均值后的样本（m），保留线性趋势 | `risP-mean(risP,1)` | 赋值 | [sea_wave/summarize_RIS_motion_existing.m:21](<C:/Users/Spidey/Desktop/code/sea_wave/summarize_RIS_motion_existing.m:21>) |
| `risPhi` | roll/pitch/yaw 三列欧拉角（rad），按变量顺序对应 | `risAngles(:,1)` | 赋值 | [sea_wave/summarize_RIS_motion_existing.m:23](<C:/Users/Spidey/Desktop/code/sea_wave/summarize_RIS_motion_existing.m:23>) |
| `risTheta` | roll/pitch/yaw 三列欧拉角（rad），按变量顺序对应 | `risAngles(:,2)` | 赋值 | [sea_wave/summarize_RIS_motion_existing.m:23](<C:/Users/Spidey/Desktop/code/sea_wave/summarize_RIS_motion_existing.m:23>) |
| `risPsi` | roll/pitch/yaw 三列欧拉角（rad），按变量顺序对应 | `risAngles(:,3)` | 赋值 | [sea_wave/summarize_RIS_motion_existing.m:23](<C:/Users/Spidey/Desktop/code/sea_wave/summarize_RIS_motion_existing.m:23>) |
| `risNormals` | 由 ZYX 欧拉角重建的单位板法向量样本（Nt×3） | `[cos(risPsi).*sin(risTheta).*cos(risPhi)+sin(risPsi).*sin(risPhi), sin(risPsi).*sin(risTheta).*cos(risPhi)-cos(risPsi).*sin(risPhi), cos(risTheta).*cos(risPhi)]` | 赋值 | [sea_wave/summarize_RIS_motion_existing.m:24](<C:/Users/Spidey/Desktop/code/sea_wave/summarize_RIS_motion_existing.m:24>) |
| `risMeanNormal` | 各样本法向平均后再次单位化的参考法向 | `mean(risNormals,1)`<br>`risMeanNormal/norm(risMeanNormal)` | 赋值 | [sea_wave/summarize_RIS_motion_existing.m:28](<C:/Users/Spidey/Desktop/code/sea_wave/summarize_RIS_motion_existing.m:28>) |
| `risNormalAngles` | 瞬时法向与平均法向的夹角（rad） | `acos(max(-1,min(1,risNormals*risMeanNormal.')))` | 赋值 | [sea_wave/summarize_RIS_motion_existing.m:29](<C:/Users/Spidey/Desktop/code/sea_wave/summarize_RIS_motion_existing.m:29>) |
| `risEulerStd` | 按 N 归一化的 roll/pitch/yaw 总体标准差（°） | `std(risAngles,1,1)*180/pi` | 赋值 | [sea_wave/summarize_RIS_motion_existing.m:30](<C:/Users/Spidey/Desktop/code/sea_wave/summarize_RIS_motion_existing.m:30>) |
| `risTrend` | 三维位置对 [1,time] 的线性回归系数（第一行 m，第二行 m/s） | `[ones(size(risTime)),risTime]\risP` | 赋值 | [sea_wave/summarize_RIS_motion_existing.m:31](<C:/Users/Spidey/Desktop/code/sea_wave/summarize_RIS_motion_existing.m:31>) |
| `risResidual` | 去除线性回归趋势后的三维位置残差（m） | `risP-[ones(size(risTime)),risTime]*risTrend` | 赋值 | [sea_wave/summarize_RIS_motion_existing.m:32](<C:/Users/Spidey/Desktop/code/sea_wave/summarize_RIS_motion_existing.m:32>) |
| `risRows.windMps` | 风速（m/s） | `str2double(risToken{1})` | 赋值 | [sea_wave/summarize_RIS_motion_existing.m:33](<C:/Users/Spidey/Desktop/code/sea_wave/summarize_RIS_motion_existing.m:33>) |
| `risRows.directionDeg` | 波浪传播方向角（°），从 +Y 朝 +X 增大 | `str2double(risToken{2})` | 赋值 | [sea_wave/summarize_RIS_motion_existing.m:34](<C:/Users/Spidey/Desktop/code/sea_wave/summarize_RIS_motion_existing.m:34>) |
| `risRows.sampleCount` | EGG 信道快照样本数或运动记录样本数 | `height(risT)` | 赋值 | [sea_wave/summarize_RIS_motion_existing.m:35](<C:/Users/Spidey/Desktop/code/sea_wave/summarize_RIS_motion_existing.m:35>) |
| `risRows.varXM2` | 去均值、保留趋势、按 N 归一化的质心 x/y/z 位移方差（m²） | `risC(1,1)` | 赋值 | [sea_wave/summarize_RIS_motion_existing.m:36](<C:/Users/Spidey/Desktop/code/sea_wave/summarize_RIS_motion_existing.m:36>) |
| `risRows.varYM2` | 去均值、保留趋势、按 N 归一化的质心 x/y/z 位移方差（m²） | `risC(2,2)` | 赋值 | [sea_wave/summarize_RIS_motion_existing.m:37](<C:/Users/Spidey/Desktop/code/sea_wave/summarize_RIS_motion_existing.m:37>) |
| `risRows.varZM2` | 去均值、保留趋势、按 N 归一化的质心 x/y/z 位移方差（m²） | `risC(3,3)` | 赋值 | [sea_wave/summarize_RIS_motion_existing.m:38](<C:/Users/Spidey/Desktop/code/sea_wave/summarize_RIS_motion_existing.m:38>) |
| `risRows.covXYM2` | 质心相应两坐标的协方差（m²），按 N 归一化 | `risC(1,2)` | 赋值 | [sea_wave/summarize_RIS_motion_existing.m:39](<C:/Users/Spidey/Desktop/code/sea_wave/summarize_RIS_motion_existing.m:39>) |
| `risRows.covXZM2` | 质心相应两坐标的协方差（m²），按 N 归一化 | `risC(1,3)` | 赋值 | [sea_wave/summarize_RIS_motion_existing.m:40](<C:/Users/Spidey/Desktop/code/sea_wave/summarize_RIS_motion_existing.m:40>) |
| `risRows.covYZM2` | 质心相应两坐标的协方差（m²），按 N 归一化 | `risC(2,3)` | 赋值 | [sea_wave/summarize_RIS_motion_existing.m:41](<C:/Users/Spidey/Desktop/code/sea_wave/summarize_RIS_motion_existing.m:41>) |
| `risRows.sigma3DM` | 三维质心离散度 sqrt(trace(covariance))（m） | `sqrt(trace(risC))` | 赋值 | [sea_wave/summarize_RIS_motion_existing.m:42](<C:/Users/Spidey/Desktop/code/sea_wave/summarize_RIS_motion_existing.m:42>) |
| `risRows.sigmaHorizontalM` | 水平质心离散度 sqrt(varX+varY)（m） | `sqrt(risC(1,1)+risC(2,2))` | 赋值 | [sea_wave/summarize_RIS_motion_existing.m:43](<C:/Users/Spidey/Desktop/code/sea_wave/summarize_RIS_motion_existing.m:43>) |
| `risRows.stdRollDeg` | 欧拉角总体标准差（°），按 N 归一化 | `risEulerStd(1)` | 赋值 | [sea_wave/summarize_RIS_motion_existing.m:44](<C:/Users/Spidey/Desktop/code/sea_wave/summarize_RIS_motion_existing.m:44>) |
| `risRows.stdPitchDeg` | 欧拉角总体标准差（°），按 N 归一化 | `risEulerStd(2)` | 赋值 | [sea_wave/summarize_RIS_motion_existing.m:45](<C:/Users/Spidey/Desktop/code/sea_wave/summarize_RIS_motion_existing.m:45>) |
| `risRows.stdYawDeg` | 欧拉角总体标准差（°），按 N 归一化 | `risEulerStd(3)` | 赋值 | [sea_wave/summarize_RIS_motion_existing.m:46](<C:/Users/Spidey/Desktop/code/sea_wave/summarize_RIS_motion_existing.m:46>) |
| `risRows.normalAngularRmsDeg` | 瞬时法向与平均单位法向夹角的 RMS（°） | `sqrt(mean(risNormalAngles.^2))*180/pi` | 赋值 | [sea_wave/summarize_RIS_motion_existing.m:47](<C:/Users/Spidey/Desktop/code/sea_wave/summarize_RIS_motion_existing.m:47>) |
| `risRows.linearDriftSpeedMps` | 位置对时间线性回归的三维斜率模长（m/s） | `norm(risTrend(2,:))` | 赋值 | [sea_wave/summarize_RIS_motion_existing.m:48](<C:/Users/Spidey/Desktop/code/sea_wave/summarize_RIS_motion_existing.m:48>) |
| `risRows.sigma3DLinearDetrendedM` | 去线性趋势后三维质心残差 RMS（m） | `sqrt(mean(sum(risResidual.^2,2)))` | 赋值 | [sea_wave/summarize_RIS_motion_existing.m:49](<C:/Users/Spidey/Desktop/code/sea_wave/summarize_RIS_motion_existing.m:49>) |
| `risRows.meanXM` | 三维质心相应坐标的时间样本均值（m） | `mean(risP(:,1))` | 赋值 | [sea_wave/summarize_RIS_motion_existing.m:50](<C:/Users/Spidey/Desktop/code/sea_wave/summarize_RIS_motion_existing.m:50>) |
| `risRows.meanYM` | 三维质心相应坐标的时间样本均值（m） | `mean(risP(:,2))` | 赋值 | [sea_wave/summarize_RIS_motion_existing.m:51](<C:/Users/Spidey/Desktop/code/sea_wave/summarize_RIS_motion_existing.m:51>) |
| `risRows.meanZM` | 三维质心相应坐标的时间样本均值（m） | `mean(risP(:,3))` | 赋值 | [sea_wave/summarize_RIS_motion_existing.m:52](<C:/Users/Spidey/Desktop/code/sea_wave/summarize_RIS_motion_existing.m:52>) |
| `risMotionSummary` | 按 directionDeg、windMps 排序的全部运动统计表 | `sortrows(struct2table(risRows),{'directionDeg','windMps'})` | 赋值 | [sea_wave/summarize_RIS_motion_existing.m:54](<C:/Users/Spidey/Desktop/code/sea_wave/summarize_RIS_motion_existing.m:54>) |
| `risShow` | 仅保留风向 30° 的用于绘图的统计子表 | `risMotionSummary(risMotionSummary.directionDeg==30,:)` | 赋值 | [sea_wave/summarize_RIS_motion_existing.m:56](<C:/Users/Spidey/Desktop/code/sea_wave/summarize_RIS_motion_existing.m:56>) |
| `risFig` | 图窗句柄；频率等重名按所在函数解释 | `figure('Color','w','Visible','off','Position',[80,80,1450,580])` | 赋值 | [sea_wave/summarize_RIS_motion_existing.m:57](<C:/Users/Spidey/Desktop/code/sea_wave/summarize_RIS_motion_existing.m:57>) |
| `risLayout` | 分面绘图布局句柄 | `tiledlayout(risFig,1,2,'TileSpacing','compact','Padding','compact')` | 赋值 | [sea_wave/summarize_RIS_motion_existing.m:58](<C:/Users/Spidey/Desktop/code/sea_wave/summarize_RIS_motion_existing.m:58>) |
| `risAx` | 绘图坐标轴句柄 | `nexttile(risLayout)` | 赋值 | [sea_wave/summarize_RIS_motion_existing.m:59](<C:/Users/Spidey/Desktop/code/sea_wave/summarize_RIS_motion_existing.m:59>) |


### 5.16 W_RIS/run_wris.m

**作用域：`script`**

| 源码符号 | 含义与单位/维度 | 源码取值/计算式（各次赋值） | 角色 | 源码位置 |
| --- | --- | --- | --- | --- |
| `wris_results` | production 档完整流水线的输出结果结构 | `wris_run('production')` | 赋值 | [W_RIS/run_wris.m:3](<C:/Users/Spidey/Desktop/code/W_RIS/run_wris.m:3>) |


### 5.17 W_RIS/wris_cached.m

**作用域：`wris_cached`**

| 源码符号 | 含义与单位/维度 | 源码取值/计算式（各次赋值） | 角色 | 源码位置 |
| --- | --- | --- | --- | --- |
| `c` | 当前 W_RIS 总配置结构 | 函数输入/输出；值随调用或上游结果确定 | 输入 | [W_RIS/wris_cached.m:1](<C:/Users/Spidey/Desktop/code/W_RIS/wris_cached.m:1>) |
| `kind` | 物理缓存类别：motion 或 optics | 函数输入/输出；值随调用或上游结果确定 | 输入 | [W_RIS/wris_cached.m:1](<C:/Users/Spidey/Desktop/code/W_RIS/wris_cached.m:1>) |
| `key` | 生成当前物理结果的所有相关输入结构 | 函数输入/输出；值随调用或上游结果确定 | 输入 | [W_RIS/wris_cached.m:1](<C:/Users/Spidey/Desktop/code/W_RIS/wris_cached.m:1>) |
| `producer` | 缓存未命中时执行的计算函数句柄 | 函数输入/输出；值随调用或上游结果确定 | 输入 | [W_RIS/wris_cached.m:1](<C:/Users/Spidey/Desktop/code/W_RIS/wris_cached.m:1>) |
| `value` | 缓存加载或新计算得到的物理结果 | `producer()`<br>`data.value` | 输出/赋值 | [W_RIS/wris_cached.m:1](<C:/Users/Spidey/Desktop/code/W_RIS/wris_cached.m:1>) |
| `parent` | W_RIS 的父目录，即工程目录 | `fileparts(c.root)` | 赋值 | [W_RIS/wris_cached.m:7](<C:/Users/Spidey/Desktop/code/W_RIS/wris_cached.m:7>) |
| `names` | motion 缓存依赖的八个 sea_wave 源码文件名 | `{'RIS_float.m','RIS_wave_components.m','RIS_wave_field.m','RIS_rotation.m', 'RIS_radiation.m','RIS_sidewall_load.m','get_default_RIS_params.m','get_default_sea_params.m'}` | 赋值 | [W_RIS/wris_cached.m:10](<C:/Users/Spidey/Desktop/code/W_RIS/wris_cached.m:10>) |
| `files` | 当前缓存应纳入哈希的物理源码完整路径集合 | `cellfun(@(n) fullfile(parent,'sea_wave',n),names,'UniformOutput',false)`<br>`{fullfile(c.root,'wris_mc.m'),fullfile(c.root,'wris_intersect.m'), fullfile(c.root,'wris_metrics.m'),fullfile(parent,'fig3-7','Hen_Green.m')}` | 赋值 | [W_RIS/wris_cached.m:12](<C:/Users/Spidey/Desktop/code/W_RIS/wris_cached.m:12>) |
| `n` | cellfun 匿名函数中的源码文件名，不是采样数 | `files=cellfun(@(n) fullfile(parent,'sea_wave',n),names,'UniformOutput',false)` | 匿名函数输入 | [W_RIS/wris_cached.m:12](<C:/Users/Spidey/Desktop/code/W_RIS/wris_cached.m:12>) |
| `digest` | Java SHA-256 摘要对象，对输入 JSON 和源码 UTF-8 内容计算哈希 | `java.security.MessageDigest.getInstance('SHA-256')` | 赋值 | [W_RIS/wris_cached.m:19](<C:/Users/Spidey/Desktop/code/W_RIS/wris_cached.m:19>) |
| `j` | 当前循环、分量或样本索引（无量纲，MATLAB 从 1 开始） | `1:numel(files)` | 循环索引 | [W_RIS/wris_cached.m:21](<C:/Users/Spidey/Desktop/code/W_RIS/wris_cached.m:21>) |
| `hash` | 小写十六进制 SHA-256 摘要字符串 | `lower(reshape(dec2hex(typecast(digest.digest(),'uint8'),2).',1,[]))` | 赋值 | [W_RIS/wris_cached.m:24](<C:/Users/Spidey/Desktop/code/W_RIS/wris_cached.m:24>) |
| `file` | 生成或读取的文件完整路径 | `fullfile(c.cache.folder,[kind,'_',hash,'.mat'])` | 赋值 | [W_RIS/wris_cached.m:26](<C:/Users/Spidey/Desktop/code/W_RIS/wris_cached.m:26>) |
| `data` | 从缓存 MAT 文件读入的结构，实际结果位于 data.value | `load(file,'value')` | 赋值 | [W_RIS/wris_cached.m:28](<C:/Users/Spidey/Desktop/code/W_RIS/wris_cached.m:28>) |


**作用域：`serializable`**

| 源码符号 | 含义与单位/维度 | 源码取值/计算式（各次赋值） | 角色 | 源码位置 |
| --- | --- | --- | --- | --- |
| `value` | 递归转换中的缓存键对象；复数改为 realPart/imaginaryPart 结构以保留相位 | `serializable(value.(names{j}))`<br>`cellfun(@serializable,value,'UniformOutput',false)`<br>`struct('realPart',real(value),'imaginaryPart',imag(value))` | 输入/输出/赋值 | [W_RIS/wris_cached.m:34](<C:/Users/Spidey/Desktop/code/W_RIS/wris_cached.m:34>) |
| `names` | 当前结构所有字段名 | `fieldnames(value)` | 赋值 | [W_RIS/wris_cached.m:38](<C:/Users/Spidey/Desktop/code/W_RIS/wris_cached.m:38>) |
| `j` | 当前循环、分量或样本索引（无量纲，MATLAB 从 1 开始） | `1:numel(names)` | 循环索引 | [W_RIS/wris_cached.m:39](<C:/Users/Spidey/Desktop/code/W_RIS/wris_cached.m:39>) |
| `value.realPart` | 缓存 JSON 保留复数波幅的实部/虚部数组 | `real(value)` | 结构字段 | [W_RIS/wris_cached.m:43](<C:/Users/Spidey/Desktop/code/W_RIS/wris_cached.m:43>) |
| `value.imaginaryPart` | 缓存 JSON 保留复数波幅的实部/虚部数组 | `imag(value)` | 结构字段 | [W_RIS/wris_cached.m:43](<C:/Users/Spidey/Desktop/code/W_RIS/wris_cached.m:43>) |


### 5.18 W_RIS/wris_config.m

**作用域：`wris_config`**

| 源码符号 | 含义与单位/维度 | 源码取值/计算式（各次赋值） | 角色 | 源码位置 |
| --- | --- | --- | --- | --- |
| `profile` | 运行档位：production（默认）或 smoke | `'production'` | 输入/赋值 | [W_RIS/wris_config.m:1](<C:/Users/Spidey/Desktop/code/W_RIS/wris_config.m:1>) |
| `c` | 全部 W_RIS 场景、物理、数值、显示及输出配置结构 | 函数输入/输出；值随调用或上游结果确定 | 输出 | [W_RIS/wris_config.m:1](<C:/Users/Spidey/Desktop/code/W_RIS/wris_config.m:1>) |
| `c.version` | 配置版本标记字符串 | `'wris-3-rotation-table'` | 赋值 | [W_RIS/wris_config.m:4](<C:/Users/Spidey/Desktop/code/W_RIS/wris_config.m:4>) |
| `c.root` | 源码 W_RIS 目录 | `fileparts(mfilename('fullpath'))` | 赋值 | [W_RIS/wris_config.m:5](<C:/Users/Spidey/Desktop/code/W_RIS/wris_config.m:5>) |
| `c.geometry.distance` | TX 至 RIS 光学中心的空间直线距离（m） | `10` | 赋值 | [W_RIS/wris_config.m:6](<C:/Users/Spidey/Desktop/code/W_RIS/wris_config.m:6>) |
| `c.geometry.incidenceDeg` | 静水初始光轴相对 RIS 法线的入射角（°） | `60` | 赋值 | [W_RIS/wris_config.m:7](<C:/Users/Spidey/Desktop/code/W_RIS/wris_config.m:7>) |
| `c.geometry.azimuthDeg` | TX 光轴水平投影相对世界 +X 的方位角（°） | `90` | 赋值 | [W_RIS/wris_config.m:8](<C:/Users/Spidey/Desktop/code/W_RIS/wris_config.m:8>) |
| `c.ris` | 继承 get_default_RIS_params 后由 wris_config 覆盖的 RIS/海浪/动力学配置 | `get_default_RIS_params()` | 赋值 | [W_RIS/wris_config.m:9](<C:/Users/Spidey/Desktop/code/W_RIS/wris_config.m:9>) |
| `c.ris.duration` | 保存的运动持续时间（s） | `20` | 赋值 | [W_RIS/wris_config.m:10](<C:/Users/Spidey/Desktop/code/W_RIS/wris_config.m:10>) |
| `c.ris.dt` | 运动输出采样间隔（s），不等于实际 RK4 内部步长 | `0.05` | 赋值 | [W_RIS/wris_config.m:10](<C:/Users/Spidey/Desktop/code/W_RIS/wris_config.m:10>) |
| `c.ris.startTime` | 绝对海浪起算时刻（s） | `25` | 赋值 | [W_RIS/wris_config.m:10](<C:/Users/Spidey/Desktop/code/W_RIS/wris_config.m:10>) |
| `c.ris.fMinHz` | 海浪最低频率（Hz） | `0.02` | 赋值 | [W_RIS/wris_config.m:11](<C:/Users/Spidey/Desktop/code/W_RIS/wris_config.m:11>) |
| `c.ris.fMaxHz` | 海浪最高频率（Hz） | `3` | 赋值 | [W_RIS/wris_config.m:11](<C:/Users/Spidey/Desktop/code/W_RIS/wris_config.m:11>) |
| `c.ris.frequencyBins` | 对数频率离散点数（不是点数减一的区间数） | `200` | 赋值 | [W_RIS/wris_config.m:12](<C:/Users/Spidey/Desktop/code/W_RIS/wris_config.m:12>) |
| `c.ris.directionBins` | 方向离散点数，包含 ±π/2 端点 | `36` | 赋值 | [W_RIS/wris_config.m:12](<C:/Users/Spidey/Desktop/code/W_RIS/wris_config.m:12>) |
| `c.ris.seed` | 随机数种子；所属结构决定是海浪、TX、MC 或湍流随机流 | `600` | 赋值 | [W_RIS/wris_config.m:13](<C:/Users/Spidey/Desktop/code/W_RIS/wris_config.m:13>) |
| `c.ris.plot` | 是否自动绘制 RIS 单次运动；RIS_float 无输出且未覆盖时设为 true | `false` | 赋值 | [W_RIS/wris_config.m:13](<C:/Users/Spidey/Desktop/code/W_RIS/wris_config.m:13>) |
| `c.comparison.incidenceDeg` | 静水初始光轴相对 RIS 法线的入射角（°） | `[60 45 30]` | 赋值 | [W_RIS/wris_config.m:15](<C:/Users/Spidey/Desktop/code/W_RIS/wris_config.m:15>) |
| `c.comparison.conditions` | 四个运动列的 [风速(m/s),传播风向(°)] 工况矩阵 | `[3 30;3 50;5 50;9 30]` | 赋值 | [W_RIS/wris_config.m:16](<C:/Users/Spidey/Desktop/code/W_RIS/wris_config.m:16>) |
| `c.comparison.phaseIndex` | 3×4 站位选择索引：1峰、2谷、3空间上升坡、4空间下降坡 | `[1 1 1 1;2 2 2 2;3 4 3 4]` | 赋值 | [W_RIS/wris_config.m:17](<C:/Users/Spidey/Desktop/code/W_RIS/wris_config.m:17>) |
| `c.sampling.comparisonTime` | 所有面板统一的相对采样时刻（s） | `20` | 赋值 | [W_RIS/wris_config.m:18](<C:/Users/Spidey/Desktop/code/W_RIS/wris_config.m:18>) |
| `c.sea.quadrature` | 频率/方向积分规则，当前仅支持 trapezoid | `'trapezoid'` | 赋值 | [W_RIS/wris_config.m:19](<C:/Users/Spidey/Desktop/code/W_RIS/wris_config.m:19>) |
| `c.sea.siteSamples` | 沿主传播方向的空间剖面采样点数 | `801` | 赋值 | [W_RIS/wris_config.m:19](<C:/Users/Spidey/Desktop/code/W_RIS/wris_config.m:19>) |
| `seaDefaults` | 由 get_default_sea_params 读取的基础海面配置 | `get_default_sea_params()` | 赋值 | [W_RIS/wris_config.m:20](<C:/Users/Spidey/Desktop/code/W_RIS/wris_config.m:20>) |
| `c.sea.neumannConstant` | 源码 Neumann 频谱经验常数，默认 3.05；维度随源码谱表达式确定 | `seaDefaults.neumannConstant` | 赋值 | [W_RIS/wris_config.m:21](<C:/Users/Spidey/Desktop/code/W_RIS/wris_config.m:21>) |
| `c.tx.seed` | 随机数种子；所属结构决定是海浪、TX、MC 或湍流随机流 | `111` | 赋值 | [W_RIS/wris_config.m:22](<C:/Users/Spidey/Desktop/code/W_RIS/wris_config.m:22>) |
| `c.tx.initialYaw` | TX 平台初始航向角（rad） | `pi/2` | 赋值 | [W_RIS/wris_config.m:22](<C:/Users/Spidey/Desktop/code/W_RIS/wris_config.m:22>) |
| `c.tx.meanUVR` | TX 体坐标 [u;v;r] 的 OU 平稳均值（m/s,m/s,rad/s） | `[0.008;0;deg2rad(0.03)]` | 赋值 | [W_RIS/wris_config.m:23](<C:/Users/Spidey/Desktop/code/W_RIS/wris_config.m:23>) |
| `c.tx.stdUVR` | TX 体坐标 [u;v;r] 的 OU 平稳标准差（m/s,m/s,rad/s） | `[0.012;0.008;deg2rad(0.12)]` | 赋值 | [W_RIS/wris_config.m:24](<C:/Users/Spidey/Desktop/code/W_RIS/wris_config.m:24>) |
| `c.tx.correlationTime` | TX 三个 OU 分量的相关时间（s） | `[2;2;2]` | 赋值 | [W_RIS/wris_config.m:25](<C:/Users/Spidey/Desktop/code/W_RIS/wris_config.m:25>) |
| `c.tx.pointing` | 指向方式，当前仅支持 body-fixed（光轴随平台偏航） | `'body-fixed'` | 赋值 | [W_RIS/wris_config.m:26](<C:/Users/Spidey/Desktop/code/W_RIS/wris_config.m:26>) |
| `c.source.power` | 发射光功率（W） | `1` | 赋值 | [W_RIS/wris_config.m:27](<C:/Users/Spidey/Desktop/code/W_RIS/wris_config.m:27>) |
| `c.source.wavelength` | 真空/模型使用的光波长（m），默认 532e-9 | `532e-9` | 赋值 | [W_RIS/wris_config.m:27](<C:/Users/Spidey/Desktop/code/W_RIS/wris_config.m:27>) |
| `c.source.divergence` | 源码高斯径向发射律的角尺度（rad） | `0.01` | 赋值 | [W_RIS/wris_config.m:28](<C:/Users/Spidey/Desktop/code/W_RIS/wris_config.m:28>) |
| `c.source.waist` | 光源发射面高斯束腰尺度 wavelength/(pi*divergence)（m） | `c.source.wavelength/(pi*c.source.divergence)` | 赋值 | [W_RIS/wris_config.m:29](<C:/Users/Spidey/Desktop/code/W_RIS/wris_config.m:29>) |
| `c.water.type` | Table_abc 水体档位字符串 | `'c0.151'` | 赋值 | [W_RIS/wris_config.m:30](<C:/Users/Spidey/Desktop/code/W_RIS/wris_config.m:30>) |
| `c.water.absorption` | 水体吸收系数（m⁻¹），Table_abc 返回 a=c-b | `Table_abc(c.water.type)` | 赋值 | [W_RIS/wris_config.m:31](<C:/Users/Spidey/Desktop/code/W_RIS/wris_config.m:31>) |
| `c.water.scattering` | 水体散射系数（m⁻¹），Table_abc 返回 b | `Table_abc(c.water.type)` | 赋值 | [W_RIS/wris_config.m:31](<C:/Users/Spidey/Desktop/code/W_RIS/wris_config.m:31>) |
| `c.water.attenuation` | 水体总衰减系数（m⁻¹），Table_abc 返回 c=a+b；MC 内重新用 a+b | `Table_abc(c.water.type)` | 赋值 | [W_RIS/wris_config.m:31](<C:/Users/Spidey/Desktop/code/W_RIS/wris_config.m:31>) |
| `c.water.anisotropy` | HG 散射相函数各向异性因子，-1<值<1 | `0.924` | 赋值 | [W_RIS/wris_config.m:32](<C:/Users/Spidey/Desktop/code/W_RIS/wris_config.m:32>) |
| `c.turbulence.kind` | EGG 档位 weak/strong/none；缓存函数中为 motion/optics 类型 | `'weak'` | 赋值 | [W_RIS/wris_config.m:33](<C:/Users/Spidey/Desktop/code/W_RIS/wris_config.m:33>) |
| `c.turbulence.seed` | 随机数种子；所属结构决定是海浪、TX、MC 或湍流随机流 | `260617` | 赋值 | [W_RIS/wris_config.m:33](<C:/Users/Spidey/Desktop/code/W_RIS/wris_config.m:33>) |
| `c.mc.photons` | 单帧发射光子包总数，或验证采用的光子数档位 | `1e6`<br>`10000` | 赋值 | [W_RIS/wris_config.m:34](<C:/Users/Spidey/Desktop/code/W_RIS/wris_config.m:34>) |
| `c.mc.chunkSize` | MC 每批最大光子包数 | `100000` | 赋值 | [W_RIS/wris_config.m:34](<C:/Users/Spidey/Desktop/code/W_RIS/wris_config.m:34>) |
| `c.mc.maxScatter` | 每个光子包允许的最大散射次数 | `1000` | 赋值 | [W_RIS/wris_config.m:35](<C:/Users/Spidey/Desktop/code/W_RIS/wris_config.m:35>) |
| `c.mc.maxPath` | 每个光子包允许的最大累计路径长度（m） | `200` | 赋值 | [W_RIS/wris_config.m:35](<C:/Users/Spidey/Desktop/code/W_RIS/wris_config.m:35>) |
| `c.mc.rouletteThreshold` | 俄罗斯轮盘赌启动权重门限（无量纲） | `1e-4` | 赋值 | [W_RIS/wris_config.m:36](<C:/Users/Spidey/Desktop/code/W_RIS/wris_config.m:36>) |
| `c.mc.rouletteSurvival` | 低权重包轮盘赌存活概率（无量纲） | `0.1` | 赋值 | [W_RIS/wris_config.m:36](<C:/Users/Spidey/Desktop/code/W_RIS/wris_config.m:36>) |
| `c.mc.seed` | 随机数种子；所属结构决定是海浪、TX、MC 或湍流随机流 | `20261014` | 赋值 | [W_RIS/wris_config.m:37](<C:/Users/Spidey/Desktop/code/W_RIS/wris_config.m:37>) |
| `c.grid.bins` | 板面 x/y 方向分箱数 [nx,ny] | `[80 80]` | 赋值 | [W_RIS/wris_config.m:38](<C:/Users/Spidey/Desktop/code/W_RIS/wris_config.m:38>) |
| `c.grid.planeScale` | 延展评分窗口相对板面的边长倍数；保持单元面积不变 | `6` | 赋值 | [W_RIS/wris_config.m:38](<C:/Users/Spidey/Desktop/code/W_RIS/wris_config.m:38>) |
| `c.validation.photons` | 单帧发射光子包总数，或验证采用的光子数档位 | `[10000 100000 1000000]`<br>`[2000 10000 30000]` | 赋值 | [W_RIS/wris_config.m:39](<C:/Users/Spidey/Desktop/code/W_RIS/wris_config.m:39>) |
| `c.validation.repeats` | 每个光子数档位独立重复次数 | `4`<br>`3` | 赋值 | [W_RIS/wris_config.m:40](<C:/Users/Spidey/Desktop/code/W_RIS/wris_config.m:40>) |
| `c.validation.seed` | 随机数种子；所属结构决定是海浪、TX、MC 或湍流随机流 | `9500` | 赋值 | [W_RIS/wris_config.m:40](<C:/Users/Spidey/Desktop/code/W_RIS/wris_config.m:40>) |
| `c.output` | 输出目录 | `fullfile(c.root,'results',profile)` | 赋值 | [W_RIS/wris_config.m:41](<C:/Users/Spidey/Desktop/code/W_RIS/wris_config.m:41>) |
| `c.cache.folder` | 缓存目录 | `fullfile(c.root,'cache')` | 赋值 | [W_RIS/wris_config.m:42](<C:/Users/Spidey/Desktop/code/W_RIS/wris_config.m:42>) |
| `c.cache.enabled` | 是否复用物理计算缓存 | `true` | 赋值 | [W_RIS/wris_config.m:42](<C:/Users/Spidey/Desktop/code/W_RIS/wris_config.m:42>) |
| `c.plot.resolution` | 导出 PNG 分辨率（dpi） | `300` | 赋值 | [W_RIS/wris_config.m:43](<C:/Users/Spidey/Desktop/code/W_RIS/wris_config.m:43>) |
| `c.plot.dynamicRange` | 共享对数色标的最大/最小辐照度比 | `1e4` | 赋值 | [W_RIS/wris_config.m:43](<C:/Users/Spidey/Desktop/code/W_RIS/wris_config.m:43>) |
| `c.plot.padding` | 板面边界外绘图区留白相对板边长的比例 | `0.15` | 赋值 | [W_RIS/wris_config.m:43](<C:/Users/Spidey/Desktop/code/W_RIS/wris_config.m:43>) |


### 5.19 W_RIS/wris_egg.m

**作用域：`wris_egg`**

| 源码符号 | 含义与单位/维度 | 源码取值/计算式（各次赋值） | 角色 | 源码位置 |
| --- | --- | --- | --- | --- |
| `kind` | EGG 混合分布档位 weak（默认）、strong 或 none | `'weak' n (1,1) {mustBeNumeric,mustBeReal,mustBeFinite,mustBeNonnegative,mustBeInteger} = 1`<br>`validatestring(kind,{'weak','strong','none'},mfilename,'kind');`<br>`'weak'`<br>`validatestring(kind,{'weak','strong','none'},mfilename,'kind')` | 输入/声明默认值/赋值 | [W_RIS/wris_egg.m:1](<C:/Users/Spidey/Desktop/code/W_RIS/wris_egg.m:1>) |
| `n` | 信道快照数量，非每帧光子数量；默认 1，可为 0 | 函数输入/输出；值随调用或上游结果确定 | 输入 | [W_RIS/wris_egg.m:1](<C:/Users/Spidey/Desktop/code/W_RIS/wris_egg.m:1>) |
| `seed` | 可选随机种子，默认 []，为空时使用调用方随机流 | `[]` | 输入/声明默认值 | [W_RIS/wris_egg.m:1](<C:/Users/Spidey/Desktop/code/W_RIS/wris_egg.m:1>) |
| `h` | 均值归一化后的信道强度增益列向量（每快照一个，无量纲） | `ones(n,1)`<br>`x/mu` | 输出/赋值 | [W_RIS/wris_egg.m:1](<C:/Users/Spidey/Desktop/code/W_RIS/wris_egg.m:1>) |
| `info` | 分布参数、矩、闪烁指数及表格元数据结构 | 函数输入/输出；值随调用或上游结果确定 | 输出 | [W_RIS/wris_egg.m:1](<C:/Users/Spidey/Desktop/code/W_RIS/wris_egg.m:1>) |
| `oldRng` | 设置局部随机种子前保存的调用方随机数状态 | `rng` | 赋值 | [W_RIS/wris_egg.m:18](<C:/Users/Spidey/Desktop/code/W_RIS/wris_egg.m:18>) |
| `cleanup` | 函数结束时自动恢复随机数状态、关闭图窗或文件的清理对象 | `onCleanup(@() rng(oldRng))` | 赋值 | [W_RIS/wris_egg.m:19](<C:/Users/Spidey/Desktop/code/W_RIS/wris_egg.m:19>) |
| `info.kind` | EGG 档位 weak/strong/none；缓存函数中为 motion/optics 类型 | `kind` | 赋值 | [W_RIS/wris_egg.m:22](<C:/Users/Spidey/Desktop/code/W_RIS/wris_egg.m:22>) |
| `info.sampleCount` | EGG 信道快照样本数或运动记录样本数 | `n` | 赋值 | [W_RIS/wris_egg.m:23](<C:/Users/Spidey/Desktop/code/W_RIS/wris_egg.m:23>) |
| `info.seed` | 随机数种子；所属结构决定是海浪、TX、MC 或湍流随机流 | `seed` | 赋值 | [W_RIS/wris_egg.m:24](<C:/Users/Spidey/Desktop/code/W_RIS/wris_egg.m:24>) |
| `info.application` | EGG 增益施加方式说明文本 | `'One gain per entire channel snapshot, shared by all photons.'` | 赋值 | [W_RIS/wris_egg.m:25](<C:/Users/Spidey/Desktop/code/W_RIS/wris_egg.m:25>) |
| `info.momentOrders` | 统计矩阶数 0:4 | `0:4` | 赋值 | [W_RIS/wris_egg.m:26](<C:/Users/Spidey/Desktop/code/W_RIS/wris_egg.m:26>) |
| `info.parameters` | 当前计算实际使用的参数结构；EGG 中为混合分布参数 | `struct('omega',0,'lambda_e',NaN,'a_g',NaN,'b_g',NaN,'c_g',NaN)`<br>`struct('omega',w,'lambda_e',lambda,'a_g',a,'b_g',b,'c_g',c)` | 赋值 | [W_RIS/wris_egg.m:29](<C:/Users/Spidey/Desktop/code/W_RIS/wris_egg.m:29>) |
| `info.parameters.omega` | EGG 指数分量混合概率（无量纲），不是海浪角频率 | `0`<br>`w` | 结构字段 | [W_RIS/wris_egg.m:29](<C:/Users/Spidey/Desktop/code/W_RIS/wris_egg.m:29>) |
| `info.parameters.lambda_e` | EGG 指数分布尺度参数 | `NaN`<br>`lambda` | 结构字段 | [W_RIS/wris_egg.m:29](<C:/Users/Spidey/Desktop/code/W_RIS/wris_egg.m:29>) |
| `info.parameters.a_g` | EGG 广义 Gamma 的形状 a、尺度 b、幂指数 c，按字段分别对应 | `NaN`<br>`a` | 结构字段 | [W_RIS/wris_egg.m:29](<C:/Users/Spidey/Desktop/code/W_RIS/wris_egg.m:29>) |
| `info.parameters.b_g` | EGG 广义 Gamma 的形状 a、尺度 b、幂指数 c，按字段分别对应 | `NaN`<br>`b` | 结构字段 | [W_RIS/wris_egg.m:29](<C:/Users/Spidey/Desktop/code/W_RIS/wris_egg.m:29>) |
| `info.parameters.c_g` | EGG 广义 Gamma 的形状 a、尺度 b、幂指数 c，按字段分别对应 | `NaN`<br>`c` | 结构字段 | [W_RIS/wris_egg.m:29](<C:/Users/Spidey/Desktop/code/W_RIS/wris_egg.m:29>) |
| `info.rawMean` | 原始 EGG 混合变量 X 的均值 E[X] | `1`<br>`mu` | 赋值 | [W_RIS/wris_egg.m:30](<C:/Users/Spidey/Desktop/code/W_RIS/wris_egg.m:30>) |
| `info.rawMoments` | 原始混合变量的 0..4 阶原点矩 | `ones(1,5)`<br>`rawMoments` | 赋值 | [W_RIS/wris_egg.m:31](<C:/Users/Spidey/Desktop/code/W_RIS/wris_egg.m:31>) |
| `info.normalizedMoments` | 均值归一化增益 H 的 0..4 阶原点矩 | `ones(1,5)`<br>`rawMoments./mu.^k` | 赋值 | [W_RIS/wris_egg.m:32](<C:/Users/Spidey/Desktop/code/W_RIS/wris_egg.m:32>) |
| `info.theoreticalScintillationIndex` | 由源码分布参数计算的闪烁指数 E[X²]/E[X]²-1 | `0`<br>`rawMoments(3)/mu^2-1` | 赋值 | [W_RIS/wris_egg.m:33](<C:/Users/Spidey/Desktop/code/W_RIS/wris_egg.m:33>) |
| `info.tabulatedScintillationIndex` | 所给参数表的闪烁指数；不强行替代理论值 | `0`<br>`tabSI` | 赋值 | [W_RIS/wris_egg.m:34](<C:/Users/Spidey/Desktop/code/W_RIS/wris_egg.m:34>) |
| `info.bubbleFlow_Lmin` | EGG 表中对应的气泡流量（L/min），仅元数据 | `0`<br>`flow` | 赋值 | [W_RIS/wris_egg.m:35](<C:/Users/Spidey/Desktop/code/W_RIS/wris_egg.m:35>) |
| `info.temperatureGradient_CperCm` | EGG 表中温度梯度（°C/cm），仅元数据 | `0`<br>`gradient` | 赋值 | [W_RIS/wris_egg.m:36](<C:/Users/Spidey/Desktop/code/W_RIS/wris_egg.m:36>) |
| `p` | 五个 EGG 参数向量 [omega,lambda_e,a_g,b_g,c_g]，此处不是 RIS 配置 | `[0.2130 0.3291 1.4299 1.1817 17.1984]`<br>`[0.4539 0.2744 0.3008 1.7053 54.1422]` | 赋值 | [W_RIS/wris_egg.m:40](<C:/Users/Spidey/Desktop/code/W_RIS/wris_egg.m:40>) |
| `tabSI` | 原表提供的闪烁指数，weak=0.1484、strong=0.4769 | `0.1484`<br>`0.4769` | 赋值 | [W_RIS/wris_egg.m:41](<C:/Users/Spidey/Desktop/code/W_RIS/wris_egg.m:41>) |
| `flow` | 原表气泡流量，weak=2.4、strong=4.7（L/min） | `2.4`<br>`4.7` | 赋值 | [W_RIS/wris_egg.m:41](<C:/Users/Spidey/Desktop/code/W_RIS/wris_egg.m:41>) |
| `gradient` | 原表温度梯度，weak=0.05、strong=0.1（°C/cm） | `0.05`<br>`0.1` | 赋值 | [W_RIS/wris_egg.m:41](<C:/Users/Spidey/Desktop/code/W_RIS/wris_egg.m:41>) |
| `w` | 指数分布混合概率 omega（无量纲），此处不是海浪结构 | `p(1)` | 赋值 | [W_RIS/wris_egg.m:46](<C:/Users/Spidey/Desktop/code/W_RIS/wris_egg.m:46>) |
| `lambda` | 指数分布尺度 lambda_e | `p(2)` | 赋值 | [W_RIS/wris_egg.m:46](<C:/Users/Spidey/Desktop/code/W_RIS/wris_egg.m:46>) |
| `a` | 广义 Gamma 形状参数 a_g | `p(3)` | 赋值 | [W_RIS/wris_egg.m:46](<C:/Users/Spidey/Desktop/code/W_RIS/wris_egg.m:46>) |
| `b` | 广义 Gamma 尺度参数 b_g | `p(4)` | 赋值 | [W_RIS/wris_egg.m:46](<C:/Users/Spidey/Desktop/code/W_RIS/wris_egg.m:46>) |
| `c` | 广义 Gamma 幂指数 c_g；此处不是 W_RIS 配置 | `p(5)` | 赋值 | [W_RIS/wris_egg.m:46](<C:/Users/Spidey/Desktop/code/W_RIS/wris_egg.m:46>) |
| `k` | 所求原点矩阶数 0:4 | `info.momentOrders` | 赋值 | [W_RIS/wris_egg.m:47](<C:/Users/Spidey/Desktop/code/W_RIS/wris_egg.m:47>) |
| `rawMoments` | 混合变量 X 的 0..4 阶解析原点矩 | `w*exp(k*log(lambda)+gammaln(k+1)) + (1-w)*exp(k*log(b)+gammaln(a+k/c)-gammaln(a))` | 赋值 | [W_RIS/wris_egg.m:48](<C:/Users/Spidey/Desktop/code/W_RIS/wris_egg.m:48>) |
| `mu` | 原始混合变量 X 的解析均值 rawMoments(2) | `rawMoments(2)` | 赋值 | [W_RIS/wris_egg.m:50](<C:/Users/Spidey/Desktop/code/W_RIS/wris_egg.m:50>) |
| `isExp` | 各快照是否由指数分量抽样的伯努利掩码 | `rand(n,1)<w` | 赋值 | [W_RIS/wris_egg.m:51](<C:/Users/Spidey/Desktop/code/W_RIS/wris_egg.m:51>) |
| `x` | 尚未除以均值的原始混合随机增益 X | `zeros(n,1)`<br>`-lambda*log(max(rand(ne,1),realmin))`<br>`b*randg(a,ng,1).^(1/c)` | 赋值 | [W_RIS/wris_egg.m:52](<C:/Users/Spidey/Desktop/code/W_RIS/wris_egg.m:52>) |
| `ne` | 分别由指数/广义 Gamma 分量抽样的样本数 | `nnz(isExp)` | 赋值 | [W_RIS/wris_egg.m:53](<C:/Users/Spidey/Desktop/code/W_RIS/wris_egg.m:53>) |
| `ng` | 分别由指数/广义 Gamma 分量抽样的样本数 | `n-ne` | 赋值 | [W_RIS/wris_egg.m:55](<C:/Users/Spidey/Desktop/code/W_RIS/wris_egg.m:55>) |
| `info.source` | EGG 参数来源说明文本 | `'Supplied image-1 Eqs. (2-10),(2-11); image-2 Eqs. (2-12),(2-13), Table 2-1.'` | 赋值 | [W_RIS/wris_egg.m:68](<C:/Users/Spidey/Desktop/code/W_RIS/wris_egg.m:68>) |


**作用域：`mustBeSeed`**

| 源码符号 | 含义与单位/维度 | 源码取值/计算式（各次赋值） | 角色 | 源码位置 |
| --- | --- | --- | --- | --- |
| `seed` | 可选随机种子验证器输入；允许 [] 或 0..2^32-1 的有限实整数标量 | 函数输入/输出；值随调用或上游结果确定 | 输入 | [W_RIS/wris_egg.m:70](<C:/Users/Spidey/Desktop/code/W_RIS/wris_egg.m:70>) |


### 5.20 W_RIS/wris_geometry.m

**作用域：`wris_geometry`**

| 源码符号 | 含义与单位/维度 | 源码取值/计算式（各次赋值） | 角色 | 源码位置 |
| --- | --- | --- | --- | --- |
| `c` | 当前 W_RIS 总配置结构 | 函数输入/输出；值随调用或上游结果确定 | 输入 | [W_RIS/wris_geometry.m:1](<C:/Users/Spidey/Desktop/code/W_RIS/wris_geometry.m:1>) |
| `m` | RIS 动力学结果结构（static 模式仅需 t） | 函数输入/输出；值随调用或上游结果确定 | 输入 | [W_RIS/wris_geometry.m:1](<C:/Users/Spidey/Desktop/code/W_RIS/wris_geometry.m:1>) |
| `tx` | 发射平台运动时序结构 | 函数输入/输出；值随调用或上游结果确定 | 输入 | [W_RIS/wris_geometry.m:1](<C:/Users/Spidey/Desktop/code/W_RIS/wris_geometry.m:1>) |
| `k` | 冻结光学快照所用的运动样本索引 | 函数输入/输出；值随调用或上游结果确定 | 输入 | [W_RIS/wris_geometry.m:1](<C:/Users/Spidey/Desktop/code/W_RIS/wris_geometry.m:1>) |
| `mode` | 快照运动组合；默认 joint，支持 ris-only、tx-only、static | `'joint'` | 输入/赋值 | [W_RIS/wris_geometry.m:1](<C:/Users/Spidey/Desktop/code/W_RIS/wris_geometry.m:1>) |
| `g` | 当前光学几何结构 | `wris_reference(c)` | 输出/赋值 | [W_RIS/wris_geometry.m:1](<C:/Users/Spidey/Desktop/code/W_RIS/wris_geometry.m:1>) |
| `g.RRis` | RIS 体坐标至世界坐标旋转矩阵 | `m.rotationMatrix(:,:,k)` | 赋值 | [W_RIS/wris_geometry.m:6](<C:/Users/Spidey/Desktop/code/W_RIS/wris_geometry.m:6>) |
| `g.pRis` | RIS 下表面光学中心世界坐标（m） | `m.centerOfMass(:,k)-c.ris.thickness/2*g.RRis(:,3)` | 赋值 | [W_RIS/wris_geometry.m:8](<C:/Users/Spidey/Desktop/code/W_RIS/wris_geometry.m:8>) |
| `g.pTx` | TX 发射面中心世界坐标（m） | `tx.position(:,k)` | 赋值 | [W_RIS/wris_geometry.m:11](<C:/Users/Spidey/Desktop/code/W_RIS/wris_geometry.m:11>) |
| `g.RTx` | 以光轴为局部 z 的 TX 光学系至世界系旋转矩阵 | `tx.rotation(:,:,k)` | 赋值 | [W_RIS/wris_geometry.m:11](<C:/Users/Spidey/Desktop/code/W_RIS/wris_geometry.m:11>) |
| `g.dTx` | TX 光轴在世界系的单位方向 | `tx.direction(:,k)` | 赋值 | [W_RIS/wris_geometry.m:11](<C:/Users/Spidey/Desktop/code/W_RIS/wris_geometry.m:11>) |
| `g.normal` | 板体 +z 轴的世界单位法向量（3×Nt） | `g.RRis(:,3)` | 赋值 | [W_RIS/wris_geometry.m:13](<C:/Users/Spidey/Desktop/code/W_RIS/wris_geometry.m:13>) |
| `g.distance` | TX 至 RIS 光学中心的空间直线距离（m） | `norm(g.pRis-g.pTx)` | 赋值 | [W_RIS/wris_geometry.m:13](<C:/Users/Spidey/Desktop/code/W_RIS/wris_geometry.m:13>) |
| `g.thetaDeg` | 当前光轴与板法线的锐角入射角 acosd(abs(dot(dTx,normal)))（°） | `acosd(min(1,abs(dot(g.dTx,g.normal))))` | 赋值 | [W_RIS/wris_geometry.m:14](<C:/Users/Spidey/Desktop/code/W_RIS/wris_geometry.m:14>) |
| `g.frontSide` | TX 是否位于接收面下方的几何标志 dot(pRis-pTx,normal)>0 | `dot(g.pRis-g.pTx,g.normal)>0` | 赋值 | [W_RIS/wris_geometry.m:15](<C:/Users/Spidey/Desktop/code/W_RIS/wris_geometry.m:15>) |
| `g.axisLocal` | 光轴与 RIS 无限平面交点的板局部坐标 [x,y,z]（m） | `wris_intersect(g.pTx.',g.dTx.',g.pRis,g.RRis,[c.ris.Lx c.ris.Ly])` | 赋值 | [W_RIS/wris_geometry.m:16](<C:/Users/Spidey/Desktop/code/W_RIS/wris_geometry.m:16>) |
| `g.axisRange` | 光轴从 TX 到平面交点的前向参数距离（m）；无效为 Inf | `wris_intersect(g.pTx.',g.dTx.',g.pRis,g.RRis,[c.ris.Lx c.ris.Ly])` | 赋值 | [W_RIS/wris_geometry.m:16](<C:/Users/Spidey/Desktop/code/W_RIS/wris_geometry.m:16>) |
| `g.axisInside` | 光轴交点是否在有限矩形板孔径内 | `wris_intersect(g.pTx.',g.dTx.',g.pRis,g.RRis,[c.ris.Lx c.ris.Ly])` | 赋值 | [W_RIS/wris_geometry.m:16](<C:/Users/Spidey/Desktop/code/W_RIS/wris_geometry.m:16>) |
| `g.time` | 几何快照的相对时间（s） | `m.t(k)` | 赋值 | [W_RIS/wris_geometry.m:17](<C:/Users/Spidey/Desktop/code/W_RIS/wris_geometry.m:17>) |
| `g.mode` | 运动组合：static/joint/ris-only/tx-only | `mode` | 赋值 | [W_RIS/wris_geometry.m:17](<C:/Users/Spidey/Desktop/code/W_RIS/wris_geometry.m:17>) |
| `localDirection` | TX 光轴变换到当前 RIS 局部体轴的方向分量 | `g.RRis.'*g.dTx` | 赋值 | [W_RIS/wris_geometry.m:20](<C:/Users/Spidey/Desktop/code/W_RIS/wris_geometry.m:20>) |
| `g.directionLocal` | TX 光轴在当前 RIS 体坐标中的分量 | `localDirection` | 赋值 | [W_RIS/wris_geometry.m:21](<C:/Users/Spidey/Desktop/code/W_RIS/wris_geometry.m:21>) |
| `g.yawRisDeg` | 由 RRis 提取的实际 RIS yaw（°） | `atan2d(g.RRis(2,1),g.RRis(1,1))` | 赋值 | [W_RIS/wris_geometry.m:22](<C:/Users/Spidey/Desktop/code/W_RIS/wris_geometry.m:22>) |
| `g.projectedOrientationDeg` | 光轴板面切向投影对应的几何主轴角（°，[-90,90)，法向入射时 NaN） | `NaN`<br>`mod(atan2d(localDirection(2),localDirection(1))+90,180)-90` | 赋值 | [W_RIS/wris_geometry.m:23](<C:/Users/Spidey/Desktop/code/W_RIS/wris_geometry.m:23>) |
| `noYaw` | 移除 RIS 世界 yaw 后的旋转矩阵，仅供对照诊断 | `RIS_rotation([0;0;-deg2rad(g.yawRisDeg)])*g.RRis` | 赋值 | [W_RIS/wris_geometry.m:28](<C:/Users/Spidey/Desktop/code/W_RIS/wris_geometry.m:28>) |
| `directionNoYaw` | 去除 RIS yaw 后的 TX 光轴局部方向分量 | `noYaw.'*g.dTx` | 赋值 | [W_RIS/wris_geometry.m:29](<C:/Users/Spidey/Desktop/code/W_RIS/wris_geometry.m:29>) |
| `g.orientationWithoutRisYawDeg` | 仅消去 RIS 世界 yaw、保留其余姿态和 TX 后的几何主轴角（°） | `NaN`<br>`mod(atan2d(directionNoYaw(2),directionNoYaw(1))+90,180)-90` | 赋值 | [W_RIS/wris_geometry.m:30](<C:/Users/Spidey/Desktop/code/W_RIS/wris_geometry.m:30>) |


### 5.21 W_RIS/wris_intersect.m

**作用域：`wris_intersect`**

| 源码符号 | 含义与单位/维度 | 源码取值/计算式（各次赋值） | 角色 | 源码位置 |
| --- | --- | --- | --- | --- |
| `p` | 每条射线世界起点（N×3，m） | 函数输入/输出；值随调用或上游结果确定 | 输入 | [W_RIS/wris_intersect.m:1](<C:/Users/Spidey/Desktop/code/W_RIS/wris_intersect.m:1>) |
| `d` | 每条射线世界单位方向（N×3） | 函数输入/输出；值随调用或上游结果确定 | 输入 | [W_RIS/wris_intersect.m:1](<C:/Users/Spidey/Desktop/code/W_RIS/wris_intersect.m:1>) |
| `center` | 接收无限平面的世界中心（3×1，m） | 函数输入/输出；值随调用或上游结果确定 | 输入 | [W_RIS/wris_intersect.m:1](<C:/Users/Spidey/Desktop/code/W_RIS/wris_intersect.m:1>) |
| `R` | 体坐标到世界坐标的主动旋转矩阵（无量纲） | 函数输入/输出；值随调用或上游结果确定 | 输入 | [W_RIS/wris_intersect.m:1](<C:/Users/Spidey/Desktop/code/W_RIS/wris_intersect.m:1>) |
| `boardSize` | 有限矩形接收板尺寸 [Lx,Ly]（m） | 函数输入/输出；值随调用或上游结果确定 | 输入 | [W_RIS/wris_intersect.m:1](<C:/Users/Spidey/Desktop/code/W_RIS/wris_intersect.m:1>) |
| `local` | 各射线与无限平面交点的板局部坐标（N×3，m）；无效为 NaN | `(p+d.*range-center.')*R`<br>`NaN` | 输出/赋值 | [W_RIS/wris_intersect.m:1](<C:/Users/Spidey/Desktop/code/W_RIS/wris_intersect.m:1>) |
| `range` | 前向交点距离参数（m）；平行或反向为 Inf | `((center.'-p)*normal)./den`<br>`Inf` | 输出/赋值 | [W_RIS/wris_intersect.m:1](<C:/Users/Spidey/Desktop/code/W_RIS/wris_intersect.m:1>) |
| `inside` | 有效前向交点位于有限矩形孔径内的布尔掩码 | `valid & abs(local(:,1))<=boardSize(1)/2 & abs(local(:,2))<=boardSize(2)/2` | 输出/赋值 | [W_RIS/wris_intersect.m:1](<C:/Users/Spidey/Desktop/code/W_RIS/wris_intersect.m:1>) |
| `normal` | 当前板面或平面单位法向量 | `R(:,3)` | 赋值 | [W_RIS/wris_intersect.m:4](<C:/Users/Spidey/Desktop/code/W_RIS/wris_intersect.m:4>) |
| `den` | 射线单位方向与平面法向点积 | `d*normal` | 赋值 | [W_RIS/wris_intersect.m:5](<C:/Users/Spidey/Desktop/code/W_RIS/wris_intersect.m:5>) |
| `valid` | abs(den)>1e-12 且 range>=0 的交点有效掩码 | `abs(den)>1e-12 & range>=0` | 赋值 | [W_RIS/wris_intersect.m:7](<C:/Users/Spidey/Desktop/code/W_RIS/wris_intersect.m:7>) |


### 5.22 W_RIS/wris_mc.m

**作用域：`wris_mc`**

| 源码符号 | 含义与单位/维度 | 源码取值/计算式（各次赋值） | 角色 | 源码位置 |
| --- | --- | --- | --- | --- |
| `c` | 当前 W_RIS 总配置结构 | 函数输入/输出；值随调用或上游结果确定 | 输入 | [W_RIS/wris_mc.m:1](<C:/Users/Spidey/Desktop/code/W_RIS/wris_mc.m:1>) |
| `g` | 当前光学几何结构 | 函数输入/输出；值随调用或上游结果确定 | 输入 | [W_RIS/wris_mc.m:1](<C:/Users/Spidey/Desktop/code/W_RIS/wris_mc.m:1>) |
| `seed` | 随机数种子（非负整数） | 函数输入/输出；值随调用或上游结果确定 | 输入 | [W_RIS/wris_mc.m:1](<C:/Users/Spidey/Desktop/code/W_RIS/wris_mc.m:1>) |
| `gain` | 单个信道快照共享的均值为 1 的 EGG 强度增益（无量纲） | `1` | 输入/赋值 | [W_RIS/wris_mc.m:1](<C:/Users/Spidey/Desktop/code/W_RIS/wris_mc.m:1>) |
| `out` | 当前函数输出结果结构 | 函数输入/输出；值随调用或上游结果确定 | 输出 | [W_RIS/wris_mc.m:1](<C:/Users/Spidey/Desktop/code/W_RIS/wris_mc.m:1>) |
| `old` | 设置局部随机种子前保存的调用方随机数状态 | `rng` | 赋值 | [W_RIS/wris_mc.m:11](<C:/Users/Spidey/Desktop/code/W_RIS/wris_mc.m:11>) |
| `guard` | 函数结束时自动恢复随机数状态、关闭图窗或文件的清理对象 | `onCleanup(@() rng(old))` | 赋值 | [W_RIS/wris_mc.m:11](<C:/Users/Spidey/Desktop/code/W_RIS/wris_mc.m:11>) |
| `N` | 本帧发射光子包总数 c.mc.photons | `c.mc.photons` | 赋值 | [W_RIS/wris_mc.m:13](<C:/Users/Spidey/Desktop/code/W_RIS/wris_mc.m:13>) |
| `a` | 水体吸收系数（m⁻¹），不是 EGG a_g | `c.water.absorption` | 赋值 | [W_RIS/wris_mc.m:13](<C:/Users/Spidey/Desktop/code/W_RIS/wris_mc.m:13>) |
| `b` | 水体散射系数（m⁻¹），不是 EGG b_g | `c.water.scattering` | 赋值 | [W_RIS/wris_mc.m:13](<C:/Users/Spidey/Desktop/code/W_RIS/wris_mc.m:13>) |
| `ext` | 水体总衰减系数 a+b（m⁻¹） | `a+b` | 赋值 | [W_RIS/wris_mc.m:13](<C:/Users/Spidey/Desktop/code/W_RIS/wris_mc.m:13>) |
| `nx` | 有限接收板 x/y 分箱数 | `c.grid.bins(1)` | 赋值 | [W_RIS/wris_mc.m:18](<C:/Users/Spidey/Desktop/code/W_RIS/wris_mc.m:18>) |
| `ny` | 有限接收板 x/y 分箱数 | `c.grid.bins(2)` | 赋值 | [W_RIS/wris_mc.m:18](<C:/Users/Spidey/Desktop/code/W_RIS/wris_mc.m:18>) |
| `xe` | 有限接收板 x/y 分箱边界（m） | `linspace(-c.ris.Lx/2,c.ris.Lx/2,nx+1)` | 赋值 | [W_RIS/wris_mc.m:19](<C:/Users/Spidey/Desktop/code/W_RIS/wris_mc.m:19>) |
| `ye` | 有限接收板 x/y 分箱边界（m） | `linspace(-c.ris.Ly/2,c.ris.Ly/2,ny+1)` | 赋值 | [W_RIS/wris_mc.m:19](<C:/Users/Spidey/Desktop/code/W_RIS/wris_mc.m:19>) |
| `area` | 单个分箱面积（m²），不是整板面积 | `(xe(2)-xe(1))*(ye(2)-ye(1))` | 赋值 | [W_RIS/wris_mc.m:20](<C:/Users/Spidey/Desktop/code/W_RIS/wris_mc.m:20>) |
| `hist` | 板内各格首次到达光子包权重和（无量纲） | `zeros(nx,ny)`<br>`hist+accumarray([ix iy],hitW(valid),[nx ny])` | 赋值 | [W_RIS/wris_mc.m:20](<C:/Users/Spidey/Desktop/code/W_RIS/wris_mc.m:20>) |
| `hist2` | 板内各格首次到达光子包权重平方和（用于标准误） | `hist`<br>`hist2+accumarray([ix iy],hitW(valid).^2,[nx ny])` | 赋值 | [W_RIS/wris_mc.m:20](<C:/Users/Spidey/Desktop/code/W_RIS/wris_mc.m:20>) |
| `planeScale` | 延展评分窗口相对板边长的整数倍数 | `c.grid.planeScale` | 赋值 | [W_RIS/wris_mc.m:23](<C:/Users/Spidey/Desktop/code/W_RIS/wris_mc.m:23>) |
| `pxe` | 延展评分窗口 x/y 分箱边界（m），格距与有限板相同 | `linspace(-planeScale*c.ris.Lx/2,planeScale*c.ris.Lx/2,planeScale*nx+1)` | 赋值 | [W_RIS/wris_mc.m:25](<C:/Users/Spidey/Desktop/code/W_RIS/wris_mc.m:25>) |
| `pye` | 延展评分窗口 x/y 分箱边界（m），格距与有限板相同 | `linspace(-planeScale*c.ris.Ly/2,planeScale*c.ris.Ly/2,planeScale*ny+1)` | 赋值 | [W_RIS/wris_mc.m:26](<C:/Users/Spidey/Desktop/code/W_RIS/wris_mc.m:26>) |
| `planeHist` | 延展评分窗口各格首次到达权重和 | `zeros(planeScale*nx,planeScale*ny)`<br>`planeHist+accumarray([pix(shown) piy(shown)],landedW(shown),size(planeHist))` | 赋值 | [W_RIS/wris_mc.m:27](<C:/Users/Spidey/Desktop/code/W_RIS/wris_mc.m:27>) |
| `ledger` | 终止路径与吸收/轮盘赌的无量纲光子包权重账本 | `struct('inside',0,'outside',0,'backEscape',0,'absorbed',0,'truncated',0,'rouletteDelta',0)` | 赋值 | [W_RIS/wris_mc.m:28](<C:/Users/Spidey/Desktop/code/W_RIS/wris_mc.m:28>) |
| `ledger.inside` | 落在有限板孔径内的权重和或光子计数；父结构区分 | `0`<br>`ledger.inside+sum(wr(inside))` | 结构字段/赋值 | [W_RIS/wris_mc.m:28](<C:/Users/Spidey/Desktop/code/W_RIS/wris_mc.m:28>) |
| `ledger.outside` | 首次到达平面但在有限孔径外的权重和或光子计数 | `0`<br>`ledger.outside+sum(wr(~inside))` | 结构字段/赋值 | [W_RIS/wris_mc.m:28](<C:/Users/Spidey/Desktop/code/W_RIS/wris_mc.m:28>) |
| `ledger.backEscape` | 反向穿出源平面的权重和或光子计数 | `0`<br>`ledger.backEscape+sum(ww(escape))` | 结构字段/赋值 | [W_RIS/wris_mc.m:28](<C:/Users/Spidey/Desktop/code/W_RIS/wris_mc.m:28>) |
| `ledger.absorbed` | 全部碰撞累计被吸收的无量纲光子包权重 | `0`<br>`ledger.absorbed+sum(w(collide))*a/ext` | 结构字段/赋值 | [W_RIS/wris_mc.m:28](<C:/Users/Spidey/Desktop/code/W_RIS/wris_mc.m:28>) |
| `ledger.truncated` | 达到最大路径或最大散射次数而终止的权重和或计数 | `0`<br>`ledger.truncated+sum(w(cut))` | 结构字段/赋值 | [W_RIS/wris_mc.m:28](<C:/Users/Spidey/Desktop/code/W_RIS/wris_mc.m:28>) |
| `ledger.rouletteDelta` | 轮盘赌引入的净权重变化，能量账本单独扣除 | `0`<br>`ledger.rouletteDelta+sum(w(low)-before)` | 结构字段/赋值 | [W_RIS/wris_mc.m:28](<C:/Users/Spidey/Desktop/code/W_RIS/wris_mc.m:28>) |
| `counts` | 各类终止事件计数结构 | `struct('inside',0,'outside',0,'backEscape',0,'rouletteKilled',0,'truncated',0)` | 赋值 | [W_RIS/wris_mc.m:29](<C:/Users/Spidey/Desktop/code/W_RIS/wris_mc.m:29>) |
| `counts.inside` | 落在有限板孔径内的权重和或光子计数；父结构区分 | `0`<br>`counts.inside+nnz(inside)` | 结构字段/赋值 | [W_RIS/wris_mc.m:29](<C:/Users/Spidey/Desktop/code/W_RIS/wris_mc.m:29>) |
| `counts.outside` | 首次到达平面但在有限孔径外的权重和或光子计数 | `0`<br>`counts.outside+nnz(~inside)` | 结构字段/赋值 | [W_RIS/wris_mc.m:29](<C:/Users/Spidey/Desktop/code/W_RIS/wris_mc.m:29>) |
| `counts.backEscape` | 反向穿出源平面的权重和或光子计数 | `0`<br>`counts.backEscape+nnz(escape)` | 结构字段/赋值 | [W_RIS/wris_mc.m:29](<C:/Users/Spidey/Desktop/code/W_RIS/wris_mc.m:29>) |
| `counts.rouletteKilled` | 轮盘赌死亡的光子包总数 | `0`<br>`counts.rouletteKilled+nnz(~survive)` | 结构字段/赋值 | [W_RIS/wris_mc.m:29](<C:/Users/Spidey/Desktop/code/W_RIS/wris_mc.m:29>) |
| `counts.truncated` | 达到最大路径或最大散射次数而终止的权重和或计数 | `0`<br>`counts.truncated+numel(cut)` | 结构字段/赋值 | [W_RIS/wris_mc.m:29](<C:/Users/Spidey/Desktop/code/W_RIS/wris_mc.m:29>) |
| `xyAll` | 按批保存的板内所有首次到达坐标（m）/包权重 | `cell(ceil(N/c.mc.chunkSize),1)`<br>`hitXY(valid,:)` | 赋值 | [W_RIS/wris_mc.m:30](<C:/Users/Spidey/Desktop/code/W_RIS/wris_mc.m:30>) |
| `wAll` | 按批保存的板内所有首次到达坐标（m）/包权重 | `xyAll`<br>`hitW(valid)` | 赋值 | [W_RIS/wris_mc.m:30](<C:/Users/Spidey/Desktop/code/W_RIS/wris_mc.m:30>) |
| `xyDirect` | 按批保存的板内未散射核心坐标（m）/包权重 | `xyAll`<br>`hitXY(core,:)` | 赋值 | [W_RIS/wris_mc.m:31](<C:/Users/Spidey/Desktop/code/W_RIS/wris_mc.m:31>) |
| `wDirect` | 按批保存的板内未散射核心坐标（m）/包权重 | `xyAll`<br>`hitW(core)` | 赋值 | [W_RIS/wris_mc.m:31](<C:/Users/Spidey/Desktop/code/W_RIS/wris_mc.m:31>) |
| `localDirectionSum` | 板内入射方向局部分量的权重加权和 | `zeros(1,3)`<br>`localDirectionSum+sum((dd(reach,:)*g.RRis).* (wr.*inside),1)` | 赋值 | [W_RIS/wris_mc.m:31](<C:/Users/Spidey/Desktop/code/W_RIS/wris_mc.m:31>) |
| `xyCorePlane` | 按批保存的完整无限面未散射核心坐标（m）/包权重 | `xyAll`<br>`hitXY(direct,:)` | 赋值 | [W_RIS/wris_mc.m:32](<C:/Users/Spidey/Desktop/code/W_RIS/wris_mc.m:32>) |
| `wCorePlane` | 按批保存的完整无限面未散射核心坐标（m）/包权重 | `xyAll`<br>`hitW(direct)` | 赋值 | [W_RIS/wris_mc.m:32](<C:/Users/Spidey/Desktop/code/W_RIS/wris_mc.m:32>) |
| `first` | 当前分块/积分子步的起始索引 | `1:c.mc.chunkSize:N` | 循环索引 | [W_RIS/wris_mc.m:33](<C:/Users/Spidey/Desktop/code/W_RIS/wris_mc.m:33>) |
| `chunk` | 当前光子批次编号 | `ceil(first/c.mc.chunkSize)` | 赋值 | [W_RIS/wris_mc.m:34](<C:/Users/Spidey/Desktop/code/W_RIS/wris_mc.m:34>) |
| `n` | 当前批次实际发射光子数 | `min(c.mc.chunkSize,N-first+1)` | 赋值 | [W_RIS/wris_mc.m:34](<C:/Users/Spidey/Desktop/code/W_RIS/wris_mc.m:34>) |
| `uv` | 发射面光学局部横向高斯位置（n×2，m），每轴标准差 waist/2 | `c.source.waist/2*randn(n,2)` | 赋值 | [W_RIS/wris_mc.m:36](<C:/Users/Spidey/Desktop/code/W_RIS/wris_mc.m:36>) |
| `radius` | 发射横向位置距光轴的径向距离（m） | `hypot(uv(:,1),uv(:,2))` | 赋值 | [W_RIS/wris_mc.m:36](<C:/Users/Spidey/Desktop/code/W_RIS/wris_mc.m:36>) |
| `angle` | 当前光子相对光轴的发射极角 radius/waist*divergence（rad） | `radius/c.source.waist*c.source.divergence` | 赋值 | [W_RIS/wris_mc.m:37](<C:/Users/Spidey/Desktop/code/W_RIS/wris_mc.m:37>) |
| `radial` | 发射横向位置的二维径向单位向量 | `uv./max(radius,realmin)` | 赋值 | [W_RIS/wris_mc.m:38](<C:/Users/Spidey/Desktop/code/W_RIS/wris_mc.m:38>) |
| `d` | 每个光子包当前的世界单位传播方向（n×3） | `[sin(angle).*radial,cos(angle)]*g.RTx.'`<br>`newDir./vecnorm(newDir,2,2)` | 赋值 | [W_RIS/wris_mc.m:39](<C:/Users/Spidey/Desktop/code/W_RIS/wris_mc.m:39>) |
| `p` | 每个光子包当前世界位置（n×3，m） | `[uv,zeros(n,1)]*g.RTx.'+g.pTx.'`<br>`pp(move,:)+dd(move,:).*travel(move)` | 赋值 | [W_RIS/wris_mc.m:40](<C:/Users/Spidey/Desktop/code/W_RIS/wris_mc.m:40>) |
| `w` | 每个光子包的当前无量纲存活权重，初始为 1 | `ones(n,1)`<br>`w(collide)*b/ext`<br>`before.*survive/c.mc.rouletteSurvival` | 赋值 | [W_RIS/wris_mc.m:41](<C:/Users/Spidey/Desktop/code/W_RIS/wris_mc.m:41>) |
| `path` | 每个光子包已累计的路径长度（m） | `zeros(n,1)`<br>`path(mid)+travel(move)` | 赋值 | [W_RIS/wris_mc.m:41](<C:/Users/Spidey/Desktop/code/W_RIS/wris_mc.m:41>) |
| `scatter` | 每个光子包已发生散射次数 | `zeros(n,1)`<br>`scatter(collide)+1` | 赋值 | [W_RIS/wris_mc.m:41](<C:/Users/Spidey/Desktop/code/W_RIS/wris_mc.m:41>) |
| `active` | 光子包是否仍在输运中的布尔掩码 | `true(n,1)`<br>`false` | 赋值 | [W_RIS/wris_mc.m:42](<C:/Users/Spidey/Desktop/code/W_RIS/wris_mc.m:42>) |
| `hitXY` | 每个光子包首次到达评分平面的局部二维位置（m） | `zeros(n,2)`<br>`at` | 赋值 | [W_RIS/wris_mc.m:42](<C:/Users/Spidey/Desktop/code/W_RIS/wris_mc.m:42>) |
| `hitW` | 光子包首次到达评分平面时的权重；未到达为 0 | `zeros(n,1)`<br>`wr` | 赋值 | [W_RIS/wris_mc.m:42](<C:/Users/Spidey/Desktop/code/W_RIS/wris_mc.m:42>) |
| `direct` | 首次到达评分面且未发生散射的核心标志 | `false(n,1)`<br>`scatter(rid)==0` | 赋值 | [W_RIS/wris_mc.m:42](<C:/Users/Spidey/Desktop/code/W_RIS/wris_mc.m:42>) |
| `accepted` | 首次到达点落在有限板孔径内的接受标志 | `false(n,1)`<br>`true` | 赋值 | [W_RIS/wris_mc.m:42](<C:/Users/Spidey/Desktop/code/W_RIS/wris_mc.m:42>) |
| `id` | 本轮仍存活的光子包索引 | `find(active)` | 赋值 | [W_RIS/wris_mc.m:44](<C:/Users/Spidey/Desktop/code/W_RIS/wris_mc.m:44>) |
| `pp` | 本轮存活光子的当前位置（m）/世界方向/无量纲权重 | `p(id,:)` | 赋值 | [W_RIS/wris_mc.m:44](<C:/Users/Spidey/Desktop/code/W_RIS/wris_mc.m:44>) |
| `dd` | 本轮存活光子的当前位置（m）/世界方向/无量纲权重 | `d(id,:)` | 赋值 | [W_RIS/wris_mc.m:44](<C:/Users/Spidey/Desktop/code/W_RIS/wris_mc.m:44>) |
| `ww` | 本轮存活光子的当前位置（m）/世界方向/无量纲权重 | `w(id)` | 赋值 | [W_RIS/wris_mc.m:44](<C:/Users/Spidey/Desktop/code/W_RIS/wris_mc.m:44>) |
| `free` | 指数自由程 -log(U)/ext（m），透明水体取 Inf | `-log(max(rand(numel(id),1),realmin))/ext`<br>`inf(numel(id),1)` | 赋值 | [W_RIS/wris_mc.m:45](<C:/Users/Spidey/Desktop/code/W_RIS/wris_mc.m:45>) |
| `remaining` | 当前光子剩余路径额度 maxPath-path（m） | `c.mc.maxPath-path(id)` | 赋值 | [W_RIS/wris_mc.m:47](<C:/Users/Spidey/Desktop/code/W_RIS/wris_mc.m:47>) |
| `travel` | 本轮可走距离 min(free,remaining)（m） | `min(free,remaining)` | 赋值 | [W_RIS/wris_mc.m:47](<C:/Users/Spidey/Desktop/code/W_RIS/wris_mc.m:47>) |
| `local` | 当前射线与无限 RIS 平面交点的局部三维坐标（m） | `wris_intersect(pp,dd,g.pRis,g.RRis,[c.ris.Lx c.ris.Ly])` | 赋值 | [W_RIS/wris_mc.m:48](<C:/Users/Spidey/Desktop/code/W_RIS/wris_mc.m:48>) |
| `dist` | 当前射线到 RIS 平面的前向距离（m），非朝向接收面时 Inf | `wris_intersect(pp,dd,g.pRis,g.RRis,[c.ris.Lx c.ris.Ly])`<br>`Inf` | 赋值 | [W_RIS/wris_mc.m:48](<C:/Users/Spidey/Desktop/code/W_RIS/wris_mc.m:48>) |
| `toward` | 当前射线沿板法线正方向传播的掩码 dd*normal>1e-12 | `dd*g.normal>1e-12` | 赋值 | [W_RIS/wris_mc.m:49](<C:/Users/Spidey/Desktop/code/W_RIS/wris_mc.m:49>) |
| `denom` | 当前方向与 TX 光轴的点积，用于源平面反向逃逸判断 | `dd*g.dTx` | 赋值 | [W_RIS/wris_mc.m:51](<C:/Users/Spidey/Desktop/code/W_RIS/wris_mc.m:51>) |
| `back` | 到 TX 源平面的反向逃逸前向距离（m），无效为 Inf | `((g.pTx.'-pp)*g.dTx)./denom`<br>`Inf` | 赋值 | [W_RIS/wris_mc.m:52](<C:/Users/Spidey/Desktop/code/W_RIS/wris_mc.m:52>) |
| `reach` | 在碰撞/路径上限/反向逃逸之前首次到达接收平面的掩码 | `dist<=travel & dist<=back & isfinite(dist)` | 赋值 | [W_RIS/wris_mc.m:54](<C:/Users/Spidey/Desktop/code/W_RIS/wris_mc.m:54>) |
| `escape` | 在到达接收平面/碰撞/路径上限前反向穿过源平面的掩码 | `back<=travel & back<dist & isfinite(back)` | 赋值 | [W_RIS/wris_mc.m:55](<C:/Users/Spidey/Desktop/code/W_RIS/wris_mc.m:55>) |
| `at` | 本轮到达评分面光子的二维局部坐标（m） | `local(reach,1:2)` | 赋值 | [W_RIS/wris_mc.m:56](<C:/Users/Spidey/Desktop/code/W_RIS/wris_mc.m:56>) |
| `wr` | 本轮到达评分面的光子权重 | `ww(reach)` | 赋值 | [W_RIS/wris_mc.m:56](<C:/Users/Spidey/Desktop/code/W_RIS/wris_mc.m:56>) |
| `inside` | 本轮到达点是否位于有限矩形板孔径内的掩码 | `abs(at(:,1))<=c.ris.Lx/2 & abs(at(:,2))<=c.ris.Ly/2` | 赋值 | [W_RIS/wris_mc.m:57](<C:/Users/Spidey/Desktop/code/W_RIS/wris_mc.m:57>) |
| `rid` | 本轮首次到达评分面光子的原始批内索引 | `id(reach)` | 赋值 | [W_RIS/wris_mc.m:58](<C:/Users/Spidey/Desktop/code/W_RIS/wris_mc.m:58>) |
| `captured` | 本轮板内接受光子的原始批内索引 | `rid(inside)` | 赋值 | [W_RIS/wris_mc.m:58](<C:/Users/Spidey/Desktop/code/W_RIS/wris_mc.m:58>) |
| `move` | 尚未到达接收面且未反向逃逸、继续走本轮距离的掩码 | `~reach & ~escape` | 赋值 | [W_RIS/wris_mc.m:71](<C:/Users/Spidey/Desktop/code/W_RIS/wris_mc.m:71>) |
| `mid` | 继续移动光子的原始批内索引 | `id(move)` | 赋值 | [W_RIS/wris_mc.m:73](<C:/Users/Spidey/Desktop/code/W_RIS/wris_mc.m:73>) |
| `cutoff` | 自由程大于剩余路径预算、达到路径上限的掩码 | `free(move)>=remaining(move)` | 赋值 | [W_RIS/wris_mc.m:75](<C:/Users/Spidey/Desktop/code/W_RIS/wris_mc.m:75>) |
| `cut` | 达到路径上限或散射次数上限的光子索引 | `mid(cutoff)`<br>`collide(scatter(collide)>=c.mc.maxScatter & active(collide))` | 赋值 | [W_RIS/wris_mc.m:75](<C:/Users/Spidey/Desktop/code/W_RIS/wris_mc.m:75>) |
| `collide` | 本轮真实碰撞点处的光子索引 | `mid(~cutoff)` | 赋值 | [W_RIS/wris_mc.m:77](<C:/Users/Spidey/Desktop/code/W_RIS/wris_mc.m:77>) |
| `low` | 正权重但低于轮盘赌门限的光子索引 | `collide(w(collide)<c.mc.rouletteThreshold & w(collide)>0)` | 赋值 | [W_RIS/wris_mc.m:81](<C:/Users/Spidey/Desktop/code/W_RIS/wris_mc.m:81>) |
| `before` | 轮盘赌前的包权重 | `w(low)` | 赋值 | [W_RIS/wris_mc.m:82](<C:/Users/Spidey/Desktop/code/W_RIS/wris_mc.m:82>) |
| `survive` | 轮盘赌存活伯努利随机掩码 | `rand(numel(low),1)<c.mc.rouletteSurvival` | 赋值 | [W_RIS/wris_mc.m:82](<C:/Users/Spidey/Desktop/code/W_RIS/wris_mc.m:82>) |
| `sid` | 发生碰撞后仍存活、需更新散射方向的光子索引 | `collide(active(collide))` | 赋值 | [W_RIS/wris_mc.m:89](<C:/Users/Spidey/Desktop/code/W_RIS/wris_mc.m:89>) |
| `ns` | 当前需抽样散射方向的光子数 | `numel(sid)` | 赋值 | [W_RIS/wris_mc.m:89](<C:/Users/Spidey/Desktop/code/W_RIS/wris_mc.m:89>) |
| `theta` | HG 抽样散射极角（rad） | `Hen_Green(ns,2).'`<br>`acos(max(-1,min(1,mu)))` | 赋值 | [W_RIS/wris_mc.m:92](<C:/Users/Spidey/Desktop/code/W_RIS/wris_mc.m:92>) |
| `u` | 非默认 HG 分支的均匀随机变量 [0,1) | `rand(ns,1)` | 赋值 | [W_RIS/wris_mc.m:94](<C:/Users/Spidey/Desktop/code/W_RIS/wris_mc.m:94>) |
| `gg` | 非默认 HG 分支使用的各向异性因子 | `c.water.anisotropy` | 赋值 | [W_RIS/wris_mc.m:94](<C:/Users/Spidey/Desktop/code/W_RIS/wris_mc.m:94>) |
| `mu` | 非默认 HG 逆变换得到的 cos(theta) | `2*u-1`<br>`(1+gg^2-((1-gg^2)./(1-gg+2*gg*u)).^2)/(2*gg)` | 赋值 | [W_RIS/wris_mc.m:95](<C:/Users/Spidey/Desktop/code/W_RIS/wris_mc.m:95>) |
| `phi` | 散射方位角，均匀 [0,2*pi)（rad） | `2*pi*rand(ns,1)` | 赋值 | [W_RIS/wris_mc.m:99](<C:/Users/Spidey/Desktop/code/W_RIS/wris_mc.m:99>) |
| `oldDir` | 散射前的世界单位方向 | `d(sid,:)` | 赋值 | [W_RIS/wris_mc.m:99](<C:/Users/Spidey/Desktop/code/W_RIS/wris_mc.m:99>) |
| `helper` | 构造局部散射基的辅助轴；靠近极点时由 z 轴改用 x 轴 | `repmat([0 0 1],ns,1)`<br>`repmat([1 0 0],nnz(abs(oldDir(:,3))>0.9),1)` | 赋值 | [W_RIS/wris_mc.m:101](<C:/Users/Spidey/Desktop/code/W_RIS/wris_mc.m:101>) |
| `e1` | 与散射前方向正交的两个单位基方向 | `cross(helper,oldDir,2)`<br>`e1./vecnorm(e1,2,2)` | 赋值 | [W_RIS/wris_mc.m:102](<C:/Users/Spidey/Desktop/code/W_RIS/wris_mc.m:102>) |
| `e2` | 与散射前方向正交的两个单位基方向 | `cross(oldDir,e1,2)` | 赋值 | [W_RIS/wris_mc.m:102](<C:/Users/Spidey/Desktop/code/W_RIS/wris_mc.m:102>) |
| `newDir` | 抽样散射角变换得到的世界新传播方向，随后再次归一化 | `cos(theta).*oldDir+sin(theta).*(cos(phi).*e1+sin(phi).*e2)` | 赋值 | [W_RIS/wris_mc.m:103](<C:/Users/Spidey/Desktop/code/W_RIS/wris_mc.m:103>) |
| `valid` | 当前批次板内接受且首次到达权重>0 的掩码 | `accepted & hitW>0` | 赋值 | [W_RIS/wris_mc.m:106](<C:/Users/Spidey/Desktop/code/W_RIS/wris_mc.m:106>) |
| `ix` | 板内命中位置在 x/y 方向的分箱索引 | `discretize(hitXY(valid,1),xe)` | 赋值 | [W_RIS/wris_mc.m:110](<C:/Users/Spidey/Desktop/code/W_RIS/wris_mc.m:110>) |
| `iy` | 板内命中位置在 x/y 方向的分箱索引 | `discretize(hitXY(valid,2),ye)` | 赋值 | [W_RIS/wris_mc.m:110](<C:/Users/Spidey/Desktop/code/W_RIS/wris_mc.m:110>) |
| `landed` | 首次到达评分面且权重>0 的掩码 | `hitW>0` | 赋值 | [W_RIS/wris_mc.m:114](<C:/Users/Spidey/Desktop/code/W_RIS/wris_mc.m:114>) |
| `pix` | 到达点在延展评分窗 x/y 方向的分箱索引 | `discretize(hitXY(landed,1),pxe)` | 赋值 | [W_RIS/wris_mc.m:115](<C:/Users/Spidey/Desktop/code/W_RIS/wris_mc.m:115>) |
| `piy` | 到达点在延展评分窗 x/y 方向的分箱索引 | `discretize(hitXY(landed,2),pye)` | 赋值 | [W_RIS/wris_mc.m:115](<C:/Users/Spidey/Desktop/code/W_RIS/wris_mc.m:115>) |
| `shown` | 落在有限延展绘图窗口内的到达点掩码 | `isfinite(pix) & isfinite(piy)` | 赋值 | [W_RIS/wris_mc.m:116](<C:/Users/Spidey/Desktop/code/W_RIS/wris_mc.m:116>) |
| `landedW` | 首次到达无限评分面光子的权重 | `hitW(landed)` | 赋值 | [W_RIS/wris_mc.m:116](<C:/Users/Spidey/Desktop/code/W_RIS/wris_mc.m:116>) |
| `core` | 板内接受并且没有散射的核心掩码 | `accepted & direct` | 赋值 | [W_RIS/wris_mc.m:120](<C:/Users/Spidey/Desktop/code/W_RIS/wris_mc.m:120>) |
| `xy` | 全部板内到达二维坐标（m） | `vertcat(xyAll{:})` | 赋值 | [W_RIS/wris_mc.m:124](<C:/Users/Spidey/Desktop/code/W_RIS/wris_mc.m:124>) |
| `weights` | 全部板内到达点对应包权重 | `vertcat(wAll{:})` | 赋值 | [W_RIS/wris_mc.m:124](<C:/Users/Spidey/Desktop/code/W_RIS/wris_mc.m:124>) |
| `scale` | 固定发射归一化功率 source.power/N（W/初始包权重） | `c.source.power/N` | 赋值 | [W_RIS/wris_mc.m:125](<C:/Users/Spidey/Desktop/code/W_RIS/wris_mc.m:125>) |
| `out.irradiance` | 板内含 EGG 公共增益的辐照度网格（W/m²） | `hist*scale/area*gain` | 赋值 | [W_RIS/wris_mc.m:126](<C:/Users/Spidey/Desktop/code/W_RIS/wris_mc.m:126>) |
| `out.noTurbulence` | 除去 EGG 后的辐照度网格 I/H（W/m²），主图实际显示此量 | `hist*scale/area` | 赋值 | [W_RIS/wris_mc.m:126](<C:/Users/Spidey/Desktop/code/W_RIS/wris_mc.m:126>) |
| `out.standardError` | 逐格含 EGG 的辐照度蒙特卡洛标准误（W/m²） | `sqrt(max(0,hist2-hist.^2/N)/(N*(N-1)))*c.source.power/area*gain` | 赋值 | [W_RIS/wris_mc.m:127](<C:/Users/Spidey/Desktop/code/W_RIS/wris_mc.m:127>) |
| `out.x` | 板面/延展面分箱中心 x/y 坐标（m），由父结构区分 | `(xe(1:end-1)+xe(2:end))/2` | 赋值 | [W_RIS/wris_mc.m:128](<C:/Users/Spidey/Desktop/code/W_RIS/wris_mc.m:128>) |
| `out.y` | 板面/延展面分箱中心 x/y 坐标（m），由父结构区分 | `(ye(1:end-1)+ye(2:end))/2` | 赋值 | [W_RIS/wris_mc.m:128](<C:/Users/Spidey/Desktop/code/W_RIS/wris_mc.m:128>) |
| `out.cellArea` | 板面单格面积（m²） | `area` | 赋值 | [W_RIS/wris_mc.m:129](<C:/Users/Spidey/Desktop/code/W_RIS/wris_mc.m:129>) |
| `out.capturedPower` | 板内接收功率，含公共 EGG 增益（W） | `ledger.inside*scale*gain` | 赋值 | [W_RIS/wris_mc.m:129](<C:/Users/Spidey/Desktop/code/W_RIS/wris_mc.m:129>) |
| `out.plane` | 延展评分面结构；不是实际接收孔径 | `struct(...)（所有字段在本节逐行列明）` | 赋值 | [W_RIS/wris_mc.m:130](<C:/Users/Spidey/Desktop/code/W_RIS/wris_mc.m:130>) |
| `out.plane.x` | 板面/延展面分箱中心 x/y 坐标（m），由父结构区分 | `(pxe(1:end-1)+pxe(2:end))/2` | 结构字段 | [W_RIS/wris_mc.m:130](<C:/Users/Spidey/Desktop/code/W_RIS/wris_mc.m:130>) |
| `out.plane.y` | 板面/延展面分箱中心 x/y 坐标（m），由父结构区分 | `(pye(1:end-1)+pye(2:end))/2` | 结构字段 | [W_RIS/wris_mc.m:131](<C:/Users/Spidey/Desktop/code/W_RIS/wris_mc.m:131>) |
| `out.plane.noTurbulence` | 除去 EGG 后的辐照度网格 I/H（W/m²），主图实际显示此量 | `planeHist*scale/area` | 结构字段 | [W_RIS/wris_mc.m:131](<C:/Users/Spidey/Desktop/code/W_RIS/wris_mc.m:131>) |
| `out.plane.powerWithoutTurbulence` | 不含 EGG 的接收功率（W）；plane 下为整个无限评分面的 inside+outside 功率 | `(ledger.inside+ledger.outside)*scale` | 结构字段 | [W_RIS/wris_mc.m:132](<C:/Users/Spidey/Desktop/code/W_RIS/wris_mc.m:132>) |
| `out.powerWithoutTurbulence` | 不含 EGG 的接收功率（W）；plane 下为整个无限评分面的 inside+outside 功率 | `ledger.inside*scale` | 赋值 | [W_RIS/wris_mc.m:133](<C:/Users/Spidey/Desktop/code/W_RIS/wris_mc.m:133>) |
| `out.powerSE` | 板内不含 EGG 的 MC 接收功率标准误（W） | `c.source.power*sqrt(max(0,sum(weights.^2)-sum(weights)^2/N)/(N*(N-1)))` | 赋值 | [W_RIS/wris_mc.m:134](<C:/Users/Spidey/Desktop/code/W_RIS/wris_mc.m:134>) |
| `out.peak` | 板内含 EGG 的最大辐照度（W/m²） | `max(out.irradiance,[],'all')` | 赋值 | [W_RIS/wris_mc.m:135](<C:/Users/Spidey/Desktop/code/W_RIS/wris_mc.m:135>) |
| `out.gain` | 本帧所有光子共享的 EGG 增益（无量纲） | `gain` | 赋值 | [W_RIS/wris_mc.m:135](<C:/Users/Spidey/Desktop/code/W_RIS/wris_mc.m:135>) |
| `out.metrics` | 板内所有命中光子的加权二阶矩光斑统计 | `wris_metrics(xy,weights)` | 赋值 | [W_RIS/wris_mc.m:136](<C:/Users/Spidey/Desktop/code/W_RIS/wris_mc.m:136>) |
| `out.directMetrics` | 板内未散射核心的加权二阶矩统计 | `wris_metrics(vertcat(xyDirect{:}),vertcat(wDirect{:}))` | 赋值 | [W_RIS/wris_mc.m:137](<C:/Users/Spidey/Desktop/code/W_RIS/wris_mc.m:137>) |
| `out.directPlaneMetrics` | 完整未散射核心的加权二阶矩统计，包括板外及绘图窗口外 | `wris_metrics(vertcat(xyCorePlane{:}),vertcat(wCorePlane{:}))` | 赋值 | [W_RIS/wris_mc.m:138](<C:/Users/Spidey/Desktop/code/W_RIS/wris_mc.m:138>) |
| `out.directPowerWithoutTurbulence` | 板内未散射核心功率（W，不含 EGG） | `sum(vertcat(wDirect{:}))*scale` | 赋值 | [W_RIS/wris_mc.m:139](<C:/Users/Spidey/Desktop/code/W_RIS/wris_mc.m:139>) |
| `out.directPlanePowerWithoutTurbulence` | 完整首次到达评分平面的未散射核心功率（W，不含 EGG） | `sum(vertcat(wCorePlane{:}))*scale` | 赋值 | [W_RIS/wris_mc.m:140](<C:/Users/Spidey/Desktop/code/W_RIS/wris_mc.m:140>) |
| `out.coreClippedFraction` | 1-板内核心功率/完整核心功率（无量纲） | `1-out.directPowerWithoutTurbulence/max(out.directPlanePowerWithoutTurbulence,realmin)` | 赋值 | [W_RIS/wris_mc.m:141](<C:/Users/Spidey/Desktop/code/W_RIS/wris_mc.m:141>) |
| `out.meanIncidentDirectionLocal` | 板内到达方向的权重平均，在 RIS 局部系表示；未再单位化 | `localDirectionSum/max(ledger.inside,realmin)` | 赋值 | [W_RIS/wris_mc.m:142](<C:/Users/Spidey/Desktop/code/W_RIS/wris_mc.m:142>) |
| `out.ledger` | 光子包无量纲权重能量账本；乘 source.power/N 转为 W | `ledger` | 赋值 | [W_RIS/wris_mc.m:143](<C:/Users/Spidey/Desktop/code/W_RIS/wris_mc.m:143>) |
| `out.counts` | 各终止类别的光子包个数 | `counts` | 赋值 | [W_RIS/wris_mc.m:143](<C:/Users/Spidey/Desktop/code/W_RIS/wris_mc.m:143>) |
| `out.seed` | 随机数种子；所属结构决定是海浪、TX、MC 或湍流随机流 | `seed` | 赋值 | [W_RIS/wris_mc.m:143](<C:/Users/Spidey/Desktop/code/W_RIS/wris_mc.m:143>) |
| `out.photons` | 单帧发射光子包总数，或验证采用的光子数档位 | `N` | 赋值 | [W_RIS/wris_mc.m:143](<C:/Users/Spidey/Desktop/code/W_RIS/wris_mc.m:143>) |
| `out.energyResidual` | 扣除轮盘赌权重变化后的相对能量守恒残差 | `(ledger.inside+ledger.outside+ledger.backEscape+ledger.absorbed+ledger.truncated-N-ledger.rouletteDelta)/N` | 赋值 | [W_RIS/wris_mc.m:144](<C:/Users/Spidey/Desktop/code/W_RIS/wris_mc.m:144>) |
| `out.geometry` | 该帧使用的光学几何结构 | `g` | 赋值 | [W_RIS/wris_mc.m:145](<C:/Users/Spidey/Desktop/code/W_RIS/wris_mc.m:145>) |


### 5.23 W_RIS/wris_metrics.m

**作用域：`wris_metrics`**

| 源码符号 | 含义与单位/维度 | 源码取值/计算式（各次赋值） | 角色 | 源码位置 |
| --- | --- | --- | --- | --- |
| `xy` | 被统计光斑的二维命中位置（N×2，m） | 函数输入/输出；值随调用或上游结果确定 | 输入 | [W_RIS/wris_metrics.m:1](<C:/Users/Spidey/Desktop/code/W_RIS/wris_metrics.m:1>) |
| `w` | 与 xy 对应的非负光子包权重 | 函数输入/输出；值随调用或上游结果确定 | 输入 | [W_RIS/wris_metrics.m:1](<C:/Users/Spidey/Desktop/code/W_RIS/wris_metrics.m:1>) |
| `s` | 光斑中心、D4σ 长短轴、轴比和主轴角结构 | `struct('center',[NaN NaN],'axes',[NaN NaN],'ratio',NaN,'orientationDeg',NaN)` | 输出/赋值 | [W_RIS/wris_metrics.m:1](<C:/Users/Spidey/Desktop/code/W_RIS/wris_metrics.m:1>) |
| `s.center` | RIS 的水平参考中心 [x;y]（m），或光斑加权中心 [x,y]（m）；以父结构为准 | `[NaN NaN]`<br>`center` | 结构字段/赋值 | [W_RIS/wris_metrics.m:4](<C:/Users/Spidey/Desktop/code/W_RIS/wris_metrics.m:4>) |
| `s.axes` | 二阶矩 D4σ 长/短轴全直径 [major,minor]（m） | `[NaN NaN]`<br>`diam.'` | 结构字段/赋值 | [W_RIS/wris_metrics.m:4](<C:/Users/Spidey/Desktop/code/W_RIS/wris_metrics.m:4>) |
| `s.ratio` | 长短轴全直径之比（无量纲） | `NaN`<br>`diam(1)/max(diam(2),realmin)` | 结构字段/赋值 | [W_RIS/wris_metrics.m:4](<C:/Users/Spidey/Desktop/code/W_RIS/wris_metrics.m:4>) |
| `s.orientationDeg` | 加权协方差最大特征向量的主轴角（°，[-90,90)） | `NaN`<br>`mod(atan2d(v(2),v(1))+90,180)-90` | 结构字段/赋值 | [W_RIS/wris_metrics.m:4](<C:/Users/Spidey/Desktop/code/W_RIS/wris_metrics.m:4>) |
| `sw` | 总包权重 sum(w) | `sum(w)` | 赋值 | [W_RIS/wris_metrics.m:6](<C:/Users/Spidey/Desktop/code/W_RIS/wris_metrics.m:6>) |
| `center` | 二维光斑权重均值中心（m） | `sum(xy.*w,1)/sw` | 赋值 | [W_RIS/wris_metrics.m:6](<C:/Users/Spidey/Desktop/code/W_RIS/wris_metrics.m:6>) |
| `q` | 命中位置减加权中心后的二维偏差（m） | `xy-center` | 赋值 | [W_RIS/wris_metrics.m:6](<C:/Users/Spidey/Desktop/code/W_RIS/wris_metrics.m:6>) |
| `covariance` | 按总权重归一化的二维加权协方差矩阵（m²） | `q.'*(q.*w)/sw` | 赋值 | [W_RIS/wris_metrics.m:7](<C:/Users/Spidey/Desktop/code/W_RIS/wris_metrics.m:7>) |
| `V` | 协方差的特征向量矩阵/特征值列向量 | `eig(covariance,'vector')` | 赋值 | [W_RIS/wris_metrics.m:7](<C:/Users/Spidey/Desktop/code/W_RIS/wris_metrics.m:7>) |
| `D` | 协方差的特征向量矩阵/特征值列向量 | `eig(covariance,'vector')` | 赋值 | [W_RIS/wris_metrics.m:7](<C:/Users/Spidey/Desktop/code/W_RIS/wris_metrics.m:7>) |
| `e` | 降序排列且截去负舍入误差的协方差特征值（m²） | `sort(max(0,D),'descend')` | 赋值 | [W_RIS/wris_metrics.m:8](<C:/Users/Spidey/Desktop/code/W_RIS/wris_metrics.m:8>) |
| `order` | 特征值由大到小的排列索引 | `sort(max(0,D),'descend')` | 赋值 | [W_RIS/wris_metrics.m:8](<C:/Users/Spidey/Desktop/code/W_RIS/wris_metrics.m:8>) |
| `diam` | D4σ 长短轴全直径 4*sqrt(e)（m） | `4*sqrt(e)` | 赋值 | [W_RIS/wris_metrics.m:8](<C:/Users/Spidey/Desktop/code/W_RIS/wris_metrics.m:8>) |
| `v` | 最大特征值对应的二维单位主轴方向 | `V(:,order(1))` | 赋值 | [W_RIS/wris_metrics.m:10](<C:/Users/Spidey/Desktop/code/W_RIS/wris_metrics.m:10>) |


### 5.24 W_RIS/wris_plot.m

**作用域：`wris_plot`**

| 源码符号 | 含义与单位/维度 | 源码取值/计算式（各次赋值） | 角色 | 源码位置 |
| --- | --- | --- | --- | --- |
| `s` | 完整仿真结果结构 | 函数输入/输出；值随调用或上游结果确定 | 输入 | [W_RIS/wris_plot.m:1](<C:/Users/Spidey/Desktop/code/W_RIS/wris_plot.m:1>) |
| `c` | 当前 W_RIS 总配置结构 | `s.config` | 赋值 | [W_RIS/wris_plot.m:3](<C:/Users/Spidey/Desktop/code/W_RIS/wris_plot.m:3>) |
| `peak` | 所有 15 个延展评分窗 I/H 的最大值（W/m²），用于共享色标 | `max(cellfun(@(p) max(p.out.plane.noTurbulence,[],'all'),s.panels),[],'all')` | 赋值 | [W_RIS/wris_plot.m:3](<C:/Users/Spidey/Desktop/code/W_RIS/wris_plot.m:3>) |
| `p` | 当前面板结构；cellfun 中也是面板输入 | `peak=max(cellfun(@(p) max(p.out.plane.noTurbulence,[],'all'),s.panels),[],'all')`<br>`s.panels{row,col}` | 匿名函数输入/赋值 | [W_RIS/wris_plot.m:3](<C:/Users/Spidey/Desktop/code/W_RIS/wris_plot.m:3>) |
| `limits` | 所有光斑图共用的对数辐照度色限（W/m²） | `[max(peak,eps)/c.plot.dynamicRange max(peak,eps)]` | 赋值 | [W_RIS/wris_plot.m:4](<C:/Users/Spidey/Desktop/code/W_RIS/wris_plot.m:4>) |
| `f` | 图窗句柄；频率等重名按所在函数解释 | `figure('Name','RIS rotation and aperture clipping','Visible','off','Color','w')` | 赋值 | [W_RIS/wris_plot.m:5](<C:/Users/Spidey/Desktop/code/W_RIS/wris_plot.m:5>) |
| `f.Position` | 图窗像素位置 [left,bottom,width,height] | `[5 5 1500 980]` | 赋值 | [W_RIS/wris_plot.m:6](<C:/Users/Spidey/Desktop/code/W_RIS/wris_plot.m:6>) |
| `guard` | 函数结束时自动恢复随机数状态、关闭图窗或文件的清理对象 | `onCleanup(@() close(f))` | 赋值 | [W_RIS/wris_plot.m:6](<C:/Users/Spidey/Desktop/code/W_RIS/wris_plot.m:6>) |
| `layout` | 分面绘图布局句柄 | `tiledlayout(f,3,5,'TileSpacing','compact','Padding','compact')` | 赋值 | [W_RIS/wris_plot.m:7](<C:/Users/Spidey/Desktop/code/W_RIS/wris_plot.m:7>) |
| `missAxes` | 光轴完全移出主显示范围的坐标轴句柄/面板结果列表 | `cell(15,1)`<br>`ax` | 赋值 | [W_RIS/wris_plot.m:8](<C:/Users/Spidey/Desktop/code/W_RIS/wris_plot.m:8>) |
| `missPanels` | 光轴完全移出主显示范围的坐标轴句柄/面板结果列表 | `cell(15,1)`<br>`p` | 赋值 | [W_RIS/wris_plot.m:8](<C:/Users/Spidey/Desktop/code/W_RIS/wris_plot.m:8>) |
| `nMiss` | 需要附加虚拟评分面小图的面板数量 | `0`<br>`nMiss+1` | 赋值 | [W_RIS/wris_plot.m:8](<C:/Users/Spidey/Desktop/code/W_RIS/wris_plot.m:8>) |
| `row` | 工况表的行或列索引 | `1:3` | 循环索引 | [W_RIS/wris_plot.m:9](<C:/Users/Spidey/Desktop/code/W_RIS/wris_plot.m:9>) |
| `reference` | 同一行静水参考工况的投影主轴角（°） | `s.panels{row,1}.out.geometry.projectedOrientationDeg` | 赋值 | [W_RIS/wris_plot.m:10](<C:/Users/Spidey/Desktop/code/W_RIS/wris_plot.m:10>) |
| `col` | 工况表的行或列索引 | `1:5` | 循环索引 | [W_RIS/wris_plot.m:11](<C:/Users/Spidey/Desktop/code/W_RIS/wris_plot.m:11>) |
| `o` | 当前面板 MC 输出结构 | `p.out` | 赋值 | [W_RIS/wris_plot.m:12](<C:/Users/Spidey/Desktop/code/W_RIS/wris_plot.m:12>) |
| `ax` | 绘图坐标轴句柄 | `nexttile(layout)` | 赋值 | [W_RIS/wris_plot.m:12](<C:/Users/Spidey/Desktop/code/W_RIS/wris_plot.m:12>) |
| `extent` | (0.5+plot.padding)*[Lx,Ly] 主图水平半宽（m） | `(0.5+c.plot.padding)*[c.ris.Lx c.ris.Ly]` | 赋值 | [W_RIS/wris_plot.m:15](<C:/Users/Spidey/Desktop/code/W_RIS/wris_plot.m:15>) |
| `delta` | 当前几何主轴相对本行静水参考的模 180° 旋转角（°） | `mod(o.geometry.projectedOrientationDeg-reference+90,180)-90` | 赋值 | [W_RIS/wris_plot.m:22](<C:/Users/Spidey/Desktop/code/W_RIS/wris_plot.m:22>) |
| `label` | 工况、相位或图形标签文本 | `'Calm \| TX/RIS fixed'`<br>`sprintf('U=%g m/s \| wind=%g deg \| %s',p.windMps,p.windDirectionDeg,shortPhase)` | 赋值 | [W_RIS/wris_plot.m:23](<C:/Users/Spidey/Desktop/code/W_RIS/wris_plot.m:23>) |
| `shortPhase` | 删去 ' site' 后的简短站位标签 | `strrep(p.phase,' site','')` | 赋值 | [W_RIS/wris_plot.m:25](<C:/Users/Spidey/Desktop/code/W_RIS/wris_plot.m:25>) |
| `cb` | 颜色条句柄 | `colorbar(ax)` | 赋值 | [W_RIS/wris_plot.m:37](<C:/Users/Spidey/Desktop/code/W_RIS/wris_plot.m:37>) |
| `cb.Layout.Tile` | 颜色条在 tiledlayout 中的位置，取 east | `'east'` | 赋值 | [W_RIS/wris_plot.m:37](<C:/Users/Spidey/Desktop/code/W_RIS/wris_plot.m:37>) |
| `cb.Label.String` | 颜色条标题文本 | `'I/H (W/m^2), shared logarithmic scale'` | 赋值 | [W_RIS/wris_plot.m:37](<C:/Users/Spidey/Desktop/code/W_RIS/wris_plot.m:37>) |
| `j` | 当前循环、分量或样本索引（无量纲，MATLAB 从 1 开始） | `1:nMiss` | 循环索引 | [W_RIS/wris_plot.m:44](<C:/Users/Spidey/Desktop/code/W_RIS/wris_plot.m:44>) |


**作用域：`drawCore`**

| 源码符号 | 含义与单位/维度 | 源码取值/计算式（各次赋值） | 角色 | 源码位置 |
| --- | --- | --- | --- | --- |
| `ax` | 绘图坐标轴句柄 | 函数输入/输出；值随调用或上游结果确定 | 输入 | [W_RIS/wris_plot.m:49](<C:/Users/Spidey/Desktop/code/W_RIS/wris_plot.m:49>) |
| `o` | 当前 RIS 运动结果或动力学模型；具体类型由所在函数输入确定 | 函数输入/输出；值随调用或上游结果确定 | 输入 | [W_RIS/wris_plot.m:49](<C:/Users/Spidey/Desktop/code/W_RIS/wris_plot.m:49>) |
| `reference` | 同一行静水参考工况的投影主轴角（°） | 函数输入/输出；值随调用或上游结果确定 | 输入 | [W_RIS/wris_plot.m:49](<C:/Users/Spidey/Desktop/code/W_RIS/wris_plot.m:49>) |
| `m` | 包含板外的完整未散射核心二阶矩统计 | `o.directPlaneMetrics` | 赋值 | [W_RIS/wris_plot.m:51](<C:/Users/Spidey/Desktop/code/W_RIS/wris_plot.m:51>) |
| `angle` | 完整核心 MC 主轴角转为 rad | `deg2rad(m.orientationDeg)` | 赋值 | [W_RIS/wris_plot.m:53](<C:/Users/Spidey/Desktop/code/W_RIS/wris_plot.m:53>) |
| `R` | 二维椭圆主轴旋转矩阵 | `[cos(angle) -sin(angle);sin(angle) cos(angle)]` | 赋值 | [W_RIS/wris_plot.m:53](<C:/Users/Spidey/Desktop/code/W_RIS/wris_plot.m:53>) |
| `parameter` | 椭圆轮廓参数角 0..2*pi 的 121 个采样点（rad） | `linspace(0,2*pi,121)` | 赋值 | [W_RIS/wris_plot.m:54](<C:/Users/Spidey/Desktop/code/W_RIS/wris_plot.m:54>) |
| `ellipse` | 完整核心 D4σ 椭圆轮廓坐标（121×2，m） | `(R*(m.axes(:)/2.*[cos(parameter);sin(parameter)])).'+m.center` | 赋值 | [W_RIS/wris_plot.m:55](<C:/Users/Spidey/Desktop/code/W_RIS/wris_plot.m:55>) |
| `half` | 核心长轴半直径（m） | `m.axes(1)/2` | 赋值 | [W_RIS/wris_plot.m:57](<C:/Users/Spidey/Desktop/code/W_RIS/wris_plot.m:57>) |
| `baseline` | 同一行静水参考方向的长轴半径向量（m） | `half*[cosd(reference);sind(reference)]` | 赋值 | [W_RIS/wris_plot.m:58](<C:/Users/Spidey/Desktop/code/W_RIS/wris_plot.m:58>) |
| `axisDirection` | 当前几何投影主轴方向的长轴半径向量（m） | `half*[cosd(o.geometry.projectedOrientationDeg);sind(o.geometry.projectedOrientationDeg)]` | 赋值 | [W_RIS/wris_plot.m:61](<C:/Users/Spidey/Desktop/code/W_RIS/wris_plot.m:61>) |


**作用域：`incomingInset`**

| 源码符号 | 含义与单位/维度 | 源码取值/计算式（各次赋值） | 角色 | 源码位置 |
| --- | --- | --- | --- | --- |
| `f` | 图窗句柄；频率等重名按所在函数解释 | 函数输入/输出；值随调用或上游结果确定 | 输入 | [W_RIS/wris_plot.m:66](<C:/Users/Spidey/Desktop/code/W_RIS/wris_plot.m:66>) |
| `parent` | 需要附加小图的父面板坐标轴句柄 | 函数输入/输出；值随调用或上游结果确定 | 输入 | [W_RIS/wris_plot.m:66](<C:/Users/Spidey/Desktop/code/W_RIS/wris_plot.m:66>) |
| `o` | 当前 RIS 运动结果或动力学模型；具体类型由所在函数输入确定 | 函数输入/输出；值随调用或上游结果确定 | 输入 | [W_RIS/wris_plot.m:66](<C:/Users/Spidey/Desktop/code/W_RIS/wris_plot.m:66>) |
| `c` | 当前 W_RIS 总配置结构 | 函数输入/输出；值随调用或上游结果确定 | 输入 | [W_RIS/wris_plot.m:66](<C:/Users/Spidey/Desktop/code/W_RIS/wris_plot.m:66>) |
| `limits` | 所有光斑图共用的对数辐照度色限（W/m²） | 函数输入/输出；值随调用或上游结果确定 | 输入 | [W_RIS/wris_plot.m:66](<C:/Users/Spidey/Desktop/code/W_RIS/wris_plot.m:66>) |
| `parent.Units` | 父面板坐标轴的位置单位，取 normalized | `'normalized'` | 赋值 | [W_RIS/wris_plot.m:70](<C:/Users/Spidey/Desktop/code/W_RIS/wris_plot.m:70>) |
| `q` | 父坐标轴在图窗中的归一化位置 [left,bottom,width,height] | `parent.Position` | 赋值 | [W_RIS/wris_plot.m:70](<C:/Users/Spidey/Desktop/code/W_RIS/wris_plot.m:70>) |
| `ax` | 绘图坐标轴句柄 | `axes(f,'Units','normalized','Position',[q(1)+0.14*q(3) q(2)+0.13*q(4) 0.74*q(3) 0.68*q(4)])` | 赋值 | [W_RIS/wris_plot.m:71](<C:/Users/Spidey/Desktop/code/W_RIS/wris_plot.m:71>) |
| `m` | 完整未散射核心统计结果 | `o.directPlaneMetrics` | 赋值 | [W_RIS/wris_plot.m:74](<C:/Users/Spidey/Desktop/code/W_RIS/wris_plot.m:74>) |
| `radius` | 完整核心最大 D4σ 半轴长（m） | `max(m.axes)/2` | 赋值 | [W_RIS/wris_plot.m:74](<C:/Users/Spidey/Desktop/code/W_RIS/wris_plot.m:74>) |
| `low` | 同时容纳有限板和完整核心范围的 x/y 下界/上界（m） | `min([-c.ris.Lx/2 -c.ris.Ly/2],m.center-radius)` | 赋值 | [W_RIS/wris_plot.m:75](<C:/Users/Spidey/Desktop/code/W_RIS/wris_plot.m:75>) |
| `high` | 同时容纳有限板和完整核心范围的 x/y 下界/上界（m） | `max([c.ris.Lx/2 c.ris.Ly/2],m.center+radius)` | 赋值 | [W_RIS/wris_plot.m:76](<C:/Users/Spidey/Desktop/code/W_RIS/wris_plot.m:76>) |
| `center` | 虚拟评分面小图显示窗口中心（m） | `(low+high)/2` | 赋值 | [W_RIS/wris_plot.m:77](<C:/Users/Spidey/Desktop/code/W_RIS/wris_plot.m:77>) |
| `half` | 正方形小图的显示半宽 0.55*max(high-low)（m） | `0.55*max(high-low)` | 赋值 | [W_RIS/wris_plot.m:77](<C:/Users/Spidey/Desktop/code/W_RIS/wris_plot.m:77>) |
| `ax.Title.BackgroundColor` | 小图标题背景色，取白色 w | `'w'` | 赋值 | [W_RIS/wris_plot.m:83](<C:/Users/Spidey/Desktop/code/W_RIS/wris_plot.m:83>) |


### 5.25 W_RIS/wris_reference.m

**作用域：`wris_reference`**

| 源码符号 | 含义与单位/维度 | 源码取值/计算式（各次赋值） | 角色 | 源码位置 |
| --- | --- | --- | --- | --- |
| `c` | 当前 W_RIS 总配置结构 | 函数输入/输出；值随调用或上游结果确定 | 输入 | [W_RIS/wris_reference.m:1](<C:/Users/Spidey/Desktop/code/W_RIS/wris_reference.m:1>) |
| `g` | 静水固定 RIS/TX 参考光学几何结构 | 函数输入/输出；值随调用或上游结果确定 | 输出 | [W_RIS/wris_reference.m:1](<C:/Users/Spidey/Desktop/code/W_RIS/wris_reference.m:1>) |
| `draft` | 静水平衡吃水 mass/(rho*Lx*Ly)（m） | `c.ris.mass/(c.ris.rho*c.ris.Lx*c.ris.Ly)` | 赋值 | [W_RIS/wris_reference.m:3](<C:/Users/Spidey/Desktop/code/W_RIS/wris_reference.m:3>) |
| `g.pRis` | RIS 下表面光学中心世界坐标（m） | `[c.ris.center;-draft]` | 赋值 | [W_RIS/wris_reference.m:4](<C:/Users/Spidey/Desktop/code/W_RIS/wris_reference.m:4>) |
| `th` | 参考入射角 incidenceDeg 转成 rad | `deg2rad(c.geometry.incidenceDeg)` | 赋值 | [W_RIS/wris_reference.m:5](<C:/Users/Spidey/Desktop/code/W_RIS/wris_reference.m:5>) |
| `az` | 参考光轴方位角 azimuthDeg 转成 rad | `deg2rad(c.geometry.azimuthDeg)` | 赋值 | [W_RIS/wris_reference.m:5](<C:/Users/Spidey/Desktop/code/W_RIS/wris_reference.m:5>) |
| `g.dTx` | TX 光轴在世界系的单位方向 | `[sin(th)*cos(az);sin(th)*sin(az);cos(th)]` | 赋值 | [W_RIS/wris_reference.m:6](<C:/Users/Spidey/Desktop/code/W_RIS/wris_reference.m:6>) |
| `g.pTx` | TX 发射面中心世界坐标（m） | `g.pRis-c.geometry.distance*g.dTx` | 赋值 | [W_RIS/wris_reference.m:7](<C:/Users/Spidey/Desktop/code/W_RIS/wris_reference.m:7>) |
| `g.RRis` | RIS 体坐标至世界坐标旋转矩阵 | `eye(3)` | 赋值 | [W_RIS/wris_reference.m:8](<C:/Users/Spidey/Desktop/code/W_RIS/wris_reference.m:8>) |
| `g.normal` | 板体 +z 轴的世界单位法向量（3×Nt） | `[0;0;1]` | 赋值 | [W_RIS/wris_reference.m:8](<C:/Users/Spidey/Desktop/code/W_RIS/wris_reference.m:8>) |
| `g.RTx` | 以光轴为局部 z 的 TX 光学系至世界系旋转矩阵 | `beamFrame(g.dTx)` | 赋值 | [W_RIS/wris_reference.m:10](<C:/Users/Spidey/Desktop/code/W_RIS/wris_reference.m:10>) |
| `g.distance` | TX 至 RIS 光学中心的空间直线距离（m） | `norm(g.pRis-g.pTx)` | 赋值 | [W_RIS/wris_reference.m:11](<C:/Users/Spidey/Desktop/code/W_RIS/wris_reference.m:11>) |
| `g.thetaDeg` | 当前光轴与板法线的锐角入射角 acosd(abs(dot(dTx,normal)))（°） | `acosd(abs(dot(g.dTx,g.normal)))` | 赋值 | [W_RIS/wris_reference.m:12](<C:/Users/Spidey/Desktop/code/W_RIS/wris_reference.m:12>) |
| `g.axisLocal` | 光轴与 RIS 无限平面交点的板局部坐标 [x,y,z]（m） | `wris_intersect(g.pTx.',g.dTx.',g.pRis,g.RRis,[c.ris.Lx c.ris.Ly])` | 赋值 | [W_RIS/wris_reference.m:13](<C:/Users/Spidey/Desktop/code/W_RIS/wris_reference.m:13>) |
| `g.axisRange` | 光轴从 TX 到平面交点的前向参数距离（m）；无效为 Inf | `wris_intersect(g.pTx.',g.dTx.',g.pRis,g.RRis,[c.ris.Lx c.ris.Ly])` | 赋值 | [W_RIS/wris_reference.m:13](<C:/Users/Spidey/Desktop/code/W_RIS/wris_reference.m:13>) |
| `g.axisInside` | 光轴交点是否在有限矩形板孔径内 | `wris_intersect(g.pTx.',g.dTx.',g.pRis,g.RRis,[c.ris.Lx c.ris.Ly])` | 赋值 | [W_RIS/wris_reference.m:13](<C:/Users/Spidey/Desktop/code/W_RIS/wris_reference.m:13>) |


**作用域：`beamFrame`**

| 源码符号 | 含义与单位/维度 | 源码取值/计算式（各次赋值） | 角色 | 源码位置 |
| --- | --- | --- | --- | --- |
| `d` | TX 世界光轴单位方向，作为光学局部 z 轴 | 函数输入/输出；值随调用或上游结果确定 | 输入 | [W_RIS/wris_reference.m:15](<C:/Users/Spidey/Desktop/code/W_RIS/wris_reference.m:15>) |
| `R` | 以 [u,v,d] 为列的光学系到世界系旋转矩阵 | `[u,v,d]` | 输出/赋值 | [W_RIS/wris_reference.m:15](<C:/Users/Spidey/Desktop/code/W_RIS/wris_reference.m:15>) |
| `u` | 与 d 正交的光学局部 x 轴世界单位方向 cross([1;0;0],d) | `cross([1;0;0],d)`<br>`u/norm(u)` | 赋值 | [W_RIS/wris_reference.m:17](<C:/Users/Spidey/Desktop/code/W_RIS/wris_reference.m:17>) |
| `v` | 光学局部 y 轴世界单位方向 cross(d,u) | `cross(d,u)` | 赋值 | [W_RIS/wris_reference.m:17](<C:/Users/Spidey/Desktop/code/W_RIS/wris_reference.m:17>) |


### 5.26 W_RIS/wris_report.m

**作用域：`wris_report`**

| 源码符号 | 含义与单位/维度 | 源码取值/计算式（各次赋值） | 角色 | 源码位置 |
| --- | --- | --- | --- | --- |
| `s` | 包含 15 个面板及验证结果的总结果结构 | 函数输入/输出；值随调用或上游结果确定 | 输入 | [W_RIS/wris_report.m:1](<C:/Users/Spidey/Desktop/code/W_RIS/wris_report.m:1>) |
| `c` | 当前 W_RIS 总配置结构 | `s.config` | 赋值 | [W_RIS/wris_report.m:3](<C:/Users/Spidey/Desktop/code/W_RIS/wris_report.m:3>) |
| `records` | 15 个面板快照数值记录的元胞数组 | `cell(15,1)`<br>`struct(...)（所有字段在本节逐行列明）` | 赋值 | [W_RIS/wris_report.m:3](<C:/Users/Spidey/Desktop/code/W_RIS/wris_report.m:3>) |
| `motionRows` | 12 个运动面板完整运动 CSV 数值块的元胞数组 | `cell(12,1)`<br>`[repmat([row col p.windMps p.windDirectionDeg p.site.center.'],n,1), m.t.',m.eta.',m.etaCenter.',p.tx.position.',p.tx.yaw.']` | 赋值 | [W_RIS/wris_report.m:3](<C:/Users/Spidey/Desktop/code/W_RIS/wris_report.m:3>) |
| `count` | 当前快照记录序号 1..15 | `0`<br>`count+1` | 赋值 | [W_RIS/wris_report.m:3](<C:/Users/Spidey/Desktop/code/W_RIS/wris_report.m:3>) |
| `nm` | 当前非静水运动记录序号 1..12 | `0`<br>`nm+1` | 赋值 | [W_RIS/wris_report.m:3](<C:/Users/Spidey/Desktop/code/W_RIS/wris_report.m:3>) |
| `row` | 工况表的行或列索引 | `1:3` | 循环索引 | [W_RIS/wris_report.m:4](<C:/Users/Spidey/Desktop/code/W_RIS/wris_report.m:4>) |
| `reference` | 同一行静水参考工况的投影主轴角（°） | `s.panels{row,1}.out.geometry.projectedOrientationDeg` | 赋值 | [W_RIS/wris_report.m:5](<C:/Users/Spidey/Desktop/code/W_RIS/wris_report.m:5>) |
| `col` | 工况表的行或列索引 | `1:5` | 循环索引 | [W_RIS/wris_report.m:6](<C:/Users/Spidey/Desktop/code/W_RIS/wris_report.m:6>) |
| `p` | 当前面板元数据结构 | `s.panels{row,col}` | 赋值 | [W_RIS/wris_report.m:7](<C:/Users/Spidey/Desktop/code/W_RIS/wris_report.m:7>) |
| `o` | 当前面板的 MC 输出结构 | `p.out` | 赋值 | [W_RIS/wris_report.m:7](<C:/Users/Spidey/Desktop/code/W_RIS/wris_report.m:7>) |
| `g` | 当前光学几何结构 | `o.geometry` | 赋值 | [W_RIS/wris_report.m:7](<C:/Users/Spidey/Desktop/code/W_RIS/wris_report.m:7>) |
| `rotation` | 当前几何主轴相对同一行静水参考的模 180° 旋转（°） | `mod(g.projectedOrientationDeg-reference+90,180)-90` | 赋值 | [W_RIS/wris_report.m:8](<C:/Users/Spidey/Desktop/code/W_RIS/wris_report.m:8>) |
| `records.row` | 主图面板行号/列号（从 1 开始） | `row` | 结构字段 | [W_RIS/wris_report.m:9](<C:/Users/Spidey/Desktop/code/W_RIS/wris_report.m:9>) |
| `records.column` | 主图面板行号/列号（从 1 开始） | `col` | 结构字段 | [W_RIS/wris_report.m:9](<C:/Users/Spidey/Desktop/code/W_RIS/wris_report.m:9>) |
| `records.initialIncidenceDeg` | 该行静水初始入射角（°） | `p.initialIncidenceDeg` | 结构字段 | [W_RIS/wris_report.m:9](<C:/Users/Spidey/Desktop/code/W_RIS/wris_report.m:9>) |
| `records.windMps` | 风速（m/s） | `p.windMps` | 结构字段 | [W_RIS/wris_report.m:10](<C:/Users/Spidey/Desktop/code/W_RIS/wris_report.m:10>) |
| `records.windDirectionDeg` | 波浪传播方向角（°），从 +Y 朝 +X 增大 | `p.windDirectionDeg` | 结构字段 | [W_RIS/wris_report.m:10](<C:/Users/Spidey/Desktop/code/W_RIS/wris_report.m:10>) |
| `records.phase` | 空间海浪站位标签；Calm 表示静水固定组 | `string(p.phase)` | 结构字段 | [W_RIS/wris_report.m:10](<C:/Users/Spidey/Desktop/code/W_RIS/wris_report.m:10>) |
| `records.timeS` | 相对采样时刻（s） | `g.time` | 结构字段 | [W_RIS/wris_report.m:11](<C:/Users/Spidey/Desktop/code/W_RIS/wris_report.m:11>) |
| `records.siteXM` | 所选初始站位的世界水平坐标（m） | `p.site.center(1)` | 结构字段 | [W_RIS/wris_report.m:11](<C:/Users/Spidey/Desktop/code/W_RIS/wris_report.m:11>) |
| `records.siteYM` | 所选初始站位的世界水平坐标（m） | `p.site.center(2)` | 结构字段 | [W_RIS/wris_report.m:11](<C:/Users/Spidey/Desktop/code/W_RIS/wris_report.m:11>) |
| `records.yawRisDeg` | 由 RRis 提取的实际 RIS yaw（°） | `g.yawRisDeg` | 结构字段 | [W_RIS/wris_report.m:12](<C:/Users/Spidey/Desktop/code/W_RIS/wris_report.m:12>) |
| `records.actualIncidenceDeg` | 运动后的实际光轴入射角（°） | `g.thetaDeg` | 结构字段 | [W_RIS/wris_report.m:12](<C:/Users/Spidey/Desktop/code/W_RIS/wris_report.m:12>) |
| `records.projectedOrientationDeg` | 光轴板面切向投影对应的几何主轴角（°，[-90,90)，法向入射时 NaN） | `g.projectedOrientationDeg` | 结构字段 | [W_RIS/wris_report.m:13](<C:/Users/Spidey/Desktop/code/W_RIS/wris_report.m:13>) |
| `records.orientationWithoutRisYawDeg` | 仅消去 RIS 世界 yaw、保留其余姿态和 TX 后的几何主轴角（°） | `g.orientationWithoutRisYawDeg` | 结构字段 | [W_RIS/wris_report.m:14](<C:/Users/Spidey/Desktop/code/W_RIS/wris_report.m:14>) |
| `records.completeCoreOrientationDeg` | 完整未散射核心的 MC 二阶矩主轴角（°） | `o.directPlaneMetrics.orientationDeg` | 结构字段 | [W_RIS/wris_report.m:15](<C:/Users/Spidey/Desktop/code/W_RIS/wris_report.m:15>) |
| `records.capturedCoreOrientationDeg` | 板内未散射核心的 MC 二阶矩主轴角（°） | `o.directMetrics.orientationDeg` | 结构字段 | [W_RIS/wris_report.m:16](<C:/Users/Spidey/Desktop/code/W_RIS/wris_report.m:16>) |
| `records.rotationVsStaticDeg` | 几何主轴相对同一行静水主轴的模 180° 旋转（°） | `rotation` | 结构字段 | [W_RIS/wris_report.m:16](<C:/Users/Spidey/Desktop/code/W_RIS/wris_report.m:16>) |
| `records.coreClippedFraction` | 1-板内核心功率/完整核心功率（无量纲） | `o.coreClippedFraction` | 结构字段 | [W_RIS/wris_report.m:17](<C:/Users/Spidey/Desktop/code/W_RIS/wris_report.m:17>) |
| `records.capturedWithoutEggW` | 板内不含/含 EGG 的接收功率（W），按字段区分 | `o.powerWithoutTurbulence` | 结构字段 | [W_RIS/wris_report.m:18](<C:/Users/Spidey/Desktop/code/W_RIS/wris_report.m:18>) |
| `records.capturedWithEggW` | 板内不含/含 EGG 的接收功率（W），按字段区分 | `o.capturedPower` | 结构字段 | [W_RIS/wris_report.m:18](<C:/Users/Spidey/Desktop/code/W_RIS/wris_report.m:18>) |
| `records.completeCorePowerWithoutEggW` | 完整/板内未散射核心功率（W，不含 EGG），按字段区分 | `o.directPlanePowerWithoutTurbulence` | 结构字段 | [W_RIS/wris_report.m:19](<C:/Users/Spidey/Desktop/code/W_RIS/wris_report.m:19>) |
| `records.capturedCorePowerWithoutEggW` | 完整/板内未散射核心功率（W，不含 EGG），按字段区分 | `o.directPowerWithoutTurbulence` | 结构字段 | [W_RIS/wris_report.m:20](<C:/Users/Spidey/Desktop/code/W_RIS/wris_report.m:20>) |
| `records.axisHitXM` | 光轴与 RIS 平面的局部 x/y 交点坐标（m） | `g.axisLocal(1)` | 结构字段 | [W_RIS/wris_report.m:21](<C:/Users/Spidey/Desktop/code/W_RIS/wris_report.m:21>) |
| `records.axisHitYM` | 光轴与 RIS 平面的局部 x/y 交点坐标（m） | `g.axisLocal(2)` | 结构字段 | [W_RIS/wris_report.m:21](<C:/Users/Spidey/Desktop/code/W_RIS/wris_report.m:21>) |
| `records.energyResidual` | 扣除轮盘赌权重变化后的相对能量守恒残差 | `o.energyResidual` | 结构字段 | [W_RIS/wris_report.m:21](<C:/Users/Spidey/Desktop/code/W_RIS/wris_report.m:21>) |
| `m` | 当前面板完整 RIS 运动记录 | `p.motion` | 赋值 | [W_RIS/wris_report.m:23](<C:/Users/Spidey/Desktop/code/W_RIS/wris_report.m:23>) |
| `n` | 当前运动记录时间点总数 | `numel(m.t)` | 赋值 | [W_RIS/wris_report.m:23](<C:/Users/Spidey/Desktop/code/W_RIS/wris_report.m:23>) |
| `summary` | 工况或运动统计汇总表 | `struct2table(vertcat(records{:}))` | 赋值 | [W_RIS/wris_report.m:29](<C:/Users/Spidey/Desktop/code/W_RIS/wris_report.m:29>) |
| `names` | motion.csv 的 18 个列名 | `{'row','column','windMps','windDirectionDeg','siteXM','siteYM','timeS', 'surgeM','swayM','heaveM','rollRad','pitchRad','yawRad','localSeaM','txXM','txYM','txZM','txYawRad'}` | 赋值 | [W_RIS/wris_report.m:31](<C:/Users/Spidey/Desktop/code/W_RIS/wris_report.m:31>) |
| `fid` | 以 UTF-8 写出 RESULTS.md 的文件句柄 | `fopen(fullfile(c.output,'RESULTS.md'),'w','n','UTF-8')` | 赋值 | [W_RIS/wris_report.m:34](<C:/Users/Spidey/Desktop/code/W_RIS/wris_report.m:34>) |
| `guard` | 函数结束时自动恢复随机数状态、关闭图窗或文件的清理对象 | `onCleanup(@() fclose(fid))` | 赋值 | [W_RIS/wris_report.m:34](<C:/Users/Spidey/Desktop/code/W_RIS/wris_report.m:34>) |
| `k` | 当前循环、分量或样本索引（无量纲，MATLAB 从 1 开始） | `1:height(summary)` | 循环索引 | [W_RIS/wris_report.m:43](<C:/Users/Spidey/Desktop/code/W_RIS/wris_report.m:43>) |
| `q` | 当前输出快照表中的一行记录；不是四元数 | `summary(k,:)` | 赋值 | [W_RIS/wris_report.m:44](<C:/Users/Spidey/Desktop/code/W_RIS/wris_report.m:44>) |


### 5.27 W_RIS/wris_run.m

**作用域：`wris_run`**

| 源码符号 | 含义与单位/维度 | 源码取值/计算式（各次赋值） | 角色 | 源码位置 |
| --- | --- | --- | --- | --- |
| `profile` | 运行档位：production（默认）或 smoke | `'production'` | 输入/赋值 | [W_RIS/wris_run.m:1](<C:/Users/Spidey/Desktop/code/W_RIS/wris_run.m:1>) |
| `s` | 保存实际配置、全部面板、验证结果与耗时的总输出结构 | `wris_simulate(c)` | 输出/赋值 | [W_RIS/wris_run.m:1](<C:/Users/Spidey/Desktop/code/W_RIS/wris_run.m:1>) |
| `c` | 当前 W_RIS 总配置结构 | `wris_config(profile)` | 赋值 | [W_RIS/wris_run.m:4](<C:/Users/Spidey/Desktop/code/W_RIS/wris_run.m:4>) |
| `tests` | tests/WrisTest.m 的 MATLAB 测试运行结果对象数组；当前源码快照未提供该文件 | `runtests(fullfile(c.root,'tests','WrisTest.m'))` | 赋值 | [W_RIS/wris_run.m:6](<C:/Users/Spidey/Desktop/code/W_RIS/wris_run.m:6>) |
| `s.validation` | 投影、能量、主轴与收敛验证结果结构 | `wris_validate(c,s)` | 赋值 | [W_RIS/wris_run.m:9](<C:/Users/Spidey/Desktop/code/W_RIS/wris_run.m:9>) |
| `s.elapsedSeconds` | 仿真与验证计算耗时（s），不含随后绘图/报告 | `toc` | 赋值 | [W_RIS/wris_run.m:10](<C:/Users/Spidey/Desktop/code/W_RIS/wris_run.m:10>) |


### 5.28 W_RIS/wris_setup.m

**作用域：`wris_setup`**

| 源码符号 | 含义与单位/维度 | 源码取值/计算式（各次赋值） | 角色 | 源码位置 |
| --- | --- | --- | --- | --- |
| `root` | 当前源码所在目录 | `fileparts(mfilename('fullpath'))` | 赋值 | [W_RIS/wris_setup.m:3](<C:/Users/Spidey/Desktop/code/W_RIS/wris_setup.m:3>) |


### 5.29 W_RIS/wris_simulate.m

**作用域：`wris_simulate`**

| 源码符号 | 含义与单位/维度 | 源码取值/计算式（各次赋值） | 角色 | 源码位置 |
| --- | --- | --- | --- | --- |
| `c` | 当前 W_RIS 总配置结构 | 函数输入/输出；值随调用或上游结果确定 | 输入 | [W_RIS/wris_simulate.m:1](<C:/Users/Spidey/Desktop/code/W_RIS/wris_simulate.m:1>) |
| `s` | 包含实际配置和 3×5 全部面板的仿真结果 | 函数输入/输出；值随调用或上游结果确定 | 输出 | [W_RIS/wris_simulate.m:1](<C:/Users/Spidey/Desktop/code/W_RIS/wris_simulate.m:1>) |
| `t` | 所有面板共同使用的相对快照时间 c.sampling.comparisonTime（s） | `c.sampling.comparisonTime` | 赋值 | [W_RIS/wris_simulate.m:3](<C:/Users/Spidey/Desktop/code/W_RIS/wris_simulate.m:3>) |
| `s.config` | 完整实际运行配置 | `c` | 赋值 | [W_RIS/wris_simulate.m:5](<C:/Users/Spidey/Desktop/code/W_RIS/wris_simulate.m:5>) |
| `s.panels` | 3×5 面板结果元胞数组 | `cell(3,5)`<br>`panel(cc,g,'Calm',0,0,[],[],struct('center',cc.ris.center))`<br>`panel(cc,g,site.label,cc.ris.windSpeed,cc.ris.windAngleDeg,m,tx,site)` | 赋值 | [W_RIS/wris_simulate.m:5](<C:/Users/Spidey/Desktop/code/W_RIS/wris_simulate.m:5>) |
| `row` | 3 个初始入射角的行索引/4 个运动工况的列索引 | `1:3` | 循环索引 | [W_RIS/wris_simulate.m:6](<C:/Users/Spidey/Desktop/code/W_RIS/wris_simulate.m:6>) |
| `cc` | 针对当前工况复制并覆盖后的 W_RIS 配置 | `c` | 赋值 | [W_RIS/wris_simulate.m:7](<C:/Users/Spidey/Desktop/code/W_RIS/wris_simulate.m:7>) |
| `cc.geometry.incidenceDeg` | 静水初始光轴相对 RIS 法线的入射角（°） | `c.comparison.incidenceDeg(row)` | 赋值 | [W_RIS/wris_simulate.m:7](<C:/Users/Spidey/Desktop/code/W_RIS/wris_simulate.m:7>) |
| `g` | 当前光学几何结构 | `wris_geometry(cc,struct('t',t),struct(),1,'static')`<br>`wris_geometry(cc,m,tx,k,'joint')` | 赋值 | [W_RIS/wris_simulate.m:8](<C:/Users/Spidey/Desktop/code/W_RIS/wris_simulate.m:8>) |
| `col` | 3 个初始入射角的行索引/4 个运动工况的列索引 | `1:4` | 循环索引 | [W_RIS/wris_simulate.m:11](<C:/Users/Spidey/Desktop/code/W_RIS/wris_simulate.m:11>) |
| `cc.ris.windSpeed` | 风速（m/s） | `c.comparison.conditions(col,1)` | 赋值 | [W_RIS/wris_simulate.m:12](<C:/Users/Spidey/Desktop/code/W_RIS/wris_simulate.m:12>) |
| `cc.ris.windAngleDeg` | 波浪传播方向角（°），从 +Y 朝 +X 增大 | `c.comparison.conditions(col,2)` | 赋值 | [W_RIS/wris_simulate.m:13](<C:/Users/Spidey/Desktop/code/W_RIS/wris_simulate.m:13>) |
| `wave` | 海浪分量结构，含 omega、kx、ky、coefficient；具体重名按下方作用域解释 | `RIS_wave_components(cc.ris,c.sea.quadrature)` | 赋值 | [W_RIS/wris_simulate.m:14](<C:/Users/Spidey/Desktop/code/W_RIS/wris_simulate.m:14>) |
| `defaults` | 由 get_default_sea_params 读取的基础海面配置 | `get_default_sea_params()` | 赋值 | [W_RIS/wris_simulate.m:15](<C:/Users/Spidey/Desktop/code/W_RIS/wris_simulate.m:15>) |
| `wave.coefficient` | 海浪复振幅 amplitude*exp(1i*phase)（m），参考世界原点与绝对 t=0 | `wave.coefficient*sqrt(c.sea.neumannConstant/defaults.neumannConstant)` | 赋值 | [W_RIS/wris_simulate.m:16](<C:/Users/Spidey/Desktop/code/W_RIS/wris_simulate.m:16>) |
| `sites` | 当前风速/风向同一时刻四类空间海浪站位 | `wris_sites(cc,wave)` | 赋值 | [W_RIS/wris_simulate.m:17](<C:/Users/Spidey/Desktop/code/W_RIS/wris_simulate.m:17>) |
| `site` | 由 comparison.phaseIndex 选中的当前行列空间站位 | `sites(c.comparison.phaseIndex(row,col))` | 赋值 | [W_RIS/wris_simulate.m:19](<C:/Users/Spidey/Desktop/code/W_RIS/wris_simulate.m:19>) |
| `cc.ris.center` | RIS 的水平参考中心 [x;y]（m），或光斑加权中心 [x,y]（m）；以父结构为准 | `site.center` | 赋值 | [W_RIS/wris_simulate.m:21](<C:/Users/Spidey/Desktop/code/W_RIS/wris_simulate.m:21>) |
| `cc.ris.wave` | 外部共享的海浪分量结构；[] 时按参数生成 | `wave` | 赋值 | [W_RIS/wris_simulate.m:21](<C:/Users/Spidey/Desktop/code/W_RIS/wris_simulate.m:21>) |
| `m` | 当前工况完整 RIS 动力学结果 | `wris_cached(cc,'motion',cc.ris,@() RIS_float(cc.ris))` | 赋值 | [W_RIS/wris_simulate.m:24](<C:/Users/Spidey/Desktop/code/W_RIS/wris_simulate.m:24>) |
| `k` | RIS 已保存时间中最接近共同快照时刻的样本索引，要求误差<1e-10 s | `min(abs(m.t-t))` | 赋值 | [W_RIS/wris_simulate.m:25](<C:/Users/Spidey/Desktop/code/W_RIS/wris_simulate.m:25>) |
| `tx` | 发射平台运动时序结构 | `wris_tx_motion(cc,m.t)` | 赋值 | [W_RIS/wris_simulate.m:27](<C:/Users/Spidey/Desktop/code/W_RIS/wris_simulate.m:27>) |


**作用域：`panel`**

| 源码符号 | 含义与单位/维度 | 源码取值/计算式（各次赋值） | 角色 | 源码位置 |
| --- | --- | --- | --- | --- |
| `c` | 当前 W_RIS 总配置结构 | 函数输入/输出；值随调用或上游结果确定 | 输入 | [W_RIS/wris_simulate.m:32](<C:/Users/Spidey/Desktop/code/W_RIS/wris_simulate.m:32>) |
| `g` | 当前光学几何结构 | 函数输入/输出；值随调用或上游结果确定 | 输入 | [W_RIS/wris_simulate.m:32](<C:/Users/Spidey/Desktop/code/W_RIS/wris_simulate.m:32>) |
| `phase` | 当前面板站位标签或 Calm | 函数输入/输出；值随调用或上游结果确定 | 输入 | [W_RIS/wris_simulate.m:32](<C:/Users/Spidey/Desktop/code/W_RIS/wris_simulate.m:32>) |
| `wind` | 当前工况风速（m/s）和传播风向（°） | 函数输入/输出；值随调用或上游结果确定 | 输入 | [W_RIS/wris_simulate.m:32](<C:/Users/Spidey/Desktop/code/W_RIS/wris_simulate.m:32>) |
| `direction` | 当前工况风速（m/s）和传播风向（°） | 函数输入/输出；值随调用或上游结果确定 | 输入 | [W_RIS/wris_simulate.m:32](<C:/Users/Spidey/Desktop/code/W_RIS/wris_simulate.m:32>) |
| `motion` | 当前 RIS 运动结果；静水固定面板为 [] | 函数输入/输出；值随调用或上游结果确定 | 输入 | [W_RIS/wris_simulate.m:32](<C:/Users/Spidey/Desktop/code/W_RIS/wris_simulate.m:32>) |
| `tx` | 发射平台运动时序结构 | 函数输入/输出；值随调用或上游结果确定 | 输入 | [W_RIS/wris_simulate.m:32](<C:/Users/Spidey/Desktop/code/W_RIS/wris_simulate.m:32>) |
| `site` | 当前空间海浪站位结构 | 函数输入/输出；值随调用或上游结果确定 | 输入 | [W_RIS/wris_simulate.m:32](<C:/Users/Spidey/Desktop/code/W_RIS/wris_simulate.m:32>) |
| `p` | 保存当前面板元数据、RIS/TX 轨迹及光学输出的结构 | `struct(...)（所有字段在本节逐行列明）` | 输出/赋值 | [W_RIS/wris_simulate.m:32](<C:/Users/Spidey/Desktop/code/W_RIS/wris_simulate.m:32>) |
| `gain` | 当前快照共享 EGG 增益；静水固定面板为 1 | `1`<br>`gains(k)` | 赋值 | [W_RIS/wris_simulate.m:34](<C:/Users/Spidey/Desktop/code/W_RIS/wris_simulate.m:34>) |
| `gains` | 与整个 RIS 时间序列一一对应的 EGG 独立增益向量 | `wris_egg(c.turbulence.kind,numel(motion.t),c.turbulence.seed)` | 赋值 | [W_RIS/wris_simulate.m:36](<C:/Users/Spidey/Desktop/code/W_RIS/wris_simulate.m:36>) |
| `k` | 与当前光学快照时间匹配的 RIS 时间索引 | `min(abs(motion.t-g.time))` | 赋值 | [W_RIS/wris_simulate.m:37](<C:/Users/Spidey/Desktop/code/W_RIS/wris_simulate.m:37>) |
| `key` | 缓存的物理输入标识结构 | `struct('source',c.source,'water',c.water,'mc',c.mc,'grid',c.grid, 'size',[c.ris.Lx c.ris.Ly],'geometry',g,'gain',gain)` | 赋值 | [W_RIS/wris_simulate.m:39](<C:/Users/Spidey/Desktop/code/W_RIS/wris_simulate.m:39>) |
| `key.source` | EGG 参数来源说明文本 | `c.source` | 结构字段 | [W_RIS/wris_simulate.m:39](<C:/Users/Spidey/Desktop/code/W_RIS/wris_simulate.m:39>) |
| `key.water` | 缓存键中对应水体/MC/评分网格配置 | `c.water` | 结构字段 | [W_RIS/wris_simulate.m:39](<C:/Users/Spidey/Desktop/code/W_RIS/wris_simulate.m:39>) |
| `key.mc` | 缓存键中对应水体/MC/评分网格配置 | `c.mc` | 结构字段 | [W_RIS/wris_simulate.m:39](<C:/Users/Spidey/Desktop/code/W_RIS/wris_simulate.m:39>) |
| `key.grid` | 缓存键中对应水体/MC/评分网格配置 | `c.grid` | 结构字段 | [W_RIS/wris_simulate.m:39](<C:/Users/Spidey/Desktop/code/W_RIS/wris_simulate.m:39>) |
| `key.size` | 缓存键中的 RIS 尺寸 [Lx,Ly]（m） | `[c.ris.Lx c.ris.Ly]` | 结构字段 | [W_RIS/wris_simulate.m:40](<C:/Users/Spidey/Desktop/code/W_RIS/wris_simulate.m:40>) |
| `key.geometry` | 该帧使用的光学几何结构 | `g` | 结构字段 | [W_RIS/wris_simulate.m:40](<C:/Users/Spidey/Desktop/code/W_RIS/wris_simulate.m:40>) |
| `key.gain` | 本帧所有光子共享的 EGG 增益（无量纲） | `gain` | 结构字段 | [W_RIS/wris_simulate.m:40](<C:/Users/Spidey/Desktop/code/W_RIS/wris_simulate.m:40>) |
| `out` | 当前函数输出结果结构 | `wris_cached(c,'optics',key,@() wris_mc(c,g,c.mc.seed,gain))` | 赋值 | [W_RIS/wris_simulate.m:41](<C:/Users/Spidey/Desktop/code/W_RIS/wris_simulate.m:41>) |
| `p.initialIncidenceDeg` | 该行静水初始入射角（°） | `c.geometry.incidenceDeg` | 结构字段 | [W_RIS/wris_simulate.m:42](<C:/Users/Spidey/Desktop/code/W_RIS/wris_simulate.m:42>) |
| `p.windMps` | 风速（m/s） | `wind` | 结构字段 | [W_RIS/wris_simulate.m:42](<C:/Users/Spidey/Desktop/code/W_RIS/wris_simulate.m:42>) |
| `p.windDirectionDeg` | 波浪传播方向角（°），从 +Y 朝 +X 增大 | `direction` | 结构字段 | [W_RIS/wris_simulate.m:43](<C:/Users/Spidey/Desktop/code/W_RIS/wris_simulate.m:43>) |
| `p.phase` | 空间海浪站位标签；Calm 表示静水固定组 | `phase` | 结构字段 | [W_RIS/wris_simulate.m:43](<C:/Users/Spidey/Desktop/code/W_RIS/wris_simulate.m:43>) |
| `p.site` | 初始空间站位结构（center/elevation/slope/label） | `site` | 结构字段 | [W_RIS/wris_simulate.m:43](<C:/Users/Spidey/Desktop/code/W_RIS/wris_simulate.m:43>) |
| `p.motion` | 当前工况 RIS 完整运动记录，静水面板为 [] | `motion` | 结构字段 | [W_RIS/wris_simulate.m:43](<C:/Users/Spidey/Desktop/code/W_RIS/wris_simulate.m:43>) |
| `p.tx` | 当前工况 TX 完整运动记录，静水面板为 [] | `tx` | 结构字段 | [W_RIS/wris_simulate.m:43](<C:/Users/Spidey/Desktop/code/W_RIS/wris_simulate.m:43>) |
| `p.out` | 当前工况 MC 结果结构 | `out` | 结构字段 | [W_RIS/wris_simulate.m:43](<C:/Users/Spidey/Desktop/code/W_RIS/wris_simulate.m:43>) |


### 5.30 W_RIS/wris_sites.m

**作用域：`wris_sites`**

| 源码符号 | 含义与单位/维度 | 源码取值/计算式（各次赋值） | 角色 | 源码位置 |
| --- | --- | --- | --- | --- |
| `c` | 当前 W_RIS 总配置结构 | 函数输入/输出；值随调用或上游结果确定 | 输入 | [W_RIS/wris_sites.m:1](<C:/Users/Spidey/Desktop/code/W_RIS/wris_sites.m:1>) |
| `wave` | 海浪分量结构，含 omega、kx、ky、coefficient；具体重名按下方作用域解释 | 函数输入/输出；值随调用或上游结果确定 | 输入 | [W_RIS/wris_sites.m:1](<C:/Users/Spidey/Desktop/code/W_RIS/wris_sites.m:1>) |
| `sites` | 同一海浪实现、同一绝对时刻的峰/谷/上升坡/下降坡四个空间站位 | `struct('label',labels,'center',cell(1,4),'elevation',cell(1,4),'slope',cell(1,4))` | 输出/赋值 | [W_RIS/wris_sites.m:1](<C:/Users/Spidey/Desktop/code/W_RIS/wris_sites.m:1>) |
| `direction` | [sind(windAngleDeg);cosd(windAngleDeg)] 世界水平主传播单位方向 | `[sind(c.ris.windAngleDeg);cosd(c.ris.windAngleDeg)]` | 赋值 | [W_RIS/wris_sites.m:4](<C:/Users/Spidey/Desktop/code/W_RIS/wris_sites.m:4>) |
| `omega` | 去重后的海浪角频率（rad/s） | `unique(wave.omega)` | 赋值 | [W_RIS/wris_sites.m:5](<C:/Users/Spidey/Desktop/code/W_RIS/wris_sites.m:5>) |
| `group` | 每个波浪分量在去重频率中的分组编号 | `unique(wave.omega)` | 赋值 | [W_RIS/wris_sites.m:5](<C:/Users/Spidey/Desktop/code/W_RIS/wris_sites.m:5>) |
| `energy` | 同一频率全部方向的复波幅平方和（m²） | `accumarray(group,abs(wave.coefficient).^2)` | 赋值 | [W_RIS/wris_sites.m:6](<C:/Users/Spidey/Desktop/code/W_RIS/wris_sites.m:6>) |
| `peak` | 离散波幅平方和最大的频率索引 | `max(energy)` | 赋值 | [W_RIS/wris_sites.m:7](<C:/Users/Spidey/Desktop/code/W_RIS/wris_sites.m:7>) |
| `wavelength` | 峰能频率的深水波长 2*pi*g/omega(peak)²（m） | `2*pi*c.ris.g/omega(peak)^2` | 赋值 | [W_RIS/wris_sites.m:8](<C:/Users/Spidey/Desktop/code/W_RIS/wris_sites.m:8>) |
| `distance` | 沿主传播方向、范围 ±wavelength 的剖面坐标（m） | `linspace(-wavelength,wavelength,c.sea.siteSamples).'` | 赋值 | [W_RIS/wris_sites.m:9](<C:/Users/Spidey/Desktop/code/W_RIS/wris_sites.m:9>) |
| `points` | 剖面采样点世界水平坐标（Nsite×2，m） | `c.ris.center.'+distance*direction.'` | 赋值 | [W_RIS/wris_sites.m:10](<C:/Users/Spidey/Desktop/code/W_RIS/wris_sites.m:10>) |
| `field` | 提供给海面重建函数的 wave/parameters 封装，或当前块复波场 | `struct('wave',wave,'parameters',c.ris)` | 赋值 | [W_RIS/wris_sites.m:11](<C:/Users/Spidey/Desktop/code/W_RIS/wris_sites.m:11>) |
| `field.wave` | 外部共享的海浪分量结构；[] 时按参数生成 | `wave` | 结构字段 | [W_RIS/wris_sites.m:11](<C:/Users/Spidey/Desktop/code/W_RIS/wris_sites.m:11>) |
| `field.parameters` | 当前计算实际使用的参数结构；EGG 中为混合分布参数 | `c.ris` | 结构字段 | [W_RIS/wris_sites.m:11](<C:/Users/Spidey/Desktop/code/W_RIS/wris_sites.m:11>) |
| `z` | 统一比较时刻的剖面海面高程（m） | `RIS_wave_field(field,points,c.sampling.comparisonTime)` | 赋值 | [W_RIS/wris_sites.m:12](<C:/Users/Spidey/Desktop/code/W_RIS/wris_sites.m:12>) |
| `gradient` | 海面 x/y 方向解析斜率（无量纲） | `RIS_wave_field(field,points,c.sampling.comparisonTime)` | 赋值 | [W_RIS/wris_sites.m:12](<C:/Users/Spidey/Desktop/code/W_RIS/wris_sites.m:12>) |
| `slope` | 海面梯度沿主传播方向的投影（无量纲） | `gradient*direction` | 赋值 | [W_RIS/wris_sites.m:13](<C:/Users/Spidey/Desktop/code/W_RIS/wris_sites.m:13>) |
| `crest` | 剖面局部波峰/波谷索引集合 | `find(z(2:end-1)>=z(1:end-2) & z(2:end-1)>z(3:end))+1` | 赋值 | [W_RIS/wris_sites.m:14](<C:/Users/Spidey/Desktop/code/W_RIS/wris_sites.m:14>) |
| `trough` | 剖面局部波峰/波谷索引集合 | `find(z(2:end-1)<=z(1:end-2) & z(2:end-1)<z(3:end))+1` | 赋值 | [W_RIS/wris_sites.m:15](<C:/Users/Spidey/Desktop/code/W_RIS/wris_sites.m:15>) |
| `rising` | 剖面空间上升/下降过零的相邻区间左端索引；不是时间变化方向 | `find(z(1:end-1)<=0 & z(2:end)>0)` | 赋值 | [W_RIS/wris_sites.m:16](<C:/Users/Spidey/Desktop/code/W_RIS/wris_sites.m:16>) |
| `falling` | 剖面空间上升/下降过零的相邻区间左端索引；不是时间变化方向 | `find(z(1:end-1)>=0 & z(2:end)<0)` | 赋值 | [W_RIS/wris_sites.m:17](<C:/Users/Spidey/Desktop/code/W_RIS/wris_sites.m:17>) |
| `a` | 最高局部峰/最低局部谷在候选集合中的位置 | `max(z(crest))` | 赋值 | [W_RIS/wris_sites.m:20](<C:/Users/Spidey/Desktop/code/W_RIS/wris_sites.m:20>) |
| `b` | 最高局部峰/最低局部谷在候选集合中的位置 | `min(z(trough))` | 赋值 | [W_RIS/wris_sites.m:20](<C:/Users/Spidey/Desktop/code/W_RIS/wris_sites.m:20>) |
| `d` | 上升坡最大/下降坡最小方向斜率在候选集合中的位置 | `max(slope(rising))` | 赋值 | [W_RIS/wris_sites.m:21](<C:/Users/Spidey/Desktop/code/W_RIS/wris_sites.m:21>) |
| `e` | 上升坡最大/下降坡最小方向斜率在候选集合中的位置 | `min(slope(falling))` | 赋值 | [W_RIS/wris_sites.m:21](<C:/Users/Spidey/Desktop/code/W_RIS/wris_sites.m:21>) |
| `ids` | 最终四类站位对应的剖面样本索引 | `[crest(a),trough(b),rising(d),falling(e)]` | 赋值 | [W_RIS/wris_sites.m:22](<C:/Users/Spidey/Desktop/code/W_RIS/wris_sites.m:22>) |
| `labels` | 工况、相位或图形标签文本 | `{'Crest site','Trough site','Rising-slope site','Falling-slope site'}` | 赋值 | [W_RIS/wris_sites.m:23](<C:/Users/Spidey/Desktop/code/W_RIS/wris_sites.m:23>) |
| `sites.label` | 空间海浪站位名称 | `labels` | 结构字段 | [W_RIS/wris_sites.m:24](<C:/Users/Spidey/Desktop/code/W_RIS/wris_sites.m:24>) |
| `sites.center` | RIS 的水平参考中心 [x;y]（m），或光斑加权中心 [x,y]（m）；以父结构为准 | `cell(1,4)`<br>`c.ris.center+coordinate*direction` | 结构字段/赋值 | [W_RIS/wris_sites.m:24](<C:/Users/Spidey/Desktop/code/W_RIS/wris_sites.m:24>) |
| `sites.elevation` | 选定初始站位在统一时刻的参考海面高度（m） | `cell(1,4)`<br>`RIS_wave_field(field,sites(j).center.',c.sampling.comparisonTime)` | 结构字段/赋值 | [W_RIS/wris_sites.m:24](<C:/Users/Spidey/Desktop/code/W_RIS/wris_sites.m:24>) |
| `sites.slope` | 沿主传播方向的海面空间斜率（无量纲） | `cell(1,4)`<br>`g*direction` | 结构字段/赋值 | [W_RIS/wris_sites.m:24](<C:/Users/Spidey/Desktop/code/W_RIS/wris_sites.m:24>) |
| `j` | 当前循环、分量或样本索引（无量纲，MATLAB 从 1 开始） | `1:4` | 循环索引 | [W_RIS/wris_sites.m:25](<C:/Users/Spidey/Desktop/code/W_RIS/wris_sites.m:25>) |
| `k` | 当前站位在剖面中的样本索引 | `ids(j)` | 赋值 | [W_RIS/wris_sites.m:26](<C:/Users/Spidey/Desktop/code/W_RIS/wris_sites.m:26>) |
| `coordinate` | 当前选定站位沿传播方向的标量坐标（m）；坡面以 fzero 细化 | `distance(k)`<br>`fzero(@(q) RIS_wave_field(field,(c.ris.center+q*direction).', c.sampling.comparisonTime),distance(k:k+1))` | 赋值 | [W_RIS/wris_sites.m:26](<C:/Users/Spidey/Desktop/code/W_RIS/wris_sites.m:26>) |
| `q` | fzero 匿名函数的沿传播方向位置自变量（m） | `coordinate=fzero(@(q) RIS_wave_field(field,(c.ris.center+q*direction).', c.sampling.comparisonTime),distance(k:k+1))` | 匿名函数输入 | [W_RIS/wris_sites.m:28](<C:/Users/Spidey/Desktop/code/W_RIS/wris_sites.m:28>) |
| `g` | 所选站位的海面 [dx,dy] 梯度；此处不是光学几何结构 | `RIS_wave_field(field,sites(j).center.',c.sampling.comparisonTime)` | 赋值 | [W_RIS/wris_sites.m:32](<C:/Users/Spidey/Desktop/code/W_RIS/wris_sites.m:32>) |


### 5.31 W_RIS/wris_tx_motion.m

**作用域：`wris_tx_motion`**

| 源码符号 | 含义与单位/维度 | 源码取值/计算式（各次赋值） | 角色 | 源码位置 |
| --- | --- | --- | --- | --- |
| `c` | 当前 W_RIS 总配置结构 | 函数输入/输出；值随调用或上游结果确定 | 输入 | [W_RIS/wris_tx_motion.m:1](<C:/Users/Spidey/Desktop/code/W_RIS/wris_tx_motion.m:1>) |
| `t` | 与 RIS 相同的相对运动时间采样向量（s） | 函数输入/输出；值随调用或上游结果确定 | 输入 | [W_RIS/wris_tx_motion.m:1](<C:/Users/Spidey/Desktop/code/W_RIS/wris_tx_motion.m:1>) |
| `tx` | 发射平台运动时序结构 | `struct(...)（所有字段在本节逐行列明）` | 输出/赋值 | [W_RIS/wris_tx_motion.m:1](<C:/Users/Spidey/Desktop/code/W_RIS/wris_tx_motion.m:1>) |
| `ref` | 静水参考光学几何 | `wris_reference(c)` | 赋值 | [W_RIS/wris_tx_motion.m:5](<C:/Users/Spidey/Desktop/code/W_RIS/wris_tx_motion.m:5>) |
| `n` | 运动输出时间点总数 | `numel(t)` | 赋值 | [W_RIS/wris_tx_motion.m:5](<C:/Users/Spidey/Desktop/code/W_RIS/wris_tx_motion.m:5>) |
| `stream` | 局部 Mersenne Twister 随机数流 | `RandStream('mt19937ar','Seed',c.tx.seed)` | 赋值 | [W_RIS/wris_tx_motion.m:6](<C:/Users/Spidey/Desktop/code/W_RIS/wris_tx_motion.m:6>) |
| `uvr` | OU 体速度输入 [u;v;r] 时序（m/s,m/s,rad/s） | `zeros(3,n)`<br>`c.tx.meanUVR+c.tx.stdUVR.*randn(stream,3,1)`<br>`c.tx.meanUVR+q.*(uvr(:,k-1)-c.tx.meanUVR)+ c.tx.stdUVR.*sqrt(1-q.^2).*randn(stream,3,1)` | 赋值 | [W_RIS/wris_tx_motion.m:7](<C:/Users/Spidey/Desktop/code/W_RIS/wris_tx_motion.m:7>) |
| `p` | 恒深 TX 世界位置时序（3×Nt，m）；此处不是 RIS 参数结构 | `repmat(ref.pTx,1,n)`<br>`p(1:2,k-1)+dt*[u*cos(psi)-v*sin(psi);u*sin(psi)+v*cos(psi)]` | 赋值 | [W_RIS/wris_tx_motion.m:8](<C:/Users/Spidey/Desktop/code/W_RIS/wris_tx_motion.m:8>) |
| `yaw` | TX 平台航向时序（rad），从 initialYaw 递推 | `zeros(1,n)`<br>`c.tx.initialYaw`<br>`psi+dt*uvr(3,k-1)` | 赋值 | [W_RIS/wris_tx_motion.m:8](<C:/Users/Spidey/Desktop/code/W_RIS/wris_tx_motion.m:8>) |
| `k` | 当前循环、分量或样本索引（无量纲，MATLAB 从 1 开始） | `2:n`<br>`1:n` | 循环索引 | [W_RIS/wris_tx_motion.m:9](<C:/Users/Spidey/Desktop/code/W_RIS/wris_tx_motion.m:9>) |
| `dt` | 相邻运动采样时间差（s） | `t(k)-t(k-1)` | 赋值 | [W_RIS/wris_tx_motion.m:10](<C:/Users/Spidey/Desktop/code/W_RIS/wris_tx_motion.m:10>) |
| `q` | 三个 OU 分量的离散相关系数 exp(-dt/correlationTime)，不是四元数 | `exp(-dt./c.tx.correlationTime)` | 赋值 | [W_RIS/wris_tx_motion.m:11](<C:/Users/Spidey/Desktop/code/W_RIS/wris_tx_motion.m:11>) |
| `psi` | 前一个样本的 TX 平台航向（rad） | `yaw(k-1)` | 赋值 | [W_RIS/wris_tx_motion.m:14](<C:/Users/Spidey/Desktop/code/W_RIS/wris_tx_motion.m:14>) |
| `u` | 前一个样本的 TX 体坐标纵向/横向速度（m/s） | `uvr(1,k-1)` | 赋值 | [W_RIS/wris_tx_motion.m:14](<C:/Users/Spidey/Desktop/code/W_RIS/wris_tx_motion.m:14>) |
| `v` | 前一个样本的 TX 体坐标纵向/横向速度（m/s） | `uvr(2,k-1)` | 赋值 | [W_RIS/wris_tx_motion.m:14](<C:/Users/Spidey/Desktop/code/W_RIS/wris_tx_motion.m:14>) |
| `rotation` | 光学坐标系到世界坐标系旋转时序（3×3×Nt） | `zeros(3,3,n)`<br>`R*ref.RTx` | 赋值 | [W_RIS/wris_tx_motion.m:18](<C:/Users/Spidey/Desktop/code/W_RIS/wris_tx_motion.m:18>) |
| `d` | 当前 TX 光轴的世界单位方向时序（3×Nt） | `zeros(3,n)`<br>`rotation(:,3,k)` | 赋值 | [W_RIS/wris_tx_motion.m:18](<C:/Users/Spidey/Desktop/code/W_RIS/wris_tx_motion.m:18>) |
| `R` | 相对初始平台 yaw 的世界竖直轴旋转矩阵 | `RIS_rotation([0;0;yaw(k)-yaw(1)])` | 赋值 | [W_RIS/wris_tx_motion.m:21](<C:/Users/Spidey/Desktop/code/W_RIS/wris_tx_motion.m:21>) |
| `tx.t` | 相对时间采样向量（s） | `t` | 结构字段 | [W_RIS/wris_tx_motion.m:24](<C:/Users/Spidey/Desktop/code/W_RIS/wris_tx_motion.m:24>) |
| `tx.position` | TX 位置世界坐标时序（m），z 恒定 | `p` | 结构字段 | [W_RIS/wris_tx_motion.m:24](<C:/Users/Spidey/Desktop/code/W_RIS/wris_tx_motion.m:24>) |
| `tx.yaw` | 右手 ZYX 横摇、纵摇、偏航（rad）；yaw 时序已 unwrap | `yaw` | 结构字段 | [W_RIS/wris_tx_motion.m:24](<C:/Users/Spidey/Desktop/code/W_RIS/wris_tx_motion.m:24>) |
| `tx.uvr` | TX 体坐标速度 [u;v;r] 时序（m/s,m/s,rad/s） | `uvr` | 结构字段 | [W_RIS/wris_tx_motion.m:24](<C:/Users/Spidey/Desktop/code/W_RIS/wris_tx_motion.m:24>) |
| `tx.rotation` | TX 光学系到世界系旋转时序 | `rotation` | 结构字段 | [W_RIS/wris_tx_motion.m:24](<C:/Users/Spidey/Desktop/code/W_RIS/wris_tx_motion.m:24>) |
| `tx.direction` | TX 世界光轴单位方向时序 | `d` | 结构字段 | [W_RIS/wris_tx_motion.m:24](<C:/Users/Spidey/Desktop/code/W_RIS/wris_tx_motion.m:24>) |
| `tx.model` | TX 运动模型说明文本 | `'111.pdf Eqs (5)-(7), constant depth, assumed OU velocity inputs'` | 结构字段 | [W_RIS/wris_tx_motion.m:25](<C:/Users/Spidey/Desktop/code/W_RIS/wris_tx_motion.m:25>) |


### 5.32 W_RIS/wris_validate.m

**作用域：`wris_validate`**

| 源码符号 | 含义与单位/维度 | 源码取值/计算式（各次赋值） | 角色 | 源码位置 |
| --- | --- | --- | --- | --- |
| `c` | 当前 W_RIS 总配置结构 | 函数输入/输出；值随调用或上游结果确定 | 输入 | [W_RIS/wris_validate.m:1](<C:/Users/Spidey/Desktop/code/W_RIS/wris_validate.m:1>) |
| `s` | 已完成的 15 面板仿真结果 | 函数输入/输出；值随调用或上游结果确定 | 输入 | [W_RIS/wris_validate.m:1](<C:/Users/Spidey/Desktop/code/W_RIS/wris_validate.m:1>) |
| `v` | 全部物理与统计验证结果结构 | 函数输入/输出；值随调用或上游结果确定 | 输出 | [W_RIS/wris_validate.m:1](<C:/Users/Spidey/Desktop/code/W_RIS/wris_validate.m:1>) |
| `ns` | 收敛检查使用的光子数档位向量 | `c.validation.photons` | 赋值 | [W_RIS/wris_validate.m:3](<C:/Users/Spidey/Desktop/code/W_RIS/wris_validate.m:3>) |
| `reps` | 每个光子数档位的独立重复次数 | `c.validation.repeats` | 赋值 | [W_RIS/wris_validate.m:3](<C:/Users/Spidey/Desktop/code/W_RIS/wris_validate.m:3>) |
| `rows` | 每档的 [光子数,平均功率,重复标准差,平均报告标准误] | `zeros(numel(ns),4)`<br>`[ns(j),mean(powers(j,:)),std(powers(j,:)),mean(reported)]` | 赋值 | [W_RIS/wris_validate.m:4](<C:/Users/Spidey/Desktop/code/W_RIS/wris_validate.m:4>) |
| `powers` | 各光子数档位各重复的不含 EGG 板内功率（W） | `zeros(numel(ns),reps)`<br>`out.powerWithoutTurbulence` | 赋值 | [W_RIS/wris_validate.m:4](<C:/Users/Spidey/Desktop/code/W_RIS/wris_validate.m:4>) |
| `g` | 收敛测试选取面板 {1,2} 的几何，或透明测试的各行静水几何 | `s.panels{1,2}.out.geometry`<br>`s.panels{row,1}.out.geometry` | 赋值 | [W_RIS/wris_validate.m:5](<C:/Users/Spidey/Desktop/code/W_RIS/wris_validate.m:5>) |
| `j` | 当前循环、分量或样本索引（无量纲，MATLAB 从 1 开始） | `1:numel(ns)` | 循环索引 | [W_RIS/wris_validate.m:6](<C:/Users/Spidey/Desktop/code/W_RIS/wris_validate.m:6>) |
| `reported` | 当前光子数档位每次重复的 MC 功率标准误（W） | `zeros(1,reps)`<br>`out.powerSE` | 赋值 | [W_RIS/wris_validate.m:8](<C:/Users/Spidey/Desktop/code/W_RIS/wris_validate.m:8>) |
| `k` | 当前循环、分量或样本索引（无量纲，MATLAB 从 1 开始） | `1:reps` | 循环索引 | [W_RIS/wris_validate.m:9](<C:/Users/Spidey/Desktop/code/W_RIS/wris_validate.m:9>) |
| `cc` | 针对当前工况复制并覆盖后的 W_RIS 配置 | `c` | 赋值 | [W_RIS/wris_validate.m:10](<C:/Users/Spidey/Desktop/code/W_RIS/wris_validate.m:10>) |
| `cc.mc.photons` | 单帧发射光子包总数，或验证采用的光子数档位 | `ns(j)`<br>`max(ns)` | 赋值 | [W_RIS/wris_validate.m:10](<C:/Users/Spidey/Desktop/code/W_RIS/wris_validate.m:10>) |
| `seed` | 独立重复种子 validation.seed+100*j+k | `c.validation.seed+100*j+k` | 赋值 | [W_RIS/wris_validate.m:10](<C:/Users/Spidey/Desktop/code/W_RIS/wris_validate.m:10>) |
| `out` | 当前函数输出结果结构 | `optical(cc,g,seed)` | 赋值 | [W_RIS/wris_validate.m:11](<C:/Users/Spidey/Desktop/code/W_RIS/wris_validate.m:11>) |
| `v.convergence` | 光子数、均值功率、重复标准差与报告标准误的收敛表 | `array2table(rows,'VariableNames',{'photons','meanPowerW','powerStdW','reportedPowerSEW'})` | 赋值 | [W_RIS/wris_validate.m:16](<C:/Users/Spidey/Desktop/code/W_RIS/wris_validate.m:16>) |
| `v.powers` | 各光子数档位各独立重复的板内无 EGG 功率矩阵（W） | `powers` | 赋值 | [W_RIS/wris_validate.m:17](<C:/Users/Spidey/Desktop/code/W_RIS/wris_validate.m:17>) |
| `v.convergenceAgreement` | 最高两档功率均值差是否小于 5 倍合成标准误 | `abs(rows(end,2)-rows(end-1,2))< 5*sqrt((rows(end,3)^2+rows(end-1,3)^2)/reps)` | 赋值 | [W_RIS/wris_validate.m:18](<C:/Users/Spidey/Desktop/code/W_RIS/wris_validate.m:18>) |
| `v.transparent` | 三种入射角的透明水体 MC 结果 | `cell(3,1)`<br>`optical(cc,g,901+row)` | 赋值 | [W_RIS/wris_validate.m:20](<C:/Users/Spidey/Desktop/code/W_RIS/wris_validate.m:20>) |
| `ratioError` | 透明水体光斑 MC 轴比与理论 sec(alpha) 的绝对误差 | `zeros(3,1)`<br>`abs(v.transparent{row}.metrics.ratio-1/cosd(c.comparison.incidenceDeg(row)))` | 赋值 | [W_RIS/wris_validate.m:20](<C:/Users/Spidey/Desktop/code/W_RIS/wris_validate.m:20>) |
| `row` | 工况表的行或列索引 | `1:3` | 循环索引 | [W_RIS/wris_validate.m:21](<C:/Users/Spidey/Desktop/code/W_RIS/wris_validate.m:21>) |
| `cc.water.absorption` | 水体吸收系数（m⁻¹），Table_abc 返回 a=c-b | `0` | 赋值 | [W_RIS/wris_validate.m:22](<C:/Users/Spidey/Desktop/code/W_RIS/wris_validate.m:22>) |
| `cc.water.scattering` | 水体散射系数（m⁻¹），Table_abc 返回 b | `0` | 赋值 | [W_RIS/wris_validate.m:22](<C:/Users/Spidey/Desktop/code/W_RIS/wris_validate.m:22>) |
| `frames` | 15 个面板 MC 输出的元胞数组 | `cellfun(@(p) p.out,s.panels,'UniformOutput',false)` | 赋值 | [W_RIS/wris_validate.m:27](<C:/Users/Spidey/Desktop/code/W_RIS/wris_validate.m:27>) |
| `p` | cellfun 匿名函数中的面板结构/MC 输出结构 | `frames=cellfun(@(p) p.out,s.panels,'UniformOutput',false)` | 匿名函数输入 | [W_RIS/wris_validate.m:27](<C:/Users/Spidey/Desktop/code/W_RIS/wris_validate.m:27>) |
| `v.maxEnergyResidual` | 全部 15 面板最大绝对相对能量残差 | `max(cellfun(@(o) abs(o.energyResidual),frames),[],'all')` | 赋值 | [W_RIS/wris_validate.m:28](<C:/Users/Spidey/Desktop/code/W_RIS/wris_validate.m:28>) |
| `o` | cellfun 匿名函数中的面板结构/MC 输出结构 | `v.maxEnergyResidual=max(cellfun(@(o) abs(o.energyResidual),frames),[],'all')`<br>`v.maxPowerIntegralErrorW=max(cellfun(@(o) abs(sum(o.noTurbulence,'all')*o.cellArea-o.powerWithoutTurbulence),frames),[],'all')`<br>`v.orientationErrorDeg=cellfun(@(o) abs(mod(o.directPlaneMetrics.orientationDeg- o.geometry.projectedOrientationDeg+90,180)-90),frames)`<br>`v.orientationToleranceDeg=cellfun(@(o) 0.25+rad2deg(4/(sqrt(max(1, o.directPlanePowerWithoutTurbulence*o.photons/c.source.power))* max(o.directPlaneMetrics.ratio-1/o.directPlaneMetrics.ratio,0.02))),frames)` | 匿名函数输入 | [W_RIS/wris_validate.m:28](<C:/Users/Spidey/Desktop/code/W_RIS/wris_validate.m:28>) |
| `v.maxPowerIntegralErrorW` | 网格辐照度面积积分与板内功率的最大绝对差（W） | `max(cellfun(@(o) abs(sum(o.noTurbulence,'all')*o.cellArea-o.powerWithoutTurbulence),frames),[],'all')` | 赋值 | [W_RIS/wris_validate.m:29](<C:/Users/Spidey/Desktop/code/W_RIS/wris_validate.m:29>) |
| `v.orientationErrorDeg` | 完整核心 MC 主轴与几何主轴的模 180° 绝对差（°） | `cellfun(@(o) abs(mod(o.directPlaneMetrics.orientationDeg- o.geometry.projectedOrientationDeg+90,180)-90),frames)` | 赋值 | [W_RIS/wris_validate.m:31](<C:/Users/Spidey/Desktop/code/W_RIS/wris_validate.m:31>) |
| `v.orientationToleranceDeg` | 由有效核心光子数和轴比计算的 4 标准误加 0.25° 容差 | `cellfun(@(o) 0.25+rad2deg(4/(sqrt(max(1, o.directPlanePowerWithoutTurbulence*o.photons/c.source.power))* max(o.directPlaneMetrics.ratio-1/o.directPlaneMetrics.ratio,0.02))),frames)` | 赋值 | [W_RIS/wris_validate.m:34](<C:/Users/Spidey/Desktop/code/W_RIS/wris_validate.m:34>) |
| `v.maxProjectionRatioError` | 透明水体 MC 轴比与 sec(初始入射角) 的最大绝对误差 | `max(ratioError)` | 赋值 | [W_RIS/wris_validate.m:37](<C:/Users/Spidey/Desktop/code/W_RIS/wris_validate.m:37>) |
| `v.physicalAssertionsPassed` | 所有物理与统计验证判据是否同时通过 | `all(ratioError<0.05) && v.maxEnergyResidual<1e-10 && v.maxPowerIntegralErrorW<1e-10 && all(v.orientationErrorDeg<v.orientationToleranceDeg,'all') && v.convergenceAgreement` | 赋值 | [W_RIS/wris_validate.m:38](<C:/Users/Spidey/Desktop/code/W_RIS/wris_validate.m:38>) |


**作用域：`optical`**

| 源码符号 | 含义与单位/维度 | 源码取值/计算式（各次赋值） | 角色 | 源码位置 |
| --- | --- | --- | --- | --- |
| `c` | 当前 W_RIS 总配置结构 | 函数输入/输出；值随调用或上游结果确定 | 输入 | [W_RIS/wris_validate.m:44](<C:/Users/Spidey/Desktop/code/W_RIS/wris_validate.m:44>) |
| `g` | 当前光学几何结构 | 函数输入/输出；值随调用或上游结果确定 | 输入 | [W_RIS/wris_validate.m:44](<C:/Users/Spidey/Desktop/code/W_RIS/wris_validate.m:44>) |
| `seed` | 随机数种子（非负整数） | 函数输入/输出；值随调用或上游结果确定 | 输入 | [W_RIS/wris_validate.m:44](<C:/Users/Spidey/Desktop/code/W_RIS/wris_validate.m:44>) |
| `out` | 固定几何和随机种子、gain=1 的缓存 MC 验证结果 | `wris_cached(c,'optics',key,@() wris_mc(c,g,seed,1))` | 输出/赋值 | [W_RIS/wris_validate.m:44](<C:/Users/Spidey/Desktop/code/W_RIS/wris_validate.m:44>) |
| `key` | 缓存的物理输入标识结构 | `struct('source',c.source,'water',c.water,'mc',c.mc,'grid',c.grid, 'size',[c.ris.Lx c.ris.Ly],'geometry',g,'seed',seed,'gain',1)` | 赋值 | [W_RIS/wris_validate.m:46](<C:/Users/Spidey/Desktop/code/W_RIS/wris_validate.m:46>) |
| `key.source` | EGG 参数来源说明文本 | `c.source` | 结构字段 | [W_RIS/wris_validate.m:46](<C:/Users/Spidey/Desktop/code/W_RIS/wris_validate.m:46>) |
| `key.water` | 缓存键中对应水体/MC/评分网格配置 | `c.water` | 结构字段 | [W_RIS/wris_validate.m:46](<C:/Users/Spidey/Desktop/code/W_RIS/wris_validate.m:46>) |
| `key.mc` | 缓存键中对应水体/MC/评分网格配置 | `c.mc` | 结构字段 | [W_RIS/wris_validate.m:46](<C:/Users/Spidey/Desktop/code/W_RIS/wris_validate.m:46>) |
| `key.grid` | 缓存键中对应水体/MC/评分网格配置 | `c.grid` | 结构字段 | [W_RIS/wris_validate.m:46](<C:/Users/Spidey/Desktop/code/W_RIS/wris_validate.m:46>) |
| `key.size` | 缓存键中的 RIS 尺寸 [Lx,Ly]（m） | `[c.ris.Lx c.ris.Ly]` | 结构字段 | [W_RIS/wris_validate.m:47](<C:/Users/Spidey/Desktop/code/W_RIS/wris_validate.m:47>) |
| `key.geometry` | 该帧使用的光学几何结构 | `g` | 结构字段 | [W_RIS/wris_validate.m:47](<C:/Users/Spidey/Desktop/code/W_RIS/wris_validate.m:47>) |
| `key.seed` | 随机数种子；所属结构决定是海浪、TX、MC 或湍流随机流 | `seed` | 结构字段 | [W_RIS/wris_validate.m:47](<C:/Users/Spidey/Desktop/code/W_RIS/wris_validate.m:47>) |
| `key.gain` | 本帧所有光子共享的 EGG 增益（无量纲） | `1` | 结构字段 | [W_RIS/wris_validate.m:47](<C:/Users/Spidey/Desktop/code/W_RIS/wris_validate.m:47>) |


### 5.33 fig3-7/Table_abc.m

**作用域：`Table_abc`**

| 源码符号 | 含义与单位/维度 | 源码取值/计算式（各次赋值） | 角色 | 源码位置 |
| --- | --- | --- | --- | --- |
| `water_type` | 水体档位字符串，当前 W_RIS 默认 'c0.151' | 函数输入/输出；值随调用或上游结果确定 | 输入 | [fig3-7/Table_abc.m:1](<C:/Users/Spidey/Desktop/code/fig3-7/Table_abc.m:1>) |
| `a` | 吸收系数 a=c-b（m⁻¹） | `c - b` | 输出/赋值 | [fig3-7/Table_abc.m:1](<C:/Users/Spidey/Desktop/code/fig3-7/Table_abc.m:1>) |
| `b` | 散射系数（m⁻¹） | `paras(1)` | 输出/赋值 | [fig3-7/Table_abc.m:1](<C:/Users/Spidey/Desktop/code/fig3-7/Table_abc.m:1>) |
| `c` | 总衰减系数 a+b（m⁻¹）；此处不是配置结构 | `paras(2)` | 输出/赋值 | [fig3-7/Table_abc.m:1](<C:/Users/Spidey/Desktop/code/fig3-7/Table_abc.m:1>) |
| `c3_2` | 对应水体档位的 [散射系数 b,总衰减系数 c]（m⁻¹）；变量名小数点用下划线表示 | `[2.6, 3.2]` | 赋值 | [fig3-7/Table_abc.m:3](<C:/Users/Spidey/Desktop/code/fig3-7/Table_abc.m:3>) |
| `c2_675` | 对应水体档位的 [散射系数 b,总衰减系数 c]（m⁻¹）；变量名小数点用下划线表示 | `[2.1734 , 2.675]` | 赋值 | [fig3-7/Table_abc.m:4](<C:/Users/Spidey/Desktop/code/fig3-7/Table_abc.m:4>) |
| `c1_92` | 对应水体档位的 [散射系数 b,总衰减系数 c]（m⁻¹）；变量名小数点用下划线表示 | `[1.56, 1.92]` | 赋值 | [fig3-7/Table_abc.m:5](<C:/Users/Spidey/Desktop/code/fig3-7/Table_abc.m:5>) |
| `c1_475` | 对应水体档位的 [散射系数 b,总衰减系数 c]（m⁻¹）；变量名小数点用下划线表示 | `[1.1984, 1.475]` | 赋值 | [fig3-7/Table_abc.m:6](<C:/Users/Spidey/Desktop/code/fig3-7/Table_abc.m:6>) |
| `c1_0` | 对应水体档位的 [散射系数 b,总衰减系数 c]（m⁻¹）；变量名小数点用下划线表示 | `[0.8125, 1]` | 赋值 | [fig3-7/Table_abc.m:7](<C:/Users/Spidey/Desktop/code/fig3-7/Table_abc.m:7>) |
| `c0_5` | 对应水体档位的 [散射系数 b,总衰减系数 c]（m⁻¹）；变量名小数点用下划线表示 | `[0.406, 0.5]` | 赋值 | [fig3-7/Table_abc.m:8](<C:/Users/Spidey/Desktop/code/fig3-7/Table_abc.m:8>) |
| `c0_275` | 对应水体档位的 [散射系数 b,总衰减系数 c]（m⁻¹）；变量名小数点用下划线表示 | `[0.2234, 0.275]` | 赋值 | [fig3-7/Table_abc.m:9](<C:/Users/Spidey/Desktop/code/fig3-7/Table_abc.m:9>) |
| `c0_15` | 对应水体档位的 [散射系数 b,总衰减系数 c]（m⁻¹）；变量名小数点用下划线表示 | `[0.08,0.15]` | 赋值 | [fig3-7/Table_abc.m:10](<C:/Users/Spidey/Desktop/code/fig3-7/Table_abc.m:10>) |
| `c0_305` | 对应水体档位的 [散射系数 b,总衰减系数 c]（m⁻¹）；变量名小数点用下划线表示 | `[0.216,0.305]` | 赋值 | [fig3-7/Table_abc.m:11](<C:/Users/Spidey/Desktop/code/fig3-7/Table_abc.m:11>) |
| `c0_151` | 对应水体档位的 [散射系数 b,总衰减系数 c]（m⁻¹）；变量名小数点用下划线表示 | `[0.037,0.151]` | 赋值 | [fig3-7/Table_abc.m:12](<C:/Users/Spidey/Desktop/code/fig3-7/Table_abc.m:12>) |
| `c0_398` | 对应水体档位的 [散射系数 b,总衰减系数 c]（m⁻¹）；变量名小数点用下划线表示 | `[0.219,0.398]` | 赋值 | [fig3-7/Table_abc.m:13](<C:/Users/Spidey/Desktop/code/fig3-7/Table_abc.m:13>) |
| `c2_190` | 对应水体档位的 [散射系数 b,总衰减系数 c]（m⁻¹）；变量名小数点用下划线表示 | `[1.824, 2.190]` | 赋值 | [fig3-7/Table_abc.m:14](<C:/Users/Spidey/Desktop/code/fig3-7/Table_abc.m:14>) |
| `c2_17` | 对应水体档位的 [散射系数 b,总衰减系数 c]（m⁻¹）；变量名小数点用下划线表示 | `[1.875,2.17]` | 赋值 | [fig3-7/Table_abc.m:15](<C:/Users/Spidey/Desktop/code/fig3-7/Table_abc.m:15>) |
| `paras` | 按 water_type 选中的 [b,c] 水体系数向量（m⁻¹） | `c3_2`<br>`c2_675`<br>`c1_92`<br>`c1_475`<br>`c1_0`<br>`c0_5`<br>`c0_275`<br>`c0_15`<br>`c0_305`<br>`c2_17`<br>`c0_151`<br>`c0_398`<br>`c2_190` | 赋值 | [fig3-7/Table_abc.m:17](<C:/Users/Spidey/Desktop/code/fig3-7/Table_abc.m:17>) |


### 5.34 fig3-7/Hen_Green.m

**作用域：`Hen_Green`**

| 源码符号 | 含义与单位/维度 | 源码取值/计算式（各次赋值） | 角色 | 源码位置 |
| --- | --- | --- | --- | --- |
| `length` | 本次生成的 HG 散射极角样本数；变量遮蔽了 MATLAB 同名函数 | 函数输入/输出；值随调用或上游结果确定 | 输入 | [fig3-7/Hen_Green.m:1](<C:/Users/Spidey/Desktop/code/fig3-7/Hen_Green.m:1>) |
| `judge` | HG 参数选择标志：1 对应 g=0.86，其余对应 g=0.924；W_RIS 传 2 | 函数输入/输出；值随调用或上游结果确定 | 输入 | [fig3-7/Hen_Green.m:1](<C:/Users/Spidey/Desktop/code/fig3-7/Hen_Green.m:1>) |
| `angle` | HG 散射极角样本，acos(val)（rad，返回行向量） | `acos(val)` | 输出/赋值 | [fig3-7/Hen_Green.m:1](<C:/Users/Spidey/Desktop/code/fig3-7/Hen_Green.m:1>) |
| `g` | HG 各向异性因子，judge==1 取 0.86，否则 0.924 | `0.86`<br>`0.924` | 赋值 | [fig3-7/Hen_Green.m:3](<C:/Users/Spidey/Desktop/code/fig3-7/Hen_Green.m:3>) |
| `val` | HG 逆累积分布采样得到的散射极角余弦 cos(theta) | `(1+g^2-((1-g^2)./(1-g+2*g*rand(1,length))).^2)*(1/(2*g))` | 赋值 | [fig3-7/Hen_Green.m:7](<C:/Users/Spidey/Desktop/code/fig3-7/Hen_Green.m:7>) |


## 6. 无独立变量名的调用参数、检查门限与绘图设置

第 5 节的赋值表达式已经保留内部数值系数、字符串选项和每次覆盖值。下表补充不出现在赋值右侧、而直接传入函数调用中的参数；绘图文字、RGB、线宽、采样数、图窗尺寸、输出格式等也在覆盖范围内。错误消息、报告正文格式串与 `disp` 本身不构成模型/显示参数，不重复长篇正文。

**sea_wave/audit_sea_spectrum.m**

| 作用域 | 含义 | 完整调用及内联参数 | 源码位置 |
| --- | --- | --- | --- |
| `audit_sea_spectrum` | 数值/工况检查判据及容差，原式中的比较值为实际门限 | `if nargin == 0` | [sea_wave/audit_sea_spectrum.m:5](<C:/Users/Spidey/Desktop/code/sea_wave/audit_sea_spectrum.m:5>) |
| `audit_sea_spectrum` | 数值/工况检查判据及容差，原式中的比较值为实际门限 | `assert(abs(variance/exact-1)<0.002,'Coarse spectral energy error exceeds 0.2%%.')` | [sea_wave/audit_sea_spectrum.m:20](<C:/Users/Spidey/Desktop/code/sea_wave/audit_sea_spectrum.m:20>) |
| `audit_sea_spectrum` | 数值/工况检查判据及容差，原式中的比较值为实际门限 | `assert(abs(denseVariance/exact-1)<0.00003,'Refined spectral quadrature did not converge.')` | [sea_wave/audit_sea_spectrum.m:21](<C:/Users/Spidey/Desktop/code/sea_wave/audit_sea_spectrum.m:21>) |
| `audit_sea_spectrum` | 数值/工况检查判据及容差，原式中的比较值为实际门限 | `assert(abs(sum(spectrum.frequencyWeightsHz)-(p.fMaxHz-p.fMinHz))<1e-12)` | [sea_wave/audit_sea_spectrum.m:22](<C:/Users/Spidey/Desktop/code/sea_wave/audit_sea_spectrum.m:22>) |
| `audit_sea_spectrum` | 数值/工况检查判据及容差，原式中的比较值为实际门限 | `assert(abs(fixedError)<0.002 && formerError>0.05)` | [sea_wave/audit_sea_spectrum.m:33](<C:/Users/Spidey/Desktop/code/sea_wave/audit_sea_spectrum.m:33>) |
| `audit_sea_spectrum` | 输出文件名、保存格式、编码或图像分辨率选项 | `writetable(cases,fullfile(outputFolder,'corrected_spectrum_audit.csv'))` | [sea_wave/audit_sea_spectrum.m:38](<C:/Users/Spidey/Desktop/code/sea_wave/audit_sea_spectrum.m:38>) |
| `audit_sea_spectrum` | 输出文件名、保存格式、编码或图像分辨率选项 | `save(fullfile(outputFolder,'corrected_spectrum_audit.mat'),'report')` | [sea_wave/audit_sea_spectrum.m:39](<C:/Users/Spidey/Desktop/code/sea_wave/audit_sea_spectrum.m:39>) |


**sea_wave/main3_windseasurface_0326.m**

| 作用域 | 含义 | 完整调用及内联参数 | 源码位置 |
| --- | --- | --- | --- |
| `script` | 绘图数据与颜色、线宽、透明度、标记或曲线样式参数 | `surf(ax,X,Y,Z,'EdgeColor','none')` | [sea_wave/main3_windseasurface_0326.m:22](<C:/Users/Spidey/Desktop/code/sea_wave/main3_windseasurface_0326.m:22>) |
| `script` | 坐标范围、色轴、刻度、图框比例、视角或等比例显示参数 | `zlim(ax,[-0.5,0.5])` | [sea_wave/main3_windseasurface_0326.m:24](<C:/Users/Spidey/Desktop/code/sea_wave/main3_windseasurface_0326.m:24>) |
| `script` | 坐标范围、色轴、刻度、图框比例、视角或等比例显示参数 | `clim(ax,[-0.5,0.5])` | [sea_wave/main3_windseasurface_0326.m:24](<C:/Users/Spidey/Desktop/code/sea_wave/main3_windseasurface_0326.m:24>) |
| `script` | 坐标范围、色轴、刻度、图框比例、视角或等比例显示参数 | `zticks(ax,[-0.5,0,0.5])` | [sea_wave/main3_windseasurface_0326.m:24](<C:/Users/Spidey/Desktop/code/sea_wave/main3_windseasurface_0326.m:24>) |
| `script` | 坐标范围、色轴、刻度、图框比例、视角或等比例显示参数 | `pbaspect(ax,[1,1,0.7])` | [sea_wave/main3_windseasurface_0326.m:25](<C:/Users/Spidey/Desktop/code/sea_wave/main3_windseasurface_0326.m:25>) |
| `script` | 坐标范围、色轴、刻度、图框比例、视角或等比例显示参数 | `view(ax,-37.5,30)` | [sea_wave/main3_windseasurface_0326.m:25](<C:/Users/Spidey/Desktop/code/sea_wave/main3_windseasurface_0326.m:25>) |
| `script` | 图形标签、字体、颜色、网格与布局等显示参数 | `grid(ax,'on')` | [sea_wave/main3_windseasurface_0326.m:25](<C:/Users/Spidey/Desktop/code/sea_wave/main3_windseasurface_0326.m:25>) |
| `script` | 图形标签、字体、颜色、网格与布局等显示参数 | `box(ax,'on')` | [sea_wave/main3_windseasurface_0326.m:25](<C:/Users/Spidey/Desktop/code/sea_wave/main3_windseasurface_0326.m:25>) |
| `script` | 图形标签、字体、颜色、网格与布局等显示参数 | `colormap(ax,turbo(256))` | [sea_wave/main3_windseasurface_0326.m:26](<C:/Users/Spidey/Desktop/code/sea_wave/main3_windseasurface_0326.m:26>) |
| `script` | 图形标签、字体、颜色、网格与布局等显示参数 | `xlabel(ax,'x (m)')` | [sea_wave/main3_windseasurface_0326.m:27](<C:/Users/Spidey/Desktop/code/sea_wave/main3_windseasurface_0326.m:27>) |
| `script` | 图形标签、字体、颜色、网格与布局等显示参数 | `ylabel(ax,'y (m)')` | [sea_wave/main3_windseasurface_0326.m:27](<C:/Users/Spidey/Desktop/code/sea_wave/main3_windseasurface_0326.m:27>) |
| `script` | 图形标签、字体、颜色、网格与布局等显示参数 | `zlabel(ax,'z (m)')` | [sea_wave/main3_windseasurface_0326.m:27](<C:/Users/Spidey/Desktop/code/sea_wave/main3_windseasurface_0326.m:27>) |
| `script` | 图形标签、字体、颜色、网格与布局等显示参数 | `title(ax,sprintf('U=%g m/s \| wind %g deg \| t_{sea}=%g s',exp_cases(i,:),params.t0))` | [sea_wave/main3_windseasurface_0326.m:28](<C:/Users/Spidey/Desktop/code/sea_wave/main3_windseasurface_0326.m:28>) |
| `script` | 输出文件名、保存格式、编码或图像分辨率选项 | `savefig(fig,fullfile(output_dir,[name,'.fig']))` | [sea_wave/main3_windseasurface_0326.m:30](<C:/Users/Spidey/Desktop/code/sea_wave/main3_windseasurface_0326.m:30>) |
| `script` | 输出文件名、保存格式、编码或图像分辨率选项 | `exportgraphics(fig,fullfile(output_dir,[name,'.png']),'Resolution',160)` | [sea_wave/main3_windseasurface_0326.m:31](<C:/Users/Spidey/Desktop/code/sea_wave/main3_windseasurface_0326.m:31>) |
| `script` | 输出文件名、保存格式、编码或图像分辨率选项 | `save(fullfile(output_dir,'sea_surface_comparison.mat'),'sim_data','params','exp_cases')` | [sea_wave/main3_windseasurface_0326.m:34](<C:/Users/Spidey/Desktop/code/sea_wave/main3_windseasurface_0326.m:34>) |


**sea_wave/oceanSurfaceModel_phiwind_260315.m**

| 作用域 | 含义 | 完整调用及内联参数 | 源码位置 |
| --- | --- | --- | --- |
| `oceanSurfaceModel_phiwind_260315` | 数值/工况检查判据及容差，原式中的比较值为实际门限 | `if nargin < 11` | [sea_wave/oceanSurfaceModel_phiwind_260315.m:25](<C:/Users/Spidey/Desktop/code/sea_wave/oceanSurfaceModel_phiwind_260315.m:25>) |
| `oceanSurfaceModel_phiwind_260315` | 数值/工况检查判据及容差，原式中的比较值为实际门限 | `if nargin < 12` | [sea_wave/oceanSurfaceModel_phiwind_260315.m:28](<C:/Users/Spidey/Desktop/code/sea_wave/oceanSurfaceModel_phiwind_260315.m:28>) |
| `oceanSurfaceModel_phiwind_260315` | 数值/工况检查判据及容差，原式中的比较值为实际门限 | `if nargin < 13` | [sea_wave/oceanSurfaceModel_phiwind_260315.m:31](<C:/Users/Spidey/Desktop/code/sea_wave/oceanSurfaceModel_phiwind_260315.m:31>) |
| `oceanSurfaceModel_phiwind_260315` | 数值/工况检查判据及容差，原式中的比较值为实际门限 | `if ~ismember(output_mode, {'all', 'single'})` | [sea_wave/oceanSurfaceModel_phiwind_260315.m:36](<C:/Users/Spidey/Desktop/code/sea_wave/oceanSurfaceModel_phiwind_260315.m:36>) |
| `oceanSurfaceModel_phiwind_260315` | 数值/工况检查判据及容差，原式中的比较值为实际门限 | `if strcmp(output_mode,'single')` | [sea_wave/oceanSurfaceModel_phiwind_260315.m:50](<C:/Users/Spidey/Desktop/code/sea_wave/oceanSurfaceModel_phiwind_260315.m:50>) |
| `fitPlaneLS` | 数值/工况检查判据及容差，原式中的比较值为实际门限 | `if normal_length > 0` | [sea_wave/oceanSurfaceModel_phiwind_260315.m:110](<C:/Users/Spidey/Desktop/code/sea_wave/oceanSurfaceModel_phiwind_260315.m:110>) |
| `fitPlaneLS` | 数值/工况检查判据及容差，原式中的比较值为实际门限 | `if normal(3) < 0` | [sea_wave/oceanSurfaceModel_phiwind_260315.m:115](<C:/Users/Spidey/Desktop/code/sea_wave/oceanSurfaceModel_phiwind_260315.m:115>) |


**sea_wave/RIS_draw_state.m**

| 作用域 | 含义 | 完整调用及内联参数 | 源码位置 |
| --- | --- | --- | --- |
| `RIS_draw_state` | 此源码调用的内联选项/固定常数；具体实参及位置完整列于右侧 | `ax (1,1) matlab.graphics.axis.Axes` | [sea_wave/RIS_draw_state.m:5](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_draw_state.m:5>) |
| `RIS_draw_state` | 此源码调用的内联选项/固定常数；具体实参及位置完整列于右侧 | `motion (1,1) struct` | [sea_wave/RIS_draw_state.m:6](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_draw_state.m:6>) |
| `RIS_draw_state` | 输入类型、维度、允许值及数值范围约束；原式列出全部界限 | `time (1,1) double {mustBeReal,mustBeFinite,mustBeNonnegative}` | [sea_wave/RIS_draw_state.m:7](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_draw_state.m:7>) |
| `RIS_draw_state` | 此源码调用的内联选项/固定常数；具体实参及位置完整列于右侧 | `scene (1,1) struct` | [sea_wave/RIS_draw_state.m:8](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_draw_state.m:8>) |
| `RIS_draw_state` | 绘图数据与颜色、线宽、透明度、标记或曲线样式参数 | `surf(ax,scene.X,scene.Y,scene.Z,'EdgeColor','none')` | [sea_wave/RIS_draw_state.m:15](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_draw_state.m:15>) |
| `RIS_draw_state` | 此源码调用的内联选项/固定常数；具体实参及位置完整列于右侧 | `hold(ax,'on')` | [sea_wave/RIS_draw_state.m:16](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_draw_state.m:16>) |
| `RIS_draw_state` | 绘图数据与颜色、线宽、透明度、标记或曲线样式参数 | `patch(ax,'Vertices',motion.corners(:,:,index),'Faces',[1,2,3,4], 'FaceColor',[1,0.48,0.08],'FaceAlpha',1,'EdgeColor','k','LineWidth',1.3)` | [sea_wave/RIS_draw_state.m:17](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_draw_state.m:17>) |
| `RIS_draw_state` | 绘图数据与颜色、线宽、透明度、标记或曲线样式参数 | `quiver3(ax,center(1),center(2),center(3),0.5*normal(1),0.5*normal(2),0.5*normal(3), 0,'k','LineWidth',1.4)` | [sea_wave/RIS_draw_state.m:20](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_draw_state.m:20>) |
| `RIS_draw_state` | 数值/工况检查判据及容差，原式中的比较值为实际门限 | `if isfield(motion,'rotationMatrix')` | [sea_wave/RIS_draw_state.m:22](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_draw_state.m:22>) |
| `RIS_draw_state` | 绘图数据与颜色、线宽、透明度、标记或曲线样式参数 | `quiver3(ax,center(1),center(2),center(3),direction(1),direction(2),direction(3), 0,'Color',colors(bodyAxis,:),'LineWidth',1.5)` | [sea_wave/RIS_draw_state.m:33](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_draw_state.m:33>) |
| `RIS_draw_state` | 坐标范围、色轴、刻度、图框比例、视角或等比例显示参数 | `pbaspect(ax,[1,1,0.7])` | [sea_wave/RIS_draw_state.m:40](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_draw_state.m:40>) |
| `RIS_draw_state` | 坐标范围、色轴、刻度、图框比例、视角或等比例显示参数 | `view(ax,-37.5,30)` | [sea_wave/RIS_draw_state.m:40](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_draw_state.m:40>) |
| `RIS_draw_state` | 图形标签、字体、颜色、网格与布局等显示参数 | `grid(ax,'on')` | [sea_wave/RIS_draw_state.m:40](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_draw_state.m:40>) |
| `RIS_draw_state` | 图形标签、字体、颜色、网格与布局等显示参数 | `box(ax,'on')` | [sea_wave/RIS_draw_state.m:40](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_draw_state.m:40>) |
| `RIS_draw_state` | 坐标范围、色轴、刻度、图框比例、视角或等比例显示参数 | `zticks(ax,[-0.5,0,0.5])` | [sea_wave/RIS_draw_state.m:41](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_draw_state.m:41>) |
| `RIS_draw_state` | 图形标签、字体、颜色、网格与布局等显示参数 | `colormap(ax,turbo(256))` | [sea_wave/RIS_draw_state.m:42](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_draw_state.m:42>) |
| `RIS_draw_state` | 图形标签、字体、颜色、网格与布局等显示参数 | `xlabel(ax,'x (m)')` | [sea_wave/RIS_draw_state.m:43](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_draw_state.m:43>) |
| `RIS_draw_state` | 图形标签、字体、颜色、网格与布局等显示参数 | `ylabel(ax,'y (m)')` | [sea_wave/RIS_draw_state.m:43](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_draw_state.m:43>) |
| `RIS_draw_state` | 图形标签、字体、颜色、网格与布局等显示参数 | `zlabel(ax,'z (m)')` | [sea_wave/RIS_draw_state.m:43](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_draw_state.m:43>) |
| `RIS_draw_state` | 图形标签、字体、颜色、网格与布局等显示参数 | `set(ax,'FontSize',9)` | [sea_wave/RIS_draw_state.m:43](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_draw_state.m:43>) |


**sea_wave/RIS_float.m**

| 作用域 | 含义 | 完整调用及内联参数 | 源码位置 |
| --- | --- | --- | --- |
| `RIS_float` | 数值/工况检查判据及容差，原式中的比较值为实际门限 | `if any(restoring<=0)` | [sea_wave/RIS_float.m:40](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:40>) |
| `initialState` | 数值/工况检查判据及容差，原式中的比较值为实际门限 | `if strcmp(p.initialization,'waveSteady')` | [sea_wave/RIS_float.m:109](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:109>) |
| `derivative` | 数值/工况检查判据及容差，原式中的比较值为实际门限 | `if nargout>1` | [sea_wave/RIS_float.m:193](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:193>) |
| `plotMotion` | 此源码调用的内联选项/固定常数；具体实参及位置完整列于右侧 | `figure('Name','Free-floating RIS six DOFs','Color','w')` | [sea_wave/RIS_float.m:314](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:314>) |
| `plotMotion` | 此源码调用的内联选项/固定常数；具体实参及位置完整列于右侧 | `tiledlayout(4,1)` | [sea_wave/RIS_float.m:315](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:315>) |
| `plotMotion` | 绘图数据与颜色、线宽、透明度、标记或曲线样式参数 | `plot(o.t,[o.etaReferenceCenter;o.heave].')` | [sea_wave/RIS_float.m:316](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:316>) |
| `plotMotion` | 图形标签、字体、颜色、网格与布局等显示参数 | `ylabel('z (m)')` | [sea_wave/RIS_float.m:316](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:316>) |
| `plotMotion` | 图形标签、字体、颜色、网格与布局等显示参数 | `legend('Reference sea','Heave')` | [sea_wave/RIS_float.m:317](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:317>) |
| `plotMotion` | 绘图数据与颜色、线宽、透明度、标记或曲线样式参数 | `plot(o.t,[o.rollDeg;o.pitchDeg;o.yawDeg].')` | [sea_wave/RIS_float.m:318](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:318>) |
| `plotMotion` | 图形标签、字体、颜色、网格与布局等显示参数 | `ylabel('Angle (deg)')` | [sea_wave/RIS_float.m:318](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:318>) |
| `plotMotion` | 图形标签、字体、颜色、网格与布局等显示参数 | `legend('Roll','Pitch','Yaw')` | [sea_wave/RIS_float.m:319](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:319>) |
| `plotMotion` | 绘图数据与颜色、线宽、透明度、标记或曲线样式参数 | `plot(o.t,[o.surge;o.sway].')` | [sea_wave/RIS_float.m:320](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:320>) |
| `plotMotion` | 图形标签、字体、颜色、网格与布局等显示参数 | `ylabel('Position (m)')` | [sea_wave/RIS_float.m:320](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:320>) |
| `plotMotion` | 图形标签、字体、颜色、网格与布局等显示参数 | `legend('Surge','Sway')` | [sea_wave/RIS_float.m:320](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:320>) |
| `plotMotion` | 图形标签、字体、颜色、网格与布局等显示参数 | `ylabel('Yaw moment (N m)')` | [sea_wave/RIS_float.m:321](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:321>) |
| `plotMotion` | 图形标签、字体、颜色、网格与布局等显示参数 | `xlabel('Time (s)')` | [sea_wave/RIS_float.m:321](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:321>) |
| `plotMotion` | 图形标签、字体、颜色、网格与布局等显示参数 | `sgtitle('Six-DOF rigid-body response with approximate linear hydrodynamics')` | [sea_wave/RIS_float.m:322](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:322>) |
| `plotMotion` | 数值/工况检查判据及容差，原式中的比较值为实际门限 | `if ~o.validity.smallAngleUnder10Deg` | [sea_wave/RIS_float.m:324](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:324>) |
| `validateParameters` | 输入类型、维度、允许值及数值范围约束；原式列出全部界限 | `validateattributes(p.(positive{j}),{'numeric'},{'scalar','real','finite','positive'}, 'RIS_float',positive{j})` | [sea_wave/RIS_float.m:332](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:332>) |
| `validateParameters` | 输入类型、维度、允许值及数值范围约束；原式列出全部界限 | `validateattributes(p.windSpeed,{'numeric'},{'scalar','real','finite','>=',3,'<=',25})` | [sea_wave/RIS_float.m:335](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:335>) |
| `validateParameters` | 输入类型、维度、允许值及数值范围约束；原式列出全部界限 | `validateattributes(p.windAngleDeg,{'numeric'},{'scalar','real','finite'})` | [sea_wave/RIS_float.m:336](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:336>) |
| `validateParameters` | 输入类型、维度、允许值及数值范围约束；原式列出全部界限 | `validateattributes(p.startTime,{'numeric'},{'scalar','real','finite','nonnegative'})` | [sea_wave/RIS_float.m:337](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:337>) |
| `validateParameters` | 输入类型、维度、允许值及数值范围约束；原式列出全部界限 | `validateattributes(p.frequencyBins,{'numeric'},{'scalar','finite','integer','>=',2})` | [sea_wave/RIS_float.m:338](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:338>) |
| `validateParameters` | 输入类型、维度、允许值及数值范围约束；原式列出全部界限 | `validateattributes(p.directionBins,{'numeric'},{'scalar','finite','integer','>=',2})` | [sea_wave/RIS_float.m:339](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:339>) |
| `validateParameters` | 输入类型、维度、允许值及数值范围约束；原式列出全部界限 | `validateattributes(p.panelCount,{'numeric'},{'scalar','finite','integer','>=',4,'<=',80})` | [sea_wave/RIS_float.m:340](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:340>) |
| `validateParameters` | 输入类型、维度、允许值及数值范围约束；原式列出全部界限 | `validateattributes(p.radiationModes,{'numeric'},{'scalar','integer','>=',4})` | [sea_wave/RIS_float.m:341](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:341>) |
| `validateParameters` | 输入类型、维度、允许值及数值范围约束；原式列出全部界限 | `validateattributes(p.radiationFrequencyBins,{'numeric'},{'scalar','integer','>=',8})` | [sea_wave/RIS_float.m:342](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:342>) |
| `validateParameters` | 输入类型、维度、允许值及数值范围约束；原式列出全部界限 | `validatestring(p.initialization,{'waveSteady','rest'})` | [sea_wave/RIS_float.m:343](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:343>) |
| `validateParameters` | 输入类型、维度、允许值及数值范围约束；原式列出全部界限 | `validateattributes(p.(names{j}),{'numeric'},{'vector','numel',sizes(j),'real','finite'})` | [sea_wave/RIS_float.m:347](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:347>) |
| `validateParameters` | 数值/工况检查判据及容差，原式中的比较值为实际门限 | `if mod(p.panelCount,2)~=0` | [sea_wave/RIS_float.m:351](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:351>) |
| `validateParameters` | 输入类型、维度、允许值及数值范围约束；原式列出全部界限 | `validateattributes(p.seed,{'numeric'},{'scalar','real','finite','integer','>=',0,'<=',2^32-1})` | [sea_wave/RIS_float.m:352](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:352>) |
| `validateParameters` | 数值/工况检查判据及容差，原式中的比较值为实际门限 | `if j==3 && isempty(p.(pairs{j}))` | [sea_wave/RIS_float.m:355](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:355>) |
| `validateParameters` | 输入类型、维度、允许值及数值范围约束；原式列出全部界限 | `validateattributes(p.(pairs{j}),{'numeric'},{'vector','numel',2,'real','finite'})` | [sea_wave/RIS_float.m:356](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:356>) |
| `validateParameters` | 输入类型、维度、允许值及数值范围约束；原式列出全部界限 | `validateattributes(p.plot,{'logical','numeric'},{'scalar','real','finite'})` | [sea_wave/RIS_float.m:359](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:359>) |
| `validateParameters` | 数值/工况检查判据及容差，原式中的比较值为实际门限 | `if ~ismember(p.plot,[0,1])` | [sea_wave/RIS_float.m:360](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:360>) |
| `makeWave` | 输入类型、维度、允许值及数值范围约束；原式列出全部界限 | `validateattributes(w.(required{j}),{'numeric'},{'vector','nonempty','finite'})` | [sea_wave/RIS_float.m:372](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:372>) |
| `makeWave` | 输入类型、维度、允许值及数值范围约束；原式列出全部界限 | `validateattributes(w.omega,{'numeric'},{'real','positive'})` | [sea_wave/RIS_float.m:375](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:375>) |
| `makeWave` | 输入类型、维度、允许值及数值范围约束；原式列出全部界限 | `validateattributes(w.kx,{'numeric'},{'real'})` | [sea_wave/RIS_float.m:376](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:376>) |
| `makeWave` | 输入类型、维度、允许值及数值范围约束；原式列出全部界限 | `validateattributes(w.ky,{'numeric'},{'real'})` | [sea_wave/RIS_float.m:377](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:377>) |
| `makeWave` | 数值/工况检查判据及容差，原式中的比较值为实际门限 | `if any([numel(w.kx),numel(w.ky),numel(w.coefficient)]~=numel(w.omega)) \|\| any(abs(hypot(w.kx,w.ky)-w.omega.^2/p.g)>1e-9*max(1,w.omega.^2/p.g))` | [sea_wave/RIS_float.m:378](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:378>) |
| `plateResponse` | 数值/工况检查判据及容差，原式中的比较值为实际门限 | `if numel(keys)>=4` | [sea_wave/RIS_float.m:402](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:402>) |
| `solvePotentials` | 数值/工况检查判据及容差，原式中的比较值为实际门限 | `if any(B(:)<-1e-6) \|\| any(~isfinite(R(:)))` | [sea_wave/RIS_float.m:479](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:479>) |
| `plotSea` | 此源码调用的内联选项/固定常数；具体实参及位置完整列于右侧 | `figure('Name','RIS and the same incident sea','Color','w')` | [sea_wave/RIS_float.m:508](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:508>) |
| `plotSea` | 图形标签、字体、颜色、网格与布局等显示参数 | `title(ax,sprintf('%s \| sea time %.3f s',labels(j),o.absoluteTime(index)))` | [sea_wave/RIS_float.m:513](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:513>) |
| `plotSea` | 图形标签、字体、颜色、网格与布局等显示参数 | `sgtitle(layout,'Four states from one trajectory; z display range +/-0.5 m')` | [sea_wave/RIS_float.m:515](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_float.m:515>) |


**sea_wave/RIS_radiation.m**

| 作用域 | 含义 | 完整调用及内联参数 | 源码位置 |
| --- | --- | --- | --- |
| `RIS_radiation` | 输入类型、维度、允许值及数值范围约束；原式列出全部界限 | `omega (:,1) double {mustBePositive,mustBeFinite}` | [sea_wave/RIS_radiation.m:6](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_radiation.m:6>) |
| `RIS_radiation` | 输入类型、维度、允许值及数值范围约束；原式列出全部界限 | `addedMass (3,:) double {mustBeFinite}` | [sea_wave/RIS_radiation.m:7](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_radiation.m:7>) |
| `RIS_radiation` | 输入类型、维度、允许值及数值范围约束；原式列出全部界限 | `damping (3,:) double {mustBeFinite}` | [sea_wave/RIS_radiation.m:8](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_radiation.m:8>) |
| `RIS_radiation` | 输入类型、维度、允许值及数值范围约束；原式列出全部界限 | `modeCount (1,1) double {mustBeInteger,mustBeGreaterThanOrEqual(modeCount,4)}` | [sea_wave/RIS_radiation.m:9](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_radiation.m:9>) |


**sea_wave/RIS_rotation.m**

| 作用域 | 含义 | 完整调用及内联参数 | 源码位置 |
| --- | --- | --- | --- |
| `RIS_rotation` | 输入类型、维度、允许值及数值范围约束；原式列出全部界限 | `value (:,1) double {mustBeReal,mustBeFinite}` | [sea_wave/RIS_rotation.m:4](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_rotation.m:4>) |
| `RIS_rotation` | 数值/工况检查判据及容差，原式中的比较值为实际门限 | `if numel(value)==3` | [sea_wave/RIS_rotation.m:6](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_rotation.m:6>) |


**sea_wave/RIS_sea_surface.m**

| 作用域 | 含义 | 完整调用及内联参数 | 源码位置 |
| --- | --- | --- | --- |
| `RIS_sea_surface` | 输入类型、维度、允许值及数值范围约束；原式列出全部界限 | `validateattributes(x,{'numeric'},{'vector','real','finite','nonempty'})` | [sea_wave/RIS_sea_surface.m:4](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_sea_surface.m:4>) |
| `RIS_sea_surface` | 输入类型、维度、允许值及数值范围约束；原式列出全部界限 | `validateattributes(y,{'numeric'},{'vector','real','finite','nonempty'})` | [sea_wave/RIS_sea_surface.m:5](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_sea_surface.m:5>) |
| `RIS_sea_surface` | 输入类型、维度、允许值及数值范围约束；原式列出全部界限 | `validateattributes(time,{'numeric'},{'scalar','real','finite','nonnegative'})` | [sea_wave/RIS_sea_surface.m:6](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_sea_surface.m:6>) |


**sea_wave/RIS_sidewall_load.m**

| 作用域 | 含义 | 完整调用及内联参数 | 源码位置 |
| --- | --- | --- | --- |
| `RIS_sidewall_load` | 此源码调用的内联选项/固定常数；具体实参及位置完整列于右侧 | `p (1,1) struct` | [sea_wave/RIS_sidewall_load.m:6](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_sidewall_load.m:6>) |
| `RIS_sidewall_load` | 此源码调用的内联选项/固定常数；具体实参及位置完整列于右侧 | `wave (1,1) struct` | [sea_wave/RIS_sidewall_load.m:7](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_sidewall_load.m:7>) |
| `RIS_sidewall_load` | 输入类型、维度、允许值及数值范围约束；原式列出全部界限 | `draft (1,1) double {mustBePositive}` | [sea_wave/RIS_sidewall_load.m:8](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_sidewall_load.m:8>) |
| `RIS_sidewall_load` | 输入类型、维度、允许值及数值范围约束；原式列出全部界限 | `zG (1,1) double {mustBeReal,mustBeFinite}` | [sea_wave/RIS_sidewall_load.m:9](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_sidewall_load.m:9>) |


**sea_wave/RIS_wave_components.m**

| 作用域 | 含义 | 完整调用及内联参数 | 源码位置 |
| --- | --- | --- | --- |
| `RIS_wave_components` | 此源码调用的内联选项/固定常数；具体实参及位置完整列于右侧 | `p (1,1) struct` | [sea_wave/RIS_wave_components.m:6](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_wave_components.m:6>) |


**sea_wave/RIS_wave_field.m**

| 作用域 | 含义 | 完整调用及内联参数 | 源码位置 |
| --- | --- | --- | --- |
| `RIS_wave_field` | 此源码调用的内联选项/固定常数；具体实参及位置完整列于右侧 | `motion (1,1) struct` | [sea_wave/RIS_wave_field.m:6](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_wave_field.m:6>) |
| `RIS_wave_field` | 输入类型、维度、允许值及数值范围约束；原式列出全部界限 | `points (:,2) double {mustBeReal,mustBeFinite,mustBeNonempty}` | [sea_wave/RIS_wave_field.m:7](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_wave_field.m:7>) |
| `RIS_wave_field` | 数值/工况检查判据及容差，原式中的比较值为实际门限 | `if ~isscalar(time) && numel(time)~=size(points,1)` | [sea_wave/RIS_wave_field.m:10](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_wave_field.m:10>) |
| `RIS_wave_field` | 数值/工况检查判据及容差，原式中的比较值为实际门限 | `if nargout>1` | [sea_wave/RIS_wave_field.m:24](<C:/Users/Spidey/Desktop/code/sea_wave/RIS_wave_field.m:24>) |


**sea_wave/run_RIS_wind_comparison.m**

| 作用域 | 含义 | 完整调用及内联参数 | 源码位置 |
| --- | --- | --- | --- |
| `run_RIS_wind_comparison` | 数值/工况检查判据及容差，原式中的比较值为实际门限 | `if nargin == 0` | [sea_wave/run_RIS_wind_comparison.m:6](<C:/Users/Spidey/Desktop/code/sea_wave/run_RIS_wind_comparison.m:6>) |
| `run_RIS_wind_comparison` | 输出文件名、保存格式、编码或图像分辨率选项 | `writetable(summary,fullfile(outputFolder,'RIS_case_summary.csv'))` | [sea_wave/run_RIS_wind_comparison.m:43](<C:/Users/Spidey/Desktop/code/sea_wave/run_RIS_wind_comparison.m:43>) |
| `run_RIS_wind_comparison` | 数值/工况检查判据及容差，原式中的比较值为实际门限 | `if k<=3` | [sea_wave/run_RIS_wind_comparison.m:51](<C:/Users/Spidey/Desktop/code/sea_wave/run_RIS_wind_comparison.m:51>) |
| `run_RIS_wind_comparison` | 图形标签、字体、颜色、网格与布局等显示参数 | `set(ax,'FontSize',9)` | [sea_wave/run_RIS_wind_comparison.m:67](<C:/Users/Spidey/Desktop/code/sea_wave/run_RIS_wind_comparison.m:67>) |
| `run_RIS_wind_comparison` | 图形标签、字体、颜色、网格与布局等显示参数 | `sgtitle(layout,{'Sea surface and RIS: winds 3, 5, 9, 12 m/s (rows)', 'Columns: wind 30 deg \| wind 60 deg \| wind 90 deg \| 30-deg centre trajectory', 'Snapshots: sea t=27 s \| paths: sea t=25-35 s \| z limits +/-0.5 m (display crop)'}, 'FontWeight','normal','FontSize',15)` | [sea_wave/run_RIS_wind_comparison.m:70](<C:/Users/Spidey/Desktop/code/sea_wave/run_RIS_wind_comparison.m:70>) |
| `run_RIS_wind_comparison` | 输出文件名、保存格式、编码或图像分辨率选项 | `save(fullfile(outputFolder,'RIS_wind_direction_results.mat'), 'motions','scenes','summary','meta','panelFiles','overviewFile')` | [sea_wave/run_RIS_wind_comparison.m:83](<C:/Users/Spidey/Desktop/code/sea_wave/run_RIS_wind_comparison.m:83>) |
| `drawCell` | 数值/工况检查判据及容差，原式中的比较值为实际门限 | `if column<=3` | [sea_wave/run_RIS_wind_comparison.m:90](<C:/Users/Spidey/Desktop/code/sea_wave/run_RIS_wind_comparison.m:90>) |
| `drawCell` | 图形标签、字体、颜色、网格与布局等显示参数 | `title(ax,sprintf('U=%g m/s \| wind %g deg \| t_{sea}=27 s', o.parameters.windSpeed,o.parameters.windAngleDeg),'FontWeight','normal')` | [sea_wave/run_RIS_wind_comparison.m:93](<C:/Users/Spidey/Desktop/code/sea_wave/run_RIS_wind_comparison.m:93>) |
| `drawCell` | 绘图数据与颜色、线宽、透明度、标记或曲线样式参数 | `plot3(ax,p(1,:),p(2,:),p(3,:),'Color',[0.15,0.3,0.65],'LineWidth',1.4)` | [sea_wave/run_RIS_wind_comparison.m:98](<C:/Users/Spidey/Desktop/code/sea_wave/run_RIS_wind_comparison.m:98>) |
| `drawCell` | 此源码调用的内联选项/固定常数；具体实参及位置完整列于右侧 | `hold(ax,'on')` | [sea_wave/run_RIS_wind_comparison.m:99](<C:/Users/Spidey/Desktop/code/sea_wave/run_RIS_wind_comparison.m:99>) |
| `drawCell` | 绘图数据与颜色、线宽、透明度、标记或曲线样式参数 | `plot3(ax,p(1,:),p(2,:),repmat(zLimits(1),size(o.t)),':', 'Color',[0.5,0.5,0.5],'LineWidth',1)` | [sea_wave/run_RIS_wind_comparison.m:101](<C:/Users/Spidey/Desktop/code/sea_wave/run_RIS_wind_comparison.m:101>) |
| `drawCell` | 绘图数据与颜色、线宽、透明度、标记或曲线样式参数 | `scatter3(ax,p(1,inRange(1)),p(2,inRange(1)),p(3,inRange(1)),35, [0.1,0.55,0.25],'filled')` | [sea_wave/run_RIS_wind_comparison.m:105](<C:/Users/Spidey/Desktop/code/sea_wave/run_RIS_wind_comparison.m:105>) |
| `drawCell` | 绘图数据与颜色、线宽、透明度、标记或曲线样式参数 | `scatter3(ax,p(1,inRange(end)),p(2,inRange(end)),p(3,inRange(end)),35, [0.8,0.15,0.1],'filled')` | [sea_wave/run_RIS_wind_comparison.m:107](<C:/Users/Spidey/Desktop/code/sea_wave/run_RIS_wind_comparison.m:107>) |
| `drawCell` | 图形标签、字体、颜色、网格与布局等显示参数 | `xlabel(ax,'x_G (m)')` | [sea_wave/run_RIS_wind_comparison.m:110](<C:/Users/Spidey/Desktop/code/sea_wave/run_RIS_wind_comparison.m:110>) |
| `drawCell` | 图形标签、字体、颜色、网格与布局等显示参数 | `ylabel(ax,'y_G (m)')` | [sea_wave/run_RIS_wind_comparison.m:110](<C:/Users/Spidey/Desktop/code/sea_wave/run_RIS_wind_comparison.m:110>) |
| `drawCell` | 图形标签、字体、颜色、网格与布局等显示参数 | `zlabel(ax,'z_G (m)')` | [sea_wave/run_RIS_wind_comparison.m:110](<C:/Users/Spidey/Desktop/code/sea_wave/run_RIS_wind_comparison.m:110>) |
| `drawCell` | 坐标范围、色轴、刻度、图框比例、视角或等比例显示参数 | `zticks(ax,[-0.5,0,0.5])` | [sea_wave/run_RIS_wind_comparison.m:111](<C:/Users/Spidey/Desktop/code/sea_wave/run_RIS_wind_comparison.m:111>) |
| `drawCell` | 坐标范围、色轴、刻度、图框比例、视角或等比例显示参数 | `pbaspect(ax,[1,1,0.7])` | [sea_wave/run_RIS_wind_comparison.m:112](<C:/Users/Spidey/Desktop/code/sea_wave/run_RIS_wind_comparison.m:112>) |
| `drawCell` | 坐标范围、色轴、刻度、图框比例、视角或等比例显示参数 | `view(ax,-37.5,30)` | [sea_wave/run_RIS_wind_comparison.m:112](<C:/Users/Spidey/Desktop/code/sea_wave/run_RIS_wind_comparison.m:112>) |
| `drawCell` | 图形标签、字体、颜色、网格与布局等显示参数 | `grid(ax,'on')` | [sea_wave/run_RIS_wind_comparison.m:112](<C:/Users/Spidey/Desktop/code/sea_wave/run_RIS_wind_comparison.m:112>) |
| `drawCell` | 图形标签、字体、颜色、网格与布局等显示参数 | `box(ax,'on')` | [sea_wave/run_RIS_wind_comparison.m:112](<C:/Users/Spidey/Desktop/code/sea_wave/run_RIS_wind_comparison.m:112>) |
| `drawCell` | 图形标签、字体、颜色、网格与布局等显示参数 | `title(ax,{sprintf('U=%g m/s \| wind 30 deg \| centre path',o.parameters.windSpeed), sprintf('t_{sea}=25-35 s \| %.0f%% outside z view',100*outside)},'FontWeight','normal')` | [sea_wave/run_RIS_wind_comparison.m:114](<C:/Users/Spidey/Desktop/code/sea_wave/run_RIS_wind_comparison.m:114>) |
| `saveFigure` | 输出文件名、保存格式、编码或图像分辨率选项 | `exportgraphics(f,fullfile(folder,[name,'.png']),'Resolution',160)` | [sea_wave/run_RIS_wind_comparison.m:126](<C:/Users/Spidey/Desktop/code/sea_wave/run_RIS_wind_comparison.m:126>) |
| `exportTrajectory` | 输出文件名、保存格式、编码或图像分辨率选项 | `writetable(array2table(values,'VariableNames',names), fullfile(folder,sprintf('RIS_trajectory_U%02d_Phi%02d.csv',wind,direction)))` | [sea_wave/run_RIS_wind_comparison.m:133](<C:/Users/Spidey/Desktop/code/sea_wave/run_RIS_wind_comparison.m:133>) |


**sea_wave/summarize_RIS_motion_existing.m**

| 作用域 | 含义 | 完整调用及内联参数 | 源码位置 |
| --- | --- | --- | --- |
| `script` | 数值/工况检查判据及容差，原式中的比较值为实际门限 | `assert(numel(risFiles)==12,'Expected 12 existing trajectories.')` | [sea_wave/summarize_RIS_motion_existing.m:9](<C:/Users/Spidey/Desktop/code/sea_wave/summarize_RIS_motion_existing.m:9>) |
| `script` | 数值/工况检查判据及容差，原式中的比较值为实际门限 | `assert(height(risT)==401 && all(isfinite([risTime,risP,risAngles]),'all'))` | [sea_wave/summarize_RIS_motion_existing.m:17](<C:/Users/Spidey/Desktop/code/sea_wave/summarize_RIS_motion_existing.m:17>) |
| `script` | 数值/工况检查判据及容差，原式中的比较值为实际门限 | `assert(all(diff(risTime)>0) && max(abs(diff(risTime)-0.025))<1e-9)` | [sea_wave/summarize_RIS_motion_existing.m:18](<C:/Users/Spidey/Desktop/code/sea_wave/summarize_RIS_motion_existing.m:18>) |
| `script` | 数值/工况检查判据及容差，原式中的比较值为实际门限 | `assert(abs(risT.seaTimeS(1)-25)<1e-9 && abs(risT.seaTimeS(end)-35)<1e-9)` | [sea_wave/summarize_RIS_motion_existing.m:19](<C:/Users/Spidey/Desktop/code/sea_wave/summarize_RIS_motion_existing.m:19>) |
| `script` | 数值/工况检查判据及容差，原式中的比较值为实际门限 | `assert(abs(trace(risC)-mean(sum(risCentered.^2,2)))<1e-10)` | [sea_wave/summarize_RIS_motion_existing.m:22](<C:/Users/Spidey/Desktop/code/sea_wave/summarize_RIS_motion_existing.m:22>) |
| `script` | 数值/工况检查判据及容差，原式中的比较值为实际门限 | `assert(max(abs(sum(risNormals.^2,2)-1))<1e-10)` | [sea_wave/summarize_RIS_motion_existing.m:27](<C:/Users/Spidey/Desktop/code/sea_wave/summarize_RIS_motion_existing.m:27>) |
| `script` | 输出文件名、保存格式、编码或图像分辨率选项 | `writetable(risMotionSummary,fullfile(risOutputFolder,'existing_trajectory_motion_statistics.csv'))` | [sea_wave/summarize_RIS_motion_existing.m:55](<C:/Users/Spidey/Desktop/code/sea_wave/summarize_RIS_motion_existing.m:55>) |
| `script` | 绘图数据与颜色、线宽、透明度、标记或曲线样式参数 | `semilogy(risAx,risShow.windMps,risShow{:,{'varXM2','varYM2','varZM2'}},'-o','LineWidth',1.6)` | [sea_wave/summarize_RIS_motion_existing.m:60](<C:/Users/Spidey/Desktop/code/sea_wave/summarize_RIS_motion_existing.m:60>) |
| `script` | 图形标签、字体、颜色、网格与布局等显示参数 | `xlabel(risAx,'Wind speed U (m/s)')` | [sea_wave/summarize_RIS_motion_existing.m:61](<C:/Users/Spidey/Desktop/code/sea_wave/summarize_RIS_motion_existing.m:61>) |
| `script` | 图形标签、字体、颜色、网格与布局等显示参数 | `ylabel(risAx,'Centroid displacement variance (m^2)')` | [sea_wave/summarize_RIS_motion_existing.m:61](<C:/Users/Spidey/Desktop/code/sea_wave/summarize_RIS_motion_existing.m:61>) |
| `script` | 坐标范围、色轴、刻度、图框比例、视角或等比例显示参数 | `xticks(risAx,[3,5,9,12])` | [sea_wave/summarize_RIS_motion_existing.m:62](<C:/Users/Spidey/Desktop/code/sea_wave/summarize_RIS_motion_existing.m:62>) |
| `script` | 图形标签、字体、颜色、网格与布局等显示参数 | `grid(risAx,'on')` | [sea_wave/summarize_RIS_motion_existing.m:62](<C:/Users/Spidey/Desktop/code/sea_wave/summarize_RIS_motion_existing.m:62>) |
| `script` | 图形标签、字体、颜色、网格与布局等显示参数 | `legend(risAx,{'x_G','y_G','z_G'},'Location','northwest')` | [sea_wave/summarize_RIS_motion_existing.m:62](<C:/Users/Spidey/Desktop/code/sea_wave/summarize_RIS_motion_existing.m:62>) |
| `script` | 图形标签、字体、颜色、网格与布局等显示参数 | `title(risAx,'Translation: mean removed, drift retained')` | [sea_wave/summarize_RIS_motion_existing.m:63](<C:/Users/Spidey/Desktop/code/sea_wave/summarize_RIS_motion_existing.m:63>) |
| `script` | 绘图数据与颜色、线宽、透明度、标记或曲线样式参数 | `plot(risAx,risShow.windMps,risShow.normalAngularRmsDeg,'-o','LineWidth',1.6)` | [sea_wave/summarize_RIS_motion_existing.m:65](<C:/Users/Spidey/Desktop/code/sea_wave/summarize_RIS_motion_existing.m:65>) |
| `script` | 图形标签、字体、颜色、网格与布局等显示参数 | `xlabel(risAx,'Wind speed U (m/s)')` | [sea_wave/summarize_RIS_motion_existing.m:66](<C:/Users/Spidey/Desktop/code/sea_wave/summarize_RIS_motion_existing.m:66>) |
| `script` | 图形标签、字体、颜色、网格与布局等显示参数 | `ylabel(risAx,'Normal angular RMS (deg)')` | [sea_wave/summarize_RIS_motion_existing.m:66](<C:/Users/Spidey/Desktop/code/sea_wave/summarize_RIS_motion_existing.m:66>) |
| `script` | 坐标范围、色轴、刻度、图框比例、视角或等比例显示参数 | `xticks(risAx,[3,5,9,12])` | [sea_wave/summarize_RIS_motion_existing.m:67](<C:/Users/Spidey/Desktop/code/sea_wave/summarize_RIS_motion_existing.m:67>) |
| `script` | 图形标签、字体、颜色、网格与布局等显示参数 | `grid(risAx,'on')` | [sea_wave/summarize_RIS_motion_existing.m:67](<C:/Users/Spidey/Desktop/code/sea_wave/summarize_RIS_motion_existing.m:67>) |
| `script` | 图形标签、字体、颜色、网格与布局等显示参数 | `title(risAx,'Orientation: angle about mean board normal')` | [sea_wave/summarize_RIS_motion_existing.m:68](<C:/Users/Spidey/Desktop/code/sea_wave/summarize_RIS_motion_existing.m:68>) |
| `script` | 图形标签、字体、颜色、网格与布局等显示参数 | `title(risLayout,{'RIS motion \| wind direction parameter 30 deg', 'Existing simulation records: sea time 25-35 s, seed 600, 401 samples; normalization N'},'FontSize',15)` | [sea_wave/summarize_RIS_motion_existing.m:69](<C:/Users/Spidey/Desktop/code/sea_wave/summarize_RIS_motion_existing.m:69>) |
| `script` | 输出文件名、保存格式、编码或图像分辨率选项 | `exportgraphics(risFig,fullfile(risOutputFolder,'RIS_motion_strength_vs_wind.png'),'Resolution',170)` | [sea_wave/summarize_RIS_motion_existing.m:71](<C:/Users/Spidey/Desktop/code/sea_wave/summarize_RIS_motion_existing.m:71>) |
| `script` | 输出文件名、保存格式、编码或图像分辨率选项 | `savefig(risFig,fullfile(risOutputFolder,'RIS_motion_strength_vs_wind.fig'))` | [sea_wave/summarize_RIS_motion_existing.m:72](<C:/Users/Spidey/Desktop/code/sea_wave/summarize_RIS_motion_existing.m:72>) |


**W_RIS/wris_cached.m**

| 作用域 | 含义 | 完整调用及内联参数 | 源码位置 |
| --- | --- | --- | --- |
| `wris_cached` | 此源码调用的内联选项/固定常数；具体实参及位置完整列于右侧 | `digest.update(unicode2native(jsonencode(serializable(key)),'UTF-8'))` | [W_RIS/wris_cached.m:20](<C:/Users/Spidey/Desktop/code/W_RIS/wris_cached.m:20>) |
| `wris_cached` | 此源码调用的内联选项/固定常数；具体实参及位置完整列于右侧 | `digest.update(unicode2native(fileread(files{j}),'UTF-8'))` | [W_RIS/wris_cached.m:22](<C:/Users/Spidey/Desktop/code/W_RIS/wris_cached.m:22>) |
| `wris_cached` | 输出文件名、保存格式、编码或图像分辨率选项 | `save(file,'value','key','-v7.3')` | [W_RIS/wris_cached.m:31](<C:/Users/Spidey/Desktop/code/W_RIS/wris_cached.m:31>) |


**W_RIS/wris_config.m**

| 作用域 | 含义 | 完整调用及内联参数 | 源码位置 |
| --- | --- | --- | --- |
| `wris_config` | 数值/工况检查判据及容差，原式中的比较值为实际门限 | `if nargin==0` | [W_RIS/wris_config.m:3](<C:/Users/Spidey/Desktop/code/W_RIS/wris_config.m:3>) |
| `wris_config` | 数值/工况检查判据及容差，原式中的比较值为实际门限 | `if strcmp(profile,'smoke')` | [W_RIS/wris_config.m:44](<C:/Users/Spidey/Desktop/code/W_RIS/wris_config.m:44>) |


**W_RIS/wris_egg.m**

| 作用域 | 含义 | 完整调用及内联参数 | 源码位置 |
| --- | --- | --- | --- |
| `wris_egg` | 随机数发生器类型和随机状态设置 | `rng(seed,'twister')` | [W_RIS/wris_egg.m:20](<C:/Users/Spidey/Desktop/code/W_RIS/wris_egg.m:20>) |
| `wris_egg` | 数值/工况检查判据及容差，原式中的比较值为实际门限 | `if strcmp(kind,'none')` | [W_RIS/wris_egg.m:27](<C:/Users/Spidey/Desktop/code/W_RIS/wris_egg.m:27>) |
| `wris_egg` | 数值/工况检查判据及容差，原式中的比较值为实际门限 | `if strcmp(kind,'weak')` | [W_RIS/wris_egg.m:39](<C:/Users/Spidey/Desktop/code/W_RIS/wris_egg.m:39>) |
| `wris_egg` | 数值/工况检查判据及容差，原式中的比较值为实际门限 | `if ng > 0` | [W_RIS/wris_egg.m:56](<C:/Users/Spidey/Desktop/code/W_RIS/wris_egg.m:56>) |
| `mustBeSeed` | 输入类型、维度、允许值及数值范围约束；原式列出全部界限 | `validateattributes(seed,{'numeric'}, {'scalar','real','finite','nonnegative','integer','<=',2^32-1})` | [W_RIS/wris_egg.m:75](<C:/Users/Spidey/Desktop/code/W_RIS/wris_egg.m:75>) |


**W_RIS/wris_geometry.m**

| 作用域 | 含义 | 完整调用及内联参数 | 源码位置 |
| --- | --- | --- | --- |
| `wris_geometry` | 数值/工况检查判据及容差，原式中的比较值为实际门限 | `if nargin<5` | [W_RIS/wris_geometry.m:3](<C:/Users/Spidey/Desktop/code/W_RIS/wris_geometry.m:3>) |
| `wris_geometry` | 数值/工况检查判据及容差，原式中的比较值为实际门限 | `if ismember(mode,{'joint','ris-only'})` | [W_RIS/wris_geometry.m:5](<C:/Users/Spidey/Desktop/code/W_RIS/wris_geometry.m:5>) |
| `wris_geometry` | 数值/工况检查判据及容差，原式中的比较值为实际门限 | `if ismember(mode,{'joint','tx-only'})` | [W_RIS/wris_geometry.m:10](<C:/Users/Spidey/Desktop/code/W_RIS/wris_geometry.m:10>) |
| `wris_geometry` | 数值/工况检查判据及容差，原式中的比较值为实际门限 | `if norm(localDirection(1:2))>1e-10` | [W_RIS/wris_geometry.m:24](<C:/Users/Spidey/Desktop/code/W_RIS/wris_geometry.m:24>) |
| `wris_geometry` | 数值/工况检查判据及容差，原式中的比较值为实际门限 | `if norm(directionNoYaw(1:2))>1e-10` | [W_RIS/wris_geometry.m:31](<C:/Users/Spidey/Desktop/code/W_RIS/wris_geometry.m:31>) |


**W_RIS/wris_mc.m**

| 作用域 | 含义 | 完整调用及内联参数 | 源码位置 |
| --- | --- | --- | --- |
| `wris_mc` | 数值/工况检查判据及容差，原式中的比较值为实际门限 | `if nargin<4` | [W_RIS/wris_mc.m:10](<C:/Users/Spidey/Desktop/code/W_RIS/wris_mc.m:10>) |
| `wris_mc` | 随机数发生器类型和随机状态设置 | `rng(seed,'twister')` | [W_RIS/wris_mc.m:12](<C:/Users/Spidey/Desktop/code/W_RIS/wris_mc.m:12>) |
| `wris_mc` | 输入类型、维度、允许值及数值范围约束；原式列出全部界限 | `validateattributes(N,{'numeric'},{'scalar','integer','positive'})` | [W_RIS/wris_mc.m:14](<C:/Users/Spidey/Desktop/code/W_RIS/wris_mc.m:14>) |
| `wris_mc` | 输入类型、维度、允许值及数值范围约束；原式列出全部界限 | `validateattributes(c.mc.chunkSize,{'numeric'},{'scalar','integer','positive'})` | [W_RIS/wris_mc.m:15](<C:/Users/Spidey/Desktop/code/W_RIS/wris_mc.m:15>) |
| `wris_mc` | 输入类型、维度、允许值及数值范围约束；原式列出全部界限 | `validateattributes([a b],{'numeric'},{'real','finite','nonnegative'})` | [W_RIS/wris_mc.m:16](<C:/Users/Spidey/Desktop/code/W_RIS/wris_mc.m:16>) |
| `wris_mc` | 输入类型、维度、允许值及数值范围约束；原式列出全部界限 | `validateattributes(c.water.anisotropy,{'numeric'},{'scalar','>',-1,'<',1})` | [W_RIS/wris_mc.m:17](<C:/Users/Spidey/Desktop/code/W_RIS/wris_mc.m:17>) |
| `wris_mc` | 输入类型、维度、允许值及数值范围约束；原式列出全部界限 | `validateattributes(planeScale,{'numeric'},{'scalar','integer','>=',1})` | [W_RIS/wris_mc.m:24](<C:/Users/Spidey/Desktop/code/W_RIS/wris_mc.m:24>) |
| `wris_mc` | 数值/工况检查判据及容差，原式中的比较值为实际门限 | `if ext>0` | [W_RIS/wris_mc.m:45](<C:/Users/Spidey/Desktop/code/W_RIS/wris_mc.m:45>) |
| `wris_mc` | 数值/工况检查判据及容差，原式中的比较值为实际门限 | `if ns==0` | [W_RIS/wris_mc.m:90](<C:/Users/Spidey/Desktop/code/W_RIS/wris_mc.m:90>) |
| `wris_mc` | 数值/工况检查判据及容差，原式中的比较值为实际门限 | `if c.water.anisotropy==0.924` | [W_RIS/wris_mc.m:91](<C:/Users/Spidey/Desktop/code/W_RIS/wris_mc.m:91>) |
| `wris_mc` | 数值/工况检查判据及容差，原式中的比较值为实际门限 | `if abs(gg)<1e-12` | [W_RIS/wris_mc.m:95](<C:/Users/Spidey/Desktop/code/W_RIS/wris_mc.m:95>) |


**W_RIS/wris_metrics.m**

| 作用域 | 含义 | 完整调用及内联参数 | 源码位置 |
| --- | --- | --- | --- |
| `wris_metrics` | 数值/工况检查判据及容差，原式中的比较值为实际门限 | `if isempty(w) \|\| sum(w)<=0` | [W_RIS/wris_metrics.m:5](<C:/Users/Spidey/Desktop/code/W_RIS/wris_metrics.m:5>) |


**W_RIS/wris_plot.m**

| 作用域 | 含义 | 完整调用及内联参数 | 源码位置 |
| --- | --- | --- | --- |
| `wris_plot` | 绘图数据与颜色、线宽、透明度、标记或曲线样式参数 | `imagesc(ax,o.x,o.y,max(o.noTurbulence,limits(1)).')` | [W_RIS/wris_plot.m:13](<C:/Users/Spidey/Desktop/code/W_RIS/wris_plot.m:13>) |
| `wris_plot` | 坐标范围、色轴、刻度、图框比例、视角或等比例显示参数 | `axis(ax,'xy','equal')` | [W_RIS/wris_plot.m:14](<C:/Users/Spidey/Desktop/code/W_RIS/wris_plot.m:14>) |
| `wris_plot` | 坐标范围、色轴、刻度、图框比例、视角或等比例显示参数 | `xlim(ax,[-extent(1) extent(1)])` | [W_RIS/wris_plot.m:16](<C:/Users/Spidey/Desktop/code/W_RIS/wris_plot.m:16>) |
| `wris_plot` | 坐标范围、色轴、刻度、图框比例、视角或等比例显示参数 | `ylim(ax,[-extent(2) extent(2)])` | [W_RIS/wris_plot.m:16](<C:/Users/Spidey/Desktop/code/W_RIS/wris_plot.m:16>) |
| `wris_plot` | 图形标签、字体、颜色、网格与布局等显示参数 | `set(ax,'ColorScale','log','Color',[0.78 0.78 0.78],'FontSize',8)` | [W_RIS/wris_plot.m:17](<C:/Users/Spidey/Desktop/code/W_RIS/wris_plot.m:17>) |
| `wris_plot` | 图形标签、字体、颜色、网格与布局等显示参数 | `colormap(ax,turbo(256))` | [W_RIS/wris_plot.m:18](<C:/Users/Spidey/Desktop/code/W_RIS/wris_plot.m:18>) |
| `wris_plot` | 此源码调用的内联选项/固定常数；具体实参及位置完整列于右侧 | `hold(ax,'on')` | [W_RIS/wris_plot.m:18](<C:/Users/Spidey/Desktop/code/W_RIS/wris_plot.m:18>) |
| `wris_plot` | 绘图数据与颜色、线宽、透明度、标记或曲线样式参数 | `rectangle(ax,'Position',[-c.ris.Lx/2 -c.ris.Ly/2 c.ris.Lx c.ris.Ly], 'EdgeColor','w','LineWidth',1.2)` | [W_RIS/wris_plot.m:19](<C:/Users/Spidey/Desktop/code/W_RIS/wris_plot.m:19>) |
| `wris_plot` | 数值/工况检查判据及容差，原式中的比较值为实际门限 | `if col>1` | [W_RIS/wris_plot.m:24](<C:/Users/Spidey/Desktop/code/W_RIS/wris_plot.m:24>) |
| `wris_plot` | 图形标签、字体、颜色、网格与布局等显示参数 | `title(ax,{label,sprintf('yaw=%+.1f deg \| rotation=%+.1f deg',o.geometry.yawRisDeg,delta), sprintf('incidence=%.1f deg \| core clipped=%.1f%%', o.geometry.thetaDeg,100*o.coreClippedFraction)},'FontSize',8,'FontWeight','normal')` | [W_RIS/wris_plot.m:28](<C:/Users/Spidey/Desktop/code/W_RIS/wris_plot.m:28>) |
| `wris_plot` | 数值/工况检查判据及容差，原式中的比较值为实际门限 | `if col==1` | [W_RIS/wris_plot.m:31](<C:/Users/Spidey/Desktop/code/W_RIS/wris_plot.m:31>) |
| `wris_plot` | 图形标签、字体、颜色、网格与布局等显示参数 | `ylabel(ax,sprintf('alpha_0 = %g deg\nRIS local y (m)',p.initialIncidenceDeg))` | [W_RIS/wris_plot.m:31](<C:/Users/Spidey/Desktop/code/W_RIS/wris_plot.m:31>) |
| `wris_plot` | 数值/工况检查判据及容差，原式中的比较值为实际门限 | `if ~o.geometry.axisInside && any(abs(o.geometry.axisLocal(1:2))>extent)` | [W_RIS/wris_plot.m:32](<C:/Users/Spidey/Desktop/code/W_RIS/wris_plot.m:32>) |
| `wris_plot` | 图形标签、字体、颜色、网格与布局等显示参数 | `title(layout,{sprintf('RIS rotation and clipping \| col 1 fixed; cols 2-5 joint TX/RIS motion \| elapsed t=%.1f s, sea t=%.1f s', c.sampling.comparisonTime,c.ris.startTime+c.sampling.comparisonTime), 'White dashed ellipse: complete incident core (MC D4sigma); dotted: static axis; white frame: RIS; grey: outside aperture'}, 'FontSize',11)` | [W_RIS/wris_plot.m:38](<C:/Users/Spidey/Desktop/code/W_RIS/wris_plot.m:38>) |
| `wris_plot` | 图形标签、字体、颜色、网格与布局等显示参数 | `xlabel(layout,'RIS local x (m); rotation is measured relative to the static reference in the same row')` | [W_RIS/wris_plot.m:42](<C:/Users/Spidey/Desktop/code/W_RIS/wris_plot.m:42>) |
| `wris_plot` | 输出文件名、保存格式、编码或图像分辨率选项 | `savefig(f,fullfile(c.output,'Fig_RIS_rotation_clipping.fig'))` | [W_RIS/wris_plot.m:46](<C:/Users/Spidey/Desktop/code/W_RIS/wris_plot.m:46>) |
| `wris_plot` | 输出文件名、保存格式、编码或图像分辨率选项 | `exportgraphics(f,fullfile(c.output,'Fig_RIS_rotation_clipping.png'),'Resolution',c.plot.resolution)` | [W_RIS/wris_plot.m:47](<C:/Users/Spidey/Desktop/code/W_RIS/wris_plot.m:47>) |
| `drawCore` | 绘图数据与颜色、线宽、透明度、标记或曲线样式参数 | `plot(ax,ellipse(:,1),ellipse(:,2),'w--','LineWidth',0.9)` | [W_RIS/wris_plot.m:56](<C:/Users/Spidey/Desktop/code/W_RIS/wris_plot.m:56>) |
| `drawCore` | 绘图数据与颜色、线宽、透明度、标记或曲线样式参数 | `plot(ax,m.center(1)+[-baseline(1) baseline(1)],m.center(2)+[-baseline(2) baseline(2)], ':','Color',[0.25 0.25 0.25],'LineWidth',1)` | [W_RIS/wris_plot.m:59](<C:/Users/Spidey/Desktop/code/W_RIS/wris_plot.m:59>) |
| `drawCore` | 绘图数据与颜色、线宽、透明度、标记或曲线样式参数 | `plot(ax,m.center(1)+[-axisDirection(1) axisDirection(1)], m.center(2)+[-axisDirection(2) axisDirection(2)],'w-','LineWidth',0.8)` | [W_RIS/wris_plot.m:62](<C:/Users/Spidey/Desktop/code/W_RIS/wris_plot.m:62>) |
| `drawCore` | 绘图数据与颜色、线宽、透明度、标记或曲线样式参数 | `plot(ax,o.geometry.axisLocal(1),o.geometry.axisLocal(2),'w+','MarkerSize',6,'LineWidth',1)` | [W_RIS/wris_plot.m:64](<C:/Users/Spidey/Desktop/code/W_RIS/wris_plot.m:64>) |
| `incomingInset` | 绘图数据与颜色、线宽、透明度、标记或曲线样式参数 | `imagesc(ax,o.plane.x,o.plane.y,max(o.plane.noTurbulence,limits(1)).')` | [W_RIS/wris_plot.m:72](<C:/Users/Spidey/Desktop/code/W_RIS/wris_plot.m:72>) |
| `incomingInset` | 坐标范围、色轴、刻度、图框比例、视角或等比例显示参数 | `axis(ax,'xy','equal')` | [W_RIS/wris_plot.m:73](<C:/Users/Spidey/Desktop/code/W_RIS/wris_plot.m:73>) |
| `incomingInset` | 此源码调用的内联选项/固定常数；具体实参及位置完整列于右侧 | `hold(ax,'on')` | [W_RIS/wris_plot.m:73](<C:/Users/Spidey/Desktop/code/W_RIS/wris_plot.m:73>) |
| `incomingInset` | 坐标范围、色轴、刻度、图框比例、视角或等比例显示参数 | `xlim(ax,center(1)+[-half half])` | [W_RIS/wris_plot.m:78](<C:/Users/Spidey/Desktop/code/W_RIS/wris_plot.m:78>) |
| `incomingInset` | 坐标范围、色轴、刻度、图框比例、视角或等比例显示参数 | `ylim(ax,center(2)+[-half half])` | [W_RIS/wris_plot.m:78](<C:/Users/Spidey/Desktop/code/W_RIS/wris_plot.m:78>) |
| `incomingInset` | 图形标签、字体、颜色、网格与布局等显示参数 | `set(ax,'ColorScale','log','FontSize',6)` | [W_RIS/wris_plot.m:79](<C:/Users/Spidey/Desktop/code/W_RIS/wris_plot.m:79>) |
| `incomingInset` | 图形标签、字体、颜色、网格与布局等显示参数 | `colormap(ax,turbo(256))` | [W_RIS/wris_plot.m:79](<C:/Users/Spidey/Desktop/code/W_RIS/wris_plot.m:79>) |
| `incomingInset` | 绘图数据与颜色、线宽、透明度、标记或曲线样式参数 | `rectangle(ax,'Position',[-c.ris.Lx/2 -c.ris.Ly/2 c.ris.Lx c.ris.Ly],'EdgeColor','w','LineWidth',1)` | [W_RIS/wris_plot.m:80](<C:/Users/Spidey/Desktop/code/W_RIS/wris_plot.m:80>) |
| `incomingInset` | 绘图数据与颜色、线宽、透明度、标记或曲线样式参数 | `plot(ax,o.geometry.axisLocal(1),o.geometry.axisLocal(2),'w+','MarkerSize',5)` | [W_RIS/wris_plot.m:81](<C:/Users/Spidey/Desktop/code/W_RIS/wris_plot.m:81>) |
| `incomingInset` | 图形标签、字体、颜色、网格与布局等显示参数 | `title(ax,{'Incoming beam / virtual plane','Core misses the finite RIS'},'FontSize',7)` | [W_RIS/wris_plot.m:82](<C:/Users/Spidey/Desktop/code/W_RIS/wris_plot.m:82>) |


**W_RIS/wris_report.m**

| 作用域 | 含义 | 完整调用及内联参数 | 源码位置 |
| --- | --- | --- | --- |
| `wris_report` | 输出文件名、保存格式、编码或图像分辨率选项 | `writetable(summary,fullfile(c.output,'snapshots.csv'))` | [W_RIS/wris_report.m:30](<C:/Users/Spidey/Desktop/code/W_RIS/wris_report.m:30>) |
| `wris_report` | 输出文件名、保存格式、编码或图像分辨率选项 | `writetable(array2table(vertcat(motionRows{:}),'VariableNames',names),fullfile(c.output,'motion.csv'))` | [W_RIS/wris_report.m:33](<C:/Users/Spidey/Desktop/code/W_RIS/wris_report.m:33>) |


**W_RIS/wris_run.m**

| 作用域 | 含义 | 完整调用及内联参数 | 源码位置 |
| --- | --- | --- | --- |
| `wris_run` | 数值/工况检查判据及容差，原式中的比较值为实际门限 | `if nargin==0` | [W_RIS/wris_run.m:3](<C:/Users/Spidey/Desktop/code/W_RIS/wris_run.m:3>) |
| `wris_run` | 数值/工况检查判据及容差，原式中的比较值为实际门限 | `assert(all([tests.Passed]),'wris:Tests','Physical regression tests failed.')` | [W_RIS/wris_run.m:7](<C:/Users/Spidey/Desktop/code/W_RIS/wris_run.m:7>) |
| `wris_run` | 输出文件名、保存格式、编码或图像分辨率选项 | `writetable(table(tests),fullfile(c.output,'tests.csv'))` | [W_RIS/wris_run.m:8](<C:/Users/Spidey/Desktop/code/W_RIS/wris_run.m:8>) |
| `wris_run` | 输出文件名、保存格式、编码或图像分辨率选项 | `save(fullfile(c.output,'wris_results.mat'),'s','-v7.3')` | [W_RIS/wris_run.m:11](<C:/Users/Spidey/Desktop/code/W_RIS/wris_run.m:11>) |


**W_RIS/wris_setup.m**

| 作用域 | 含义 | 完整调用及内联参数 | 源码位置 |
| --- | --- | --- | --- |
| `wris_setup` | 源码路径或测试文件入口 | `addpath(root,fullfile(fileparts(root),'sea_wave'),fullfile(fileparts(root),'fig3-7'))` | [W_RIS/wris_setup.m:4](<C:/Users/Spidey/Desktop/code/W_RIS/wris_setup.m:4>) |


**W_RIS/wris_simulate.m**

| 作用域 | 含义 | 完整调用及内联参数 | 源码位置 |
| --- | --- | --- | --- |
| `wris_simulate` | 数值/工况检查判据及容差，原式中的比较值为实际门限 | `assert(t>=0 && t<=c.ris.duration,'Comparison time must lie inside the simulated interval.')` | [W_RIS/wris_simulate.m:4](<C:/Users/Spidey/Desktop/code/W_RIS/wris_simulate.m:4>) |
| `wris_simulate` | 数值/工况检查判据及容差，原式中的比较值为实际门限 | `assert(abs(m.t(k)-t)<1e-10,'Comparison time must coincide with a saved motion sample.')` | [W_RIS/wris_simulate.m:26](<C:/Users/Spidey/Desktop/code/W_RIS/wris_simulate.m:26>) |


**W_RIS/wris_sites.m**

| 作用域 | 含义 | 完整调用及内联参数 | 源码位置 |
| --- | --- | --- | --- |
| `wris_sites` | 数值/工况检查判据及容差，原式中的比较值为实际门限 | `assert(~isempty(crest) && ~isempty(trough) && ~isempty(rising) && ~isempty(falling), 'wris:Sites','Spatial transect lacks four sea positions; increase the transect span.')` | [W_RIS/wris_sites.m:18](<C:/Users/Spidey/Desktop/code/W_RIS/wris_sites.m:18>) |
| `wris_sites` | 数值/工况检查判据及容差，原式中的比较值为实际门限 | `if j>=3` | [W_RIS/wris_sites.m:27](<C:/Users/Spidey/Desktop/code/W_RIS/wris_sites.m:27>) |


**W_RIS/wris_tx_motion.m**

| 作用域 | 含义 | 完整调用及内联参数 | 源码位置 |
| --- | --- | --- | --- |
| `wris_tx_motion` | 数值/工况检查判据及容差，原式中的比较值为实际门限 | `if ~strcmp(c.tx.pointing,'body-fixed')` | [W_RIS/wris_tx_motion.m:20](<C:/Users/Spidey/Desktop/code/W_RIS/wris_tx_motion.m:20>) |


**W_RIS/wris_validate.m**

| 作用域 | 含义 | 完整调用及内联参数 | 源码位置 |
| --- | --- | --- | --- |
| `wris_validate` | 数值/工况检查判据及容差，原式中的比较值为实际门限 | `assert(v.physicalAssertionsPassed,'wris:Validation','Projection, rotation, energy or convergence check failed.')` | [W_RIS/wris_validate.m:41](<C:/Users/Spidey/Desktop/code/W_RIS/wris_validate.m:41>) |
| `wris_validate` | 输出文件名、保存格式、编码或图像分辨率选项 | `writetable(v.convergence,fullfile(c.output,'convergence.csv'))` | [W_RIS/wris_validate.m:42](<C:/Users/Spidey/Desktop/code/W_RIS/wris_validate.m:42>) |


**fig3-7/Table_abc.m**

| 作用域 | 含义 | 完整调用及内联参数 | 源码位置 |
| --- | --- | --- | --- |
| `Table_abc` | 数值/工况检查判据及容差，原式中的比较值为实际门限 | `if strcmp(water_type,'c3.2')` | [fig3-7/Table_abc.m:16](<C:/Users/Spidey/Desktop/code/fig3-7/Table_abc.m:16>) |


**fig3-7/Hen_Green.m**

| 作用域 | 含义 | 完整调用及内联参数 | 源码位置 |
| --- | --- | --- | --- |
| `Hen_Green` | 数值/工况检查判据及容差，原式中的比较值为实际门限 | `if judge == 1` | [fig3-7/Hen_Green.m:2](<C:/Users/Spidey/Desktop/code/fig3-7/Hen_Green.m:2>) |


## 7. 文件覆盖核对

以下校验针对当前磁盘源码的静态参数清单；本任务没有重新运行海浪动力学或百万光子仿真。当前 `wris_run.m` 引用的 `W_RIS/tests/WrisTest.m` 在目录中不存在，因此没有从缺失测试文件推定或补造任何参数。已有 `results` 下的运行结果不作为当前默认值来源。

| 源码文件 | 源码行数 | 具名符号记录 | SHA-256（前 16 位） |
| --- | --- | --- | --- |
| `sea_wave/audit_sea_spectrum.m` | 48 | 32 | `444eae9b35a88e8f` |
| `sea_wave/get_default_RIS_params.m` | 12 | 31 | `999c22b354b5e9b2` |
| `sea_wave/get_default_sea_params.m` | 38 | 18 | `8c2988d59f0f34a9` |
| `sea_wave/main3_windseasurface_0326.m` | 34 | 33 | `2b032362601d36b1` |
| `sea_wave/oceanSurfaceModel_phiwind_260315.m` | 140 | 90 | `0888cd99d15d5932` |
| `sea_wave/RIS_draw_state.m` | 44 | 15 | `ac82f8ab26dcd4a4` |
| `sea_wave/RIS_float.m` | 516 | 423 | `0a9b11f7dabc49c2` |
| `sea_wave/RIS_radiation.m` | 39 | 31 | `462b09d847c603ae` |
| `sea_wave/RIS_rotation.m` | 22 | 10 | `685db31fcde841e3` |
| `sea_wave/RIS_sea_surface.m` | 9 | 7 | `d50b9ac5286231d2` |
| `sea_wave/RIS_sidewall_load.m` | 41 | 28 | `870f5e95ad5d7ad9` |
| `sea_wave/RIS_wave_components.m` | 30 | 30 | `12b45ce3565dcfd1` |
| `sea_wave/RIS_wave_field.m` | 28 | 12 | `cdd50b9d8cb3c83b` |
| `sea_wave/run_RIS_wind_comparison.m` | 135 | 83 | `74f57157be181b26` |
| `sea_wave/summarize_RIS_motion_existing.m` | 74 | 48 | `60232f455cc8efd0` |
| `W_RIS/run_wris.m` | 3 | 1 | `0d9be0d92c44c720` |
| `W_RIS/wris_cached.m` | 45 | 19 | `2464fe278ce822d5` |
| `W_RIS/wris_config.m` | 50 | 60 | `23ab085433cbeb49` |
| `W_RIS/wris_egg.m` | 78 | 43 | `f0cbab92485977de` |
| `W_RIS/wris_geometry.m` | 34 | 27 | `850d71e0b77f33ac` |
| `W_RIS/wris_intersect.m` | 12 | 11 | `70a0f419c9e68d1e` |
| `W_RIS/wris_mc.m` | 146 | 139 | `9086eaaffe57fce7` |
| `W_RIS/wris_metrics.m` | 11 | 17 | `c465e56d9f261e51` |
| `W_RIS/wris_plot.m` | 84 | 51 | `cdeaa036e1f06671` |
| `W_RIS/wris_reference.m` | 19 | 20 | `c7514ec42208f6fb` |
| `W_RIS/wris_report.m` | 59 | 45 | `ede123bb561d6b97` |
| `W_RIS/wris_run.m` | 14 | 6 | `9d8218f430b06148` |
| `W_RIS/wris_setup.m` | 5 | 1 | `4f72a92e3f7bec66` |
| `W_RIS/wris_simulate.m` | 44 | 51 | `d7368c9865dbe457` |
| `W_RIS/wris_sites.m` | 35 | 36 | `b625782f5f000989` |
| `W_RIS/wris_tx_motion.m` | 26 | 25 | `08031fb6b40e5e9c` |
| `W_RIS/wris_validate.m` | 49 | 45 | `2b8c83083280e2af` |
| `fig3-7/Table_abc.m` | 47 | 18 | `35b5e43e76214de6` |
| `fig3-7/Hen_Green.m` | 11 | 5 | `3087566fb7dce1f8` |


本清单合计覆盖 32 个目标目录源码文件、2 个直接依赖源码文件、1511 条按函数作用域区分的具名符号记录、234 条额外内联调用记录。所有具名符号均有含义、取值或计算来源、源码定位；程序核对未发现未说明条目。源码变更后应重新核对；源码行号和校验哈希对应本次生成时的文件。

可复核生成依据：`研究执行/参数清单核验/build_parameter_inventory.py`（提取）、`parameter_meanings.py`（人工核对的含义字典）、`render_parameter_document.py`（生成）、`inventory.json`（原始清单）、`coverage_audit.json`（覆盖结果）。
