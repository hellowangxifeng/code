# RIS_dynamic：冻结粗糙海面下的浮动 RIS 与无 RIS 直达链路对比（v4）

算法、公式、每张图的仿真条件见 **[算法与图件说明.md](算法与图件说明.md)**，本文件只说明怎么运行、文件放在哪里。

## 运行

在 MATLAB（R2022a 及以上，需要 Statistics and Machine Learning Toolbox 和 Signal Processing Toolbox）中打开 `main_ris_dynamic.m`，按 F5 运行。

| 开关 | 作用 |
|---|---|
| `profile = 'smoke'` | 两个高度、一种风速、很少的光子，几分钟内跑通整条流程，结果写到 `results/smoke` |
| `profile = 'production'` | 论文用设置 |
| `doCompute = false` | 不重新计算，只用已保存的 `results/ris_dynamic.mat` 重画图 |
| `c.run.parallel = true` | 对同一条件下的多个设计开 `parfor`（需要 Parallel Computing Toolbox） |

- **断点续算**：每个“风速 × 高度”条件、每个 RIS 设计都有独立 checkpoint（`results/cache`），按输入参数和相关源码的 SHA-256 判断是否复用。中断后再按 F5 会从断点继续。
- **改动后哪些会重算**：只改湍流、功率范围、选方案规则或画图参数，不会重算光学部分。

## 目录

```
RIS_dynamic/
├─ main_ris_dynamic.m        一键入口（参数 → 计算 → 画图）
├─ algorithm/                计算（不画图）
│   rd_config / rd_check / rd_run / rd_setup
│   rd_scene                 冻结时刻：RIS 位姿、发射机重新指向、粗糙海面、MC 输入
│   rd_motion                sea_wave/RIS_float 六自由度运动（缓存）
│   rd_rx_direction          RIS 可达锥内的公共接收方向
│   rd_condition             单个“风速 × 高度”条件：公共 RX、直达 MC、全部 RIS 设计
│   rd_candidate             单个 RIS 设计：复光场 + MC 校验
│   rd_phase                 Steering-only / Flat top / Robust 相位（瞄准公共 RX）
│   rd_wave, rd_pixel_intensity   弹道复光场（RIS_phase rp_wave）+ 散射光
│   rd_trace, rd_ris_packets      RIS 链路全 MC（几何光学）
│   rd_direct_mc, rd_surface, rd_crossing, rd_to_plane   无 RIS 直达 MC
│   rd_receiver, rd_capture, rd_compress, rd_map          接收孔径与光强图
│   rd_link_metrics, rd_select, rd_rank, rd_evaluate, rd_tables   统计、选方案、排名、导出
│   rd_cached, rd_cache_valid, rd_code_hash, rd_physical, rd_transport
├─ plots/                    画图（只读 ris_dynamic.mat）
│   plot_ris_dynamic, rd_plot_options, rd_fig_*.m
├─ tests/RdPhysicsTest.m     快速物理/记账单元测试
└─ results/                  运行后生成：ris_dynamic.mat、*.csv、figures/、cache/
```

`RIS_phase`、`W_RIS`、`sea_wave`、`yimiao/fig3-7` 只被调用，从不写入。缓存和结果都放在 `RIS_dynamic/results`。

## 输出

| 文件 | 内容 |
|---|---|
| `results/ris_dynamic.mat` | 变量 `R`：每个条件的位姿、公共 RX、五条链路的曲线（四组湍流）、MC 校验曲线、选中的设计、排名、光强图 |
| `results/designs.csv` | 每个条件下每个 RIS 设计的光斑统计、LC 程序、复光场和 MC 的平均增益、数值诊断 |
| `results/validation.csv` | 设计功率下各方案选中的设计：复光场与 MC 的对比 |
| `results/ranking.csv` | 达到 ABER/OP 目标所需的发射功率、最优方案、领先幅度 |
| `results/curves.csv` | 主湍流下五条链路的 ABER、OP、平均接收功率，以及每个功率点用到的设计编号 |
| `results/figures/fig1…fig7` | 见说明文档第 4 节 |

## 测试

```matlab
cd RIS_dynamic
r = runtests('tests/RdPhysicsTest.m'); table(r)
```
