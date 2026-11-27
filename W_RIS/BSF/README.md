# RIS 下表面 BSF 与原 MC 对比

全部新增代码和产物均在本目录内，原 `W_RIS`、`sea_wave`、`3yimiao/fig3-10` 代码和原 MC 数据不改动。

先读 [FORMULAS.md](FORMULAS.md)：给出式（3-18）的显式 \(I_b+I_s\) 分解、完整 HG 输运修正、误差定义及后续调相接口的假设。结果完成后读 `results/REPORT.md`、`results/comparison.csv`，或 `BSF_note.pdf` 中的公式与验证说明。

## 运行

在 MATLAB R2024b 中：

```matlab
cd('C:/Users/Spidey/Desktop/code/W_RIS/BSF')
run_bsf
```

默认读取 `../results/wris_results.mat` 的 18 幅原始快照，在完全相同的 100 × 100 网格上计算。每批 800 万个强制散射样本，按场景自适应增加至整体二维标准误差 ≤0.1%。每批保存，可直接重新运行续算。无需重新求解海浪运动。

在原工程的单程、首次到达模型下，同一冻结快照的 RIS 相位优化可复用已计算的下表面 `Ib` 和散射入射光子集合；相位迭代不必重复水下入射输运。改变几何、光源或水质后应重算输入。

`bsf_run_all(TargetRelativeSE=0.001,BatchPhotons=8e6,MaxPhotons=1024e6)` 可修改数值目标。达到样本上限却未达标时，CSV 的 `precisionPassed=false`，不会伪报通过。若物理配置或算法改变，另建结果目录或归档已有 `results` 后重新计算，避免混用旧缓存。

快速试算：`bsf_pilot`（首幅、200 万散射样本）；物理测试：`runtests('bsfPhysicsTest')`；重画图：`bsf_plot_results`。正式计算仅需基础 MATLAB 及原工程模块；`bsf_qmc_pilot` 是可选随机 Sobol 实验，需要 Statistics and Machine Learning Toolbox，正式结果未使用它。

`run_bsf` 在主计算后自动执行 `bsf_validate_outputs`，保存独立种子复验和首场景的调相输入包。可选的独立 Python 数值审核命令如下（需要 numpy、h5py；不影响仅用 MATLAB 计算）：

```powershell
python bsf_independent_ballistic.py
```

## 代码入口

| 文件 | 功能 |
|---|---|
| `bsf_ballistic.m` | 与原高斯源一致的零散射辐照度 |
| `bsf_paper.m` / `bsf_paper_map.m` | 未拟合的论文 BSF，分离弹道/散射 |
| `bsf_scatter_transport.m` | 完整 HG 精确路径积分的独立、无偏数值估计 |
| `bsf_run_all.m` | 18 场景自适应精度、逐批保存、断点续算 |
| `bsf_compare.m` / `bsf_export_summary.m` | 原始像素、功率、空间分布及噪声比较 |
| `bsf_plot_results.m` | 相同色标热力图、未平滑残差与精度总览 |
| `bsf_phase_input.m` | 弹道振幅及假设的高斯入射相位接口 |
| `bsfPhysicsTest.m` | 真空归一化、Beer-Lambert、刚体坐标不变性等测试 |

## 单幅使用

```matlab
bsf_setup;
d=load('../results/wris_results.mat','s');
c=d.s.config;
o=d.s.panels{1,1}.out;
g=o.geometry;
Ib=bsf_cell_average(@(x,y) bsf_ballistic(c,g,x,y),o.x,o.y,6);
sc=bsf_scatter_transport(c,g,Photons=8e6,Seed=12345);
IBSF=Ib+sc.Is;
imagesc(o.x,o.y,IBSF.'); axis image; set(gca,'YDir','normal'); colorbar
```

矩阵约定与原 MC 一致：`I(ix,iy)`，画图需转置。`bsf_ballistic` 的 x、y 必须是相同形状的坐标数组；不能直接传入两个坐标向量而期望它内部形成二维网格。

## 调相所需的两个独立输入

```matlab
[X,Y]=ndgrid(o.x,o.y);
field=bsf_phase_input(c,g,X,Y,zeros(size(X))); % 点采样复场；不是单元平均强度
[sc,packets]=bsf_scatter_transport(c,g,Photons=8e6,Seed=23456,SavePackets=true);
% field 用于检查连续相位公式；packets -> RIS 散射出射方向模型/光线传播器
```

原 1 cm 网格是强度统计网格，不能直接用作 532 nm 光学复场的 FFT 传播网格。倾斜入射载波本身的 Nyquist 间距可小至约 0.2–0.4 μm；还需检查高斯曲率和 RIS 相位梯度。函数返回 `carrierNyquistPitchM`，并明确标记 `gridPhaseSamplingCertified=false`。后续应使用连续相位重新采样，或在经过验证的倾斜载波/相位补偿模型中传播包络。

复场接口已施加真实有限 RIS 的孔径掩模 `aperture`，扩展网格的板外复振幅为零；`bsf_ballistic` 本身表示无限计数平面的入射密度，正式 RIS 比较只在真实板面单元上积分。

若只想使用原 MC 的散射光子方向与权重，也可只读重跑原实现：

```matlab
[originalMC,allPackets]=wris_mc(c,g,o.seed,1);
mask=~allPackets.ballistic;
scatterPos=allPackets.pos(mask,:);
scatterDir=allPackets.dir(mask,:);
scatterPowerW=allPackets.weight(mask)*c.source.power/allPackets.photons;
```

正式高精度面强度与原 MC 光子包属于不同的独立样本集合，不要把一套的强度与另一套的样本权重强行逐像素相等。图示的平均强度相加依赖非相干散射集合假设；复场相位、RIS 出射方向规则和 RX 传播仍需各自验证。
