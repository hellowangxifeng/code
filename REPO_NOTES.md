# 仓库说明

本仓库同步自本地工作目录 `C:\Users\Spidey\Desktop\code`，用于备份和版本管理
**代码、文档与小体积图片**。远程地址：
<https://github.com/hellowangxifeng/code>（私有仓库）。

## 已同步

- MATLAB 代码（`.m`）与 Python 代码（`.py`）
- 文档（`.md`、`.tex`、`.bib`、`.txt` 等）
- 小体积图片与示意图（`.png`、`.eps`、`.jpg`、`.vsdx`）
- 配置与小型表格（`.json`、`.yml`、`.xlsx`）

## 未同步（保留在本地）

以下内容由 `.gitignore` 排除，体积约 3.6 GB：

| 类型 | 体积 | 说明 |
|---|---:|---|
| `*.mat` | ~2.7 GB | 仿真/实验数据，最大单文件 897 MB |
| `*.csv` | ~440 MB | 原始与中间数据表 |
| `*.pdf` | ~205 MB | 论文、参考文献 |
| `*.fig` / `*.pickle` / `*.npz` | ~72 MB | 中间结果 |
| LaTeX 编译产物、缓存、临时文件 | — | 可重新生成 |

这样仓库保持在约 56 MB，不会触发 GitHub 单文件 100 MB 的限制。
需要完整数据时，请使用本地目录或另行做冷备份。

## 项目一览

| 目录 | 内容 |
|---|---|
| `FOWT_RIS_eval/` | 浮式平台六自由度实测数据驱动的 RIS 链路评估与改进 |
| `RIS_phase/` | RIS 相位算法、光斑与链路仿真（含 W_RIS） |
| `sea_wave/` | 海面波场、浮框配平与 RIS 风速/风向分析 |
| `W_RIS/` | RIS 下表面 BSF 与原始 MC 对比 |
| `RIS_dynamic/` | 动态海况下浮动 RIS 与粗糙面直达链路 MC 对照 |
| `RIS_model/` | RIS 器件参数与设计说明 |
| `RIS_ML/` | 机器学习辅助的相位/调控策略方案 |
| `date_set/` | 海洋浮体运动预测数据集（SKLOE） |
| `yimiao/` | 独立子项目 |
| `对比文献/` | 参考文献（PDF 未入库） |
| `毕设论文/` | 电子科技大学学位论文 LaTeX 工程（PDF 未入库） |

## 日常同步

```powershell
cd C:\Users\Spidey\Desktop\code
git add -A
git commit -m "说明这次改了什么"
git push
```

拉取远程更新：

```powershell
git pull --rebase
```
