# 冗余代码清理

保留 `README.md` 原始要求。W_RIS、BSF的源码及原数值结果保持只读。
保留原 `rp_phase` 作为实际量化的 +z 基准，以及BSF、MC、分层光路、CZT、像素
相干积分等共用光学核心。`rp_run/rp_plot/rp_report/rp_validate` 公共入口已转向新任务。

删除18个由新流程替代、且工作区其他MATLAB源码没有引用的旧文件：

- 原海浪快照流程：`rp_audit`, `rp_simulate`, `rp_snapshot`, `rp_rx_motion`。
- 旧阈值覆盖圆/局部加密流程：`rp_footprint`, `rp_refine_receiver`。
- 无标定操作的兼容空包装：`rp_calibrate_input_field`。
- 旧二次/tanh模型比较及依赖旧缓存的检查：`rp_compare_fresnel_line`,
  `rp_compare_models`, `rp_compare_wave`, `rp_plot_phase_comparison`,
  `rp_select_stress_unit`, `rp_verify_models`, `run_ris_phase_comparison`。
- 与当前接收性能目标分离的旧全半球扫描：`rp_direction_classes`,
  `rp_reachability`, `rp_run_reachability`, `run_ris_reachability`。

这些文件的原路径、行数、函数依赖、SHA256均在首次审计清单中。
三份旧推导文档保留为历史资料，顶部标注本次新模型说明的入口。

2026-11-12公式与实现一致性复核：另删除 `rp_rx_complex_field.m`。
该函数只对调用者额外提供的模式数组做线性叠加，工作区MATLAB源码无调用，
没有下表面到LC的散射传播或互相干模式求解；保留它容易误解当前主流程能力。
相应历史文档已明确区分理论模式扩展与实际混合传播。未删除仍用于查询实际
弹道调制场的 `rp_modulated_field`，也未改动共用数值传播器。
