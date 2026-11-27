> Fig2/3的本表旧SHA与旧场景为历史记录，已由共同固定接收端修订替代；当前验收见`plots/stage2/fixed_receiver/VISUAL_QA.md`和相邻artifact_validation.json。其他图仍按其原场景解释。

> 后续修订：Fig1已由`../fig1_revision/VISUAL_QA.md`和对应manifest替代；Fig2/3后续固定接收机版本见`../fixed_receiver/`。本文件其余内容保留为初版历史记录。

# 第二阶段最终PNG视觉复核

状态：passed。以下为对最终PNG实际查看的记录；文件SHA绑定本次复核版本。
保存后重新打开11份FIG的检查通过：52个坐标区刻度10 pt、32个图例均在坐标区四角，密集线段遮挡样本数全部为0。
全部图号、轴标签、色标/曲线示例可读；无Fig7。

|图|实际查看结论|最终PNG SHA-256|
|---|---|---|
|Fig1_irradiance|12热图，红圈完整，椭圆十字清楚，图例避开圈/轴；统一色标尺度。|`3cdcc31de972991842bb4583ada520d2d1371b24f4a59be1f7c1043d4297a429`|
|Fig2_aber|0–100 dBm、10^0刻度与左上低功率平台完整；130 mrad，TNR10，2.67×2.37 inch；FOV floor可见，图例零遮挡。|`68a141d0819df9215112be875e59136254e5149c55f44009f82b12b4d404ee4f`|
|Fig3_op|0–100 dBm、10^0刻度与左上低功率平台完整；130 mrad，TNR10，2.67×2.37 inch；FOV floor可见，图例零遮挡。|`6abb3adf5cb0d7aa4b65e442e75b1cef89aae7e60df206b1c337c1dde762a130`|
|Fig4_width_tradeoff|3面板，7.5 m/1 W，ABER和OP对数轴，实际半径上限18 cm，图例右上。|`e015ac03ed4a58961acb457ede2e57523e46c365a5ccba132614402b6eea886e`|
|Fig5_sensitivity|9面板，位置/FOV/角σ范围正确，强弱EGG标记可区分，90 mrad独立基准明确，图例右上。|`d838148178677a371d8cefcc274af91c734bc9b9b754663e7a9cb43194ab4c04`|
|Fig6_objective_phase|仅10 m ABER相位，40×40单元，0–2π色标和图号清晰。|`bd6bf234b4be3af058b903f0e51be1391ab43332b2f48ac39639c675f8dd4a2e`|
|Fig8_radial_profiles|6面板，0–40 cm，环均值与累计能量、六方案和图例清晰。|`43ae513a159199a68ac6ef669ed94c5e3c17f035a458985415e8a71ce13882e2`|
|Fig9_turbulence|4面板，海洋条件分列，大气强弱标记区分，ABER/OP六方案，图例左下。|`f93aefc95b706fd981a83300a3bf9c49007b069a0f2ee4e2782bf2b8cc2edcdd`|
|Fig10_incidence|4面板，7.5/10 m，全反射跃升与临界线清晰，左上紧凑双列图例不挡直接束或红线。|`86174a44f6fb48e74ab8955e2e93ec931741d3dd5c24f6c3403406d2d35e2aef`|
|Fig11_width_jitter_maps|8面板，实际半径/每轴σ，白色选宽路径，深灰图例使白线可见，灰色无支持区明确。|`d6d35e0258531bc90c8c7a43618fb1f32a27639c9e62728d779920111a894b71`|
|Fig12_adaptive_widths|3面板，5项选宽准则和分别对应ABER/OP，低σ平顶优势与高σ鲁棒优势可辨。|`d31b185257dec928bed10626b7db40e20f44fecc22299e519eb073e71e987738`|

最新两图独立验收记录：`plots/stage2/Fig2_Fig3_validation.json`。原阶段全量数值验收对应其完成时的输入版本。
