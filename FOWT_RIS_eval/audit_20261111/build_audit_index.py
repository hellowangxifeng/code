"""Build the review index and archive validation output without editing sources."""
from pathlib import Path
import json, hashlib, subprocess, sys

HERE=Path(__file__).resolve().parent
ROOT=HERE.parent
ev=json.loads((HERE/'python_evidence.json').read_text(encoding='utf-8'))
me=json.loads((HERE/'matlab_evidence.json').read_text(encoding='utf-8'))
notes={
 'src/fx_phase_design.m':'P1-05/06/07，P2-06/08/10/11：相位、曲率、可达性和量化',
 'src/fx_motion_from_experiment.m':'P2-16：校准、中心说明与安装坐标；缺失数据处理需约束',
 'src/fx_wave_from_gauge.m':'P2-12/13：FFT正频率与输入有效性；深水/单向波是假设',
 'tools/fx_optics.py':'P2-06/07/08/09：合法模式、共同方向网格；分层几何主体合理',
 'tools/replay.py':'P1-01/02/03/08及P2-01/02/04/05：接收机、代理场、积分、单位',
 'tools/control_study.py':'P1-04，P2-03：物理灵敏度缺项与启动因果性',
 'tools/fowt_data.py':'P2-14：原始路径；全缺测填0、NaN边界与MATLAB行为不一致',
 'tools/build_dataset.py':'P2-14/15：路径和同步；频谱/衰减估计边界见报告',
 'tools/model_replay.py':'继承replay问题；0.01秒回放中断时长被按0.02秒统计',
 'tools/analyze_results.py':'P2-17：worst方向错误，功率积分继承replay问题',
 'tools/export_nominal.m':'名义姿态与当前地图混用；improvedQ硬编码；曲率锚点标为改进方案',
 'tools/validate.py':'运行通过；以硬编码结果区间作回归，未覆盖独立物理约束',
 'tools/probe_case.m':'过时字段rxAperturePowerW导致运行失败',
 'tools/probe_wave.m':'原缓存经rp_refine_receiver后含pixelSizeM，字段可用；光斑范围/接收点诊断',
 'tools/probe_err.m':'矩阵幂^2及子集加权均值错误',
 'tools/probe_model.m':'mass.total类型错误，兼容输出DOF索引错误',
 'tools/probe_motion_bridge.m':'失效原始数据路径；桥接校准问题',
 'tools/probe_phase_res.m':'粗采样原始光学相位unwrap混叠',
 'tools/probe_curv.m':'调用曲率限幅存在P1-05',
 'tools/probe_curv2.m':'未约束曲率，只能作为理想诊断',
 'tools/probe_improved.m':'P1-01、P2-07；生成方向锚点不等于固定RX聚焦',
 'tools/probe_feasible.m':'连续盒/圆盘可达性诊断，不等于离散命中；未发现另一个独立公式错误',
 'tools/probe_pw.m':'P2-07；含未完成空循环与旧网格探索',
 'tools/probe_pareto.m':'P2-07；空集合时max输出与标签需谨慎',
 'tools/probe_land.m':'hOut硬编码为半厚度，仅默认screenFraction成立',
 'tools/probe_measured_wave.m':'深水单向假设；局部重建正确不代表空间场验证',
 'tools/probe_wave_bridge.m':'继承FFT桥接；结束点<=与<不一致',
 'tools/probe_nat.m':'前三模态兼容顺序可用，频率公式只作对角/Ainf诊断',
 'tools/probe_damp.m':'同上；常系数等效阻尼不等于完整频率相关系统模态阻尼',
 'tools/probe_sea_rao.m':'六DOF标准顺序与RAO用途可辨；不是实测验证',
 'tools/probe_sea_rao2.m':'只显示数组类型/形状；未见独立数值错误',
 'tools/probe_sea_results.m':'生成多风向风速包络；下游重复风速插值需修复',
 'tools/probe_cmd.m':'饱和诊断；饱和本身不是算法求解错误',
 'tools/probe_sign.m':'正号约定对照；未发现与exp(+iPhi)不一致',
 'tools/probe_centroid.m':'参考光线/实际场质心不同；不得互相等同',
 'tools/probe_aperture.m':'直接圆孔求和可作独立积分对照，不能证明动态光斑近似',
 'tools/probe_coh.m':'相干/非相干图诊断；继承锚点相位模型假设',
 'tools/probe_power.m':'功率收支诊断；窗口功率不等于整个无限平面功率',
 'tools/probe_radial.m':'阈值覆盖率不是包围能量90%；受原地图范围限制',
 'tools/probe_speed.m':'性能/复算探针；不能单独验证波动算法',
 'tools/make_figures_all.py':'图1–4主要读取既有数据；继承回放错误；含硬编码参数/标注',
 'tools/make_figure5.py':'模型/实测比较不能升级为RIS硬件验证；依赖输出目录已存在',
 'figure_atlas_matlab/export_atlas_data.py':'P1-09及P2-18至26：合成聚焦、缩放、统计/扫参模型',
 'figure_atlas_matlab/tmp_enrich05.py':'非幂等文本替换，重复运行重复曲线',
 'figure_atlas_matlab/fig01_flow.m':'流程示意；没有增加验证证据',
 'figure_atlas_matlab/fig02_slope.m':'独立Gaussian斜率PDF假设，不能当作实测斜率分布',
 'figure_atlas_matlab/fig03_slope.m':'同上；斜率tan角与rad标签仅小角度近似等同',
 'figure_atlas_matlab/fig05_link.m':'继承合成功率/OP/ABER与功率参数不一致问题',
 'figure_atlas_matlab/fig06_spots.m':'合成聚焦；孔径图形改变了数据半径；截面额外平滑',
 'figure_atlas_matlab/fig07_motion.m':'合成聚焦时程固定增益；孔径位移方向和图形比例问题',
 'figure_atlas_matlab/fig08_distance.m':'几何缩放合成焦图，非重传播；高斯捕获近似',
 'figure_atlas_matlab/fig09_turb.m':'P2-19：对数正态均值和p10公式错误',
 'figure_atlas_matlab/fig10_dynamic.m':'P2-23：拟合外推与p10比例替代统计',
 'figure_atlas_matlab/fig11_conditions.m':'合成扫参；部分低功率曲线被坐标裁去',
 'figure_atlas_matlab/fig14_centroid.m':'来自代理模型参考光线质心，不是实测光斑位置',
 'figure_atlas_matlab/fig16_incidence.m':'P1-10：插值点标记为MC',
 'figure_atlas_matlab/fig17_aperture.m':'完全合成15mW模型，忽略输入S；系数不是能量比例',
 'figure_atlas_matlab/fig18_op.m':'继承export数据，log插值不能新增证据',
 'figure_atlas_matlab/fig19_aber.m':'继承export数据与OOK功率定义问题',
 'figure_atlas_matlab/fig20_distance_turb.m':'继承距离/湍流模型问题',
 'figure_atlas_matlab/fig21_delay.m':'继承错误J；Pnom硬编码14.8',
 'figure_atlas_matlab/fig22_update.m':'继承错误J，20–50ms推荐缺少有效依据',
 'figure_atlas_matlab/atlas_panel.m':'排版辅助；未见独立算法问题',
 'figure_atlas_matlab/atlas_pchip.m':'保形插值；重复x仅留首个值需调用端约定',
 'figure_atlas_matlab/atlas_save.m':'图形导出；依赖当前目录，未见独立数学问题',
 'figure_atlas_matlab/atlas_style.m':'全局图形样式；未见独立算法问题',
 'figure_atlas_matlab/make_all.m':'只重画已有atlas_data，不会从原始数据重建光学结果',
 '浮式平台六自由度实验数据/FOWT Platform Testing/MainCode.m':'固定输入/人工同步演示；手动lap未被新批处理继承',
}
lines=['# 源文件审查清单','',f"原始文件：{ev['counts']['m']} 个 .m、{ev['counts']['py']} 个 .py；{ev['counts']['total_lines']} 行。",'',
'所有原始MATLAB文件已运行checkcode；Python文件通过AST解析。未发现独立问题不代表所有物理条件已验证，调用方仍受主报告中的模型问题影响。','',
'| 文件 | 行数 | 审查要点 |','|---|---:|---|']
missing=[]; changed=[]
for item in ev['inventory']:
    rel=item['file'].replace('\\','/')
    note=notes.get(rel)
    if note is None: missing.append(rel); note='需核对覆盖记录'
    f=ROOT/Path(item['file'])
    if hashlib.sha256(f.read_bytes()).hexdigest()!=item['sha256']: changed.append(rel)
    lines.append(f"| {rel} | {item['lines']} | {note} |")
if missing or changed: raise RuntimeError({'missing_notes':missing,'changed_sources':changed})
(HERE/'源文件清单.md').write_text('\n'.join(lines)+'\n',encoding='utf-8')
run=subprocess.run([sys.executable,str(ROOT/'tools/validate.py')],capture_output=True,text=True,encoding='utf-8',errors='replace')
(HERE/'原验证输出.txt').write_text(run.stdout+'\nSTDERR:\n'+run.stderr+f'\nEXIT CODE: {run.returncode}\n',encoding='utf-8')
issue_count=sum(len(x['issues']) if isinstance(x['issues'],list) else int(bool(x['issues'])) for x in me['staticAnalysis'])
print(json.dumps({'source_files':len(ev['inventory']),'unmodified_hashes':not changed,'matlab_checked':len(me['staticAnalysis']),'matlab_issue_count':issue_count,'python_parsed':len(ev['python_compile']),'original_validate_exit':run.returncode,'atlas_focus_gain':ev['atlas']['focus_gain'],'atlas_focus_aperture_mw':ev['atlas']['focus_P5_mw']},ensure_ascii=False))
