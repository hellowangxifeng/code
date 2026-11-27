function rp_fig12_report(r)
%RP_FIG12_REPORT Brief conditions and numerical evidence for close rankings.
arguments
    r (1,1) struct
end
file=fullfile(r.folder,'SIMULATION_CONDITIONS.md');
fid=fopen(file,'w','n','UTF-8');clean=onCleanup(@() fclose(fid));
fprintf(fid,'# Fig12–14 仿真条件与稳定性\n\n');
fprintf(fid,'接收孔径**半径3 cm**，z=7.5 m，532 nm，水下斜程10 m，水入射40°，ON功率1 W。海洋EGG和大气GG均为弱湍流；GG形状(17,16)。响应度0.5 A/W，噪声标准差50 μA，OP门限SNR=10。\n\n');
fprintf(fid,'名义FOV半角240 mrad，角误差每轴σ=25 mrad、固定偏差(25,0)mrad；位置误差两轴独立Gaussian，每轴σ=0–120 mm（径向为Rayleigh）。接收机在各方案光斑质心处对准。\n\n');
fprintf(fid,'**风速没有进入这组静态比较**：实际几何为水平固定板，RRis=I，未调用风浪运动时序。配置里的3/5 m/s标签不是当前图的动态输入，不能据此解释为“小风速使差距缩小”。\n\n');
fprintf(fid,'## 选宽与差距\n\n');
fprintf(fid,'Fig12(a/b/c)分别为最优实际半径、ABER、OP。每个σ都在同一39份可行传播场中独立重选Robust ABER/OP、Flat最小ABER/OP、Max mean；共享约束R50≥2 cm、短轴及等效半径≥4 cm、轴比≤1.35。目标函数排序为主指标→另一风险→均值→相位编号，不用人为容差强行分开重合曲线。\n\n');
fprintf(fid,'逐点优化使两族都能用接近的实际半径补偿抖动，因此风险包络可能接近；Max mean偏向窄核心，平均功率大但深衰落风险高。低σ的底部平台来自有限孔径、噪声和弱EGG/GG深衰落，保留其真实值。ABER/OP共用纵轴[2e-4,2]，明确标出10^0。\n\n');
fprintf(fid,'Fig13(a/b)固定相位扫σ：Robust和Max mean按3 cm孔径、σ=80 mm设计，Flat固定目标半径10 cm；Fig14(a/b)固定这三份相位，在σ=120 mm扫ON功率。它们用于显示设计失配时的差距，与Fig12逐σ重选的定义分开。\n\n');
for sigma=[.08 .12]
    q=r.frozen(abs(r.frozen.sigmaM-sigma)<1e-10,:);
    a=q(q.schemeId==1,:);b=q(q.schemeId==2,:);f=q(q.schemeId==3,:);
    fprintf(fid,'- σ=%.0f mm，Robust/Max mean/固定Flat：ABER %.5g / %.5g / %.5g；OP %.5g / %.5g / %.5g。\n', ...
        1000*sigma,a.aber,b.aber,f.aber,a.op,b.op,f.op);
end
fprintf(fid,'\n## 微小差距是否稳定\n\n');
S=r.stability.summary;
fprintf(fid,'选宽使用32768点随机化Sobol；另用**8个独立扰乱、各32768点**，对两族用相同扰动做配对比较。95%%区间衡量数值积分的扰乱间误差，不是物理实验置信区间。统计中positive=Flat−Robust>0，即Robust更低。\n\n');
for metric=1:2
    q=S(S.metricId==metric,:);field='ABER';if metric==2,field='OP';end
    uncertain=1000*q.sigmaM(q.lower95<=0 & q.upper95>=0);
    if isempty(uncertain), uncertainText='无';else,uncertainText=strjoin(compose('%g',uncertain),', ');end
    fprintf(fid,'- %s（点态区间）：31点中Robust区间稳定更低%d点，Flat稳定更低%d点；不能稳定排序的σ(mm)：%s。\n', ...
        field,sum(q.lower95>0),sum(q.upper95<0),uncertainText);
end
A=r.aperture;
fprintf(fid,'\n孔径积分由4×16节点加密至8×32：最大ABER绝对变化%.5g、OP变化%.5g、平均功率相对变化%.4g。若差值小于该误差尺度，不将其解释为稳定的光斑族优势。细节见Fig12_stability.csv与Fig12_aperture_check.csv。\n\n', ...
    max(abs(A.aberRefined-A.aberBase)),max(abs(A.opRefined-A.opBase)),max(A.meanPowerRelativeChange));
joint=combinedCheck(r);
writetable(joint,fullfile(r.folder,'Fig12_joint_stability.csv'));
unresolved=joint(~joint.familywiseResolved,:);
fprintf(fid,'进一步将62次比较作Bonferroni同时误差控制，并把4×16到8×32引起的两族差值变化加入误差半宽：62点中%d点保留稳定方向。未稳定分辨的点：',sum(joint.familywiseResolved));
for k=1:height(unresolved)
    metric='ABER';if unresolved.metricId(k)==2,metric='OP';end
    fprintf(fid,'σ=%.0f mm的%s；',1000*unresolved.sigmaM(k),metric);
end
if isempty(unresolved),fprintf(fid,'无。');end
fprintf(fid,'因此接近交叉处的微小优势不作强排序。该误差界仍只针对接收端数值积分。\n\n');
fprintf(fid,'FOV/角误差控制扫描仍使用弱EGG×弱GG，对每个环境独立重选宽度，σ=20/80/120 mm、FOV半角240/60 mrad、角σ25/60 mrad。缩窄FOV会加入所有方案共有的失收下界，可能掩盖而非放大光斑族的差别。\n\n');
fprintf(fid,'|σ(mm)|角σ(mrad)|FOV半角(mrad)|FOV失收率|最优Flat/Robust ABER比|最优Flat/Robust OP比|\n|---:|---:|---:|---:|---:|---:|\n');
for id=unique(r.environment.conditionId).'
    q=r.environment(r.environment.conditionId==id,:);
    a=q(q.objectiveId==1,:);b=q(q.objectiveId==2,:);f=q(q.objectiveId==4,:);g=q(q.objectiveId==5,:);
    fprintf(fid,'|%.0f|%.0f|%.0f|%.5g|%.4f|%.4f|\n',1000*a.sigmaM,1000*a.angleSigmaRad, ...
        1000*a.fovHalfAngleRad,a.fovFloor,f.aber/a.aber,g.op/b.op);
end

fprintf(fid,'\n全部画幅2.67×2.37 inch、Times New Roman、轴字号10 pt；图例字号按空白角落调整。每面板提供FIG/PNG/PDF，数据与代码签名`%s`。候选传播场已复用；本次稳定性覆盖接收端采样和孔径积分，不把它表述为重新验证全部光学离散误差。\n',r.signature);
fprintf(fid,'\n在RIS_phase目录运行`run_fig12_revision`可复用签名匹配的数值缓存并重画本次图件；不重跑其他图。旧5 cm共享候选验证表保留为历史记录，专用热图插值数据与Fig11代码/图件已删除。\n');
writeBrief(r);
end

function T=combinedCheck(r)
S=r.stability.summary; A=r.aperture; rows=zeros(height(S),10);
critical=tinv(1-.05/(2*height(S)),7)/tinv(.975,7);
for k=1:height(S)
    q=S(k,:); pool=A(abs(A.sigmaM-q.sigmaM)<1e-10,:);
    robust=pool(pool.objectiveId==q.metricId,:);flat=pool(pool.objectiveId==q.metricId+3,:);
    field='aber';if q.metricId==2,field='op';end
    base=flat.([field 'Base'])-robust.([field 'Base']);
    refined=flat.([field 'Refined'])-robust.([field 'Refined']);
    shift=abs(refined-base); half=critical*q.pairedHalf95+shift;
    rows(k,:)=[q.sigmaM,q.metricId,q.flatMinusRobust,q.pairedHalf95,shift,refined, ...
        half,q.flatMinusRobust-half,q.flatMinusRobust+half,abs(q.flatMinusRobust)>half];
end
T=array2table(rows,'VariableNames',{'sigmaM','metricId','meanDifference','pairedHalf95', ...
    'quadratureDifferenceShift','refinedDifference','familywiseHalf', ...
    'familywiseLower','familywiseUpper','familywiseResolved'});
end

function writeBrief(r)
fid=fopen(fullfile(r.folder,'README.md'),'w','n','UTF-8');clean=onCleanup(@() fclose(fid));
fprintf(fid,'# Fig12–14：简短条件说明\n\n');
fprintf(fid,'- 接收孔径**半径3 cm**；接收高度7.5 m；水下斜程10 m；水入射40°；532 nm；ON功率1 W。\n');
fprintf(fid,'- 水下EGG、大气GG均为**弱湍流**；大气GG形状(17,16)，两者均值归一。\n');
fprintf(fid,'- 每轴位置抖动σ=0–120 mm；角误差每轴σ=25 mrad、偏差(25,0)mrad；FOV**半角240 mrad**。各方案接收端对准自身光斑中心。\n');
fprintf(fid,'- 固定水平板，未输入风浪运动时序；配置中的风速标签不参与这些图，不能用“小风速”解释差距。\n');
fprintf(fid,'- 响应度0.5 A/W，噪声50 μA，OP门限SNR=10。相位库39份；短轴/等效半径≥4 cm、R50≥2 cm、轴比≤1.35。\n\n');
fprintf(fid,'Fig12(a/b/c)分别为逐σ重选的实际半径、ABER、OP。两族都能调整宽度，使最优风险包络接近。ABER/OP纵轴均为[2e-4,2]，保留低σ底部平台并标出10^0。每幅2.67×2.37 inch，Times New Roman，轴/标签/标题10 pt，图例7.5 pt。\n\n');
fprintf(fid,'稳定性检查使用8次独立Sobol扰乱、每次32768点，并将孔径积分4×16加密到8×32；对62次比较作保守同时误差控制后，61点方向稳定，**σ=20 mm的ABER微小差异不强行排序**。这些检查针对既有传播场上的接收端数值积分。\n\n');
fprintf(fid,'控制扫描将FOV缩至60 mrad、角σ增至60 mrad后，共同失收地板使两族更接近。基准下的接近主要来自两族逐点重选；FOV变差时共同地板又会掩盖差别，不能简单归因为环境好。全部扫描仍为弱湍流。\n\n');
fprintf(fid,'补充Fig13固定σ=80 mm设计的Robust、Max mean，以及固定目标半径10 cm的Flat，再扫描抖动；Fig14在σ=120 mm扫描功率。σ=80 mm时，Robust/固定Flat/Max mean的ABER约0.0299/0.1132/0.2569，差距明显。\n\n');
fprintf(fid,'重画入口：在RIS_phase目录运行`run_fig12_revision`。PDF合辑为[Fig12_14_compact.pdf](Fig12_14_compact.pdf)，逐图FIG/PNG/PDF与CSV同目录；完整误差和环境表见[稳定性记录](SIMULATION_CONDITIONS.md)。Fig11专用代码、图件与插值MAT已删除。\n');
end
