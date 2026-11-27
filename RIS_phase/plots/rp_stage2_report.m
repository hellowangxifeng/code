function rp_stage2_report(s,extra,v)
%RP_STAGE2_REPORT Export exact selected phases, metadata, raw curves and findings.
c=s.config; candidates=[]; selected=[]; phaseDifference=[];
for j=1:3
    d=s.distances{j}; T=d.summary; T.distanceM=repmat(d.distanceM,height(T),1);
    T.family=string(cellfun(@(b) b.family,d.bank,'UniformOutput',false)).';
    T.centerXM=cellfun(@(b) b.spot.centerXY(1),d.bank).';
    T.centerYM=cellfun(@(b) b.spot.centerXY(2),d.bank).';
    T.windowHalfWidthXM=cellfun(@(b) b.phase.receiverHalfWidthM(1),d.bank).';
    T.windowHalfWidthYM=cellfun(@(b) b.phase.receiverHalfWidthM(2),d.bank).';
    T.scatteredMapFraction=cellfun(@(b) b.wave.scatteredMapPowerW/max(b.wave.mapPowerW,realmin),d.bank).';
    candidates=[candidates;T]; %#ok<AGROW>
    for k=1:3
        b=d.bank{d.objectiveIndices(k)}; p=b.phase; row=T(d.objectiveIndices(k),:);
        row.objective=string(d.objectiveNames{k}); selected=[selected;row]; %#ok<AGROW>
        phase=array2table([p.unitIds,p.referenceXY,p.commandDeg,p.gradientOverK0,p.phaseAtReferenceRad], ...
            'VariableNames',{'unitId','referenceXM','referenceYM','commandXDeg','commandYDeg','gradientXOverK0','gradientYOverK0','pistonRad'});
        writetable(phase,fullfile(c.output,sprintf('phase_z%g_%s.csv',d.distanceM,d.objectiveNames{k})));
        N=round(c.array.apertureM(1)/c.device.electrodePitchM);
        phaseLUT=struct('electrodeCoordinateM',((1:N)-(N+1)/2)*c.device.electrodePitchM, ...
            'pistonRad',p.phaseAtReferenceRad,'gradientOverK0',p.gradientOverK0, ...
            'referenceXY',p.referenceXY,'tileCenters',d.illumination.unitCenters(p.unitIds,:), ...
            'wavelengthM',c.source.wavelength,'periodRad',2*pi, ...
            'formula','mod(piston + k0*(b_x*(x-ref_x)+b_y*(y-ref_y)),2*pi)');
        save(fullfile(c.output,sprintf('phase_lut_z%g_%s.mat',d.distanceM,d.objectiveNames{k})),'phaseLUT');
    end
    pairs=[1 2;1 3;2 3]; weight=d.illumination.unitDirectPowerW; weight=weight/sum(weight);
    for pair=1:3
        a=d.bank{d.objectiveIndices(pairs(pair,1))}.phase; b=d.bank{d.objectiveIndices(pairs(pair,2))}.phase;
        phaseDifference=[phaseDifference;d.distanceM,pairs(pair,:), ...
            sqrt(sum(weight.*sum((a.commandDeg-b.commandDeg).^2,2))), ...
            sqrt(sum(weight.*angle(exp(1i*(a.phaseAtReferenceRad-b.phaseAtReferenceRad))).^2))]; %#ok<AGROW>
    end
end
writetable(candidates,fullfile(c.output,'all_candidates.csv')); writetable(selected,fullfile(c.output,'objective_comparison.csv'));
writetable(array2table(phaseDifference,'VariableNames',{'distanceM','objectiveA','objectiveB','weightedCommandRMSDeg','weightedPistonRMSRad'}),fullfile(c.output,'objective_phase_differences.csv'));
fid=fopen(fullfile(c.output,'RESULTS.md'),'w','n','UTF-8'); clean=onCleanup(@() fclose(fid));
fprintf(fid,'# README 第二阶段结果\n\n');
if isfile(fullfile(c.output,'fixed_receiver','comparison.mat'))
    fprintf(fid,'> Fig2/3已改为共同固定接收端、逐功率重新选宽。最新报告见 [fixed_receiver/RESULTS.md](fixed_receiver/RESULTS.md)。本文原三距离表和9份相位不定义新版Fig2/3的相位调度。\n\n');
end
fprintf(fid,'本报告仅对应 README 第3节。数值签名 `%s`。\n\n',c.study.signature);
fprintf(fid,'名义水下入射角40°、水下斜程10 m、每轴位置抖动80 mm、角误差每轴25 mrad、偏差(25,0) mrad、FOV半角240 mrad、接收半径5 cm、设计ON符号功率1 W。\n\n');
fprintf(fid,'所有通信计算均使用独立、均值归一的海洋 EGG 与大气 Gamma-Gamma 乘积。弱/强GG形状为(17,16)/(4.2,1.4)，闪烁指数分别0.125/1.12245；这些是路径等效情景，未把它们误称为短距离大气Cn²实测值。\n\n');
fprintf(fid,'## 三目标实际选相\n\n|z(m)|目标|实际光斑族|目标半径(cm)|实际半径(cm)|轴比|偏离法向(°)|均值(mW)|ABER|OP|\n|---:|---|---|---:|---:|---:|---:|---:|---:|---:|\n');
for k=1:height(selected)
    q=selected(k,:);
    fprintf(fid,'|%g|%s|%s|%.1f|%.3f|%.3f|%.3f|%.4f|%.6g|%.6g|\n', ...
        q.distanceM,q.objective,q.family,100*q.targetWidthM,q.actualRadiusCM,q.axisRatio,q.tiltDeg,1000*q.meanReceivedPowerW,q.aber,q.op);
end
fprintf(fid,'\nRobust ABER/OP分别在广义高斯肩部/尾部族内按相应风险选宽；平顶族单独重新选宽，用于检验其胜出的区间。参考Flat top固定目标半径10 cm。Max mean power在两族可行相位中按平均功率选宽。同一相位若同时最优，保留真实重合结果；表中实际光斑族明确给出。\n\n');
fprintf(fid,'## 大抖动与重新选宽\n\n');
T=extra.sensitivity.table;
for ocean=1:2
    flat=T(T.oceanId==ocean & T.schemeId==3 & T.parameterId==1,:);
    robust=T(T.oceanId==ocean & T.schemeId==4 & T.parameterId==1,:);
    [~,idx]=max(flat.parameterValue);
    fprintf(fid,'- 海洋条件%d、弱GG、z=5 m、σ=%.0f mm：固定Flat top / Robust ABER的ABER为 %.6g / %.6g，OP为 %.6g / %.6g。\n',ocean,1000*flat.parameterValue(idx),flat.aber(idx),robust.aber(idx),flat.op(idx),robust.op(idx));
end
fprintf(fid,'\nFig12已按2026-11-14修订改为3 cm接收孔径、弱EGG×弱GG，并补充固定设计的Fig13/14。最新条件、独立配对稳定性与FOV控制扫描见 [简短说明](fig12_revision/SIMULATION_CONDITIONS.md)。旧width_jitter/adaptive_widths表属于先前5 cm协议，不用于新图。已删除二维热图及其专用显示数据。\n\n');
fprintf(fid,'7.5 m的Robust目标3 cm候选在加密光学网格上R50=1.98349 cm，低于2 cm下限；已提升为加密光场并排除出可行库，再重新执行全部选宽。该修正有独立缓存与可重跑恢复流程。\n\n');
fprintf(fid,'## 直接传播、未调相与全反射\n\n');
fprintf(fid,'未调相RIS与Steering-only使用相同玻璃和14.4 mm通光掩模；前者真实相位和梯度全零。No RIS direct移除玻璃和器件掩模，使用连续水—空气界面的Snell几何映射及原Gaussian角密度，不人为施加RIS填充因子。该基准是直接弹道束，不收集未定义接收角分布的漫散射背景。\n\n');
fprintf(fid,'水—空气轴向临界角 $\\theta_c=\\arcsin(n_a/n_w)=%.5f^\\circ$。有限角谱在临界角附近仍含亚临界尾部，因此不在临界点把整束功率强制清零。60–65°时直接光束的可接收功率已数值不可用，OP趋近1、ABER趋近1/2；RIS可通过实际梯度恢复可传播空气模。\n\n',asind(c.optics.nAir/c.optics.nWater));
fprintf(fid,'入射角图的载波随角度重新计算并量化，宽度固定为名义三目标所选宽度；没有宣称每个角度重新完成全宽度全局优化。接收孔径分别跟随各方案自身实际光斑中心并朝向其主光线，主比较包含接收端重新对准条件。\n\n');
fprintf(fid,'## 图与数据\n\n');
fprintf(fid,'Fig2/3使用130 mrad半角、2.67 × 2.37 inch画幅，所有文字Times New Roman 10 pt；功率0–100 dBm，纵轴显示10^0、低功率平台及FOV floor。共同固定接收位置/朝向，每个功率点重新选择宽度；两图共用communication_data.mat与communication_curves.csv；重画入口与图例全称见Fig2_Fig3_notes.md。\n\n');
fprintf(fid,'Fig1三距离热图已按后续修订改为裸海面MC直接通信基准，束腰按协方差绘制真实椭圆虚线，长短轴为点虚线；主图保留完整计算图域，直接通信小图包含全部透射MC光子。验证见audit/fig1_revision。Fig2/3只含10 m与130 mrad半角。Fig4只含7.5 m、1 W。Fig5只含5 m及强/弱EGG×弱GG。Fig6仅10 m ABER相位。Fig7已删除。Fig8限0–40 cm。Fig9为海洋/大气两因素交叉，每组条件在1 W重新选相后固定该相位扫描功率。Fig10为7.5/10 m入射角。Fig11已删除；Fig12/13/14为3 cm孔径修订的自适应、固定抖动及固定功率比较。\n\n');
fprintf(fid,'曲线图刻度统一10 pt，图例限定坐标区四角；图形按紧凑物理尺寸输出300 dpi PNG及可编辑FIG。240 mrad的名义FOV地板极低，若低于实际扫描纵轴范围，不用截断后的假地板代替；可见扫描范围内的地板保留原值参考。\n\n');
fprintf(fid,'Fig5的独立基准FOV半角为90 mrad，位于30–150 mrad扫描区间内；图2/3当前局部重画为130 mrad。湍流条件用标记区分，各方案用颜色和线型区分；重合指标保留真实重合。\n\n');
fprintf(fid,'- 全部实际相位：phase_z*.csv及phase_lut_z*.mat（每份1600个单元、1 μm电极坐标）。\n- 三维度候选：all_candidates.csv；三目标：objective_comparison.csv、objective_phase_differences.csv。\n- 曲线与扫描：communication_curves.csv、sensitivity.csv、turbulence_scan.csv、turbulence_optimized_phases.csv、incidence_z*.csv、width_jitter.csv、adaptive_widths.csv。\n- 完整MAT：study.mat、stage2.mat、各扫描MAT和wave_cache/。\n- 验收：validation.json、validation_*.csv及audit/stage2/。\n\n');
fprintf(fid,'验证状态：`%s`。公式、变量与适用条件见 ../../MODEL_AND_FORMULAS.md；入口 run_ris_phase.m。绘图代码在上一级plots目录，和算法根目录分开。\n\n',v.status);
fprintf(fid,'模型保留标量弹道Fresnel/CZT、非相干HG散射、窄角主光线FOV、静态板姿态、等效玻璃和理想相位LC条件。大气GG是乘性闪烁，不另加未标定的随机相位屏、光束漂移或孔径平均修正。有限候选库的最优不代表1600单元任意相位全局最优。\n');
end
