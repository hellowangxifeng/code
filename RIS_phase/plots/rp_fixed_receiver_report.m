function rp_fixed_receiver_report(c,d,data,unmodulated,v)
%RP_FIXED_RECEIVER_REPORT Explain the comparison protocol and actual gains.
folder=fullfile(c.output,'fixed_receiver');fid=fopen(fullfile(folder,'RESULTS.md'),'w','n','UTF-8');clean=onCleanup(@()fclose(fid));
fprintf(fid,'# Fig2/3：共同固定接收端、逐功率选相\n\n');
fprintf(fid,'用户确认：所有方案同一固定位置/朝向；每个功率点重新选择最优宽度。数值签名`%s`。\n\n',data.fixedReceiverSignature);
fprintf(fid,'接收点世界坐标为(%.9f, %.9f, %.9f)m，距玻璃上表面10 m。名义来光传播方向/孔径法向为(%.9f, %.9f, %.9f)，接收视轴朝其反向。共同方向由硬件可达载波及固定余量0.012决定，与功率、方案和误码排名无关。\n\n', ...
    data.receiverCenterXY,d.geometry.rxPlaneZ,data.receiverDirection);
fprintf(fid,'沿用本次修正前的图设置：FOV半角130 mrad、每轴位置σ80 mm、角σ25 mrad、角偏差(25,0)mrad、接收半径5 cm、弱EGG×弱GG、0–100 dBm、10 pt与Times New Roman。0–100 dBm为该理想线性模型的扫描，100 dBm等于10 MW；高功率端趋于共同FOV地板不代表实际设备功率能力。\n\n');
fprintf(fid,'画幅严格保留2.67 × 2.37 inch；图例使用10 pt紧凑行距，文字和线段仍是可编辑FIG对象。两图保存后重开验证，图例遮线计数均为0。黑色虚线为精确FOV地板。\n\n');
fprintf(fid,'## 原图排序的原因\n\n');
fprintf(fid,'原图把每个接收机移动到各自光斑质心并跟随主光线。10 m时U-RIS光斑质心在(9.521,-13.464)m，约58.77°，却仍被完全追踪接收；Steer的自然宽斑也提供了较强抖动容限。Flat固定目标10 cm，实际R90约10.1 cm，σ80 mm时约45%%的位置样本越出R90，所以均值较高仍会深衰落。这是各自对准后的宽度比较，掩盖了把光送到同一目标的作用。\n\n');
fprintf(fid,'## 本次实际改动\n\n');
fprintf(fid,'- 所有接收积分使用上述同一位置和朝向，不再跟随实际质心，也不随方案/功率移动。\n- Steer使用对准共同目标的全局常梯度LC相位，不作曲率补偿或聚焦；没有故意让Steer瞄错位置。\n- 对Flat和Robust重新生成实际限幅/量化相位并波动传播，共37个候选（Steer 1份、两族各18份，目标半径1–32 cm）。约束仍为实际短轴/等效半径≥4 cm、轴比≤1.35、R50≥2 cm。\n- R-BER/R-OP在全部可行候选中逐功率最小化相应风险，含Steer incumbent。目标相等时用另一风险及均值打破平局；不会因OP=1饱和而任意选择较差BER。\n- Flat在两张图分别选族内最小ABER/OP；Max P最大化平均接收功率。曲线可对应不同功率下的不同实际相位，调度表完整保存。\n\n');
fprintf(fid,'## 数值结果\n\n|方案|1 W平均接收功率(mW，ABER图相位)|ABER|OP|\n|---|---:|---:|---:|\n');
[~,j]=min(abs(data.powerDBm-30));
for k=1:6
    m=data.curves{k};fprintf(fid,'|%s|%.6g|%.8g|%.8g|\n',data.schemeNames{k},1000*m.meanReceivedPowerW(j),m.aber(j),m.op(j));
end
fprintf(fid,'\n同一指标下功率越低越好；比较达到相同误码/中断门限的功率，避免把共同地板处的重合解释为无优势。\n\n');
penalties=[];fields={'aber','op'};thresholds={[.01 .001 .0001],[.1 .01 .001]};
fprintf(fid,'|指标|门限|Steer所需dBm|相应最优方案所需dBm|相对Steer节省dB|\n|---|---:|---:|---:|---:|\n');
for metric=1:2
    k=3+metric;
    for threshold=thresholds{metric}
        a=powerAt(data.powerDBm,data.curves{2}.(fields{metric}),threshold);
        b=powerAt(data.powerDBm,data.curves{k}.(fields{metric}),threshold);
        fprintf(fid,'|%s|%.4g|%.5g|%.5g|%.5g|\n',upper(fields{metric}),threshold,a,b,a-b);
        penalties=[penalties;metric,threshold,a,b,a-b]; %#ok<AGROW>
    end
end
writetable(array2table(penalties,'VariableNames',{'metricId','threshold','steeringDBm','optimizedDBm','savedDB'}),fullfile(folder,'power_advantage.csv'));
fprintf(fid,'\nMax P最优不等于ABER/OP最优。强位置抖动下，窄光斑能提高平均功率，却可能增加离开高强度核心的概率；因此不能预设Steer的BER/OP一定比Max P差。同样，若Steer自然宽度已接近风险最优，优势可能只有有限dB；最优相位也可能退化为Steer。重合结果保留。\n\n');
fprintf(fid,'## U-RIS星号与弱信号复核\n\n');
fprintf(fid,'U-RIS保留相同玻璃/透明器件、相位全零。固定接收区域原300k发射光子仅27条散射包；直接把稀疏像素中的空格当成物理零会制造假地板。本次改为%.0f光子，%d条区域样本，预先由pilot样本按Scott规则确定带宽%.6g m，不按曲线排名选择带宽。两个独立半运行的全功率最大差：ABER %.6g、OP %.6g；另保存带宽±25%%敏感性。\n\n', ...
    unmodulated.scatterPhotons,unmodulated.scatterPacketsInWindow,unmodulated.scatterBandwidthM, ...
    unmodulated.scatterIndependentABERDifference,unmodulated.scatterIndependentOPDifference);
fprintf(fid,'原光场窗口远离这个固定接收点。这里没有将窗口外插零当成衍射计算，而以矩形孔径两次分部积分的非驻相包络处理弹道尾部；17/33两网格收敛，最小相位振荡>1856周期，并留20%%数值余量。全相干强度包络为%.8g W/m²。U-RIS*画加入该最乐观弹道包络后的ABER/OP下包络；细灰带到仅散射估计。实际标量模型结果应位于该范围，未声称求出了尾部的精确瞬时相干场。\n\n',unmodulated.tailIrradianceEnvelopeWm2);
fprintf(fid,'## 验收与文件\n\n');
fprintf(fid,'状态`%s`；独立抖动样本%d。已选%d份不同相位，连同U-RIS共%d份1600单元程序。\n\n', ...
    v.status,c.study.validationSamples,v.selectedCandidateCount,v.programCount);
fprintf(fid,'- `comparison.mat`、`bank.mat`：共同接收数据和全部实际候选。\n- `curves_and_phase_schedule.csv`：逐功率/指标/方案对应相位编号和半径；`phase_programs/`含实际命令及1 μm LUT。\n- `candidates.csv`、`validation*.csv/json`：可行性、能量、独立统计。\n- `power_advantage.csv`：同一目标指标下相对Steer的功率差。\n- `scatter_full.mat`、`mc_batches/`：独立光子流与可恢复检查点。\n- 上一级`Fig2_aber`、`Fig3_op`为新主图；旧逐方案追踪图与旧绘图代码在`plots/archive_tracking_receiver_20261113/`。\n\n');
fprintf(fid,'本次只修正Fig2/3的比较协议；其他图仍按各自已注明的原始场景解释。传播器仍采用标量弹道波动、非相干HG散射、窄角主光线FOV和均值1的EGG×GG。有限候选库最优不等于任意1600单元相位的全局最优。\n');
end

function value=powerAt(x,y,target)
j=find(y<=target,1);value=NaN;if isempty(j),return;end
if j==1,value=x(1);return;end
value=x(j-1)+(x(j)-x(j-1))*(log(target)-log(y(j-1)))/(log(y(j))-log(y(j-1)));
end
