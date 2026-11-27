function fig10_dynamic(S)
%FIG10_DYNAMIC Dynamic received power versus wave amplitude for 0.4/0.8/1.2 Hz.
C_LIN=[0.00 0.45 0.74]; C_FOC=[0.00 0.62 0.45];
fig=figure('Position',[80 80 1050 600]);
t=tiledlayout(fig,1,2,'TileSpacing','loose','Padding','loose');
lss={'-','--',':'};
% frequency-dependent 10% ratio from the measured sea states
fs=S.wavegrid_freq; r10=zeros(1,3);
for i=1:3
    sel=abs(S.wave_freq-fs(i))<1e-6;
    r10(i)=mean(S.p_p10(sel,5)./S.p_mean(sel,5));
end
% (a) mean power vs amplitude
ax=nexttile(t); hold(ax,'on');
for i=1:3
    plotS(ax,S.wavegrid_amp,S.wavegrid_P_linear(i,:),C_LIN,lss{i});
    plotS(ax,S.wavegrid_amp,S.wavegrid_P_focus(i,:),C_FOC,lss{i});
end
xlim(ax,[10 50]); ylim(ax,[8 26]);
xlabel(ax,{'波幅 A (mm)';'(a)'},'FontSize',17); ylabel(ax,'平均接收功率 (mW)','FontSize',17);
set(ax,'FontSize',14);
lg=legend(ax,{'0.4 Hz 线性','0.4 Hz 汇聚','0.8 Hz 线性','0.8 Hz 汇聚','1.2 Hz 线性','1.2 Hz 汇聚'}, ...
    'Location','northeast'); set(lg,'Box','off','FontSize',8,'NumColumns',2);
% (b) 10% percentile power vs amplitude
ax=nexttile(t); hold(ax,'on');
for i=1:3
    plotS(ax,S.wavegrid_amp,S.wavegrid_P_linear(i,:)*r10(i),C_LIN,lss{i});
    plotS(ax,S.wavegrid_amp,S.wavegrid_P_focus(i,:)*r10(i),C_FOC,lss{i});
end
xlim(ax,[10 50]); ylim(ax,[8 26]);
xlabel(ax,{'波幅 A (mm)';'(d)'},'FontSize',17); ylabel(ax,'10% 分位功率 (mW)','FontSize',17);
set(ax,'FontSize',14);
atlas_save(fig,'fig10_dynamic_power');
end

function plotS(ax,x,y,col,ls)
[xx,yy]=atlas_pchip(x,y,250);
plot(ax,xx,yy,'Color',col,'LineStyle',ls,'LineWidth',1.9);
end
