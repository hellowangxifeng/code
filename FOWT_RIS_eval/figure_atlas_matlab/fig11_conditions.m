function fig11_conditions(S)
%FIG11_CONDITIONS Received power versus wave, incidence, wind, distance and aperture.
C_LIN=[0.00 0.45 0.74]; C_FOC=[0.00 0.62 0.45];
C398=exp(-(0.398-0.151)*10.0); C25=exp(-(0.25-0.151)*10.0);
fig=figure('Position',[60 60 1600 900]);
t=tiledlayout(fig,2,4,'TileSpacing','loose','Padding','loose');
lss={'-','--',':'};
% (a) wave amplitude at 0.8 Hz
ax=nexttile(t); hold(ax,'on');
plotS(ax,S.wavegrid_amp,S.wavegrid_P_linear(2,:),C_LIN,'-');
plotS(ax,S.wavegrid_amp,S.wavegrid_P_focus(2,:),C_FOC,'-');
plotS(ax,S.wavegrid_amp,S.wavegrid_P_linear(2,:)*C25,C_LIN,'-.');
plotS(ax,S.wavegrid_amp,S.wavegrid_P_focus(2,:)*C25,C_FOC,'-.');
plotS(ax,S.wavegrid_amp,S.wavegrid_P_linear(2,:)*C398,C_LIN,'--');
plotS(ax,S.wavegrid_amp,S.wavegrid_P_focus(2,:)*C398,C_FOC,'--');
xlim(ax,[10 50]); ylim(ax,[12 28]);
xlabel(ax,{'波幅 A (mm)';' '},'FontSize',16); ylabel(ax,'接收功率 (mW)','FontSize',16);
set(ax,'FontSize',13); atlas_panel(ax,'(a)' ,-0.17);
lg=legend(ax,{'线性 c=0.151','汇聚 c=0.151','线性 c=0.25','汇聚 c=0.25', ...
    '线性 c=0.398','汇聚 c=0.398'},'Location','northeast');
set(lg,'Box','off','FontSize',7,'NumColumns',3);
% (b) wave frequency
ax=nexttile(t); hold(ax,'on');
plotS(ax,S.freq_x,S.freq_P_linear,C_LIN,'-');
plotS(ax,S.freq_x,S.freq_P_focus,C_FOC,'-');
plotS(ax,S.freq_x,S.freq_P_linear*C25,C_LIN,'-.');
plotS(ax,S.freq_x,S.freq_P_focus*C25,C_FOC,'-.');
plotS(ax,S.freq_x,S.freq_P_linear*C398,C_LIN,'--');
plotS(ax,S.freq_x,S.freq_P_focus*C398,C_FOC,'--');
xlim(ax,[0.3 1.3]);
xlabel(ax,{'波频 f (Hz)';' '},'FontSize',16); ylabel(ax,'接收功率 (mW)','FontSize',16);
set(ax,'FontSize',13); atlas_panel(ax,'(b)' ,-0.17);
% (c) incidence
ax=nexttile(t); hold(ax,'on');
plotS(ax,S.inc_theta,S.inc_P_linear,C_LIN,'-');
plotS(ax,S.inc_theta,S.inc_P_focus,C_FOC,'-');
plotS(ax,S.inc_theta,S.inc_P_linear*C25,C_LIN,'-.');
plotS(ax,S.inc_theta,S.inc_P_focus*C25,C_FOC,'-.');
plotS(ax,S.inc_theta,S.inc_P_linear*C398,C_LIN,'--');
plotS(ax,S.inc_theta,S.inc_P_focus*C398,C_FOC,'--');
xlim(ax,[10 60]);
xlabel(ax,{'水入射角 \theta (度)';' '},'FontSize',16); ylabel(ax,'接收功率 (mW)','FontSize',16);
set(ax,'FontSize',13); atlas_panel(ax,'(c)' ,-0.17);
% (d) wind
ax=nexttile(t); hold(ax,'on');
plotS(ax,S.wind_U,S.wind_P_linear,C_LIN,'-');
plotS(ax,S.wind_U,S.wind_P_focus,C_FOC,'-');
plotS(ax,S.wind_U,S.wind_P_linear*C25,C_LIN,'-.');
plotS(ax,S.wind_U,S.wind_P_focus*C25,C_FOC,'-.');
plotS(ax,S.wind_U,S.wind_P_linear_c398,C_LIN,'--');
plotS(ax,S.wind_U,S.wind_P_focus_c398,C_FOC,'--');
xlim(ax,[3 12]);
xlabel(ax,{'风速 U (m/s)';' '},'FontSize',16); ylabel(ax,'接收功率 (mW)','FontSize',16);
set(ax,'FontSize',13); atlas_panel(ax,'(d)' ,-0.17);
% (e) air path
ax=nexttile(t); hold(ax,'on');
plotS(ax,S.dist_R,S.dist_P_linear,C_LIN,'-');
plotS(ax,S.dist_R,S.dist_P_focus,C_FOC,'-');
plotS(ax,S.dist_R,S.dist_P_linear*C25,C_LIN,'-.');
plotS(ax,S.dist_R,S.dist_P_focus*C25,C_FOC,'-.');
plotS(ax,S.dist_R,S.dist_P_linear*C398,C_LIN,'--');
plotS(ax,S.dist_R,S.dist_P_focus*C398,C_FOC,'--');
xlim(ax,[5 15]);
xlabel(ax,{'空气路径 L (m)';' '},'FontSize',16); ylabel(ax,'接收功率 (mW)','FontSize',16);
set(ax,'FontSize',13); atlas_panel(ax,'(e)' ,-0.17);
% (f) aperture
ax=nexttile(t); hold(ax,'on');
plotS(ax,S.ap_radius,S.ap_P_linear,C_LIN,'-');
plotS(ax,S.ap_radius,S.ap_P_focus,C_FOC,'-');
plotS(ax,S.ap_radius,S.ap_P_linear*C25,C_LIN,'-.');
plotS(ax,S.ap_radius,S.ap_P_focus*C25,C_FOC,'-.');
plotS(ax,S.ap_radius,S.ap_P_linear*C398,C_LIN,'--');
plotS(ax,S.ap_radius,S.ap_P_focus*C398,C_FOC,'--');
xlim(ax,[2.5 10]);
xlabel(ax,{'接收孔径半径 a_r (cm)';' '},'FontSize',16); ylabel(ax,'接收功率 (mW)','FontSize',16);
set(ax,'FontSize',13); atlas_panel(ax,'(f)' ,-0.17);
% (g) wind direction at U = 12 m/s
ax=nexttile(t); hold(ax,'on');
plotS(ax,S.winddir_x,S.winddir_P_linear,C_LIN,'-');
plotS(ax,S.winddir_x,S.winddir_P_focus,C_FOC,'-');
plotS(ax,S.winddir_x,S.winddir_P_linear*C25,C_LIN,'-.');
plotS(ax,S.winddir_x,S.winddir_P_focus*C25,C_FOC,'-.');
plotS(ax,S.winddir_x,S.winddir_P_linear*C398,C_LIN,'--');
plotS(ax,S.winddir_x,S.winddir_P_focus*C398,C_FOC,'--');
xlim(ax,[0 90]);
xlabel(ax,{'海面风向角 \phi (度)';' '},'FontSize',16); ylabel(ax,'接收功率 (mW)','FontSize',16);
set(ax,'FontSize',13); atlas_panel(ax,'(g)' ,-0.17);
% (h) wind direction at U = 6 m/s
ax=nexttile(t); hold(ax,'on');
plotS(ax,S.winddir6_x,S.winddir6_P_linear,C_LIN,'-');
plotS(ax,S.winddir6_x,S.winddir6_P_focus,C_FOC,'-');
plotS(ax,S.winddir6_x,S.winddir6_P_linear*C25,C_LIN,'-.');
plotS(ax,S.winddir6_x,S.winddir6_P_focus*C25,C_FOC,'-.');
plotS(ax,S.winddir6_x,S.winddir6_P_linear*C398,C_LIN,'--');
plotS(ax,S.winddir6_x,S.winddir6_P_focus*C398,C_FOC,'--');
xlim(ax,[0 90]);
xlabel(ax,{'海面风向角 \phi (度)';' '},'FontSize',16); ylabel(ax,'接收功率 (mW)','FontSize',16);
set(ax,'FontSize',13); atlas_panel(ax,'(h)' ,-0.17);
atlas_save(fig,'fig11_power_conditions');
end

function plotS(ax,x,y,col,ls)
[xx,yy]=atlas_pchip(x,y,250);
plot(ax,xx,yy,'Color',col,'LineStyle',ls,'LineWidth',1.8);
end
