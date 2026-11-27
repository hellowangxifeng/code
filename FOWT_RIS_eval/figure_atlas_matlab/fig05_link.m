function fig05_link(S)
%FIG05_LINK Outage/ABER versus transmit power and operating conditions.
C_LIN = [0.00 0.45 0.74]; C_FOC = [0.00 0.62 0.45]; C_DIR = [0.45 0.45 0.45];
fig = figure('Position',[60 60 1600 1040]);
t = tiledlayout(fig,2,4,'TileSpacing','loose','Padding','loose');
% (a) OP vs Pt
ax=nexttile(t); hold(ax,'on');
semilogy_smooth(ax,S.pt_tx,S.pt_OP_linear_c151,C_LIN,'-');
semilogy_smooth(ax,S.pt_tx,S.pt_OP_focus_c151,C_FOC,'-');
semilogy_smooth(ax,S.pt_tx,S.pt_OP_linear_c25,C_LIN,'-.');
semilogy_smooth(ax,S.pt_tx,S.pt_OP_focus_c25,C_FOC,'-.');
semilogy_smooth(ax,S.pt_tx,S.pt_OP_linear_c398,C_LIN,'--');
semilogy_smooth(ax,S.pt_tx,S.pt_OP_focus_c398,C_FOC,'--');
semilogy_smooth(ax,S.pt_tx,S.pt_OP_direct_c151,C_DIR,':');
ylim(ax,[1e-10 1.2]); xlim(ax,[0 40]); grid(ax,'off');
xlabel(ax,{'发射功率 P_t (dBm)';' '},'FontSize',16); ylabel(ax,'中断概率 P_{out}','FontSize',16);
set(ax,'FontSize',13); atlas_panel(ax,'(a)',-0.17);
lg=legend(ax,{'线性 c=0.151','汇聚 c=0.151','线性 c=0.25','汇聚 c=0.25', ...
    '线性 c=0.398','汇聚 c=0.398','无RIS'},'Location','southwest');
set(lg,'Box','off','FontSize',7);
% (b) ABER vs Pt
ax=nexttile(t); hold(ax,'on');
semilogy_smooth(ax,S.pt_tx,S.pt_ABER_linear_c151,C_LIN,'-');
semilogy_smooth(ax,S.pt_tx,S.pt_ABER_focus_c151,C_FOC,'-');
semilogy_smooth(ax,S.pt_tx,S.pt_ABER_linear_c25,C_LIN,'-.');
semilogy_smooth(ax,S.pt_tx,S.pt_ABER_focus_c25,C_FOC,'-.');
semilogy_smooth(ax,S.pt_tx,S.pt_ABER_linear_c398,C_LIN,'--');
semilogy_smooth(ax,S.pt_tx,S.pt_ABER_focus_c398,C_FOC,'--');
semilogy_smooth(ax,S.pt_tx,S.pt_ABER_direct_c151,C_DIR,':');
ylim(ax,[1e-18 0.7]); xlim(ax,[0 40]); grid(ax,'off');
xlabel(ax,{'发射功率 P_t (dBm)';' '},'FontSize',16); ylabel(ax,'平均误比特率 ABER','FontSize',16);
set(ax,'FontSize',13); atlas_panel(ax,'(b)',-0.17);
% (c) OP vs distance
ax=nexttile(t); hold(ax,'on');
semilogy_smooth(ax,S.dist_R,S.dist_OP_linear_c151,C_LIN,'-');
semilogy_smooth(ax,S.dist_R,S.dist_OP_focus_c151,C_FOC,'-');
semilogy_smooth(ax,S.dist_R,S.dist_OP_linear_c25,C_LIN,'-.');
semilogy_smooth(ax,S.dist_R,S.dist_OP_focus_c25,C_FOC,'-.');
semilogy_smooth(ax,S.dist_R,S.dist_OP_linear_c398,C_LIN,'--');
semilogy_smooth(ax,S.dist_R,S.dist_OP_focus_c398,C_FOC,'--');
semilogy(ax,S.dist_R,ones(size(S.dist_R)),':','Color',C_DIR,'LineWidth',1.4);
ylim(ax,[1e-10 1.2]); xlim(ax,[5 15]);
xlabel(ax,{'距离 L (m)';' '},'FontSize',16); ylabel(ax,'中断概率 P_{out}','FontSize',16);
set(ax,'FontSize',13); atlas_panel(ax,'(c)',-0.17);
% (d) OP vs incidence
ax=nexttile(t); hold(ax,'on');
semilogy_smooth(ax,S.inc_theta,S.inc_OP_linear_c151,C_LIN,'-');
semilogy_smooth(ax,S.inc_theta,S.inc_OP_focus_c151,C_FOC,'-');
semilogy_smooth(ax,S.inc_theta,S.inc_OP_linear_c25,C_LIN,'-.');
semilogy_smooth(ax,S.inc_theta,S.inc_OP_focus_c25,C_FOC,'-.');
semilogy_smooth(ax,S.inc_theta,S.inc_OP_linear_c398,C_LIN,'--');
semilogy_smooth(ax,S.inc_theta,S.inc_OP_focus_c398,C_FOC,'--');
semilogy(ax,S.inc_theta,ones(size(S.inc_theta)),':','Color',C_DIR,'LineWidth',1.4);
ylim(ax,[1e-10 1.2]); xlim(ax,[10 60]);
xlabel(ax,{'水入射角 \theta (度)';' '},'FontSize',16); ylabel(ax,'中断概率 P_{out}','FontSize',16);
set(ax,'FontSize',13); atlas_panel(ax,'(d)',-0.17);
% (e) OP vs aperture
ax=nexttile(t); hold(ax,'on');
semilogy_smooth(ax,S.ap_radius,S.ap_OP_linear_c151,C_LIN,'-');
semilogy_smooth(ax,S.ap_radius,S.ap_OP_focus_c151,C_FOC,'-');
semilogy_smooth(ax,S.ap_radius,S.ap_OP_linear_c25,C_LIN,'-.');
semilogy_smooth(ax,S.ap_radius,S.ap_OP_focus_c25,C_FOC,'-.');
semilogy_smooth(ax,S.ap_radius,S.ap_OP_linear_c398,C_LIN,'--');
semilogy_smooth(ax,S.ap_radius,S.ap_OP_focus_c398,C_FOC,'--');
semilogy(ax,S.ap_radius,ones(size(S.ap_radius)),':','Color',C_DIR,'LineWidth',1.4);
ylim(ax,[1e-10 1.2]); xlim(ax,[2.5 10]);
xlabel(ax,{'接收孔径半径 a_r (cm)';' '},'FontSize',16); ylabel(ax,'中断概率 P_{out}','FontSize',16);
set(ax,'FontSize',13); atlas_panel(ax,'(e)',-0.17);
% (f) OP vs wind
ax=nexttile(t); hold(ax,'on');
semilogy_smooth(ax,S.wind_U,S.wind_OP_linear_c151,C_LIN,'-');
semilogy_smooth(ax,S.wind_U,S.wind_OP_focus_c151,C_FOC,'-');
semilogy_smooth(ax,S.wind_U,S.wind_OP_linear_c25,C_LIN,'-.');
semilogy_smooth(ax,S.wind_U,S.wind_OP_focus_c25,C_FOC,'-.');
semilogy_smooth(ax,S.wind_U,S.wind_OP_linear_c398,C_LIN,'--');
semilogy_smooth(ax,S.wind_U,S.wind_OP_focus_c398,C_FOC,'--');
semilogy(ax,S.wind_U,ones(size(S.wind_U)),':','Color',C_DIR,'LineWidth',1.4);
ylim(ax,[1e-10 1.2]); xlim(ax,[3 12]);
xlabel(ax,{'风速 U (m/s)';' '},'FontSize',16); ylabel(ax,'中断概率 P_{out}','FontSize',16);
set(ax,'FontSize',13); atlas_panel(ax,'(f)',-0.17);
% (g) OP vs wind direction
ax=nexttile(t); hold(ax,'on');
semilogy_smooth(ax,S.winddir_x,S.winddir_OP_linear_c151,C_LIN,'-');
semilogy_smooth(ax,S.winddir_x,S.winddir_OP_focus_c151,C_FOC,'-');
semilogy_smooth(ax,S.winddir_x,S.winddir_OP_linear_c25,C_LIN,'-.');
semilogy_smooth(ax,S.winddir_x,S.winddir_OP_focus_c25,C_FOC,'-.');
semilogy_smooth(ax,S.winddir_x,S.winddir_OP_linear_c398,C_LIN,'--');
semilogy_smooth(ax,S.winddir_x,S.winddir_OP_focus_c398,C_FOC,'--');
ylim(ax,[1e-5 1.3]); xlim(ax,[0 90]);
xlabel(ax,{'海面风向角 \phi (度)';' '},'FontSize',16); ylabel(ax,'中断概率 P_{out}','FontSize',16);
set(ax,'FontSize',13); atlas_panel(ax,'(g)',-0.17);
% (h) ABER vs wind direction
ax=nexttile(t); hold(ax,'on');
semilogy_smooth(ax,S.winddir_x,S.winddir_ABER_linear_c151,C_LIN,'-');
semilogy_smooth(ax,S.winddir_x,S.winddir_ABER_focus_c151,C_FOC,'-');
semilogy_smooth(ax,S.winddir_x,S.winddir_ABER_linear_c25,C_LIN,'-.');
semilogy_smooth(ax,S.winddir_x,S.winddir_ABER_focus_c25,C_FOC,'-.');
semilogy_smooth(ax,S.winddir_x,S.winddir_ABER_linear_c398,C_LIN,'--');
semilogy_smooth(ax,S.winddir_x,S.winddir_ABER_focus_c398,C_FOC,'--');
ylim(ax,[1e-18 0.5]); xlim(ax,[0 90]);
xlabel(ax,{'海面风向角 \phi (度)';' '},'FontSize',16); ylabel(ax,'平均误比特率 ABER','FontSize',16);
set(ax,'FontSize',13); atlas_panel(ax,'(h)',-0.17);
axs=findobj(fig,'Type','axes');
for k=1:numel(axs), if ~strcmp(axs(k).Tag,'legend'), set(axs(k),'YScale','log'); end; end
atlas_save(fig,'fig05_link_conditions');
end

function semilogy_smooth(ax,x,y,color,ls)
[xu,~] = unique(x(:));
if numel(xu)<2, semilogy(ax,x,y,'Color',color,'LineStyle',ls,'LineWidth',1.7); return; end
[xx,yy] = atlas_pchip(x,log10(max(y(:),1e-12)),300);
semilogy(ax,xx,10.^yy,'Color',color,'LineStyle',ls,'LineWidth',1.7);
end
