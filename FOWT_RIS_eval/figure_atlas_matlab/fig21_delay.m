function fig21_delay(S)
%FIG21_DELAY Centroid, received power and normalized error versus delay.
C45=[0.00 0.45 0.74]; C30=[0.85 0.33 0.10];
fig=figure('Position',[60 60 1280 430]);
t=tiledlayout(fig,1,3,'TileSpacing','loose','Padding','loose');
% (a) centroid rms
ax=nexttile(t); hold(ax,'on');
plot(ax,S.delay_ms,S.delay_rms_ZOH_45,'o-','Color',C45,'LineWidth',1.7,'MarkerSize',5,'MarkerFaceColor','w');
plot(ax,S.delay_ms,S.delay_rms_lin_45,'s-','Color',C45,'LineWidth',1.7,'MarkerSize',5,'MarkerFaceColor','w');
plot(ax,S.delay_ms,S.delay_rms_ZOH_30,'o--','Color',C30,'LineWidth',1.6,'MarkerSize',5,'MarkerFaceColor','w');
plot(ax,S.delay_ms,S.delay_rms_lin_30,'s--','Color',C30,'LineWidth',1.6,'MarkerSize',5,'MarkerFaceColor','w');
xlim(ax,[0 100]);
xlabel(ax,{'处理时延 (ms)';' '},'FontSize',16); ylabel(ax,'质心 rms (mm)','FontSize',16);
set(ax,'FontSize',13); atlas_panel(ax,'(a)',-0.24);
lg=legend(ax,{'45 mm ZOH','45 mm 预测','30 mm ZOH','30 mm 预测'},'Location','northwest'); set(lg,'Box','off','FontSize',9);
% (b) received power for two water qualities and two sea states
ax=nexttile(t); hold(ax,'on');
Pnom=14.8;
P451=Pnom*10.^(S.delay_loss_lin/10);
P452=Pnom*10.^(S.delay_loss_c398/10);
P301=Pnom*10.^(S.delay_loss_lin_30/10);
P302=Pnom*10.^(S.delay_loss_c398_30/10);
plot(ax,S.delay_ms,P451,'-o','Color',C45,'LineWidth',1.8,'MarkerFaceColor','w','MarkerSize',5);
plot(ax,S.delay_ms,P452,'--s','Color',C45,'LineWidth',1.6,'MarkerFaceColor','w','MarkerSize',5);
plot(ax,S.delay_ms,P301,'-o','Color',C30,'LineWidth',1.8,'MarkerFaceColor','w','MarkerSize',5);
plot(ax,S.delay_ms,P302,'--s','Color',C30,'LineWidth',1.6,'MarkerFaceColor','w','MarkerSize',5);
xlim(ax,[0 100]);
xlabel(ax,{'处理时延 (ms)';' '},'FontSize',16); ylabel(ax,'接收功率 (mW)','FontSize',16);
set(ax,'FontSize',13); atlas_panel(ax,'(b)',-0.24);
lg=legend(ax,{'45 mm, c=0.151','45 mm, c=0.398','30 mm, c=0.151','30 mm, c=0.398'},'Location','southeast'); set(lg,'Box','off','FontSize',8);
% (c) normalized error
ax=nexttile(t); hold(ax,'on');
plot(ax,S.delay_ms,S.delay_rms_lin_45/max(S.delay_rms_lin_45),'-o','Color',C45,'LineWidth',1.8,'MarkerFaceColor','w','MarkerSize',5);
plot(ax,S.delay_ms,S.delay_rms_lin_30/max(S.delay_rms_lin_30),'--s','Color',C30,'LineWidth',1.7,'MarkerFaceColor','w','MarkerSize',5);
xlim(ax,[0 100]);
xlabel(ax,{'处理时延 (ms)';' '},'FontSize',16); ylabel(ax,'归一化误差','FontSize',16);
set(ax,'FontSize',13); atlas_panel(ax,'(c)',-0.24);
atlas_save(fig,'fig21_delay');
end
