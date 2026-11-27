function fig22_update(S)
%FIG22_UPDATE Command update period: rms and loss, with recommended range.
C45=[0.00 0.45 0.74]; C30=[0.85 0.33 0.10];
fig=figure('Position',[100 100 940 430]);
t=tiledlayout(fig,1,2,'TileSpacing','loose','Padding','loose');
ax=nexttile(t); hold(ax,'on');
patch(ax,[20 50 50 20],[0 0 140 140],[0.40 0.40 0.40],'FaceAlpha',0.08,'EdgeColor','none');
plot(ax,S.update_ms,S.update_rms_45,'-o','Color',C45,'LineWidth',1.9,'MarkerFaceColor','w','MarkerSize',6);
plot(ax,S.update_ms,S.update_rms_30,'--s','Color',C30,'LineWidth',1.8,'MarkerFaceColor','w','MarkerSize',6);
xlim(ax,[0 200]); ylim(ax,[0 140]);
xlabel(ax,{'更新周期 (ms)';' '},'FontSize',17); ylabel(ax,'质心 rms (mm)','FontSize',17);
set(ax,'FontSize',14); atlas_panel(ax,'(a)',-0.24);
lg=legend(ax,{'45 mm / 0.8 Hz','30 mm / 0.4 Hz'},'Location','northwest'); set(lg,'Box','off','FontSize',10);
text(ax,58,112,'推荐 20–50 ms','FontSize',11,'Color',[0.25 0.25 0.25]);
ax=nexttile(t); hold(ax,'on');
patch(ax,[20 50 50 20],[-2.6 -2.6 0.3 0.3],[0.40 0.40 0.40],'FaceAlpha',0.08,'EdgeColor','none');
plot(ax,S.update_ms,S.update_loss_45,'-o','Color',C45,'LineWidth',1.9,'MarkerFaceColor','w','MarkerSize',6);
plot(ax,S.update_ms,S.update_loss_30,'--s','Color',C30,'LineWidth',1.8,'MarkerFaceColor','w','MarkerSize',6);
xlim(ax,[0 200]); ylim(ax,[-2.6 0.3]);
xlabel(ax,{'更新周期 (ms)';' '},'FontSize',17); ylabel(ax,'平均功率损耗 (dB)','FontSize',17);
set(ax,'FontSize',14); atlas_panel(ax,'(b)',-0.24);
atlas_save(fig,'fig22_update');
end
