function fig14_centroid(S)
%FIG14_CENTROID Centroid scatter without and with beam tracking.
C_NOT=[0.00 0.45 0.74]; C_TRK=[0.85 0.20 0.20];
fig=figure('Position',[100 100 1120 620]);
t=tiledlayout(fig,1,2,'TileSpacing','loose','Padding','loose');
% (a) without tracking
ax=nexttile(t); hold(ax,'on');
scatter(ax,S.cent_x_notrack,S.cent_y_notrack,12,'filled','MarkerFaceColor',C_NOT, ...
    'MarkerFaceAlpha',0.28,'MarkerEdgeColor','none');
plot(ax,mean(S.cent_x_notrack),mean(S.cent_y_notrack),'+','Color','k', ...
    'MarkerSize',13,'LineWidth',1.6);
xlim(ax,[-8 8]); ylim(ax,[-18 18]);
xlabel(ax,{'x (cm)';' '},'FontSize',17); ylabel(ax,'y (cm)','FontSize',17);
set(ax,'FontSize',14,'LineWidth',1.1,'TickDir','in','Box','on');
atlas_panel(ax,'(a)',-0.14);
lg=legend(ax,{'无跟踪'},'Location','northeast'); set(lg,'Box','off','FontSize',12);
% (b) with tracking
ax=nexttile(t); hold(ax,'on');
scatter(ax,S.cent_x_track,S.cent_y_track,12,'filled','MarkerFaceColor',C_TRK, ...
    'MarkerFaceAlpha',0.35,'MarkerEdgeColor','none');
plot(ax,mean(S.cent_x_track),mean(S.cent_y_track),'+','Color','k', ...
    'MarkerSize',13,'LineWidth',1.6);
xlim(ax,[-0.8 1.0]); ylim(ax,[-2.2 0.4]);
xlabel(ax,{'x (cm)';' '},'FontSize',17); ylabel(ax,'y (cm)','FontSize',17);
set(ax,'FontSize',14,'LineWidth',1.1,'TickDir','in','Box','on');
atlas_panel(ax,'(b)',-0.14);
lg=legend(ax,{'有跟踪'},'Location','northeast'); set(lg,'Box','off','FontSize',12);
atlas_save(fig,'fig14_centroid_scatter');
end
