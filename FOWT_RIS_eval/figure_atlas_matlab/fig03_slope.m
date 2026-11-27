function fig03_slope(S)
%FIG03_SLOPE Along/cross-wind slope PDFs and slope rms versus wind speed.
fig = figure('Position',[80 80 1360 660]);
t = tiledlayout(fig,1,3,'TileSpacing','loose','Padding','loose');
cols = [0.00 0.45 0.74; 0.85 0.33 0.10; 0.00 0.62 0.45; 0.49 0.18 0.56];
Us = [3 5 9 12];
x = S.slope_x;
% (a) along-wind PDFs
ax=nexttile(t); hold(ax,'on');
for k=1:numel(Us)
    i=find(S.slope_U==Us(k),1); s=S.slope_sig(i,1);
    y=exp(-x.^2/(2*s^2))/(s*sqrt(2*pi));
    plot(ax,x,y,'Color',cols(k,:),'LineWidth',1.9);
end
xlim(ax,[-0.30 0.30]);
xlabel(ax,{'顺风向斜率 s_x';' '},'FontSize',17); ylabel(ax,'概率密度','FontSize',17);
set(ax,'FontSize',14); atlas_panel(ax,'(a)',-0.14);
lg=legend(ax,{'U = 3 m/s','U = 5 m/s','U = 9 m/s','U = 12 m/s'},'Location','northeast');
set(lg,'Box','off','FontSize',10);
% (b) cross-wind PDFs
ax=nexttile(t); hold(ax,'on');
for k=1:numel(Us)
    i=find(S.slope_U==Us(k),1); s=S.slope_sig(i,2);
    y=exp(-x.^2/(2*s^2))/(s*sqrt(2*pi));
    plot(ax,x,y,'Color',cols(k,:),'LineWidth',1.9);
end
xlim(ax,[-0.20 0.20]);
xlabel(ax,{'侧风向斜率 s_y';' '},'FontSize',17); ylabel(ax,'概率密度','FontSize',17);
set(ax,'FontSize',14); atlas_panel(ax,'(b)',-0.14);
% (c) slope rms vs U
ax=nexttile(t); hold(ax,'on');
for j=1:4
    [xx,yy]=atlas_pchip(S.slope_U,S.slope_sig(:,j),200);
    if j<=2, cc=cols(1,:); sty='-'; else, cc=cols(2,:); sty='--'; end
    plot(ax,xx,yy,sty,'Color',cc,'LineWidth',1.9);
end
xlim(ax,[3 12]);
xlabel(ax,{'风速 U (m/s)';' '},'FontSize',17); ylabel(ax,'斜率 rms (rad)','FontSize',17);
set(ax,'FontSize',14); atlas_panel(ax,'(c)',-0.14);
lg=legend(ax,{'s_x 全谱','s_y 全谱','s_x 足迹滤波','s_y 足迹滤波'},'Location','northwest');
set(lg,'Box','off','FontSize',9);
atlas_save(fig,'fig03_slope_pdf');
end
