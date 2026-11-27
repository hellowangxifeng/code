function fig02_slope(S)
%FIG02_SLOPE Joint slope PDF contours (lines only, no filled colour).
fig = figure('Position',[100 100 980 660]); ax = axes(fig); hold(ax,'on');
X = S.slope_jX; Y = S.slope_jY;
cols = [0.00 0.45 0.74; 0.85 0.33 0.10; 0.00 0.62 0.45];
Us = [3 9 12];
for k = 1:numel(Us)
    i = find(S.slope_U == Us(k),1);
    sx = S.slope_sig(i,1); sy = S.slope_sig(i,2);
    Z = exp(-0.5*((X/sx).^2 + (Y/sy).^2))/(2*pi*sx*sy);
    lv = linspace(0.12,0.92,5)*max(Z(:));
    contour(ax,X,Y,Z,lv,'LineColor',cols(k,:),'LineWidth',1.7);
    sxf = S.slope_sig(i,3); syf = S.slope_sig(i,4);
    Zf = exp(-0.5*((X/sxf).^2 + (Y/syf).^2))/(2*pi*sxf*syf);
    lvf = linspace(0.12,0.92,5)*max(Zf(:));
    contour(ax,X,Y,Zf,lvf,'LineColor',cols(k,:),'LineStyle','--','LineWidth',1.4);
end
xlim(ax,[-0.30 0.30]); ylim(ax,[-0.18 0.18]); daspect(ax,[1 1 1]);
xlabel(ax,{'顺风向斜率 s_x';' '},'FontSize',18);
ylabel(ax,'侧风向斜率 s_y','FontSize',18);
set(ax,'FontSize',15,'LineWidth',1.1,'TickDir','in','Box','on');
h = gobjects(1,6);
for k = 1:3
    h(2*k-1) = plot(ax,nan,nan,'Color',cols(k,:),'LineWidth',1.8);
    h(2*k)   = plot(ax,nan,nan,'Color',cols(k,:),'LineStyle','--','LineWidth',1.5);
end
lg = legend(ax,h,{'U = 3 m/s 全谱','U = 3 m/s 足迹滤波', ...
    'U = 9 m/s 全谱','U = 9 m/s 足迹滤波', ...
    'U = 12 m/s 全谱','U = 12 m/s 足迹滤波'},'Location','eastoutside');
set(lg,'Box','off','FontSize',10);
atlas_save(fig,'fig02_slope_jpdf');
end
