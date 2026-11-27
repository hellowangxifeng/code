function fig01_flow()
%FIG01_FLOW Method flow: measured data -> optics -> phase design -> metrics.
fig = figure('Position',[60 60 1500 720]);
ax = axes(fig,'Position',[0 0 1 1]); axis(ax,'off'); hold(ax,'on');
xlim(ax,[0 100]); ylim(ax,[0 100]);
rectangle(ax,'Position',[0 0 100 100],'EdgeColor','w','LineWidth',0.1,'FaceColor','none');
plot(ax,[0 100],[0 0],'Color','w','LineWidth',0.1);
plot(ax,[0 100],[100 100],'Color','w','LineWidth',0.1);

c_in  = [0.85 0.92 0.97];
c_mod = [0.87 0.95 0.89];
c_met = [0.99 0.94 0.84];
c_cmp = [0.94 0.89 0.96];
c_par = [0.94 0.94 0.94];

drawBox(ax, 2, 58, 17, 22, '① 实测数据', ...
    {'6-DOF  183 Hz','波高  128 Hz'}, c_in);
drawBox(ax, 22, 58, 19, 22, '② 分层光学模型', ...
    {'姿态 / 平移','Fermat 路径 + 玻璃折射','离散指令量化'}, c_mod);
drawBox(ax, 44, 58, 18, 22, '③ 相位设计', ...
    {'功率加权公共方向','点目标聚焦','±15.414° 限幅'}, c_mod);
drawBox(ax, 65, 58, 17, 22, '④ 光斑与孔径', ...
    {'全波锚点光斑','质心平移','5 cm 孔径功率'}, c_mod);
drawBox(ax, 85, 58, 13, 22, '⑤ 性能评估', ...
    {'OP / ABER','质心 rms','可用率'}, c_met);

arrow(ax, 19.2, 69, 21.6, 69, false);
arrow(ax, 41.2, 69, 43.6, 69, false);
arrow(ax, 62.2, 69, 64.6, 69, false);
arrow(ax, 82.2, 69, 84.6, 69, false);

drawBox(ax, 3, 18, 22, 17, '公共参数', ...
    {'532 nm · 40×40 RIS','c = 0.151 / 0.398 m^{-1}','10 m 水 + 5–15 m 空气'}, c_par);
drawBox(ax, 34, 18, 30, 17, '对比方案', ...
    {'无RIS直接','RIS线性调相','RIS汇聚调相'}, c_cmp);
drawBox(ax, 73, 18, 24, 17, '控制评估', ...
    {'处理时延 0–100 ms','更新周期 20–50 ms'}, c_met);

arrow(ax, 14, 35.2, 29, 57.6, true);
arrow(ax, 49, 35.2, 53, 57.6, true);
arrow(ax, 85, 35.2, 91.5, 57.6, true);

set(fig,'PaperPositionMode','auto');
if ~isfolder('figures'), mkdir('figures'); end
print(fig,fullfile('figures','fig01_flow.png'),'-dpng','-r300');
savefig(fig,fullfile('figures','fig01_flow.fig'));
close(fig);
end
function drawBox(ax,x,y,w,h,t,body,col)
rectangle(ax,'Position',[x y w h],'Curvature',0.08,'FaceColor',col, ...
    'EdgeColor',[0.30 0.30 0.30],'LineWidth',1.2);
text(ax,x+w/2,y+h-3.4,t,'HorizontalAlignment','center','VerticalAlignment','top', ...
    'FontSize',15,'FontWeight','bold');
text(ax,x+w/2,y+h/2-2.2,strjoin(body,newline),'HorizontalAlignment','center', ...
    'VerticalAlignment','middle','FontSize',11.5);
end

function arrow(ax,x1,y1,x2,y2,dashed)
a = annotation('arrow',[x1 x2]/100,[y1 y2]/100);
a.LineWidth = 1.5;
a.Color = [0.25 0.25 0.25];
a.HeadWidth = 8;
a.HeadLength = 10;
if dashed
    a.LineStyle = '--';
end
end
