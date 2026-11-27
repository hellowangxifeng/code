function fig07_motion(S)
%FIG07_MOTION Motion-induced spot offset and power time series.
C_LIN = [0.00 0.45 0.74]; C_FOC = [0.00 0.62 0.45];
fig = figure('Position',[60 60 1000 1160]);
t = tiledlayout(fig,3,2,'TileSpacing','loose','Padding','loose');
off = S.offv_collimated_static; % cm, N x 2
[~,idx] = sort(sum(off.^2,2)); idx = idx(max(1,round(0.90*numel(idx))));
items = {
    S.map_x_col,S.map_y_col,S.map_I_col,S.map_c_col,S.offv_collimated_static(idx,:),'(a)';
    S.map_x_col,S.map_y_col,S.map_I_col,S.map_c_col,S.offv_collimated_adaptive(idx,:),'(b)';
    S.map_x_foc,S.map_y_foc,S.map_I_foc,S.map_c_foc,S.offv_collimated_static(idx,:),'(c)';
    S.map_x_foc,S.map_y_foc,S.map_I_foc,S.map_c_foc,S.offv_point_adaptive(idx,:),'(d)'};
for k=1:4
    ax=nexttile(t);
    x=items{k,1}; y=items{k,2}; I=items{k,3}; c=items{k,4}; o=items{k,5};
    imagesc(ax,x,y,I); axis(ax,'xy');
    set(ax,'ColorScale','log'); clim(ax,[0.02 80]); colormap(ax,turbo);
    xlim(ax,[x(1) x(end)]); ylim(ax,[y(1) y(end)]); hold(ax,'on');
    cc=c+o/100;
    apertureCircle(ax,x,y,cc,0.05);
    text(ax,0.04,0.07,sprintf('偏移 %.1f mm',norm(o)*10),'Units','normalized','Color','w', ...
        'FontSize',11,'BackgroundColor','k');
    set(ax,'FontSize',12); atlas_panel(ax,items{k,6},-0.13);
end
cb=colorbar(ax); cb.Label.String='辐照度 (W/m^2)'; cb.Label.FontSize=12;
% (e) time series
ax=nexttile(t,[1 2]); hold(ax,'on');
plot(ax,S.ts_t,S.ts_p(:,3),'Color',C_LIN,'LineWidth',1.4);
plot(ax,S.ts_t,S.ts_p(:,4)*S.focus_gain,'Color',C_FOC,'LineWidth',1.4);
xlim(ax,[min(S.ts_t) max(S.ts_t)]);
xlabel(ax,{'时间 t (s)';'(e)'},'FontSize',16); ylabel(ax,'接收功率 (mW)','FontSize',16);
set(ax,'FontSize',13); ylim(ax,[0 24]);
lg=legend(ax,{'RIS线性调相','RIS汇聚调相'},'Location','southwest'); set(lg,'Box','off','FontSize',11);
atlas_save(fig,'fig07_motion_spots');
end

function apertureCircle(ax,x,y,c,r)
dx=x(end)-x(1); dy=y(end)-y(1);
w=2*r; h=2*r;
if dx>dy, w=2*r*dy/dx; else, h=2*r*dx/dy; end
rectangle(ax,'Position',[c(1)-w/2,c(2)-h/2,w,h],'Curvature',[1 1], ...
    'EdgeColor','w','LineWidth',1.6);
end
