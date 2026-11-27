function fig20_distance_turb(S)
%FIG20_DISTANCE_TURB OP/ABER versus distance and turbulence strength.
C_LIN=[0.00 0.45 0.74]; C_FOC=[0.00 0.62 0.45];
fig=figure('Position',[60 60 1000 920]);
t=tiledlayout(fig,2,2,'TileSpacing','loose','Padding','loose');
ax=nexttile(t); condPlot(ax,S.dist_R,S.dist_OP_linear_c151,S.dist_OP_focus_c151, ...
    S.dist_OP_linear_c25,S.dist_OP_focus_c25,S.dist_OP_linear_c398,S.dist_OP_focus_c398,C_LIN,C_FOC,1e-9,1.3);
xlabel(ax,{'距离 L (m)';' '},'FontSize',16); ylabel(ax,'中断概率 P_{out}','FontSize',16);
set(ax,'FontSize',13); atlas_panel(ax,'(a)' ,-0.155); xlim(ax,[5 15]);
lg=legend(ax,{'线性 c=0.151','汇聚 c=0.151','线性 c=0.25','汇聚 c=0.25', ...
    '线性 c=0.398','汇聚 c=0.398'}, ...
    'Location','southeast'); set(lg,'Box','off','FontSize',7.5,'NumColumns',3);
ax=nexttile(t); condPlot(ax,S.dist_R,S.dist_ABER_linear_c151,S.dist_ABER_focus_c151, ...
    S.dist_ABER_linear_c25,S.dist_ABER_focus_c25,S.dist_ABER_linear_c398,S.dist_ABER_focus_c398,C_LIN,C_FOC,1e-18,0.5);
xlabel(ax,{'距离 L (m)';' '},'FontSize',16); ylabel(ax,'平均误比特率 ABER','FontSize',16);
set(ax,'FontSize',13); atlas_panel(ax,'(b)' ,-0.155); xlim(ax,[5 15]);
ax=nexttile(t); condPlot(ax,S.turb_sigma,S.turb_OP_linear_c151,S.turb_OP_focus_c151, ...
    S.turb_OP_linear_c25,S.turb_OP_focus_c25,S.turb_OP_linear_c398,S.turb_OP_focus_c398,C_LIN,C_FOC,1e-9,1.3);
xlabel(ax,{'湍流强度 \sigma_{ln I}';' '},'FontSize',16); ylabel(ax,'中断概率 P_{out}','FontSize',16);
set(ax,'FontSize',13); atlas_panel(ax,'(c)' ,-0.155); xlim(ax,[0.4 1.0]);
ax=nexttile(t); condPlot(ax,S.turb_sigma,S.turb_ABER_linear_c151,S.turb_ABER_focus_c151, ...
    S.turb_ABER_linear_c25,S.turb_ABER_focus_c25,S.turb_ABER_linear_c398,S.turb_ABER_focus_c398,C_LIN,C_FOC,1e-18,0.5);
xlabel(ax,{'湍流强度 \sigma_{ln I}';' '},'FontSize',16); ylabel(ax,'平均误比特率 ABER','FontSize',16);
set(ax,'FontSize',13); atlas_panel(ax,'(d)' ,-0.155); xlim(ax,[0.4 1.0]);
atlas_save(fig,'fig20_distance_turb');
end

function condPlot(ax,x,y1,y2,y3,y4,y5,y6,c1,c2,lo,hi)
hold(ax,'on'); set(ax,'YScale','log');
semilogy_s(ax,x,y1,c1,'-'); semilogy_s(ax,x,y2,c2,'-');
semilogy_s(ax,x,y3,c1,'-.'); semilogy_s(ax,x,y4,c2,'-.');
semilogy_s(ax,x,y5,c1,'--'); semilogy_s(ax,x,y6,c2,'--');
ylim(ax,[lo hi]);
end

function semilogy_s(ax,x,y,c,ls)
[xx,yy]=atlas_pchip(x,log10(max(y(:),1e-16)),250);
semilogy(ax,xx,10.^yy,'Color',c,'LineStyle',ls,'LineWidth',1.8);
end
