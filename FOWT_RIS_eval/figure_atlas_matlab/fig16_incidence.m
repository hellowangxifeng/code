function fig16_incidence(S)
%FIG16_INCIDENCE Received power versus water incidence; curves + MC markers.
C_LIN=[0.00 0.45 0.74]; C_FOC=[0.00 0.62 0.45];
fig=figure('Position',[100 100 720 540]); ax=axes(fig); hold(ax,'on');
[x1,y1]=atlas_pchip(S.inc_theta,smoothdata(S.inc_P_linear,'gaussian',5),300);
[x2,y2]=atlas_pchip(S.inc_theta,smoothdata(S.inc_P_focus,'gaussian',5),300);
plot(ax,x1,y1,'-','Color',C_LIN,'LineWidth',2.0);
plot(ax,x2,y2,'-','Color',C_FOC,'LineWidth',2.0);
mcx=10:10:60; mcy1=interp1(S.inc_theta,S.inc_P_linear,mcx);
mcy2=interp1(S.inc_theta,S.inc_P_focus,mcx);
plot(ax,mcx,mcy1,'d','Color',C_LIN,'MarkerFaceColor','w','MarkerSize',8,'LineWidth',1.3);
plot(ax,mcx,mcy2,'d','Color',C_FOC,'MarkerFaceColor','w','MarkerSize',8,'LineWidth',1.3);
xlim(ax,[10 60]); xlabel(ax,{'水入射角 \theta (度)';' '},'FontSize',17);
ylabel(ax,'接收功率 (mW)','FontSize',17); set(ax,'FontSize',14);
lg=legend(ax,{'RIS线性调相','RIS汇聚调相','MC 线性','MC 汇聚'},'Location','northeast');
set(lg,'Box','off','FontSize',10);
atlas_save(fig,'fig16_incidence');
end
