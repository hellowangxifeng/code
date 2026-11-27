function fig09_turb(S)
%FIG09_TURB Outage, ABER and normalized power versus turbulence strength.
C_LIN=[0.00 0.45 0.74]; C_FOC=[0.00 0.62 0.45];
fig=figure('Position',[80 80 1420 540]);
t=tiledlayout(fig,1,3,'TileSpacing','loose','Padding','loose');
ax=nexttile(t); hold(ax,'on');
semilogy_s(ax,S.turb_sigma,S.turb_OP_linear_c151,C_LIN,'-');
semilogy_s(ax,S.turb_sigma,S.turb_OP_focus_c151,C_FOC,'-');
semilogy_s(ax,S.turb_sigma,S.turb_OP_linear_c25,C_LIN,'-.');
semilogy_s(ax,S.turb_sigma,S.turb_OP_focus_c25,C_FOC,'-.');
semilogy_s(ax,S.turb_sigma,S.turb_OP_linear_c398,C_LIN,'--');
semilogy_s(ax,S.turb_sigma,S.turb_OP_focus_c398,C_FOC,'--');
xlim(ax,[0.4 1.0]); ylim(ax,[1e-10 1.3]);
xlabel(ax,{'对数强度标准差 \sigma_{ln I}','(a)'},'FontSize',17); ylabel(ax,'中断概率 P_{out}','FontSize',17);
set(ax,'FontSize',14); axa=ax;
lg=legend(ax,{'线性 c=0.151','汇聚 c=0.151','线性 c=0.25','汇聚 c=0.25', ...
    '线性 c=0.398','汇聚 c=0.398'},'Location','southeast');
set(lg,'Box','off','FontSize',7.5,'NumColumns',3);
ax=nexttile(t); hold(ax,'on');
semilogy_s(ax,S.turb_sigma,S.turb_ABER_linear_c151,C_LIN,'-');
semilogy_s(ax,S.turb_sigma,S.turb_ABER_focus_c151,C_FOC,'-');
semilogy_s(ax,S.turb_sigma,S.turb_ABER_linear_c25,C_LIN,'-.');
semilogy_s(ax,S.turb_sigma,S.turb_ABER_focus_c25,C_FOC,'-.');
semilogy_s(ax,S.turb_sigma,S.turb_ABER_linear_c398,C_LIN,'--');
semilogy_s(ax,S.turb_sigma,S.turb_ABER_focus_c398,C_FOC,'--');
xlim(ax,[0.4 1.0]); ylim(ax,[1e-18 0.5]);
xlabel(ax,{'对数强度标准差 \sigma_{ln I}','(b)'},'FontSize',17); ylabel(ax,'平均误比特率 ABER','FontSize',17);
set(ax,'FontSize',14); axb=ax;
% (c) normalized mean / 10% power penalty of a log-normal channel
ax=nexttile(t); hold(ax,'on');
sg=S.turb_sigma;
plot(ax,sg,exp(-0.5*sg.^2),'-','Color',C_LIN,'LineWidth',2.0);
plot(ax,sg,exp(-1.2816*sg),'--','Color',[0.85 0.33 0.10],'LineWidth',1.9);
xlim(ax,[0.4 1.0]); ylim(ax,[0 1.05]);
xlabel(ax,{'对数强度标准差 \sigma_{ln I}','(c)'},'FontSize',17); ylabel(ax,'归一化接收功率','FontSize',17);
set(ax,'FontSize',14);
lg=legend(ax,{'均值','10% 分位'},'Location','southwest'); set(lg,'Box','off','FontSize',11);
set(axa,'YScale','log'); set(axb,'YScale','log');
atlas_save(fig,'fig09_turbulence');
end

function semilogy_s(ax,x,y,col,ls)
[xx,yy]=atlas_pchip(x,log10(max(y(:),1e-16)),300);
semilogy(ax,xx,10.^yy,'Color',col,'LineStyle',ls,'LineWidth',1.8);
end
