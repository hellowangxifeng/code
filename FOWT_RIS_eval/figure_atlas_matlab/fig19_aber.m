function fig19_aber(S)
%FIG19_ABER Average BER versus amplitude, frequency, incidence and wind.
C_LIN=[0.00 0.45 0.74]; C_FOC=[0.00 0.62 0.45];
fig=figure('Position',[60 60 1000 860]);
t=tiledlayout(fig,2,2,'TileSpacing','loose','Padding','loose');
ax=nexttile(t); condPlot(ax,S.wave_x,S.wave_ABER_linear_c151,S.wave_ABER_focus_c151, ...
    S.wave_ABER_linear_c25,S.wave_ABER_focus_c25,S.wave_ABER_linear_c398,S.wave_ABER_focus_c398,C_LIN,C_FOC);
xlabel(ax,{'波幅 A (mm)';' '},'FontSize',16); ylabel(ax,'平均误比特率 ABER','FontSize',16);
set(ax,'FontSize',13); atlas_panel(ax,'(a)' ,-0.17); xlim(ax,[10 50]);
lg=legend(ax,{'线性 c=0.151','汇聚 c=0.151','线性 c=0.25','汇聚 c=0.25', ...
    '线性 c=0.398','汇聚 c=0.398'}, ...
    'Location','southwest'); set(lg,'Box','off','FontSize',7.5,'NumColumns',3);
ax=nexttile(t); condPlot(ax,S.freq_x,S.freq_ABER_linear_c151,S.freq_ABER_focus_c151, ...
    S.freq_ABER_linear_c25,S.freq_ABER_focus_c25,S.freq_ABER_linear_c398,S.freq_ABER_focus_c398,C_LIN,C_FOC);
xlabel(ax,{'波频 f (Hz)';' '},'FontSize',16); ylabel(ax,'平均误比特率 ABER','FontSize',16);
set(ax,'FontSize',13); atlas_panel(ax,'(b)' ,-0.17); xlim(ax,[0.3 1.3]);
ax=nexttile(t); condPlot(ax,S.inc_theta,S.inc_ABER_linear_c151,S.inc_ABER_focus_c151, ...
    S.inc_ABER_linear_c25,S.inc_ABER_focus_c25,S.inc_ABER_linear_c398,S.inc_ABER_focus_c398,C_LIN,C_FOC);
xlabel(ax,{'水入射角 \theta (度)';' '},'FontSize',16); ylabel(ax,'平均误比特率 ABER','FontSize',16);
set(ax,'FontSize',13); atlas_panel(ax,'(c)' ,-0.17); xlim(ax,[10 60]);
ax=nexttile(t); condPlot(ax,S.wind_U,S.wind_ABER_linear_c151,S.wind_ABER_focus_c151, ...
    S.wind_ABER_linear_c25,S.wind_ABER_focus_c25,S.wind_ABER_linear_c398,S.wind_ABER_focus_c398,C_LIN,C_FOC);
xlabel(ax,{'风速 U (m/s)';' '},'FontSize',16); ylabel(ax,'平均误比特率 ABER','FontSize',16);
set(ax,'FontSize',13); atlas_panel(ax,'(d)' ,-0.17); xlim(ax,[3 12]);
atlas_save(fig,'fig19_aber_conditions');
end

function condPlot(ax,x,y1,y2,y3,y4,y5,y6,c1,c2)
hold(ax,'on'); set(ax,'YScale','log');
semilogy_s(ax,x,y1,c1,'-'); semilogy_s(ax,x,y2,c2,'-');
semilogy_s(ax,x,y3,c1,'-.'); semilogy_s(ax,x,y4,c2,'-.');
semilogy_s(ax,x,y5,c1,'--'); semilogy_s(ax,x,y6,c2,'--');
ylim(ax,[1e-18 0.5]);
end

function semilogy_s(ax,x,y,c,ls)
[xx,yy]=atlas_pchip(x,log10(max(y(:),1e-16)),250);
semilogy(ax,xx,10.^yy,'Color',c,'LineStyle',ls,'LineWidth',1.8);
end
