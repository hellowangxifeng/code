function fig08_distance(S)
%FIG08_DISTANCE Spot growth with air path and distance curves.
C_LIN=[0.00 0.45 0.74]; C_FOC=[0.00 0.62 0.45]; C_DIR=[0.45 0.45 0.45];
C3=[0.00 0.45 0.74; 0.85 0.33 0.10; 0.00 0.62 0.45];
fig=figure('Position',[60 60 1250 830]);
t=tiledlayout(fig,2,3,'TileSpacing','loose','Padding','loose');
R0=7.9; c=S.map_c_foc;
for k=1:3
    ax=nexttile(t);
    R=[5 10 15]; sc=R(k)/R0;
    xs=c(1)+(S.map_x_foc-c(1))*sc; ys=c(2)+(S.map_y_foc-c(2))*sc;
    I=S.map_I_foc/sc^2;
    imagesc(ax,xs,ys,I); axis(ax,'xy');
    set(ax,'ColorScale','log'); clim(ax,[0.02 80]); colormap(ax,turbo);
    xlim(ax,[xs(1) xs(end)]); ylim(ax,[ys(1) ys(end)]);
    hold(ax,'on'); rectangle(ax,'Position',[c(1)-0.05,c(2)-0.05,0.10,0.10], ...
        'Curvature',[1 1],'EdgeColor','w','LineWidth',1.6);
    set(ax,'FontSize',12); atlas_panel(ax,sprintf('(%c)',char('a'+k-1)),-0.13);
end
cb=colorbar(ax); cb.Label.String='辐照度 (W/m^2)'; cb.Label.FontSize=12;
% (d) power vs distance
ax=nexttile(t); hold(ax,'on');
semilogy_s(ax,S.dist_R,S.dist_P_linear,C_LIN,'-');
semilogy_s(ax,S.dist_R,S.dist_P_focus,C_FOC,'-');
semilogy(ax,S.dist_R,ones(size(S.dist_R))*1e-3,':','Color',C_DIR,'LineWidth',1.4);
xlabel(ax,{'空气路径 L (m)';' '},'FontSize',16); ylabel(ax,'接收功率 (mW)','FontSize',16);
set(ax,'FontSize',13); atlas_panel(ax,'(d)',-0.22);
lg=legend(ax,{'RIS线性调相','RIS汇聚调相','无RIS直接'},'Location','northeast');
set(lg,'Box','off','FontSize',10);
% (e) capture fraction for three initial waists
ax=nexttile(t); hold(ax,'on');
w0s=[0.10 0.144 0.20];
for k=1:3
    w=w0s(k)*S.dist_R/R0;
    frac=1-exp(-2*0.05^2./w.^2);
    plot(ax,S.dist_R,frac*100,'-','Color',C3(k,:),'LineWidth',1.9);
end
xlabel(ax,{'空气路径 L (m)';' '},'FontSize',16); ylabel(ax,'5 cm 捕获率 (%)','FontSize',16);
set(ax,'FontSize',13); atlas_panel(ax,'(e)',-0.22);
lg=legend(ax,{'w_0 = 10 cm','w_0 = 14.4 cm','w_0 = 20 cm'},'Location','northeast');
set(lg,'Box','off','FontSize',10);
% (f) capture power for three aperture radii
ax=nexttile(t); hold(ax,'on');
for k=1:3
    a=[2.5 5 10]; a=a(k)/100;
    w=0.144*S.dist_R/R0;
    ratio=(1-exp(-2*a^2./w.^2))./(1-exp(-2*0.05^2./w.^2));
    plot(ax,S.dist_R,S.dist_P_focus.*ratio,'-','Color',C3(k,:),'LineWidth',1.9);
end
xlabel(ax,{'空气路径 L (m)';' '},'FontSize',16); ylabel(ax,'接收功率 (mW)','FontSize',16);
set(ax,'FontSize',13); atlas_panel(ax,'(f)',-0.22);
lg=legend(ax,{'a_r = 2.5 cm','a_r = 5 cm','a_r = 10 cm'},'Location','northeast');
set(lg,'Box','off','FontSize',10);
atlas_save(fig,'fig08_distance');
end

function semilogy_s(ax,x,y,col,ls)
[xx,yy]=atlas_pchip(x,log10(max(y(:),1e-12)),300);
semilogy(ax,xx,10.^yy,'Color',col,'LineStyle',ls,'LineWidth',1.7);
end
