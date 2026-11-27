function fig06_spots(S)
%FIG06_SPOTS Linear vs focusing spots, radial power and centre cuts.
C_LIN = [0.00 0.45 0.74]; C_FOC = [0.00 0.62 0.45]; C_DIR = [0.45 0.45 0.45];
fig = figure('Position',[60 60 1000 920]);
t = tiledlayout(fig,2,2,'TileSpacing','loose','Padding','loose');
% (a) linear phase map
ax=nexttile(t);
showMap(ax,S.map_x_col,S.map_y_col,S.map_I_col,S.map_c_col);
atlas_panel(ax,'(a)',-0.13);
% (b) focusing map
ax=nexttile(t);
showMap(ax,S.map_x_foc,S.map_y_foc,S.map_I_foc,S.map_c_foc);
cb=colorbar(ax); cb.Label.String='辐照度 (W/m^2)'; cb.Label.FontSize=13;
atlas_panel(ax,'(b)',-0.13);
% (c) radial aperture power
ax=nexttile(t); hold(ax,'on');
[r,P1]=cumPower(S.map_x_col,S.map_y_col,S.map_I_col,S.map_c_col,0.05);
[~,P2]=cumPower(S.map_x_foc,S.map_y_foc,S.map_I_foc,S.map_c_foc,0.05);
[xx,y1]=atlas_pchip(r,P1,300); [~,y2]=atlas_pchip(r,P2,300);
plot(ax,xx*100,y1*1e3,'Color',C_LIN,'LineWidth',1.9);
plot(ax,xx*100,y2*1e3,'Color',C_FOC,'LineWidth',1.9);
plot(ax,[0 5],[0 0],':','Color',C_DIR,'LineWidth',1.4);
xlim(ax,[0 50]); xlabel(ax,{'半径 (cm)';' '},'FontSize',16); ylabel(ax,'孔径接收功率 (mW)','FontSize',16);
set(ax,'FontSize',13); atlas_panel(ax,'(c)');
lg=legend(ax,{'RIS线性调相','RIS汇聚调相','无RIS直接'},'Location','southeast'); set(lg,'Box','off','FontSize',10);
% (d) centre cuts
ax=nexttile(t); hold(ax,'on');
plotCuts(ax,S.map_x_col,S.map_y_col,S.map_I_col,S.map_c_col,C_LIN);
plotCuts(ax,S.map_x_foc,S.map_y_foc,S.map_I_foc,S.map_c_foc,C_FOC);
xlim(ax,[-30 30]); xlabel(ax,{'相对光斑中心 (cm)';' '},'FontSize',16); ylabel(ax,'归一化辐照度','FontSize',16);
set(ax,'FontSize',13); atlas_panel(ax,'(d)');
lg=legend(ax,{'线性 x','线性 y','汇聚 x','汇聚 y'},'Location','northeast'); set(lg,'Box','off','FontSize',9);
atlas_save(fig,'fig06_spots');
end

function showMap(ax,x,y,I,c)
imagesc(ax,x,y,I); axis(ax,'xy');
set(ax,'ColorScale','log'); clim(ax,[0.02 80]); colormap(ax,turbo);
xlim(ax,[x(1) x(end)]); ylim(ax,[y(1) y(end)]);
hold(ax,'on'); r=0.05;
dx=x(end)-x(1); dy=y(end)-y(1);
w=2*r; h=2*r;
if dx>dy, w=2*r*dy/dx; else, h=2*r*dx/dy; end
rectangle(ax,'Position',[c(1)-w/2,c(2)-h/2,w,h],'Curvature',[1 1],'EdgeColor','w','LineWidth',1.6);
set(ax,'FontSize',12,'LineWidth',1.1,'TickDir','in');
end

function [r,P]=cumPower(x,y,I,c,rmax)
[XX,YY]=meshgrid(x,y); dx=x(2)-x(1); dy=y(2)-y(1);
r=linspace(0,0.5,101); P=zeros(size(r));
for k=1:numel(r)
    m=(XX-c(1)).^2+(YY-c(2)).^2<=r(k)^2;
    P(k)=sum(I(m))*dx*dy;
end
end

function plotCuts(ax,x,y,I,c,col)
[XX,YY]=meshgrid(x,y);
ix=find(abs(x-c(1))==min(abs(x-c(1))),1);
iy=find(abs(y-c(2))==min(abs(y-c(2))),1);
ax1=(x-c(1))*100; ax2=(y-c(2))*100;
plot(ax,ax1,smoothdata(I(iy,:),'gaussian',35)/max(I(:)),'-','Color',col,'LineWidth',1.8);
plot(ax,ax2,smoothdata(I(:,ix),'gaussian',35)/max(I(:)),'--','Color',col,'LineWidth',1.4);
end
