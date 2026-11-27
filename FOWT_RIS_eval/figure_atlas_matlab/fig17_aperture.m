function fig17_aperture(S)
%FIG17_APERTURE Aperture power for focused circular, aligned elliptical and
%   unmodulated broad spots (synthetic beam-shape model, 15 mW total).
fig=figure('Position',[100 100 760 560]); ax=axes(fig); hold(ax,'on');
n=501; L=0.45; x=linspace(-L,L,n); [X,Y]=meshgrid(x,x); dx=x(2)-x(1);
% focused circular: 90% core + 10% ring
core=exp(-2*(X.^2+Y.^2)/0.03^2);
ring=exp(-((sqrt(X.^2+Y.^2)-0.08)/0.018).^2);
If=0.9*core+0.1*ring; If=If/sum(If(:))/dx^2;
% non-focused elliptical
Ie=exp(-2*((X/0.04).^2+(Y/0.10).^2)); Ie=Ie/sum(Ie(:))/dx^2;
% no-phase broad
Id=exp(-2*(X.^2+Y.^2)/0.15^2); Id=Id/sum(Id(:))/dx^2;
r=linspace(0.5,15,80); Pf=zeros(size(r)); Pe=Pf; Pd=Pf;
for k=1:numel(r)
    m=(X.^2+Y.^2)<=(r(k)/100)^2;
    Pf(k)=sum(If(m))*dx^2; Pe(k)=sum(Ie(m))*dx^2; Pd(k)=sum(Id(m))*dx^2;
end
Pf=Pf*15; Pe=Pe*15; Pd=Pd*15;
plot(ax,r,Pf,'-','Color',[0.00 0.62 0.45],'LineWidth',2.0);
plot(ax,r,Pe,'--','Color',[0.85 0.33 0.10],'LineWidth',1.9);
plot(ax,r,Pd,':','Color',[0.45 0.45 0.45],'LineWidth',1.9);
xlim(ax,[0 15]); ylim(ax,[0 16]);
xlabel(ax,{'接收孔径半径 a_r (cm)';' '},'FontSize',17); ylabel(ax,'接收功率 (mW)','FontSize',17);
set(ax,'FontSize',14);
lg=legend(ax,{'RIS汇聚（圆斑）','RIS线性（椭圆斑）','无调相（宽斑）'},'Location','southeast');
set(lg,'Box','off','FontSize',10);
atlas_save(fig,'fig17_aperture');
end
