function files = rd_fig_irradiance(R,o)
%RD_FIG_IRRADIANCE Fig1: receiver-plane irradiance of the six links at the main condition.
%   One frozen sea instant, one RX plane. Every panel is centred on the
%   common receiver (cm); colour = log10 irradiance on ONE scale for all
%   panels, so the gain of RIS modulation over the no-RIS link is visible.
%   Direct: MC through the rough surface; RIS designs: complex field
%   (ballistic Fresnel/CZT + scattered halo), the designs picked at the
%   design power with the main fading row. Solid circle: RX aperture;
%   dashed ellipse: second-moment spot waist (RIS_phase eq. 70).
iz = rd_fig_index(R.distancesM,o.mainDistanceM);
iw = rd_fig_index(R.winds,o.mainWindMps);
q = R.cond(iz,iw);
maps = q.maps;
center = q.targetXY;
bins = 201;
axisM = linspace(-o.mapHalfWidthM,o.mapHalfWidthM,bins);
[X,Y] = meshgrid(center(1)+axisM,center(2)+axisM);
shown = cell(1,numel(maps));
peak = 0;
for j = 1:numel(maps)
    I = interp2(maps(j).x,maps(j).y,double(maps(j).field),X,Y,'linear',0);
    shown{j} = I;
    peak = max(peak,max(I(:)));
end
files = strings(0,1);
if ~(peak>0), return; end
cmax = log10(peak);
cmin = cmax-o.mapDecades;
a = R.config.rx.apertureRadius;
theta = linspace(0,2*pi,181);
[fig,layout] = rd_fig_new(2,3,o,1.85,1.85);
for j = 1:numel(maps)
    ax = nexttile(layout);
    imagesc(ax,100*axisM,100*axisM,log10(max(shown{j},10^cmin)));
    set(ax,'YDir','normal');
    axis(ax,'equal'); axis(ax,'tight');
    clim(ax,[cmin cmax]);
    colormap(ax,o.mapColormap);
    hold(ax,'on');
    plot(ax,100*a*cos(theta),100*a*sin(theta),'w-','LineWidth',1.0);
    if ~isempty(maps(j).spotCovariance)
        [V,L] = eig((maps(j).spotCovariance+maps(j).spotCovariance.')/2);
        ring = (V*diag(2*sqrt(max(diag(L),0)))*[cos(theta);sin(theta)]).'+maps(j).spotCenter-center;
        plot(ax,100*ring(:,1),100*ring(:,2),'w--','LineWidth',0.8);
    end
    set(ax,'FontName',o.fontName,'FontSize',o.fontSize,'Layer','top', ...
        'XLim',100*o.mapHalfWidthM*[-1 1],'YLim',100*o.mapHalfWidthM*[-1 1]);
    title(ax,maps(j).label,'FontName',o.fontName,'FontSize',o.fontSize,'FontWeight','normal');
    text(ax,0.04,0.95,sprintf('\\langleG\\rangle = %.2e',maps(j).meanGainField),'Units','normalized', ...
        'Color','w','FontName',o.fontName,'FontSize',o.fontSize-1,'VerticalAlignment','top');
    if mod(j,3)==1, ylabel(ax,'y (cm)'); end
    if j>3, xlabel(ax,'x (cm)'); end
end
cb = colorbar(ax);
cb.Layout.Tile = 'east';
cb.Label.String = 'log_{10} irradiance (W m^{-2})';
cb.FontName = o.fontName;
cb.FontSize = o.fontSize;
title(layout,sprintf('z = %g m, U = %g m/s, frozen t^* = %.1f s, P_t = 1 W reference', ...
    q.distanceM,q.wind,q.pose.timeS),'FontName',o.fontName,'FontSize',o.fontSize);
files = rd_fig_save(fig,'fig1_irradiance_rx',o);
end
