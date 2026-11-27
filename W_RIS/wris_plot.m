function files = wris_plot(s)
%WRIS_PLOT Render SAVED results only; never invokes motion/optical solvers.
% Set s.config.plot fields and call again to change fonts, layout or quantity.
c=s.config; opts=c.plot;
quantity=validatestring(opts.quantity,{'noTurbulence','irradiance'});
validateattributes(opts.rowsPerFigure,{'numeric'},{'scalar','integer','positive'});
validateattributes(opts.dynamicRange,{'numeric'},{'scalar','finite','>',1});
[nAngles,nColumns]=size(s.panels);
peak=max(cellfun(@(p) max(p.out.(quantity),[],'all'),s.panels),[],'all');
limits=[max(peak,eps)/opts.dynamicRange max(peak,eps)];
nPages=ceil(nAngles/opts.rowsPerFigure); files=strings(nPages,2);
if ~isfolder(c.output), mkdir(c.output); end
for page=1:nPages
    rows=(page-1)*opts.rowsPerFigure+1:min(page*opts.rowsPerFigure,nAngles);
    f=figure('Name',sprintf('RIS lower-face irradiance %d/%d',page,nPages), ...
        'Visible','off','Color','w','Units','pixels','Position',[20 20 1250 980]);
    guard=onCleanup(@() closeIfHidden(f,opts.visible));
    layout=tiledlayout(f,numel(rows),nColumns,'TileSpacing','compact','Padding','loose');
    for row=rows
        reference=s.panels{row,1}.out.geometry.projectedOrientationDeg;
        for col=1:nColumns
            p=s.panels{row,col}; o=p.out; ax=nexttile(layout);
            imagesc(ax,o.x,o.y,max(o.(quantity),limits(1)).');
            axis(ax,'xy','equal');
            extent=(0.5+opts.padding)*[c.ris.Lx c.ris.Ly];
            xlim(ax,[-extent(1) extent(1)]); ylim(ax,[-extent(2) extent(2)]);
            set(ax,'ColorScale','log','Color',[0.78 0.78 0.78], ...
                'FontSize',opts.tickFontSize,'LineWidth',1,'TickDir','out','Tag','RISMainAxes');
            xticks(ax,[-c.ris.Lx/2 0 c.ris.Lx/2]); yticks(ax,[-c.ris.Ly/2 0 c.ris.Ly/2]);
            clim(ax,limits); colormap(ax,turbo(256)); hold(ax,'on');
            rectangle(ax,'Position',[-c.ris.Lx/2 -c.ris.Ly/2 c.ris.Lx c.ris.Ly], ...
                'EdgeColor','w','LineWidth',1.5);
            drawCore(ax,o,reference);
            delta=mod(o.geometry.projectedOrientationDeg-reference+90,180)-90;
            a=p.initialRotationXYDeg;
            label=sprintf('[\\alpha_x,\\alpha_y]=[%g^{\\circ},%g^{\\circ}] | fixed',a);
            if col>1, label=sprintf('U=%g m/s | wind=%g^{\\circ}',p.windMps,p.windDirectionDeg); end
            rotationText='n/a';
            if isfinite(delta), rotationText=sprintf('%+.1f^{\\circ}',delta); end
            title(ax,{label,sprintf('yaw=%+.1f^{\\circ} | rot=%s',o.geometry.yawRisDeg,rotationText), ...
                sprintf('\\theta=%.1f^{\\circ} | clipped=%.1f%%',o.geometry.thetaDeg,100*o.coreClippedFraction)}, ...
                'FontSize',opts.titleFontSize,'FontWeight','normal');
            if col==1, ylabel(ax,'y_{RIS} (m)','FontSize',opts.labelFontSize); end
            if row==rows(end), xlabel(ax,'x_{RIS} (m)','FontSize',opts.labelFontSize); end
            if ~o.geometry.axisInside
                % Keep the aperture map fully visible. State the true off-board
                % beam position instead of covering it with an inset heatmap.
                text(ax,0.5,0.02,sprintf('Outside (%.2f,%.2f) m',o.geometry.axisLocal(1:2)), ...
                    'Units','normalized','Color',[0.12 0.12 0.12], ...
                    'HorizontalAlignment','center','VerticalAlignment','bottom','Clipping','on', ...
                    'FontSize',max(10,opts.titleFontSize-2),'Interpreter','none');
            end
        end
    end
    cb=colorbar(ax); cb.Layout.Tile='east'; cb.FontSize=opts.tickFontSize;
    cb.Label.String='I/H (W/m^2)';
    if strcmp(quantity,'irradiance'), cb.Label.String='I (W/m^2), EGG snapshot'; end
    cb.Label.FontSize=opts.labelFontSize;
    title(layout,{sprintf('RIS lower-face irradiance | t=%.1f s | sea t=%.1f s | page %d/%d', ...
        c.sampling.comparisonTime,c.ris.startTime+c.sampling.comparisonTime,page,nPages), ...
        'TX/RIS moving | dashed: complete core | dotted: static axis | frame: RIS'}, ...
        'FontSize',opts.figureTitleFontSize,'FontWeight','normal');
    drawnow;
    base=fullfile(c.output,sprintf('Fig_RIS_rotation_clipping_page%d',page));
    files(page,:)=[string([base '.fig']),string([base '.png'])];
    savefig(f,files(page,1));
    exportgraphics(f,files(page,2),'Resolution',opts.resolution);
    if strcmp(opts.visible,'on'), f.Visible='on'; end
    clear guard
end
end

function drawCore(ax,o,reference)
m=o.directPlaneMetrics;
if ~all(isfinite([m.center m.axes m.orientationDeg])), return; end
angle=deg2rad(m.orientationDeg); R=[cos(angle) -sin(angle);sin(angle) cos(angle)];
parameter=linspace(0,2*pi,181);
ellipse=(R*(m.axes(:)/2.*[cos(parameter);sin(parameter)])).'+m.center;
plot(ax,ellipse(:,1),ellipse(:,2),'w--','LineWidth',1.1);
half=m.axes(1)/2;
if isfinite(reference)
    baseline=half*[cosd(reference);sind(reference)];
    plot(ax,m.center(1)+[-baseline(1) baseline(1)],m.center(2)+[-baseline(2) baseline(2)], ...
        ':','Color',[0.2 0.2 0.2],'LineWidth',1.3);
end
if isfinite(o.geometry.projectedOrientationDeg)
    direction=half*[cosd(o.geometry.projectedOrientationDeg);sind(o.geometry.projectedOrientationDeg)];
    plot(ax,m.center(1)+[-direction(1) direction(1)],m.center(2)+[-direction(2) direction(2)],'w-','LineWidth',1);
end
plot(ax,o.geometry.axisLocal(1),o.geometry.axisLocal(2),'w+','MarkerSize',8,'LineWidth',1.2);
end

function closeIfHidden(f,visibility)
if isgraphics(f) && strcmp(visibility,'off'), close(f); end
end
