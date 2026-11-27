function files = rd_fig_irradiance(R,o)
%RD_FIG_IRRADIANCE Fig1: receiver-plane irradiance, calm vs rough sea.
%   Rows: calm sea; one rough-sea frame; rough-sea time mean.
%   Columns: o.mapSchemes. Axes are centred on each link's fixed receiver.
%   Solid circle: receiver aperture; dashed circle: RMS position jitter
%   radius sqrt(2) sigma_p. Colour: log10 irradiance on a common scale.
files = strings(0,1);
keys = R.schemes.keys;
wanted = o.mapSchemes(ismember(o.mapSchemes,keys));
ks = cellfun(@(key) find(strcmp(keys,key),1),wanted);
if isempty(ks), return; end
iz = rd_fig_index(R.distancesM,o.mapDistanceM);
calm = find(R.winds==0,1);
windy = find(R.winds>0);
rows = struct('iw',calm,'field','mean','label','Calm sea');
if ~isempty(windy)
    iw = windy(rd_fig_index(R.winds(windy),o.mapWindMps));
    rows(2) = struct('iw',iw,'field','snapshot','label', ...
        sprintf('U = %g m/s, t = %.1f s',R.winds(iw),R.main.maps{iw}.snapshotTimeS));
    rows(3) = struct('iw',iw,'field','mean','label',sprintf('U = %g m/s, time mean',R.winds(iw)));
end
peak = -Inf;
for r = 1:numel(rows)
    maps = R.main.maps{rows(r).iw}.(rows(r).field);
    if isempty(maps), continue; end
    values = maps(:,:,iz,ks);
    peak = max(peak,max(double(values(:))));
end
if ~(peak>0), return; end
cmax = log10(peak);
cmin = cmax-o.mapDecades;
a = R.config.rx.apertureRadius;
sigma = R.config.study.positionSigmaM;
theta = linspace(0,2*pi,181);
[fig,layout] = rd_fig_new(numel(rows),numel(ks),o,1.75,1.75);
for r = 1:numel(rows)
    block = R.main.maps{rows(r).iw};
    for jc = 1:numel(ks)
        k = ks(jc);
        ax = nexttile(layout);
        if isempty(block.(rows(r).field)), axis(ax,'off'); continue; end
        I = double(block.(rows(r).field)(:,:,iz,k));
        half = block.halfWidthM(iz,k);
        bins = size(I,1);
        axisCM = 100*(((1:bins)-0.5)/bins*2*half-half);
        imagesc(ax,axisCM,axisCM,log10(max(I,10^cmin)));
        set(ax,'YDir','normal');
        axis(ax,'equal');
        axis(ax,'tight');
        clim(ax,[cmin cmax]);
        colormap(ax,o.mapColormap);
        hold(ax,'on');
        plot(ax,100*a*cos(theta),100*a*sin(theta),'w-','LineWidth',0.9);
        plot(ax,100*sqrt(2)*sigma*cos(theta),100*sqrt(2)*sigma*sin(theta),'w--','LineWidth',0.8);
        set(ax,'FontName',o.fontName,'FontSize',o.fontSize,'Layer','top');
        if r==1
            title(ax,R.schemes.names{k},'FontName',o.fontName,'FontSize',o.fontSize, ...
                'FontWeight','normal');
        end
        if jc==1
            ylabel(ax,{rows(r).label,'y (cm)'},'FontName',o.fontName,'FontSize',o.fontSize);
        end
        if r==numel(rows)
            xlabel(ax,'x (cm)','FontName',o.fontName,'FontSize',o.fontSize);
        end
    end
end
cb = colorbar(ax);
cb.Layout.Tile = 'east';
cb.Label.String = 'log_{10} irradiance (W m^{-2})';
cb.FontName = o.fontName;
cb.FontSize = o.fontSize;
title(layout,sprintf('Receiver plane z = %g m',R.distancesM(iz)),'FontName',o.fontName, ...
    'FontSize',o.fontSize+1);
files = rd_fig_save(fig,'fig1_irradiance_maps',o);
end
