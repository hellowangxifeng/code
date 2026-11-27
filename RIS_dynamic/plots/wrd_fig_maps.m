function paths = wrd_fig_maps(results,opts)
%WRD_FIG_MAPS Fig. 4: irradiance on the common detector plane at one time.
%   Rows: TX angles. Columns: links. All tiles share one colour scale. The
%   map is the irradiance BEFORE the FOV gate; the white circle is the
%   physical RX aperture (same position and size for every link).
items = results.cases(:,opts.mapCondition);
links = wrd_fig_links(results);
maps = cell(numel(items),numel(links));
times = zeros(numel(items),1);
maxValue = 0;
for row = 1:numel(items)
    f = items{row}.frames;
    [~,k] = min(abs([f.t]-opts.mapTimeS));
    times(row) = f(k).t;
    for j = 1:numel(links)
        maps{row,j} = f(k).(links(j).key);
        maxValue = max(maxValue,double(max(maps{row,j}.mapWm2,[],'all')));
    end
end
maxValue = max(maxValue,1e-7);
minimum = maxValue*1e-5;
radius = results.config.receiver.radius;
fig = wrd_fig_new([330*numel(links)+150,320*numel(items)+60],opts);
layout = tiledlayout(fig,numel(items),numel(links),'TileSpacing','compact','Padding','compact');
colormap(fig,parula(256));
angle = linspace(0,2*pi,161);
for row = 1:numel(items)
    for j = 1:numel(links)
        ax = nexttile(layout);
        link = maps{row,j};
        imagesc(ax,link.axisM,link.axisM,10*log10(max(minimum,double(link.mapWm2))));
        axis(ax,'image');
        ax.YDir = 'normal';
        clim(ax,10*log10([minimum maxValue]));
        hold(ax,'on');
        plot(ax,radius*cos(angle),radius*sin(angle),'w','LineWidth',1.2);
        xlabel(ax,'Detector u (m)');
        ylabel(ax,'Detector v (m)');
        title(ax,{wrd_fig_title(items{row}), ...
            sprintf('%s, t = %.1f s, P_{rx} = %.3g mW',links(j).label,times(row),1e3*link.rxPowerW)}, ...
            'FontWeight','normal');
    end
end
cb = colorbar(ax);
cb.Layout.Tile = 'east';
cb.Label.String = 'Irradiance before FOV gate (dB re 1 W m^{-2})';
cb.FontName = opts.fontName;
cb.FontSize = opts.fontSize;
title(layout,'Irradiance on the common detector plane (white circle: RX aperture)', ...
    'FontName',opts.fontName,'FontSize',opts.fontSize+1);
paths = wrd_fig_save(fig,'fig4_rx_maps',opts);
end
