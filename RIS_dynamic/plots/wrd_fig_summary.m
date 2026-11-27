function paths = wrd_fig_summary(results,opts)
%WRD_FIG_SUMMARY Fig. 5: time statistics vs wind speed, one column per TX angle.
%   Top: time-mean RX power, bars span the min..max frame. Bottom: outage
%   probability P(P_rx < threshold). Read from results.caseTable.
T = results.caseTable;
links = wrd_fig_links(results);
keys = struct('flat','flatDirect','direct','roughDirect','ris','ris');
angles = unique(T.nominalIncidenceDeg,'stable');
[floorW,topW] = wrd_fig_range(results,links);
fig = wrd_fig_new([max(520,380*numel(angles)),560],opts);
layout = tiledlayout(fig,2,numel(angles),'TileSpacing','compact','Padding','compact');
axesTop = gobjects(1,numel(angles));
axesBottom = gobjects(1,numel(angles));
for a = 1:numel(angles)
    axesTop(a) = nexttile(layout,a);
    axesBottom(a) = nexttile(layout,a+numel(angles));
end
handles = gobjects(1,numel(links));   % legend entries, taken from the first column
for a = 1:numel(angles)
    for j = 1:numel(links)
        rows = T(abs(T.nominalIncidenceDeg-angles(a))<1e-9 & T.link==keys.(links(j).key),:);
        rows = sortrows(rows,'windMps');
        u = rows.windMps;
        mu = max(rows.meanW,floorW);
        low = max(rows.minW,floorW);
        high = max(rows.maxW,floorW);
        ax = axesTop(a);
        hold(ax,'on');
        h = errorbar(ax,u,mu*1e3,(mu-low)*1e3,(high-mu)*1e3,'LineStyle',links(j).style, ...
            'Marker','o','MarkerSize',5,'Color',links(j).color,'LineWidth',1.3,'CapSize',4);
        if a==1, handles(j) = h; end
        ax = axesBottom(a);
        hold(ax,'on');
        plot(ax,u,rows.outageProbability,'LineStyle',links(j).style,'Marker','o', ...
            'MarkerSize',5,'Color',links(j).color,'LineWidth',1.3);
    end
end
winds = unique(T.windMps);
span = [min(winds)-0.5 max(winds)+0.5];
for a = 1:numel(angles)
    ax = axesTop(a);
    set(ax,'YScale','log');
    ylim(ax,1e3*[floorW/1.6 topW]);
    xlim(ax,span);
    ylabel(ax,'Time-mean RX power (mW)');
    title(ax,sprintf('\\theta_0 = %.1f^{\\circ}',angles(a)),'FontWeight','normal');
    grid(ax,'on');
    ax.YMinorGrid = 'off';
    ax = axesBottom(a);
    xlim(ax,span);
    ylim(ax,[-0.02 1.02]);
    xlabel(ax,'Wind speed U (m/s)');
    ylabel(ax,sprintf('Outage prob. (P_{rx} < %.3g mW)',1e3*results.config.metrics.outageThresholdW));
    grid(ax,'on');
end
lgd = legend(axesTop(1),handles,{links.label},'Orientation','horizontal');
lgd.Layout.Tile = 'south';
lgd.Box = 'off';
lgd.FontName = opts.fontName;
lgd.FontSize = opts.fontSize;
title(layout,'Time statistics over the dynamic sea (bars: min..max frame)', ...
    'FontName',opts.fontName,'FontSize',opts.fontSize+1);
paths = wrd_fig_save(fig,'fig5_summary_wind',opts);
end
