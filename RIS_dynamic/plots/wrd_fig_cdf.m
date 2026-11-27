function paths = wrd_fig_cdf(results,opts)
%WRD_FIG_CDF Fig. 2: empirical CDF of RX power over the saved frames.
%   F(P) = fraction of frames with RX power <= P. The value of a curve at the
%   dotted threshold line is the outage probability listed in
%   wris_dynamic_cases.csv. Zero-power frames are placed on the floor.
shape = size(results.cases);
links = wrd_fig_links(results);
[floorW,topW] = wrd_fig_range(results,links);
threshold = results.config.metrics.outageThresholdW;
fig = wrd_fig_new([max(560,430*shape(2)),290*shape(1)+90],opts);
layout = tiledlayout(fig,shape(1),shape(2),'TileSpacing','compact','Padding','compact');
for row = 1:shape(1)
    for col = 1:shape(2)
        item = results.cases{row,col};
        f = item.frames;
        ax = nexttile(layout);
        hold(ax,'on');
        handles = gobjects(1,numel(links));
        for j = 1:numel(links)
            link = [f.(links(j).key)];
            p = sort(max([link.rxPowerW],floorW));
            F = (1:numel(p))/numel(p);
            handles(j) = stairs(ax,[p(1),p,topW]*1e3,[0,F,1], ...
                'LineStyle',links(j).style,'Color',links(j).color,'LineWidth',1.4);
        end
        hThreshold = xline(ax,threshold*1e3,':','Color',[0.2 0.2 0.2],'LineWidth',1);
        set(ax,'XScale','log');
        xlim(ax,1e3*[floorW/1.6 topW]);
        ylim(ax,[0 1]);
        xlabel(ax,'RX power P (mW)');
        ylabel(ax,'Pr(P_{rx} \leq P)');
        title(ax,wrd_fig_title(item),'FontWeight','normal');
        grid(ax,'on');
        ax.XMinorGrid = 'off';
        if row==1 && col==1
            lgd = legend(ax,[handles,hThreshold],[{links.label},{'Outage threshold'}], ...
                'Orientation','horizontal','NumColumns',min(4,numel(links)+1));
            lgd.Layout.Tile = 'south';
            lgd.Box = 'off';
            lgd.FontName = opts.fontName;
            lgd.FontSize = opts.fontSize;
        end
    end
end
title(layout,'Distribution of RX power over time (dynamic sea)','FontName',opts.fontName, ...
    'FontSize',opts.fontSize+1);
paths = wrd_fig_save(fig,'fig2_rx_power_cdf',opts);
end
