function paths = wrd_fig_power(results,opts)
%WRD_FIG_POWER Fig. 1: RX aperture power vs time for every link.
%   One tile per (angle row, wind column). Error bars are the 95 % MC
%   sampling interval of each frame; frames without any accepted packet
%   are drawn as crosses on the floor line.
shape = size(results.cases);
links = wrd_fig_links(results);
[floorW,topW] = wrd_fig_range(results,links);
threshold = results.config.metrics.outageThresholdW;
fig = wrd_fig_new([max(560,430*shape(2)),300*shape(1)+90],opts);
layout = tiledlayout(fig,shape(1),shape(2),'TileSpacing','compact','Padding','compact');
for row = 1:shape(1)
    for col = 1:shape(2)
        item = results.cases{row,col};
        f = item.frames;
        t = [f.t];
        ax = nexttile(layout);
        hold(ax,'on');
        handles = gobjects(1,numel(links));
        for j = 1:numel(links)
            handles(j) = powerTrace(ax,t,[f.(links(j).key)],floorW,links(j));
        end
        hThreshold = yline(ax,threshold*1e3,':','Color',[0.2 0.2 0.2],'LineWidth',1);
        hMiss = plot(ax,NaN,NaN,'x','Color',[0.25 0.25 0.25],'MarkerSize',5);
        set(ax,'YScale','log');
        ylim(ax,1e3*[floorW/1.6 topW]);
        if numel(t)>1, xlim(ax,[t(1) t(end)]); end
        xlabel(ax,'Elapsed time (s)');
        ylabel(ax,'RX power (mW)');
        title(ax,wrd_fig_title(item),'FontWeight','normal');
        grid(ax,'on');
        ax.YMinorGrid = 'off';
        if row==1 && col==1
            lgd = legend(ax,[handles,hThreshold,hMiss], ...
                [{links.label},{'Outage threshold','No packet accepted'}], ...
                'Orientation','horizontal','NumColumns',min(3,numel(links)+2));
            lgd.Layout.Tile = 'south';
            lgd.Box = 'off';
            lgd.FontName = opts.fontName;
            lgd.FontSize = opts.fontSize;
        end
    end
end
title(layout,'Same sea realization, same TX motion, fixed receiver','FontName',opts.fontName, ...
    'FontSize',opts.fontSize+1);
paths = wrd_fig_save(fig,'fig1_rx_power_time',opts);
end

function h = powerTrace(ax,t,link,floorW,style)
value = [link.rxPowerW];
ci = 1.96*[link.powerSEW];
hasHit = [link.hitCount]>0;
shown = value;
shown(~hasHit) = NaN;
lower = min(ci,max(0,value-floorW));
h = errorbar(ax,t,shown*1e3,lower*1e3,ci*1e3,'LineStyle',style.style,'Marker',style.marker, ...
    'Color',style.color,'LineWidth',1.15,'MarkerSize',9,'CapSize',0);
if any(~hasHit)
    plot(ax,t(~hasHit),repmat(floorW*1e3,1,nnz(~hasHit)),'x','Color',style.color, ...
        'MarkerSize',5,'HandleVisibility','off');
end
end
