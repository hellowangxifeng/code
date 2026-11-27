function paths = wrd_fig_tir(results,opts)
%WRD_FIG_TIR Fig. 3: why the direct link fluctuates.
%   Left axis: share of the power reaching the water-air boundary that is
%   totally internally reflected (rough sea, and the still-water baseline if
%   saved). Right axis: incident sea slope at the nominal board centre.
shape = size(results.cases);
links = wrd_fig_links(results);
hasFlat = any(strcmp({links.key},'flat'));
red = [0.74 0.29 0.19];
grey = [0.45 0.45 0.45];
blue = [0.10 0.36 0.55];
fig = wrd_fig_new([max(580,430*shape(2)),290*shape(1)+80],opts);
layout = tiledlayout(fig,shape(1),shape(2),'TileSpacing','compact','Padding','compact');
for row = 1:shape(1)
    for col = 1:shape(2)
        item = results.cases{row,col};
        f = item.frames;
        t = [f.t];
        ax = nexttile(layout);
        yyaxis(ax,'left');
        hold(ax,'on');
        direct = [f.direct];
        handles = plot(ax,t,share(direct),'-','Color',red,'LineWidth',1.4);
        labels = {'TIR share, rough sea'};
        if hasFlat
            flat = [f.flat];
            handles(end+1) = plot(ax,t,share(flat),'--','Color',grey,'LineWidth',1.2); %#ok<AGROW>
            labels{end+1} = 'TIR share, still water'; %#ok<AGROW>
        end
        ylabel(ax,'TIR share of surface power (%)');
        ylim(ax,[0 100]);
        yyaxis(ax,'right');
        hold(ax,'on');
        slopes = vertcat(f.slope);
        handles(end+1) = plot(ax,t,atand(slopes(:,1)),':','Color',blue,'LineWidth',1.5);
        labels{end+1} = 'Sea slope along x at board centre';
        ylabel(ax,'Sea slope (deg)');
        ax.YAxis(1).Color = red;
        ax.YAxis(2).Color = blue;
        if numel(t)>1, xlim(ax,[t(1) t(end)]); end
        xlabel(ax,'Elapsed time (s)');
        title(ax,wrd_fig_title(item),'FontWeight','normal');
        grid(ax,'on');
        if row==1 && col==1
            lgd = legend(ax,handles,labels,'Orientation','horizontal');
            lgd.Layout.Tile = 'south';
            lgd.Box = 'off';
            lgd.FontName = opts.fontName;
            lgd.FontSize = opts.fontSize;
        end
    end
end
title(layout,'Direct link: refraction and total internal reflection at the sea surface', ...
    'FontName',opts.fontName,'FontSize',opts.fontSize+1);
paths = wrd_fig_save(fig,'fig3_sea_tir',opts);
end

function value = share(link)
value = 100*[link.tirPowerW]./max(realmin,[link.incidentPowerW]);
value([link.incidentPowerW]<=0) = NaN;
end
