function files = rd_fig_width(R,o)
%RD_FIG_WIDTH Fig8: sea-aware training, time-averaged ABER / OP versus spot size.
%   Each curve is one design family in one wind, scored on the INDEPENDENT
%   training sea at the design power; the x axis is the calm-sea equivalent
%   spot radius. Filled markers: the selected sea-aware design.
T = R.trainingTable;
winds = unique(T.windMps).';
z = R.distancesM;
colors = parula(numel(winds)+1);
families = ["robust","flat_top"];
lineStyles = {'-','--'};
metrics = {'aber','op'};
[fig,layout] = rd_fig_new(2,numel(z),o);
handles = gobjects(1,numel(winds)+2);
labels = [arrayfun(@(U) sprintf('U = %g m/s',U),winds,'UniformOutput',false), ...
    {'robust family','flat-top family'}];
for im = 1:2
    for iz = 1:numel(z)
        ax = nexttile(layout); hold(ax,'on');
        low = Inf;
        for iw = 1:numel(winds)
            for jf = 1:2
                rows = T(T.windMps==winds(iw) & abs(T.distanceM-z(iz))<1e-9 & ...
                    T.family==families(jf),:);
                if height(rows)==0, continue; end
                rows = sortrows(rows,'equivalentRadiusM');
                y = rows.(metrics{im});
                y(y<=0) = NaN;
                h = plot(ax,100*rows.equivalentRadiusM,y,'Color',colors(iw,:), ...
                    'LineStyle',lineStyles{jf},'LineWidth',1.1,'Marker','o','MarkerSize',3.5, ...
                    'MarkerFaceColor','w');
                if jf==1 || ~isgraphics(handles(iw)), handles(iw) = h; end
                chosen = rows.chosen;
                if any(chosen)
                    plot(ax,100*rows.equivalentRadiusM(chosen),y(chosen),'o','MarkerSize',6, ...
                        'MarkerFaceColor',colors(iw,:),'MarkerEdgeColor','k','HandleVisibility','off');
                end
                if any(isfinite(y)), low = min(low,min(y(isfinite(y)))); end
            end
        end
        b = 10^floor(log10(max(low,1e-12)));
        if ~isfinite(b) || b>=1, b = 1e-5; end
        set(ax,'YScale','log','YLim',[b 1]);
        grid(ax,'on');
        if im==2, xlabel(ax,'Calm-sea equivalent spot radius (cm)'); end
        if iz==1, ylabel(ax,['Training ' upper(metrics{im})]); end
        text(ax,0.97,0.96,sprintf('z = %g m',z(iz)),'Units','normalized', ...
            'HorizontalAlignment','right','VerticalAlignment','top','FontName',o.fontName, ...
            'FontSize',o.fontSize,'BackgroundColor','w','Margin',1);
    end
end
handles(end-1) = plot(ax,NaN,NaN,'k-','LineWidth',1.1);
handles(end) = plot(ax,NaN,NaN,'k--','LineWidth',1.1);
rd_fig_legend(layout,handles,labels,o,min(4,numel(labels)));
files = rd_fig_save(fig,'fig8_sea_aware_width',o);
end
