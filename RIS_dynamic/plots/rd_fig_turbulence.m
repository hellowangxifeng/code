function files = rd_fig_turbulence(R,o)
%RD_FIG_TURBULENCE Fig3: the five links under weak/strong ocean EGG and atmospheric GG.
%   Main condition. Rows: ABER, OP. Columns: ocean EGG strength. In every
%   panel, solid = weak atmospheric GG, dashed = strong GG. Designs are
%   reselected for every fading condition (rd_select). Thick: clear winner(s).
iz = rd_fig_index(R.distancesM,o.mainDistanceM);
iw = rd_fig_index(R.winds,o.mainWindMps);
q = R.cond(iz,iw);
T = R.turbulence;
P = R.powerDBm;
oceans = unique(T(:,1),'stable');
styles = {'-','--',':','-.'};
metrics = {'aber','op'};
labels = {'ABER','OP'};
floors = [q.fovOutageExact/2,q.fovOutageExact];
stack = zeros(0,numel(P));
for t = 1:size(T,1), stack = [stack;q.curves{t}.aber;q.curves{t}.op]; end %#ok<AGROW>
xLimits = rd_fig_xlimits(P,stack,o.yLimits,o.xPaddingDB);
[fig,layout] = rd_fig_new(2,numel(oceans),o);
h = gobjects(1,5);
for k = 1:2
    for jo = 1:numel(oceans)
        ax = nexttile(layout);
        hold(ax,'on');
        rows = find(strcmp(T(:,1),oceans{jo}));
        for jr = 1:numel(rows)
            t = rows(jr);
            rank = q.rank{t}.(metrics{k});
            for s = 1:5
                winner = ismember(s,rank.winners) && rank.clear;
                pattern = styles{min(jr,numel(styles))};
                if s==1, pattern = ':'; if jr>1, pattern = '-.'; end, end
                hh = rd_fig_curve(ax,P,q.curves{t}.(metrics{k})(s,:),s,o,winner,pattern);
                if jr==1, h(s) = hh; end
            end
        end
        rd_fig_axes(ax,xLimits,o,floors(k),labels{k},k==2,jo==1);
        atm = strjoin(T(rows,2).',' / ');
        rd_fig_tag(ax,sprintf('EGG %s; GG %s',oceans{jo},atm),o,'northeast');
    end
end
atmospheres = T(strcmp(T(:,1),oceans{1}),2).';
key = gobjects(1,numel(atmospheres));
for j = 1:numel(atmospheres)
    key(j) = plot(ax,NaN,NaN,['k' styles{min(j,numel(styles))}],'LineWidth',o.lineWidth);
end
rd_fig_legend(layout,[h,key],[R.schemeNames,strcat({'GG '},atmospheres)],o,4);
title(layout,sprintf('z = %g m, U = %g m/s (Direct: dotted / dash-dot)',q.distanceM,q.wind), ...
    'FontName',o.fontName,'FontSize',o.fontSize);
files = rd_fig_save(fig,'fig3_turbulence',o);
end
