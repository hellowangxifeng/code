function files = rd_fig_main(R,o)
%RD_FIG_MAIN Fig2: ABER and OP of the five links at the main condition.
%   Lines: complex field for the RIS links, MC for the direct link.
%   Hollow markers: the SAME selected designs evaluated with their full MC
%   gains (validation of the complex-field formula). Thick lines: the
%   clearly best RIS optimizer(s) of rd_rank. Dotted: exact FOV floor.
style = rd_fig_style();
iz = rd_fig_index(R.distancesM,o.mainDistanceM);
iw = rd_fig_index(R.winds,o.mainWindMps);
q = R.cond(iz,iw);
P = R.powerDBm;
metrics = {'aber','op'};
labels = {'ABER','OP'};
floors = [q.fovOutageExact/2,q.fovOutageExact];
xLimits = rd_fig_xlimits(P,[q.curves{1}.aber;q.curves{1}.op],o.yLimits,o.xPaddingDB);
step = max(1,round(o.mcMarkerStepDB/(P(2)-P(1))));
marks = find(P>=xLimits(1) & P<=xLimits(2));
marks = marks(1:step:end);
[fig,layout] = rd_fig_new(1,2,o,3.3,2.9);
for k = 1:2
    ax = nexttile(layout);
    hold(ax,'on');
    rank = q.rank{1}.(metrics{k});
    names = R.schemeNames;
    names{4} = ['Robust ' labels{k}];
    h = gobjects(1,6);
    for s = 1:5
        winner = ismember(s,rank.winners) && rank.clear;
        h(s) = rd_fig_curve(ax,P,q.curves{1}.(metrics{k})(s,:),s,o,winner);
        if winner, names{s} = [names{s} ' (best)']; end
        if s>1
            y = double(q.mcCurves.(metrics{k})(s,marks));
            y(y<=0) = NaN;
            plot(ax,P(marks),y,style.markers{s},'Color',style.colors(s,:), ...
                'MarkerSize',4.5,'LineWidth',0.8,'MarkerFaceColor','none','HandleVisibility','off');
        end
    end
    h(6) = plot(ax,NaN,NaN,'ko','MarkerSize',4.5,'LineWidth',0.8);
    rd_fig_axes(ax,xLimits,o,floors(k),labels{k},true,true);
    legend(ax,h,[names,{'MC check (markers)'}],'Location','southwest','Box','off', ...
        'FontName',o.fontName,'FontSize',o.fontSize-1);
    rd_fig_tag(ax,sprintf('z = %g m, U = %g m/s',q.distanceM,q.wind),o,'northeast');
end
title(layout,sprintf('ocean EGG %s, atmosphere GG %s',R.turbulence{1,1},R.turbulence{1,2}), ...
    'FontName',o.fontName,'FontSize',o.fontSize);
files = rd_fig_save(fig,'fig2_aber_op_main',o);
end
