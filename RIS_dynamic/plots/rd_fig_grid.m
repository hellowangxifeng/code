function files = rd_fig_grid(R,o,metric)
%RD_FIG_GRID Fig4 (ABER) / Fig5 (OP): every RX height x every sea state.
%   Rows: wind speeds; columns: RX heights above the calm RIS centre. Main
%   fading row; each panel has its own frozen sea instant and its own
%   common receiver inside the RIS-reachable cone. Thick: clear winner(s).
P = R.powerDBm;
label = upper(metric);
nz = numel(R.distancesM);
nw = numel(R.winds);
[fig,layout] = rd_fig_new(nw,nz,o,o.gridTileWidthIn,o.gridTileHeightIn);
h = gobjects(1,5);
for iw = 1:nw
    for iz = 1:nz
        q = R.cond(iz,iw);
        ax = nexttile(layout);
        hold(ax,'on');
        Y = q.curves{1}.(metric);
        rank = q.rank{1}.(metric);
        for s = 1:5
            h(s) = rd_fig_curve(ax,P,Y(s,:),s,o,ismember(s,rank.winners) && rank.clear);
        end
        floorValue = q.fovOutageExact;
        if strcmp(metric,'aber'), floorValue = floorValue/2; end
        rd_fig_axes(ax,rd_fig_xlimits(P,Y,o.yLimits,o.xPaddingDB),o,floorValue,label,iw==nw,iz==1);
        rd_fig_tag(ax,sprintf('z = %g m, U = %g m/s',q.distanceM,q.wind),o,'northeast');
    end
end
names = R.schemeNames;
names{4} = ['Robust ' label];
rd_fig_legend(layout,h,names,o,5);
number = 4+strcmp(metric,'op');
files = rd_fig_save(fig,sprintf('fig%d_%s_distance_wind',number,lower(metric)),o);
end
