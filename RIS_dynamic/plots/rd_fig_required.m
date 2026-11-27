function files = rd_fig_required(R,o)
%RD_FIG_REQUIRED Fig7: transmit power needed to reach the ABER / OP targets.
%   Rows: ABER target, OP target (c.rank.targets). Columns: wind speeds.
%   x: RX height; one line per link, main fading row. Links that never
%   reach the target inside the computed power range are left out (Direct
%   is usually among them). The vertical gap between lines is the RIS gain.
style = rd_fig_style();
metrics = {'aber','op'};
targets = R.config.rank.targets;
nw = numel(R.winds);
z = R.distancesM;
[fig,layout] = rd_fig_new(2,nw,o,o.gridTileWidthIn,o.gridTileHeightIn);
h = gobjects(1,5);
for k = 1:2
    for iw = 1:nw
        ax = nexttile(layout);
        hold(ax,'on');
        need = nan(5,numel(z));
        for iz = 1:numel(z)
            need(:,iz) = R.cond(iz,iw).rank{1}.(metrics{k}).requiredDBm(:);
        end
        need(~isfinite(need)) = NaN;
        for s = 1:5
            h(s) = plot(ax,z,need(s,:),'Color',style.colors(s,:),'LineStyle',style.lines{s}, ...
                'Marker',style.markers{s},'MarkerSize',4,'LineWidth',o.lineWidth);
        end
        set(ax,'FontName',o.fontName,'FontSize',o.fontSize,'Box','on','XLim',[z(1) z(end)], ...
            'XTick',z);
        grid(ax,'on');
        if k==2, xlabel(ax,'RX height z (m)'); end
        if iw==1
            ylabel(ax,sprintf('P_t at %s = 10^{%g} (dBm)',upper(metrics{k}),log10(targets(k))));
        end
        rd_fig_tag(ax,sprintf('U = %g m/s',R.winds(iw)),o,'northwest');
    end
end
rd_fig_legend(layout,h,R.schemeNames,o,5);
files = rd_fig_save(fig,'fig7_required_power',o);
end
