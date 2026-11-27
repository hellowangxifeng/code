function files = rd_fig_wind(R,o)
%RD_FIG_WIND Fig6: ABER and OP at o.fixedPowerDBm versus wind speed.
%   U = 0 is the calm sea (motionless board, flat water). At U = 0 the
%   sea-aware link coincides with Robust ABER by construction.
files = strings(0,1);
if numel(R.winds)<2, return; end
style = rd_fig_style();
keys = R.schemes.keys;
P = R.powerDBm;
z = R.distancesM;
metrics = {'aber','op'};
[fig,layout] = rd_fig_new(2,numel(z),o);
handles = gobjects(1,numel(keys));
for im = 1:2
    for iz = 1:numel(z)
        ax = nexttile(layout); hold(ax,'on');
        low = Inf;
        for k = 1:numel(keys)
            y = arrayfun(@(iw) rd_fig_at(R.main.(metrics{im})(iw,iz,k,o.turbulenceId,:), ...
                P,o.fixedPowerDBm),1:numel(R.winds));
            y(y<=0) = NaN;
            if ~any(isfinite(y)), continue; end
            handles(k) = plot(ax,R.winds,y,'Color',style.colors(k,:),'LineStyle',style.lines{k}, ...
                'LineWidth',style.widths(k),'Marker',style.markers{k},'MarkerSize',4, ...
                'MarkerFaceColor','w');
            low = min(low,min(y(isfinite(y))));
        end
        set(ax,'YScale','log','YLim',[bottomOf(low) 1],'XLim',[min(R.winds) max(R.winds)]);
        grid(ax,'on');
        if im==2, xlabel(ax,'Wind speed U (m/s)'); end
        if iz==1, ylabel(ax,upper(metrics{im})); end
        text(ax,0.03,0.96,sprintf('z = %g m, P_t = %g dBm',z(iz),o.fixedPowerDBm), ...
            'Units','normalized','VerticalAlignment','top','FontName',o.fontName, ...
            'FontSize',o.fontSize,'BackgroundColor','w','Margin',1);
    end
end
rd_fig_legend(layout,handles,R.schemes.names,o);
files = rd_fig_save(fig,'fig6_wind_speed',o);
end

function b = bottomOf(low)
b = 10^floor(log10(max(low,1e-12)));
if ~isfinite(b) || b>=1, b = 1e-5; end
end
