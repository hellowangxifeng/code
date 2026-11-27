function files = rd_fig_curves(R,o,metric)
%RD_FIG_CURVES Fig2 (ABER) / Fig3 (OP) versus transmit power.
%   Rows: receiver heights. Columns: sea states o.curveWinds (0 = calm sea).
%   Every curve is the time average over the frames of the sea state, with
%   jitter, pointing error, per-ray FOV and EGG x GG fading inside each frame.
%   The dotted line is the pointing-only FOV floor P(psi > Theta).
style = rd_fig_style();
keys = R.schemes.keys;
names = R.schemes.names;
P = R.powerDBm;
z = R.distancesM;
columns = unique(rd_fig_index(R.winds,o.curveWinds),'stable');
[fig,layout] = rd_fig_new(numel(z),numel(columns),o);
handles = gobjects(1,numel(keys));
label = upper(metric);
for iz = 1:numel(z)
    for jc = 1:numel(columns)
        iw = columns(jc);
        ax = nexttile(layout);
        hold(ax,'on');
        for k = 1:numel(keys)
            if R.winds(iw)==0 && strcmp(keys{k},'sea_aware'), continue; end
            y = reshape(R.main.(metric)(iw,iz,k,o.turbulenceId,:),1,[]);
            if ~any(isfinite(y)), continue; end
            y(y<=0) = NaN;
            handles(k) = plot(ax,P,y,'Color',style.colors(k,:),'LineStyle',style.lines{k}, ...
                'LineWidth',style.widths(k)*o.lineWidth/1.2);
        end
        floorValue = min(R.main.fovFloor(iz,:));
        if isfinite(floorValue) && floorValue>o.yLimits(1)
            yline(ax,floorValue,':k','LineWidth',0.8,'HandleVisibility','off');
        end
        set(ax,'YScale','log','XLim',[P(1) P(end)],'YLim',o.yLimits,'FontName',o.fontName, ...
            'FontSize',o.fontSize,'XMinorTick','on');
        grid(ax,'on');
        text(ax,0.97,0.96,sprintf('z = %g m, %s',z(iz),rd_fig_wind_label(R.winds(iw))), ...
            'Units','normalized','HorizontalAlignment','right','VerticalAlignment','top', ...
            'FontName',o.fontName,'FontSize',o.fontSize,'BackgroundColor','w','Margin',1);
        if iz==numel(z)
            xlabel(ax,'$P_t$ (dBm)','Interpreter','latex','FontSize',o.fontSize+1);
        end
        if jc==1
            ylabel(ax,label,'FontName',o.fontName,'FontSize',o.fontSize+1);
        end
    end
end
rd_fig_legend(layout,handles,names,o);
number = 2+strcmp(metric,'op');
files = rd_fig_save(fig,sprintf('fig%d_%s_vs_power',number,lower(metric)),o);
end
