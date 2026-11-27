function files = rd_fig_incidence(R,o)
%RD_FIG_INCIDENCE Fig7: ABER and OP at o.fixedPowerDBm versus incidence angle.
%   Solid: rough sea U = R.incidence.windMps; dashed: calm sea. The vertical
%   dotted line is the calm water-air critical angle asin(n_a/n_w). Beyond it
%   the calm direct link is cut off by total internal reflection, while the
%   rough surface still lets part of the beam out.
I = R.incidence;
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
            for cond = [2 1]
                y = arrayfun(@(ia) rd_fig_at(I.(metrics{im})(ia,cond,iz,k,:),P,o.fixedPowerDBm), ...
                    1:numel(I.anglesDeg));
                y(y<=0) = NaN;
                if ~any(isfinite(y)), continue; end
                if cond==2
                    handles(k) = plot(ax,I.anglesDeg,y,'Color',style.colors(k,:),'LineStyle','-', ...
                        'LineWidth',style.widths(k),'Marker',style.markers{k},'MarkerSize',4, ...
                        'MarkerFaceColor','w');
                else
                    plot(ax,I.anglesDeg,y,'Color',style.colors(k,:),'LineStyle','--', ...
                        'LineWidth',0.9,'HandleVisibility','off');
                end
                low = min(low,min(y(isfinite(y))));
            end
        end
        xline(ax,I.criticalAngleDeg,':k','LineWidth',0.9,'HandleVisibility','off');
        set(ax,'YScale','log','YLim',[bottomOf(low) 1], ...
            'XLim',[min(I.anglesDeg) max(I.anglesDeg)+eps]);
        grid(ax,'on');
        if im==2, xlabel(ax,'Incidence angle in water (deg)'); end
        if iz==1, ylabel(ax,upper(metrics{im})); end
        text(ax,0.03,0.96,{sprintf('z = %g m, P_t = %g dBm',z(iz),o.fixedPowerDBm), ...
            sprintf('solid: U = %g m/s, dashed: calm',I.windMps)},'Units','normalized', ...
            'VerticalAlignment','top','FontName',o.fontName,'FontSize',o.fontSize-1, ...
            'BackgroundColor','w','Margin',1);
    end
end
rd_fig_legend(layout,handles,R.schemes.names,o);
files = rd_fig_save(fig,'fig7_incidence_angle',o);
end

function b = bottomOf(low)
b = 10^floor(log10(max(low,1e-12)));
if ~isfinite(b) || b>=1, b = 1e-5; end
end
