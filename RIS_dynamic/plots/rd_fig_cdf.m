function files = rd_fig_cdf(R,o)
%RD_FIG_CDF Fig5: CDF of the instantaneous received power (no fading).
%   Pooled over frames and receiver jitter/pointing samples; solid = rough
%   sea (o.cdfWindMps), dashed = calm sea. P_r = gain + o.cdfPowerDBm.
files = strings(0,1);
windy = find(R.winds>0);
calm = find(R.winds==0,1);
if isempty(windy), return; end
iw = windy(rd_fig_index(R.winds(windy),o.cdfWindMps));
style = rd_fig_style();
keys = R.schemes.keys;
x = R.main.cdfEdgesDB+o.cdfPowerDBm;
z = R.distancesM;
[fig,layout] = rd_fig_new(1,numel(z),o);
handles = gobjects(1,numel(keys));
for iz = 1:numel(z)
    ax = nexttile(layout); hold(ax,'on');
    lo = Inf; hi = -Inf;
    for k = 1:numel(keys)
        y = reshape(R.main.cdf(iw,iz,k,:),1,[]);
        if ~any(isfinite(y)), continue; end
        handles(k) = plot(ax,x,max(y,o.cdfFloor/10),'Color',style.colors(k,:), ...
            'LineStyle','-','LineWidth',style.widths(k));
        [lo,hi] = visibleRange(y,x,lo,hi,o.cdfFloor);
        if ~strcmp(keys{k},'sea_aware') && ~isempty(calm)
            y0 = reshape(R.main.cdf(calm,iz,k,:),1,[]);
            if any(isfinite(y0))
                plot(ax,x,max(y0,o.cdfFloor/10),'Color',style.colors(k,:),'LineStyle','--', ...
                    'LineWidth',0.9,'HandleVisibility','off');
                [lo,hi] = visibleRange(y0,x,lo,hi,o.cdfFloor);
            end
        end
    end
    if isfinite(lo) && isfinite(hi) && hi>lo, xlim(ax,[lo hi]); end
    set(ax,'YScale','log','YLim',[o.cdfFloor 1]); grid(ax,'on');
    xlabel(ax,'Received power P_r (dBm)');
    if iz==1, ylabel(ax,'CDF'); end
    text(ax,0.03,0.96,{sprintf('z = %g m, P_t = %g dBm',z(iz),o.cdfPowerDBm), ...
        sprintf('solid: U = %g m/s, dashed: calm',R.winds(iw))},'Units','normalized', ...
        'VerticalAlignment','top','FontName',o.fontName,'FontSize',o.fontSize-1, ...
        'BackgroundColor','w','Margin',1);
end
rd_fig_legend(layout,handles,R.schemes.names,o);
files = rd_fig_save(fig,'fig5_power_cdf',o);
end

function [lo,hi] = visibleRange(y,x,lo,hi,floorValue)
first = find(y>floorValue,1);
last = find(y<1-1e-12,1,'last');
if ~isempty(first), lo = min(lo,x(first)-2); end
if ~isempty(last), hi = max(hi,x(min(last+1,numel(x)))+2); end
end
