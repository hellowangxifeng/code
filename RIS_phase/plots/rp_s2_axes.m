function rp_s2_axes(ax,xRange,isLog)
%RP_S2_AXES Tight data-driven limits with room for line caps and floor lines.
grid(ax,'on'); ax.GridAlpha=.18;
ax.XMinorGrid='off'; ax.YMinorGrid='off';
if isLog, ax.YScale='log'; end
xlim(ax,xRange); data=[];
lines=findobj(ax,'Type','line');
for k=1:numel(lines)
    x=lines(k).XData; y=lines(k).YData;
    if numel(x)==numel(y)
        use=x>=xRange(1) & x<=xRange(2) & isfinite(y);
        if isLog, use=use & y>0; end
        data=[data,y(use)]; %#ok<AGROW>
    end
end
if isempty(data), return; end
lo=min(data); hi=max(data);
if isLog
    span=max(.1,log10(hi)-log10(lo)); ylim(ax,[10^(log10(lo)-.045*span),10^(log10(hi)+.045*span)]);
    limits=ax.YLim;
    if strcmp(ax.YLabel.String,'OP'), limits(2)=min(limits(2),1.02); end
    if strcmp(ax.YLabel.String,'ABER'), limits(2)=min(limits(2),.51); end
    ylim(ax,limits);
else
    span=max(hi-lo,max(abs(hi),1e-9)*.05);
    ylim(ax,[max(0,lo-.04*span),hi+.06*span]);
end
ax.FontSize=10;
end
