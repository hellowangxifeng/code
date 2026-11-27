function rectangle=rp_s2_plot_box(ax)
%RP_S2_PLOT_BOX Actual 2-D plot box inside axes with a fixed aspect ratio.
rectangle=getpixelposition(ax,true); ratio=[];
if strcmp(ax.PlotBoxAspectRatioMode,'manual')
    ratio=ax.PlotBoxAspectRatio(1)/ax.PlotBoxAspectRatio(2);
elseif strcmp(ax.DataAspectRatioMode,'manual')
    dx=diff(ax.XLim); dy=diff(ax.YLim);
    if strcmp(ax.XScale,'log'), dx=diff(log10(ax.XLim)); end
    if strcmp(ax.YScale,'log'), dy=diff(log10(ax.YLim)); end
    ratio=(dx/ax.DataAspectRatio(1))/(dy/ax.DataAspectRatio(2));
end
if isempty(ratio), return; end
if rectangle(3)/rectangle(4)>ratio
    width=rectangle(4)*ratio; rectangle(1)=rectangle(1)+(rectangle(3)-width)/2; rectangle(3)=width;
else
    height=rectangle(3)/ratio; rectangle(2)=rectangle(2)+(rectangle(4)-height)/2; rectangle(4)=height;
end
end
