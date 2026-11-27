function h = rd_fig_curve(ax,P,y,k,o,winner,lineStyle)
%RD_FIG_CURVE One scheme curve; winners of the ranking are drawn thicker.
s = rd_fig_style();
if nargin<7 || isempty(lineStyle), lineStyle = s.lines{k}; end
width = o.lineWidth;
if winner, width = o.winnerLineWidth; end
y = double(y);
y(y<=0) = NaN;
h = plot(ax,P,y,'Color',s.colors(k,:),'LineStyle',lineStyle,'LineWidth',width);
end
