function lg=rp_s2_legend(ax,handles,names,columns,fontSize,tokenSize,appearance)
%RP_S2_LEGEND Place an actual-sized legend in the least occupied data corner.
if nargin<4, columns=1; end
if nargin<5, fontSize=8; end
if nargin<6, tokenSize=[15 8]; end
if nargin<7, appearance=struct(); end
locations={'northwest','northeast','southwest','southeast'};
lg=legend(ax,handles,names,'Location',locations{1},'FontSize',fontSize, ...
    'NumColumns',columns,'Interpreter','none','Box','off');
lg.ItemTokenSize=tokenSize;
% Font/box setters can trigger a native normalized-size layout. Apply them
% before measuring and converting the final pixel rectangle.
if ~isempty(fieldnames(appearance)), set(lg,appearance); end
if isempty(ax.Tag)
    ax.Tag=sprintf('s2_legend_axis_%d',numel(findall(ancestor(ax,'figure'),'Type','legend')));
end
sourceTags=cell(1,numel(handles));
for j=1:numel(handles)
    if isempty(handles(j).Tag), handles(j).Tag=sprintf('%s_item_%d',ax.Tag,j); end
    sourceTags{j}=handles(j).Tag;
end
xx=[]; yy=[]; objects=findobj(ax,'Type','line');
for j=1:numel(objects)
    x=objects(j).XData; y=objects(j).YData;
    if numel(x)>1 && numel(x)==numel(y)
        parameter=1:numel(x); fine=linspace(1,numel(x),max(200,10*numel(x)));
        x=interp1(parameter,x,fine,'linear'); y=interp1(parameter,y,fine,'linear');
    end
    xx=[xx,x]; yy=[yy,y]; %#ok<AGROW>
end
x=ax.XLim; y=ax.YLim; xx=(xx-x(1))/diff(x);
if strcmp(ax.YScale,'log')
    yy=(log10(max(yy,realmin))-log10(y(1)))/diff(log10(y));
else
    yy=(yy-y(1))/diff(y);
end
scores=zeros(1,4);
tight=strcmp(ax.DataAspectRatioMode,'manual') || strcmp(ax.PlotBoxAspectRatioMode,'manual');
for corner=1:4
    placeCorner(ax,lg,locations{corner},tight);
    a=rp_s2_plot_box(ax); bounds=getpixelposition(lg,true);
    px=a(1)+xx*a(3); py=a(2)+yy*a(4);
    scores(corner)=sum(px>bounds(1)-2 & px<sum(bounds([1 3]))+2 & ...
        py>bounds(2)-2 & py<sum(bounds([2 4]))+2);
end
[score,best]=min(scores); placeCorner(ax,lg,locations{best},tight);
lg.UserData=struct('cornerScores',scores,'selectedCorner',locations{best}, ...
    'overlapSamples',score,'criterion','dense line segments vs actual legend pixel bounds', ...
    'axisTag',ax.Tag,'sourceTags',{sourceTags},'tightCorner',tight);
end

function placeCorner(ax,lg,location,tight)
% MATLAB's automatic legend layout stores normalized sizes. Reset units
% before each native placement, then convert the complete rectangle once.
lg.Units='normalized'; lg.Location=location; drawnow;
if ~tight, return; end
fig=ancestor(ax,'figure');
box=rp_s2_plot_box(ax); current=hgconvertunits(fig,lg.Position,'normalized','pixels',fig); desired=current(1:2);
if contains(location,'west'), desired(1)=box(1)+1; else, desired(1)=box(1)+box(3)-current(3)-1; end
if contains(location,'south'), desired(2)=box(2)+1; else, desired(2)=box(2)+box(4)-current(4)-1; end
% Tiled layouts recompute normalized legend sizes during save/open. Keep the
% complete manual rectangle normalized so that native updates cannot mix units.
lg.Units='normalized'; lg.Position=hgconvertunits(fig,[desired current(3:4)],'pixels','normalized',fig);
drawnow;
end
