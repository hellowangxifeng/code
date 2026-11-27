function rp_fixed_receiver_legend(f,ax)
%RP_FIXED_RECEIVER_LEGEND Compact 10-point swatches without shrinking glyphs.
% Native legend row padding cannot fit six entries on a 2.67-inch canvas.
lg=findall(f,'Type','legend');assert(isscalar(lg));
names=cellstr(string(lg.String));original=lg.UserData;handles=gobjects(1,numel(names));
for k=1:numel(names),handles(k)=findobj(ax,'Type','line','Tag',original.sourceTags{k});end
box=rp_s2_plot_box(ax);widthPt=44;heightPt=52.8;
sizePx=hgconvertunits(f,[0 0 widthPt heightPt],'points','pixels',f);
insetPx=hgconvertunits(f,[0 0 .75 8],'points','pixels',f);inset=insetPx(3:4);
corners={'northwest','northeast','southwest','southeast'};rects=zeros(4,4);scores=zeros(1,4);
xx=[];yy=[];
for k=1:numel(handles)
    x=handles(k).XData;y=handles(k).YData;t=linspace(1,numel(x),10*numel(x));
    xx=[xx,interp1(1:numel(x),x,t)];yy=[yy,interp1(1:numel(y),y,t)]; %#ok<AGROW>
end
floorLine=findobj(ax,'Tag','exact_fov_floor');
if ~isempty(floorLine),xx=[xx,linspace(ax.XLim(1),ax.XLim(2),1000)];yy=[yy,repmat(floorLine.Value,1,1000)];end
px=box(1)+(xx-ax.XLim(1))/diff(ax.XLim)*box(3);
py=box(2)+(log10(yy)-log10(ax.YLim(1)))/diff(log10(ax.YLim))*box(4);
for k=1:4
    if contains(corners{k},'west'),x=box(1)+inset(1);else,x=box(1)+box(3)-sizePx(3)-inset(1);end
    if contains(corners{k},'south'),y=box(2)+inset(2);else,y=box(2)+box(4)-sizePx(4)-inset(2);end
    rects(k,:)=[x,y,sizePx(3:4)];
    scores(k)=sum(px>x-2 & px<x+sizePx(3)+2 & py>y-2 & py<y+sizePx(4)+2);
end
[score,k]=min(scores);rect=rects(k,:);delete(lg);
rectPt=hgconvertunits(f,rect,'pixels','points',f);
for row=1:numel(names)
    centerY=rectPt(2)+heightPt-(row-.5)*8.8;
    lineRect=hgconvertunits(f,[rectPt(1),centerY,8,0],'points','normalized',f);
    annotation(f,'line',lineRect(1)+[0 lineRect(3)],repmat(lineRect(2),1,2), ...
        'Color',handles(row).Color,'LineStyle',handles(row).LineStyle,'LineWidth',1, ...
        'Tag',sprintf('fixed_legend_line_%02d',row));
    textRect=hgconvertunits(f,[rectPt(1)+11,centerY-4.4,33,8.8],'points','normalized',f);
    annotation(f,'textbox',textRect,'String',names{row},'FontName','Times New Roman','FontSize',10, ...
        'Interpreter','none','Margin',0,'EdgeColor','none','VerticalAlignment','middle', ...
        'Tag',sprintf('fixed_legend_text_%02d',row));
end
info=struct('selectedCorner',corners{k},'overlapSamples',score,'cornerScores',scores, ...
    'rectangle',hgconvertunits(f,rect,'pixels','normalized',f),'insetPoints',[.75 8], ...
    'entries',{names},'axisTag',ax.Tag,'fontPt',10,'sourceTags',{original.sourceTags});
hggroup(ax,'Tag','fixed_compact_legend','UserData',info);
f.UserData.compactLegend=info;drawnow;
assert(score==0,'rp:FixedCompactLegend','Compact 10 pt legend still overlaps curves.');
end
