function file=rp_stage2_plot_link_metric(data,c,metric)
%RP_STAGE2_PLOT_LINK_METRIC Render one metric from the common Fig2/3 dataset.
assert(ismember(metric,{'aber','op'}),'rp:Metric','Expected aber or op.');
[colors,styles]=rp_s2_style(); number=2+strcmp(metric,'op');
stem=sprintf('Fig%d_%s',number,metric); sizeInches=[2.67 2.37];
powerDBm=data.powerDBm;
floorValue=data.([metric 'Floor']);
f=figure('Visible','off','Color','w','Units','inches', ...
    'Position',[.4 .4 sizeInches],'DefaultAxesFontName','Times New Roman', ...
    'DefaultTextFontName','Times New Roman','DefaultLegendFontName','Times New Roman');
ax=axes(f,'Units','normalized','Position',[.23 .21 .71 .67]); hold(ax,'on');
ax.TitleFontSizeMultiplier=1; ax.LabelFontSizeMultiplier=1;
handles=gobjects(1,6);
if isfield(data,'unmodulatedScatterOnly')
    lowerCurve=data.curves{1}.(metric); upperCurve=data.unmodulatedScatterOnly.(metric);
    fill(ax,[powerDBm fliplr(powerDBm)],[lowerCurve fliplr(upperCurve)], ...
        [.5 .5 .5],'FaceAlpha',.15,'EdgeColor','none','HandleVisibility','off');
end
for k=1:6
    m=data.curves{k};
    handles(k)=plot(ax,powerDBm,m.(metric),'Color',colors(k,:), ...
        'LineStyle',styles{k},'LineWidth',1.0);
end
xlabel(ax,'P_t (dBm)','FontSize',10); ylabel(ax,upper(metric),'FontSize',10);
title(ax,sprintf('Fig. %d | Fixed RX, %g mrad',number,1000*data.fovHalfAngleRad), ...
    'FontSize',10,'FontWeight','normal');
ax.YScale='log'; xlim(ax,[0 100]); ylim(ax,[floorValue/1.8,1.15]);
ax.XTick=0:20:100; ax.YTick=10.^(-4:0); ax.TickLength=[.025 .025];
grid(ax,'on'); ax.GridAlpha=.18; ax.XMinorGrid='off'; ax.YMinorGrid='off';
floorLine=yline(ax,floorValue,'k--','LineWidth',.8,'HandleVisibility','off');
floorLine.Tag='exact_fov_floor';
text(ax,98,1.55*floorValue,'Floor','VerticalAlignment','bottom','HorizontalAlignment','right', ...
    'FontSize',10,'FontName','Times New Roman','Clipping','on','Tag','fov_floor_label');
names={'U-RIS*','Steer','Flat','R-BER','R-OP','Max P'};
drawnow; ax.PlotBoxAspectRatioMode='manual';
rp_s2_legend(ax,handles,names,1,10,[8 6]);
set(findall(f,'-property','FontUnits'),'FontUnits','points');
set(findall(f,'-property','FontSize'),'FontSize',10);
set(findall(f,'-property','FontName'),'FontName','Times New Roman');
f.PaperUnits='inches'; f.PaperSize=sizeInches; f.PaperPosition=[0 0 sizeInches];
f.UserData=rmfield(data,{'curves','powerDBm'});
f.UserData.fovFloor=floorValue;
file=rp_s2_finish(f,c,stem);
% PRINT preserves the requested canvas rather than cropping to text extents.
f=openfig(fullfile(c.output,[stem '.fig']),'invisible'); f.PaperPositionMode='manual';
ax=findall(f,'Type','axes');rp_fixed_receiver_legend(f,ax);
savefig(f,fullfile(c.output,[stem '.fig']));
print(f,file,'-dpng',sprintf('-r%d',c.plot.resolution)); close(f);
end
