function files=rp_fig12_plot(r)
%RP_FIG12_PLOT Independent 2.67 x 2.37 inch journal panels, with matched risks.
arguments
    r (1,1) struct
end
files={}; B=r.selection; colors=[0 .48 .20;.52 .12 .68;.83 .12 .12;.94 .52 0;.94 .52 0];
styles={'-','--',':','-.','--'}; names={'Robust ABER','Robust OP','Max mean','Flat (ABER)','Flat (OP)'};
fields={'actualRadiusCM','aber','op'}; stems={'Fig12a_radius','Fig12b_aber','Fig12c_op'};
labels={'Optimal radius (cm)','ABER','OP'}; auditRows=cell(7,10); figureId=0;
for panel=1:3
    [f,ax]=newPanel(); hold(ax,'on'); ids=1:5;
    if panel==2,ids=[1 2 3 4];elseif panel==3,ids=[1 2 3 5];end
    hh=gobjects(1,numel(ids));
    for k=1:numel(ids)
        id=ids(k); q=B(B.objectiveId==id,:);
        hh(k)=plot(ax,1000*q.sigmaM,q.(fields{panel}),'Color',colors(id,:), ...
            'LineStyle',styles{id},'LineWidth',1,'Marker','.','MarkerSize',4,'MarkerIndices',1:3:height(q));
    end
    xlim(ax,[-5 125]);xlabel(ax,'Position jitter \sigma_p (mm)');ylabel(ax,labels{panel});
    if panel==1,ylim(ax,[0 1.24*max(B.actualRadiusCM)]);else,riskAxis(ax,[2e-4 2]);end
    title(ax,sprintf('Fig. 12(%c)',char('a'+panel-1)),'FontSize',10,'FontWeight','normal');
    lg=placeLegend(ax,hh,names(ids));
    files{end+1}=savePanel(f,r,stems{panel}); %#ok<AGROW>
    figureId=figureId+1;auditRows(figureId,:)=inspectPanel(f,ax,lg,stems{panel});close(f);
end
schemeNames={'Robust','Max mean','Flat, fixed'}; colors3=colors([1 3 4],:); styles3={'-',':','--'};
for group=1:2
    T=r.frozen; xlabelText='Position jitter \sigma_p (mm)'; xfield='sigmaM'; multiplier=1000;
    limits=[-5 125]; low=2e-4;
    if group==2
        T=r.power;xlabelText='ON power (dBm)';xfield='powerDBm';multiplier=1;limits=[-8 48];
        low=10^floor(log10(min([T.aber;T.op]))-.3);
    end
    for metric=1:2
        field='aber';if metric==2,field='op';end
        [f,ax]=newPanel();hold(ax,'on');hh=gobjects(1,3);
        for scheme=1:3
            q=T(T.schemeId==scheme,:);
            hh(scheme)=plot(ax,multiplier*q.(xfield),q.(field),'Color',colors3(scheme,:), ...
                'LineStyle',styles3{scheme},'LineWidth',1.1);
        end
        xlim(ax,limits);riskAxis(ax,[low 2]);xlabel(ax,xlabelText);ylabel(ax,upper(field));
        title(ax,sprintf('Fig. %d(%c)',12+group,char('a'+metric-1)),'FontSize',10,'FontWeight','normal');
        lg=placeLegend(ax,hh,schemeNames);
        stem=sprintf('Fig%d%s_%s',12+group,char('a'+metric-1),field);
        files{end+1}=savePanel(f,r,stem); %#ok<AGROW>
        figureId=figureId+1;auditRows(figureId,:)=inspectPanel(f,ax,lg,stem);close(f);
    end
end
audit=cell2table(auditRows,'VariableNames',{'figure','widthIn','heightIn','font','tickPt', ...
    'legendPt','corner','overlapSamples','yMin','yMax'});
assert(all(abs(audit.widthIn-2.67)<1e-9 & abs(audit.heightIn-2.37)<1e-9));
assert(all(audit.tickPt==10 & audit.font=="Times New Roman" & audit.overlapSamples==0));
writetable(audit,fullfile(r.folder,'figure_audit.csv'));
save(fullfile(r.folder,'figure_files.mat'),'files','audit');
script=fullfile(fileparts(mfilename('fullpath')),'rp_fig12_finalize.py');
[status,message]=system(sprintf('python "%s" "%s"',script,r.folder));
assert(status==0,'rp:Fig12PDF','PDF frame verification failed: %s',message);
end

function [f,ax]=newPanel()
f=figure('Visible','off','Units','inches','Position',[1 1 2.67 2.37], ...
    'PaperUnits','inches','PaperSize',[2.67 2.37],'PaperPosition',[0 0 2.67 2.37], ...
    'Color','white','InvertHardcopy','off');
ax=axes(f,'Units','normalized','Position',[.235 .215 .725 .70], ...
    'FontName','Times New Roman','FontSize',10,'LineWidth',.65,'Box','on', ...
    'LabelFontSizeMultiplier',1,'TitleFontSizeMultiplier',1, ...
    'TickLabelInterpreter','tex','XGrid','on','YGrid','on','XMinorGrid','off','YMinorGrid','off');
end

function riskAxis(ax,limits)
ax.YScale='log';ax.YLim=limits;
ticks=10.^(ceil(log10(limits(1))):0);ax.YTick=ticks;
ax.YTickLabel=compose('10^{%d}',round(log10(ticks)));
end

function lg=placeLegend(ax,handles,names)
for fontSize=[7.5 7 6.5 6]
    lg=rp_s2_legend(ax,handles,names,1,fontSize,[10 6], ...
        struct('FontName','Times New Roman','Color','white','Box','off'));
    if lg.UserData.overlapSamples==0,break;end
end
assert(lg.UserData.overlapSamples==0,'rp:Fig12Legend','No unobstructed corner for this panel.');
end

function path=savePanel(f,r,stem)
f.UserData=struct('signature',r.signature,'apertureRadiusM',r.config.rx.apertureRadius, ...
    'ocean',r.config.study.turbulence,'atmosphere',r.config.study.atmosphere, ...
    'distanceM',r.distanceM);
drawnow;path=fullfile(r.folder,[stem '.png']);
savefig(f,fullfile(r.folder,[stem '.fig']));
print(f,path,'-dpng','-r600');
print(f,fullfile(r.folder,[stem '.pdf']),'-dpdf','-vector');
end

function row=inspectPanel(f,ax,lg,stem)
position=f.Position;
assert(ax.XLabel.FontSize==10 && ax.YLabel.FontSize==10 && ax.Title.FontSize==10, ...
    'rp:Fig12TextSize','Axis labels and titles must be 10 pt.');
row={string(stem),position(3),position(4),string(ax.FontName),ax.FontSize,lg.FontSize, ...
    string(lg.Location),lg.UserData.overlapSamples,ax.YLim(1),ax.YLim(2)};
end
