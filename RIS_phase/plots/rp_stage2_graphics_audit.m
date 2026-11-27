function audit=rp_stage2_graphics_audit(c)
%RP_STAGE2_GRAPHICS_AUDIT Inspect saved FIG objects, rather than source promises.
stems={'Fig1_irradiance','Fig2_aber','Fig3_op','Fig4_width_tradeoff','Fig5_sensitivity', ...
    'Fig6_objective_phase','Fig8_radial_profiles','Fig9_turbulence','Fig10_incidence', ...
    'fig12_revision/Fig12a_radius','fig12_revision/Fig12b_aber','fig12_revision/Fig12c_op', ...
    'fig12_revision/Fig13a_aber','fig12_revision/Fig13b_op', ...
    'fig12_revision/Fig14a_aber','fig12_revision/Fig14b_op'};
expectedAxes=[15 1 1 3 9 1 6 4 4 1 1 1 1 1 1 1]; axesRows={}; legendRows={};
allowed={'northwest','northeast','southwest','southeast'};
assert(~isfile(fullfile(c.output,'Fig7_distance_width.fig')),'rp:Fig7','Removed Fig7 is still present.');
for j=1:numel(stems)
    file=fullfile(c.output,[stems{j},'.fig']); assert(isfile(file),'rp:FigureMissing','Missing %s.',file);
    f=openfig(file,'invisible'); clean=onCleanup(@() close(f)); drawnow;
    axesHandles=findall(f,'Type','axes');
    assert(numel(axesHandles)==expectedAxes(j),'rp:PanelCount','Unexpected panel count in %s.',stems{j});
    for k=1:numel(axesHandles)
        ax=axesHandles(k); limits=[ax.XLim ax.YLim];
        assert(abs(ax.FontSize-10)<1e-9 && strcmp(ax.FontUnits,'points'),'rp:TickFont','Ticks are not 10 pt.');
        if j==4
            assert(abs(ax.XLim(2)-18)<1e-9,'rp:Fig4Range','Fig4 must end at 18 cm.');
            if ismember(ax.YLabel.String,{'ABER','OP'}), assert(strcmp(ax.YScale,'log')); end
        elseif j==5
            labels={'Position jitter sigma (mm)','FOV half-angle (mrad)','Angular error sigma (mrad)'};
            expected=[20 120;30 150;10 60]; id=find(strcmp(labels,ax.XLabel.String),1);
            assert(~isempty(id) && max(abs(ax.XLim-expected(id,:)))<1e-9,'rp:Fig5Range','Incorrect sensitivity range.');
        elseif j==7
            assert(max(abs(ax.XLim-[0 40]))<1e-9,'rp:Fig8Range','Fig8 must show 0--40 cm.');
        end
        axesRows(end+1,:)={stems{j},k,ax.FontSize,string(ax.XScale),string(ax.YScale),limits(1),limits(2),limits(3),limits(4)}; %#ok<AGROW>
    end
    legends=findall(f,'Type','legend');
    for k=1:numel(legends)
        lg=legends(k); corner=lg.Location;
        if strcmp(corner,'none') && isstruct(lg.UserData) && isfield(lg.UserData,'tightCorner') && lg.UserData.tightCorner
            corner=lg.UserData.selectedCorner; ax=findobj(f,'Type','axes','Tag',lg.UserData.axisTag);
            box=rp_s2_plot_box(ax); bounds=getpixelposition(lg,true);
            assert(bounds(1)>=box(1)-1 && bounds(2)>=box(2)-1 && ...
                bounds(1)+bounds(3)<=box(1)+box(3)+1 && bounds(2)+bounds(4)<=box(2)+box(4)+1, ...
                'rp:LegendBounds','Manual image legend is outside the plot box.');
            if contains(corner,'east'), dx=box(1)+box(3)-bounds(1)-bounds(3); else, dx=bounds(1)-box(1); end
            if contains(corner,'north'), dy=box(2)+box(4)-bounds(2)-bounds(4); else, dy=bounds(2)-box(2); end
            assert(abs(dx-1)<2 && abs(dy-1)<2,'rp:LegendCorner', ...
                '%s legend %d corner gaps are %.4g / %.4g pixels.',stems{j},k,dx,dy);
        end
        assert(ismember(corner,allowed),'rp:LegendCorner','Legend is not in a coordinate corner.');
        names=string(lg.String);
        assert(~any(contains(lower(names),'circular') | contains(names,'Original RIS')),'rp:LegacyLabel','Legacy scheme name in a current figure.');
        overlaps=NaN; if isstruct(lg.UserData) && isfield(lg.UserData,'overlapSamples'), overlaps=lg.UserData.overlapSamples; end
        legendRows(end+1,:)={stems{j},k,string(corner),overlaps,strjoin(names,' | ')}; %#ok<AGROW>
    end
    compact=findall(f,'Tag','fixed_compact_legend');
    for k=1:numel(compact)
        info=compact(k).UserData;ax=findobj(f,'Type','axes','Tag',info.axisTag);
        box=rp_s2_plot_box(ax);bounds=hgconvertunits(f,info.rectangle,'normalized','pixels',f);
        inset=hgconvertunits(f,[0 0 info.insetPoints],'points','pixels',f);
        if contains(info.selectedCorner,'east'),dx=sum(box([1 3]))-sum(bounds([1 3]));else,dx=bounds(1)-box(1);end
        if contains(info.selectedCorner,'north'),dy=sum(box([2 4]))-sum(bounds([2 4]));else,dy=bounds(2)-box(2);end
        assert(abs(dx-inset(3))<2 && abs(dy-inset(4))<2 && ismember(info.selectedCorner,allowed), ...
            'rp:CompactCorner','Compact legend is not at its declared coordinate corner.');
        for row=1:numel(info.entries)
            txt=findall(f,'Tag',sprintf('fixed_legend_text_%02d',row));
            assert(isscalar(txt) && txt.FontSize==10 && strcmp(txt.FontName,'Times New Roman') && ...
                strcmp(txt.String,info.entries{row}),'rp:CompactFont','Compact legend label/font changed.');
        end
        legendRows(end+1,:)={stems{j},k,string(info.selectedCorner),info.overlapSamples,strjoin(string(info.entries),' | ')}; %#ok<AGROW>
    end
    clear clean
end
audit=struct('axes',cell2table(axesRows,'VariableNames',{'figure','axesId','tickPt','xScale','yScale','xMin','xMax','yMin','yMax'}), ...
    'legends',cell2table(legendRows,'VariableNames',{'figure','legendId','corner','overlapSamples','entries'}));
writetable(audit.axes,fullfile(c.root,'audit','stage2','figure_axes.csv'));
writetable(audit.legends,fullfile(c.root,'audit','stage2','figure_legends.csv'));
save(fullfile(c.root,'audit','stage2','graphics_audit.mat'),'audit');
assert(~any(audit.legends.overlapSamples>0),'rp:LegendOverlap', ...
    'At least one saved curve legend overlaps plotted line segments; inspect figure_legends.csv.');
end
