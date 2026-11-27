function files=rp_stage2_plot(s,extra)
%RP_STAGE2_PLOT Stage-two figures only; all drawing code lives in plots/.
c=s.config; [colors,~]=rp_s2_style(); files={};
files{end+1}=rp_stage2_plot_irradiance(s);
files=[files,rp_stage2_plot_links(s.distances{3},c)];
% Fig4: z=7.5 m, 1 W; only the requested 0--18 cm horizontal extent.
d=s.distances{2}; fields={'meanReceivedPowerW','aber','op'}; labels={'Mean received power (mW)','ABER','OP'};
shape=d.summary.minorRadiusM>=c.study.minimumWidthM & d.summary.actualRadiusCM>=100*c.study.minimumWidthM & ...
    d.summary.axisRatio<=c.study.maximumAxisRatio & d.summary.r50M>=.02 & ...
    ismember(d.summary.family,string(c.study.families));
lowerWidth=max(100*c.study.minimumWidthM,floor(2*min(d.summary.actualRadiusCM(shape)))/2);
f=figure('Visible','off','Color','w','Position',[40 40 1060 350]); t=tiledlayout(f,1,3,'TileSpacing','compact','Padding','compact');
for row=1:3
    ax=nexttile(t); hold(ax,'on'); hh=gobjects(1,3);
    for k=1:3
        T=d.summary(d.summary.family==string(c.study.families{k}),:);
        T=T(T.minorRadiusM>=c.study.minimumWidthM & T.actualRadiusCM>=100*c.study.minimumWidthM & ...
            T.axisRatio<=c.study.maximumAxisRatio & T.r50M>=.02,:);
        [x,idx]=sort(T.actualRadiusCM); y=T.(fields{row})(idx); if row==1, y=1000*y; end
        xx=linspace(min(x),min(18,max(x)),160);
        if row>1, yy=exp(interp1(x,log(y),xx,'pchip')); else, yy=interp1(x,y,xx,'pchip'); end
        colorIds=[2 3 4]; hh(k)=plot(ax,xx,yy,'Color',colors(colorIds(k),:),'LineWidth',1.4);
        plot(ax,x,y,'.','Color',colors(colorIds(k),:),'HandleVisibility','off','MarkerSize',7);
    end
    xlabel(ax,'Actual equivalent radius (cm)'); ylabel(ax,labels{row});
    rp_s2_axes(ax,[lowerWidth 18],row>1);
    if row>1
        ticks=[.001 .002 .003 .005 .01 .02 .03 .05 .1 .2 .3 .5 1];
        ticks=ticks(ticks>=ax.YLim(1) & ticks<=ax.YLim(2));
        ax.YTick=ticks; ax.YTickLabel=compose('%.2g',ticks);
    end
    rp_s2_legend(ax,hh,{'Gaussian','Flat top','Robust family'});
end
title(t,'Fig 4 | z=7.5 m, ON power 1 W, jitter 80 mm; feasible widths, full scan in Fig 11/12','FontSize',10);
files{end+1}=rp_s2_finish(f,c,'Fig4_width_tradeoff');
files=[files,rp_stage2_plot_sensitivity(s,extra.sensitivity)];
% Fig6: only former seventh panel, i.e. z=10 m ABER objective.
d=s.distances{3}; b=d.bank{d.objectiveIndices(1)};
f=figure('Visible','off','Color','w','Position',[40 40 430 390]); ax=axes(f);
imagesc(ax,reshape(b.phase.phaseAtReferenceRad,c.array.count).'); axis(ax,'xy','image'); clim(ax,[0 2*pi]); colormap(ax,hsv);
cb=colorbar(ax); cb.Label.String='LC piston phase (rad)';
xlabel(ax,'Tile x index'); ylabel(ax,'Tile y index'); title(ax,'Fig 6 | z=10 m, Robust ABER, actual wrapped phase','FontSize',10);
files{end+1}=rp_s2_finish(f,c,'Fig6_objective_phase');
files{end+1}=rp_stage2_plot_radial(s);
files{end+1}=rp_stage2_plot_turbulence(s,extra.turbulence);
if isfield(extra,'angles') && numel(extra.angles)==2 && all(~cellfun(@isempty,extra.angles))
    files{end+1}=rp_stage2_plot_angles(s,extra.angles);
end
[~,adaptiveFiles]=rp_fig12_run(s); files=[files,adaptiveFiles];
save(fullfile(c.output,'figure_files.mat'),'files');
end
