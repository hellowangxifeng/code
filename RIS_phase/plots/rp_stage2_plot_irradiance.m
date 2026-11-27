function file=rp_stage2_plot_irradiance(s)
%RP_STAGE2_PLOT_IRRADIANCE Bare-surface MC reference and true waist ellipses.
arguments
    s (1,1) struct
end
c=s.config; panels=cell(3,4); direct=cell(1,3); rows=zeros(12,9); peak=0;
names={'No RIS (MC)','Steering-only','Flat top','Robust ABER'};
for j=1:3
    d=s.distances{j}; direct{j}=rp_mc_water_air(c,d.geometry,d.mc);
    map=rp_fig1_mc_map(direct{j}); panels{j,1}=map;
    for k=2:4, panels{j,k}=d.bank{d.comparisonIndices(k)}; end
    for k=1:4
        b=panels{j,k}; peak=max(peak,b.spot.peakWm2);
        rows((j-1)*4+k,:)=[d.distanceM,k,b.spot.powerW,b.spot.peakWm2, ...
            b.spot.majorRadiusM,b.spot.minorRadiusM,b.spot.axisRatio,b.spot.centerXY];
    end
end
summary=array2table(rows,'VariableNames',{'distanceM','schemeId','mapPowerW','peakWm2', ...
    'majorRadiusM','minorRadiusM','axisRatio','centroidXM','centroidYM'});
summary.scheme=string(names(summary.schemeId)).';
save(fullfile(c.output,'Fig1_data.mat'),'panels','direct','summary','names','-v7.3');
writetable(summary,fullfile(c.output,'Fig1_summary.csv'));
f=figure('Visible','off','Color','w','Position',[40 40 1260 890]);
t=tiledlayout(f,3,4,'TileSpacing','compact','Padding','compact'); axesMain=gobjects(3,4);
floorI=peak*1e-4;
for j=1:3
    for k=1:4
        b=panels{j,k}; w=b.wave; e=rp_spot_ellipse(b.spot);
        ax=nexttile(t);axesMain(j,k)=ax;hold(ax,'on');
        imagesc(ax,100*(w.x-b.spot.centerXY(1)),100*(w.y-b.spot.centerXY(2)),max(w.irradianceWm2,floorI));
        axis(ax,'xy','equal');clim(ax,[floorI peak]);ax.ColorScale='log';colormap(ax,turbo);
        % Show every pixel edge, plus the complete rotated ellipse and axes.
        dx=diff(w.x(1:2));dy=diff(w.y(1:2));
        lo=[w.x(1)-dx/2,w.y(1)-dy/2]-b.spot.centerXY;
        hi=[w.x(end)+dx/2,w.y(end)+dy/2]-b.spot.centerXY;
        lo=min(lo,min(e.offsetXY,[],1));hi=max(hi,max(e.offsetXY,[],1));
        margin=.03*(hi-lo);xlim(ax,100*[lo(1)-margin(1),hi(1)+margin(1)]);
        ylim(ax,100*[lo(2)-margin(2),hi(2)+margin(2)]);
        h=plot(ax,100*e.offsetXY(:,1),100*e.offsetXY(:,2),'k--','LineWidth',.9,'Tag','waist_ellipse');
        if k<=2
            for dimension=1:2
                endpoint=100*e.axisEndsXY(:,dimension);
                plot(ax,[-endpoint(1) endpoint(1)],[-endpoint(2) endpoint(2)], ...
                    'k:','LineWidth',.7,'Tag',sprintf('waist_axis_%d',dimension));
            end
        end
        lg=rp_s2_legend(ax,h,{'RMS ellipse'},1,7,[12 6]);lg.Color='white';lg.Box='on';
        cb=colorbar(ax);cb.Label.String='W m^{-2}';
        title(ax,sprintf('z=%g m | %s\nw_{maj}/w_{min}=%.1f/%.1f cm', ...
            s.distances{j}.distanceM,names{k},100*e.semiAxesM),'FontSize',9);
        xlabel(ax,'x - centroid (cm)');ylabel(ax,'y - centroid (cm)');
        ax.UserData=struct('row',j,'column',k,'worldCentroidXY',b.spot.centerXY, ...
            'displayedMapFraction',1,'waistCovarianceM2',b.spot.covarianceM2);
    end
end
title(t,{'Fig 1 | Same 1 W launch and receiver planes; dashed waist ellipses, dotted principal axes', ...
    'No RIS: bare water-air MC; context insets use asinh(Delta x / m), asinh(Delta y / m) and include all transmitted rays'},'FontSize',10);
drawnow;
% Context insets include the complete MC landing distribution. Invertible
% asinh coordinates compress remote grazing rays without clipping packets.
for j=1:3
    parent=axesMain(j,1);box=rp_s2_plot_box(parent);pos=f.Position;
    halo=panels{j,1}.halo;
    inset=axes(f,'Units','normalized','Position', ...
        [(box(1)+.015*box(3))/pos(3),(box(2)+.015*box(4))/pos(4),.33*box(3)/pos(3),.23*box(4)/pos(4)], ...
        'Tag',sprintf('Fig1_MC_full_halo_%d',j));
    imagesc(inset,halo.u,halo.v,max(halo.irradianceWm2,floorI));axis(inset,'xy');
    inset.Tag=sprintf('Fig1_MC_full_halo_%d',j);
    inset.ColorScale='log';clim(inset,[floorI peak]);inset.XTick=[];inset.YTick=[];
    text(inset,.5,.96,'Full MC','Units','normalized','Color','white', ...
        'HorizontalAlignment','center','VerticalAlignment','top','FontSize',7);
    inset.UserData=struct('row',j,'capturedPowerFraction',halo.capturedPowerFraction);
end
% Export to a fresh file first; keep the previous PNG intact if rendering fails.
rp_s2_finish(f,c,'Fig1_render');
movefile(fullfile(c.output,'Fig1_render.fig'),fullfile(c.output,'Fig1_irradiance.fig'),'f');
file=fullfile(c.output,'Fig1_irradiance.png');movefile(fullfile(c.output,'Fig1_render.png'),file,'f');
end
