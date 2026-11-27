function file=rp_stage2_plot_radial(s)
%RP_STAGE2_PLOT_RADIAL Actual annular means and encircled energy to 40 cm.
c=s.config; [colors,styles]=rp_s2_style(); rows=[]; handles=gobjects(1,6);
labels={'Annular mean / peak','Encircled power fraction'};
f=figure('Visible','off','Color','w','Position',[40 40 1060 650]); t=tiledlayout(f,2,3,'TileSpacing','compact','Padding','compact');
for j=1:3
    d=s.distances{j};
    for scheme=1:6
        b=d.bank{d.comparisonIndices(scheme)}; w=b.wave; [X,Y]=meshgrid(w.x,w.y);
        radius=hypot(X-b.spot.centerXY(1),Y-b.spot.centerXY(2));
        edges=0:.005:1.5; rr=(edges(1:end-1)+edges(2:end))/2; n=numel(rr);
        bin=discretize(radius(:),edges); good=isfinite(bin); intensity=w.irradianceWm2(:);
        total=accumarray(bin(good),intensity(good),[n 1],@sum,0); count=accumarray(bin(good),1,[n 1],@sum,0);
        meanI=total./max(count,1); meanI(count==0)=NaN; F=cumsum(total)*w.cellAreaM2/b.spot.powerW;
        ax=nexttile(t,j); hold(ax,'on');
        h=plot(ax,100*rr,meanI/max(meanI,[],'omitnan'),'Color',colors(scheme,:),'LineStyle',styles{scheme},'LineWidth',1.2);
        if j==1, handles(scheme)=h; end
        ax=nexttile(t,j+3); hold(ax,'on'); plot(ax,100*rr,F,'Color',colors(scheme,:),'LineStyle',styles{scheme},'LineWidth',1.2);
        rows=[rows;repmat([d.distanceM scheme],n,1),rr(:),meanI,F]; %#ok<AGROW>
    end
    for row=1:2
        ax=nexttile(t,j+3*(row-1)); xlabel(ax,'Radius from centroid (cm)');
        ylabel(ax,labels{row}); title(ax,sprintf('z=%g m',d.distanceM),'FontSize',10);
        rp_s2_axes(ax,[0 40],false); if row==2, ylim(ax,[0 1.02]); end
    end
end
rp_s2_legend(nexttile(t,1),handles,s.distances{1}.comparisonNames);
title(t,'Fig 8 | Actual radial shape and energy concentration','FontSize',11);
file=rp_s2_finish(f,c,'Fig8_radial_profiles');
writetable(array2table(rows,'VariableNames',{'distanceM','schemeId','radiusM','annularMeanWm2','encircledFraction'}),fullfile(c.output,'radial_profiles.csv'));
end
