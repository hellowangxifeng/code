function file=rp_stage2_plot_angles(s,angles)
%RP_STAGE2_PLOT_ANGLES Incidence dependence including direct-path critical angle.
c=s.config; [colors,styles]=rp_s2_style();
names={'No RIS direct','Steering-only','Flat top','Robust ABER','Robust OP','Max mean power'};
f=figure('Visible','off','Color','w','Position',[40 40 850 700]); t=tiledlayout(f,2,2,'TileSpacing','compact','Padding','compact');
for row=1:2
    for distance=1:2
        result=angles{distance}; ax=nexttile(t); hold(ax,'on'); hh=gobjects(1,6);
        fields={'aber','op'}; field=fields{row};
        for scheme=1:6
            T=result.table(result.table.schemeId==scheme,:);
            hh(scheme)=plot(ax,T.incidenceDeg,T.(field),'Color',colors(scheme,:),'LineStyle',styles{scheme},'LineWidth',1.3,'Marker','.','MarkerSize',7);
        end
        xlabel(ax,'Underwater incidence to normal (deg)'); ylabel(ax,upper(field));
        title(ax,sprintf('z=%g m | ON power 1 W',result.distanceM),'FontSize',10);
        rp_s2_axes(ax,[min(c.study.angleGridDeg),max(c.study.angleGridDeg)],true);
        xline(ax,result.criticalAngleDeg,'k:','Direct critical angle','LabelOrientation','horizontal','LabelVerticalAlignment','bottom','FontSize',8);
        if row==1 && distance==1, rp_s2_legend(ax,hh,names,2,7,[12 7]); end
    end
end
title(t,'Fig 10 | Carrier follows incidence; nominal widths frozen; direct beam includes subcritical angular tails','FontSize',10);
file=rp_s2_finish(f,c,'Fig10_incidence');
end
