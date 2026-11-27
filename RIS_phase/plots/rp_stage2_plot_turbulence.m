function file=rp_stage2_plot_turbulence(s,result)
%RP_STAGE2_PLOT_TURBULENCE Two GG strengths in each ocean/metric panel.
c=s.config; [colors,styles]=rp_s2_style();
f=figure('Visible','off','Color','w','Position',[40 40 850 700]); t=tiledlayout(f,2,2,'TileSpacing','compact','Padding','compact');
for row=1:2
    for ocean=1:2
        ax=nexttile(t); hold(ax,'on'); hh=gobjects(1,8); field={'aber','op'}; field=field{row};
        for atmosphere=1:2
            marker='none'; if atmosphere==2, marker='x'; end
            for scheme=1:6
                T=result.table; T=T(T.oceanId==ocean & T.atmosphereId==atmosphere & T.schemeId==scheme,:);
                h=plot(ax,T.transmitDBm,T.(field),'Color',colors(scheme,:),'LineStyle',styles{scheme}, ...
                    'LineWidth',1.2,'Marker',marker,'MarkerIndices',1:10:height(T),'MarkerSize',3);
                if atmosphere==1, hh(scheme)=h; end
            end
        end
        hh(7)=plot(ax,NaN,NaN,'k-'); hh(8)=plot(ax,NaN,NaN,'k-x','MarkerSize',3);
        xlabel(ax,'ON-symbol transmit power (dBm)'); ylabel(ax,upper(field));
            kinds={'Weak','Strong'}; title(ax,[kinds{ocean},' ocean EGG'],'FontSize',10);
        rp_s2_axes(ax,[c.study.powerDBm(1),c.study.powerDBm(end)],true);
        if row==1 && ocean==1, rp_s2_legend(ax,hh,[result.schemeNames,{'Weak GG: line','Strong GG: x'}],2); end
    end
end
title(t,'Fig 9 | z=7.5 m: ocean EGG x atmospheric GG; condition-specific 1 W design','FontSize',11);
file=rp_s2_finish(f,c,'Fig9_turbulence');
end
