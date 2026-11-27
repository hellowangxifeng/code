function files=rp_stage2_plot_sensitivity(s,result)
%RP_STAGE2_PLOT_SENSITIVITY One compact 3x3 plot, z=5 m, weak GG and two oceans.
c=s.config; [colors,styles]=rp_s2_style(); fields={'meanReceivedPowerW','aber','op'};
labels={'Mean received power (mW)','ABER','OP'};
xlabels={'Position jitter sigma (mm)','FOV half-angle (mrad)','Angular error sigma (mrad)'};
ranges=[20 120;30 150;10 60];
titles={'Rayleigh position','Finite FOV','Angular error'};
f=figure('Visible','off','Color','w','Position',[40 40 1100 820]); t=tiledlayout(f,3,3,'TileSpacing','compact','Padding','compact');
for row=1:3
    for parameter=1:3
        ax=nexttile(t); hold(ax,'on'); hh=gobjects(1,8); names=[result.schemeNames,{'Weak EGG: line','Strong EGG: x'}];
        for ocean=1:2
            for scheme=1:6
                T=result.table; T=T(T.oceanId==ocean & T.schemeId==scheme & T.parameterId==parameter,:);
                y=T.(fields{row}); if row==1, y=1000*y; end
                marker='none'; if ocean==2, marker='x'; end
                h=plot(ax,T.parameterValue*1000,y,'Color',colors(scheme,:),'LineStyle',styles{scheme}, ...
                    'LineWidth',1.1,'Marker',marker,'MarkerIndices',1:6:height(T),'MarkerSize',3);
                if ocean==1, hh(scheme)=h; end
            end
        end
        hh(7)=plot(ax,NaN,NaN,'k-'); hh(8)=plot(ax,NaN,NaN,'k-x','MarkerSize',3);
        xlabel(ax,xlabels{parameter}); ylabel(ax,labels{row}); rp_s2_axes(ax,ranges(parameter,:),row>1);
        if row==1, title(ax,titles{parameter},'FontSize',10); end
        if row==1 && parameter==1, rp_s2_legend(ax,hh,names,1); end
    end
end
title(t,'Fig 5 | z=5 m, weak GG, baseline FOV half-angle 90 mrad; x markers: strong EGG','FontSize',10);
files={rp_s2_finish(f,c,'Fig5_sensitivity')};
end
