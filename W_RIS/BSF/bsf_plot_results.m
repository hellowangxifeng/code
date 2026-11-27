function bsf_plot_results(folder)
%BSF_PLOT_RESULTS Identical coordinates/scales and honest unsmoothed residuals.
arguments
    folder (1,1) string = fullfile(fileparts(mfilename('fullpath')),'results')
end
root = fileparts(mfilename('fullpath'));
loaded = load(fullfile(fileparts(root),'results','wris_results.mat'),'s');
summary = bsf_export_summary(folder);
for row = 1:6
    for col = 1:3
        file = fullfile(folder,sprintf('panel_r%d_c%d.mat',row,col));
        if ~isfile(file)
            continue
        end
        stored = load(file,'result');
        r = stored.result;
        if ~isfield(r,'transportMetrics')
            continue
        end
        original = loaded.s.panels{row,col}.out;
        f = figure('Visible','off','Position',[100,100,1500,900],'Color','w');
        cleanup = onCleanup(@() close(f));
        layout = tiledlayout(f,2,3,'TileSpacing','compact','Padding','compact');
        maxValue = max([r.IBSF(:);original.noTurbulence(:);r.paper.IBSF(:)]);
        maps = {original.noTurbulence,r.paper.IBSF,r.IBSF,r.Ib,r.Is, ...
            (r.IBSF-original.noTurbulence)./sqrt(r.standardError.^2+(original.standardError/original.gain).^2)};
        names = {'Original MC (1e6 photons)','Eq. (3-18): I_b + I_s', ...
            'Full HG transport: I_b + I_s','Ballistic I_b (deterministic)', ...
            'Scattered I_s (own color scale)','Residual / combined standard error'};
        for k = 1:6
            ax = nexttile(layout);
            values = maps{k};
            if k==6
                referenceSE = original.standardError/original.gain;
                populated = referenceSE>0 & original.noTurbulence>=5*referenceSE;
                values(~isfinite(values) | ~populated) = NaN;
            end
            heatmapImage = imagesc(ax,r.x,r.y,values.');
            if k==6
                set(heatmapImage,'AlphaData',isfinite(values.'));
                set(ax,'Color',[0.93,0.93,0.93]);
            end
            set(ax,'YDir','normal','FontSize',11);
            axis(ax,'image');
            xlabel(ax,'RIS x (m)');
            ylabel(ax,'RIS y (m)');
            title(ax,names{k},'Interpreter','tex');
            cb = colorbar(ax);
            if k==6
                clim(ax,[-4,4]);
                colormap(ax,divergingMap());
                cb.Label.String = 'z score; populated MC bins only';
            else
                if k==5
                    clim(ax,[0,max(r.Is,[],'all')]);
                else
                    clim(ax,[0,maxValue]);
                end
                colormap(ax,turbo(256));
                cb.Label.String = 'W m^{-2}';
            end
        end
        title(layout,sprintf(['Row %d, condition %d | incidence %.3f deg | model SE %.4f%% | ' ...
            'raw L2 %.2f%%, MC noise %.2f%%'],row,col,r.geometry.thetaDeg,100*r.relativeSE, ...
            100*r.transportMetrics.relativeL2,100*r.transportMetrics.referenceNoiseL2), ...
            'FontSize',15);
        exportgraphics(f,fullfile(folder,sprintf('comparison_r%d_c%d.png',row,col)),'Resolution',160);
        clear cleanup
    end
end
if ~isempty(summary)
    f = figure('Visible','off','Position',[100,100,1400,550],'Color','w');
    cleanup = onCleanup(@() close(f));
    layout = tiledlayout(f,1,2,'TileSpacing','compact');
    ax = nexttile(layout);
    plot(ax,1:height(summary),summary.rawL2pct,'o-', ...
        1:height(summary),summary.originalNoisePct,'s--','LineWidth',1.4);
    ylabel(ax,'Relative L2 (%)');
    xlabel(ax,'Snapshot: rows 1-6, conditions static / U=3 / U=5');
    legend(ax,'Transport vs original MC','Original MC noise estimate','Location','northwest');
    grid(ax,'on');
    title(ax,'Raw disagreement compared with the ORIGINAL noise floor');
    ax = nexttile(layout);
    plot(ax,1:height(summary),summary.numericalSEpct,'o-','LineWidth',1.4);
    hold(ax,'on');
    yline(ax,0.1,'--','0.1% numerical SE target');
    ylim(ax,[0,0.11]);
    ylabel(ax,'Numerical RMS standard error (%)');
    xlabel(ax,'Snapshot');
    grid(ax,'on');
    title(ax,'Uncertainty of the transport formula evaluation');
    exportgraphics(f,fullfile(folder,'precision_overview.png'),'Resolution',160);
end
end

function map = divergingMap()
t = linspace(0,1,128).';
blue = [t,t,ones(size(t))];
red = [ones(size(t)),flipud(t),flipud(t)];
map = [blue;red];
end
