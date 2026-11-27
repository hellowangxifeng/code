function suite = run_RIS_wind_comparison(outputFolder,mode)
%RUN_RIS_WIND_COMPARISON Four wind rows, three directions and a 30-deg path.
% Snapshots use the same absolute sea time 27 s. Paths cover sea time 25-35 s.
% Saves the 4-by-4 overview, 16 individual FIG/PNG panels and unscaled data.
% mode='replot' redraws saved motions/scenes without rerunning dynamics.
    folder = fileparts(mfilename('fullpath'));
    if nargin == 0
        outputFolder = fullfile(folder,'FigOutput_RIS_WindDirection');
    end
    if nargin<2, mode='simulate'; end
    mode = validatestring(mode,{'simulate','replot'});
    if strcmp(mode,'replot')
        resultFile = fullfile(outputFolder,'RIS_wind_direction_results.mat');
        suite = load(resultFile);
        [panelFiles,overviewFile] = plotResults(suite.motions,suite.scenes,suite.meta,outputFolder);
        meta = suite.meta; meta.risDisplay = 'glass plane only; frame retained in dynamics';
        save(resultFile,'panelFiles','overviewFile','meta','-append');
        suite.panelFiles = panelFiles; suite.overviewFile = overviewFile; suite.meta = meta;
        return
    end
    if ~isfolder(outputFolder), mkdir(outputFolder); end
    sea = get_default_sea_params();
    winds = [3,5,9,12]; directions = [30,60,90];
    snapshotTime = 2; startTime = 25; duration = 10; seed = 600;
    zLimits = [-0.5,0.5];
    x = linspace(sea.x_range(1),sea.x_range(2),sea.grid_density);
    y = linspace(sea.y_range(1),sea.y_range(2),sea.grid_density);
    motions = cell(4,3); scenes = cell(4,3);
    values = zeros(12,20); row = 0;
    for j = 1:4
        for k = 1:3
            fprintf('Sea motion: U=%g m/s, wind=%g deg\n',winds(j),directions(k));
            o = RIS_float(struct('windSpeed',winds(j),'windAngleDeg',directions(k), ...
                'startTime',startTime,'duration',duration,'seed',seed,'plot',false));
            [X,Y,Z] = RIS_sea_surface(o,x,y,snapshotTime);
            motions{j,k} = o;
            scenes{j,k} = struct('X',X,'Y',Y,'Z',Z);
            sample = find(abs(o.t-snapshotTime)<1e-10,1);
            center = o.centerOfMass(:,sample);
            row = row+1;
            values(row,:) = [winds(j),directions(k),o.spectrum.significantWaveHeightM, ...
                center.',o.rollDeg(sample),o.pitchDeg(sample),o.yawDeg(sample), ...
                min(Z,[],'all'),max(Z,[],'all'),mean(abs(Z(:))>zLimits(2)), ...
                min(o.centerOfMass(3,:)),max(o.centerOfMass(3,:)), ...
                o.validity.maxTiltDeg,max(o.validity.radiationFitRelativeRmsError), ...
                o.validity.maxIncidentWaveBoardMismatchM,o.validity.bottomGapSampleFraction, ...
                o.validity.topWettedSampleFraction,o.validity.minTopClearanceM];
            exportTrajectory(o,outputFolder,winds(j),directions(k));
            fprintf('U=%g, direction=%g: Hs=%.4f m, centre z at t=27 s: %+.4f m\n', ...
                winds(j),directions(k),o.spectrum.significantWaveHeightM,center(3));
        end
    end
    summary = array2table(values,'VariableNames',{'windMps','directionDeg','HsM', ...
        'snapshotCenterXM','snapshotCenterYM','snapshotCenterZM','snapshotRollDeg', ...
        'snapshotPitchDeg','snapshotYawDeg','seaMinM','seaMaxM','seaOutsideZFraction', ...
        'pathMinZM','pathMaxZM','maxTiltDeg','radiationFitRelativeRmsError', ...
        'maxBottomGapAbsM','bottomGapSampleFraction','topWettedSampleFraction','minTopClearanceM'});
    writetable(summary,fullfile(outputFolder,'RIS_case_summary.csv'));
    meta = struct('winds',winds,'directionsDeg',directions,'seed',seed, ...
        'startTimeS',startTime,'durationS',duration,'snapshotElapsedTimeS',snapshotTime, ...
        'snapshotSeaTimeS',startTime+snapshotTime,'zLimitsM',zLimits, ...
        'frequencyBandHz',[sea.f_min,sea.f_max],'frequencyBins',sea.Fre_segment, ...
        'directionBins',sea.Dir_segment,'quadrature','trapezoid','amplitudeScale',1, ...
        'plotBoxAspectRatio',[1,1,0.7],'trajectoryPoint','centre of mass', ...
        'displayCropOnly',true,'hydrodynamicallyValidated',false, ...
        'risDisplay','glass plane only; frame retained in dynamics');
    meta.glassSizeM=[motions{1,1}.parameters.Lx motions{1,1}.parameters.Ly];
    meta.flotation=motions{1,1}.flotation;
    [panelFiles,overviewFile] = plotResults(motions,scenes,meta,outputFolder);
    save(fullfile(outputFolder,'RIS_wind_direction_results.mat'), ...
        'motions','scenes','summary','meta','panelFiles','overviewFile');
    suite = struct('motions',{motions},'scenes',{scenes},'summary',summary, ...
        'meta',meta,'panelFiles',panelFiles,'overviewFile',overviewFile);
end

function [panelFiles,overviewFile] = plotResults(motions,scenes,meta,outputFolder)
    winds = meta.winds; directions = meta.directionsDeg;
    snapshotTime = meta.snapshotElapsedTimeS; zLimits = meta.zLimitsM;
    panelFiles = strings(4,4);
    % Save each table cell at a legible standalone size.
    for j = 1:4
        for k = 1:4
            f = newFigure(880,720);
            ax = nexttile(tiledlayout(f,1,1,'Padding','compact'));
            drawCell(ax,j,k,motions,scenes,snapshotTime,zLimits);
            if k<=3
                name = sprintf('RIS_sea_U%02d_Phi%02d_t27',winds(j),directions(k));
            else
                name = sprintf('RIS_center_path_U%02d_Phi30',winds(j));
            end
            panelFiles(j,k) = saveFigure(f,outputFolder,name);
            close(f);
        end
    end

    f = newFigure(2000,1600);
    layout = tiledlayout(f,4,4,'Padding','compact','TileSpacing','compact');
    for j = 1:4
        for k = 1:4
            ax = nexttile(layout);
            drawCell(ax,j,k,motions,scenes,snapshotTime,zLimits);
            set(ax,'FontSize',9);
        end
    end
    sgtitle(layout,{'Sea surface and RIS: winds 3, 5, 9, 12 m/s (rows)', ...
        'Columns: wind 30 deg | wind 60 deg | wind 90 deg | 30-deg centre trajectory', ...
        'Snapshots: sea t=27 s | paths: sea t=25-35 s | z limits +/-0.5 m (display crop)'}, ...
        'FontWeight','normal','FontSize',15);
    overviewFile = saveFigure(f,outputFolder,'RIS_wind_direction_table');
    close(f);
end

function drawCell(ax,row,column,motions,scenes,time,zLimits)
    if column<=3
        o = motions{row,column};
        RIS_draw_state(ax,o,time,scenes{row,column},zLimits);
        title(ax,sprintf('U=%g m/s | wind %g deg | t_{sea}=27 s', ...
            o.parameters.windSpeed,o.parameters.windAngleDeg),'FontWeight','normal');
    else
        o = motions{row,1};
        p = o.centerOfMass;
        plot3(ax,p(1,:),p(2,:),p(3,:),'Color',[0.15,0.3,0.65],'LineWidth',1.4);
        hold(ax,'on');
        % Projection keeps the horizontal path visible when z is cropped.
        plot3(ax,p(1,:),p(2,:),repmat(zLimits(1),size(o.t)),':', ...
            'Color',[0.5,0.5,0.5],'LineWidth',1);
        inRange = find(p(3,:)>=zLimits(1) & p(3,:)<=zLimits(2));
        if ~isempty(inRange)
            scatter3(ax,p(1,inRange(1)),p(2,inRange(1)),p(3,inRange(1)),35, ...
                [0.1,0.55,0.25],'filled');
            scatter3(ax,p(1,inRange(end)),p(2,inRange(end)),p(3,inRange(end)),35, ...
                [0.8,0.15,0.1],'filled');
        end
        xlabel(ax,'x_G (m)'); ylabel(ax,'y_G (m)'); zlabel(ax,'z_G (m)');
        zlim(ax,zLimits); zticks(ax,[-0.5,0,0.5]);
        pbaspect(ax,[1,1,0.7]); view(ax,-37.5,30); grid(ax,'on'); box(ax,'on');
        outside = mean(abs(p(3,:))>zLimits(2));
        title(ax,{sprintf('U=%g m/s | wind 30 deg | centre path',o.parameters.windSpeed), ...
            sprintf('t_{sea}=25-35 s | %.0f%% outside z view',100*outside)},'FontWeight','normal');
    end
end

function f = newFigure(width,height)
    f = figure('Color','w','Visible','off','Units','pixels','Position',[40,40,width,height]);
end

function path = saveFigure(f,folder,name)
    path = string(fullfile(folder,[name,'.fig']));
    savefig(f,path);
    exportgraphics(f,fullfile(folder,[name,'.png']),'Resolution',160);
end

function exportTrajectory(o,folder,wind,direction)
    values = [o.t.',o.absoluteTime.',o.centerOfMass.',o.eta(4:6,:).',o.bodyVelocity.'];
    names = {'elapsedTimeS','seaTimeS','centerXM','centerYM','centerZM', ...
        'rollRad','pitchRad','yawRad','uMps','vMps','wMps','pRadps','qRadps','rRadps'};
    writetable(array2table(values,'VariableNames',names), ...
        fullfile(folder,sprintf('RIS_trajectory_U%02d_Phi%02d.csv',wind,direction)));
end
