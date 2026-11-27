function risMotionSummary = summarize_RIS_motion_existing(risInputFolder,risOutputFolder)
% Descriptive motion statistics from existing CSV trajectories; no rerun.
% Equal sample weights, normalization N, mean removed, trend retained.
% These are finite-record statistics, not calibrated sea-state predictions.
risStatsFolder = fileparts(mfilename('fullpath'));
if nargin<1, risInputFolder=fullfile(risStatsFolder,'FigOutput_RIS_WindDirection'); end
if nargin<2, risOutputFolder=fullfile(risStatsFolder,'motion_parameterization'); end
if ~isfolder(risOutputFolder), mkdir(risOutputFolder); end
risFiles = dir(fullfile(risInputFolder,'RIS_trajectory_U*_Phi*.csv'));
assert(numel(risFiles)==12,'Expected 12 existing trajectories.');
risRows = cell(numel(risFiles),1);
for risIndex = 1:numel(risFiles)
    risToken = regexp(risFiles(risIndex).name,'U(\d+)_Phi(\d+)','tokens','once');
    risT = readtable(fullfile(risFiles(risIndex).folder,risFiles(risIndex).name),TextType="string");
    risTime = risT.elapsedTimeS;
    risP = risT{:,{'centerXM','centerYM','centerZM'}};
    risAngles = risT{:,{'rollRad','pitchRad','yawRad'}};
    assert(height(risT)==401 && all(isfinite([risTime,risP,risAngles]),'all'));
    assert(all(diff(risTime)>0) && max(abs(diff(risTime)-0.025))<1e-9);
    assert(abs(risT.seaTimeS(1)-25)<1e-9 && abs(risT.seaTimeS(end)-35)<1e-9);
    risC = cov(risP,1);
    risCentered = risP-mean(risP,1);
    assert(abs(trace(risC)-mean(sum(risCentered.^2,2)))<1e-10);
    risPhi=risAngles(:,1); risTheta=risAngles(:,2); risPsi=risAngles(:,3);
    risNormals=[cos(risPsi).*sin(risTheta).*cos(risPhi)+sin(risPsi).*sin(risPhi), ...
        sin(risPsi).*sin(risTheta).*cos(risPhi)-cos(risPsi).*sin(risPhi), ...
        cos(risTheta).*cos(risPhi)];
    assert(max(abs(sum(risNormals.^2,2)-1))<1e-10);
    risMeanNormal=mean(risNormals,1); risMeanNormal=risMeanNormal/norm(risMeanNormal);
    risNormalAngles=acos(max(-1,min(1,risNormals*risMeanNormal.')));
    risEulerStd=std(risAngles,1,1)*180/pi;
    risTrend=[ones(size(risTime)),risTime]\risP;
    risResidual=risP-[ones(size(risTime)),risTime]*risTrend;
    risRow.windMps=str2double(risToken{1});
    risRow.directionDeg=str2double(risToken{2});
    risRow.sampleCount=height(risT);
    risRow.varXM2=risC(1,1);
    risRow.varYM2=risC(2,2);
    risRow.varZM2=risC(3,3);
    risRow.covXYM2=risC(1,2);
    risRow.covXZM2=risC(1,3);
    risRow.covYZM2=risC(2,3);
    risRow.sigma3DM=sqrt(trace(risC));
    risRow.sigmaHorizontalM=sqrt(risC(1,1)+risC(2,2));
    risRow.stdRollDeg=risEulerStd(1);
    risRow.stdPitchDeg=risEulerStd(2);
    risRow.stdYawDeg=risEulerStd(3);
    risRow.normalAngularRmsDeg=sqrt(mean(risNormalAngles.^2))*180/pi;
    risRow.linearDriftSpeedMps=norm(risTrend(2,:));
    risRow.sigma3DLinearDetrendedM=sqrt(mean(sum(risResidual.^2,2)));
    risRow.meanXM=mean(risP(:,1));
    risRow.meanYM=mean(risP(:,2));
    risRow.meanZM=mean(risP(:,3));
    risRows{risIndex}=risRow;
end
risMotionSummary=sortrows(struct2table(vertcat(risRows{:})),{'directionDeg','windMps'});
writetable(risMotionSummary,fullfile(risOutputFolder,'existing_trajectory_motion_statistics.csv'));
risShow=risMotionSummary(risMotionSummary.directionDeg==30,:);
risFig=figure('Color','w','Visible','off','Position',[80,80,1450,580]);
risLayout=tiledlayout(risFig,1,2,'TileSpacing','compact','Padding','compact');
risAx=nexttile(risLayout);
semilogy(risAx,risShow.windMps,risShow{:,{'varXM2','varYM2','varZM2'}},'-o','LineWidth',1.6);
xlabel(risAx,'Wind speed U (m/s)'); ylabel(risAx,'Centroid displacement variance (m^2)');
xticks(risAx,[3,5,9,12]); grid(risAx,'on'); legend(risAx,{'x_G','y_G','z_G'},'Location','northwest');
title(risAx,'Translation: mean removed, drift retained');
risAx=nexttile(risLayout);
plot(risAx,risShow.windMps,risShow.normalAngularRmsDeg,'-o','LineWidth',1.6);
xlabel(risAx,'Wind speed U (m/s)'); ylabel(risAx,'Normal angular RMS (deg)');
xticks(risAx,[3,5,9,12]); grid(risAx,'on');
title(risAx,'Orientation: angle about mean board normal');
title(risLayout,{'RIS motion | wind direction parameter 30 deg', ...
    'Existing simulation records: sea time 25-35 s, seed 600, 401 samples; normalization N'},'FontSize',15);
exportgraphics(risFig,fullfile(risOutputFolder,'RIS_motion_strength_vs_wind.png'),'Resolution',170);
savefig(risFig,fullfile(risOutputFolder,'RIS_motion_strength_vs_wind.fig'));
close(risFig);
disp(jsonencode(table2struct(risShow)));
end
