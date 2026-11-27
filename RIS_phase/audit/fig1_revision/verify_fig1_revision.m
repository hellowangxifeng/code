function verification=verify_fig1_revision(s)
%VERIFY_FIG1_REVISION Inspect MC physics, density convergence and saved graphics.
arguments
    s (1,1) struct
end
c=s.config;data=load(fullfile(c.output,'Fig1_data.mat'));rows=zeros(3,10);
for j=1:3
    d=s.distances{j};p=data.direct{j};map=data.panels{j,1};
    balance=p.incidentPowerW-p.reflectedPowerW-p.receiverPowerW-p.nonForwardPowerW;
    assert(abs(balance)<1e-12,'rp:Fig1Energy','MC interface energy does not balance.');
    assert(max(abs(vecnorm(p.directionAir,2,2)-1))<1e-12,'rp:Fig1Snell','Air directions are not unit vectors.');
    assert(abs(map.halo.capturedPowerFraction-1)<1e-12,'rp:Fig1Halo','Some MC packet power is missing from the context map.');
    assert(map.ballisticPowerFractionInMainView>.999,'rp:Fig1MainBeam','Main beam is cropped.');
    leak=abs(map.mainMapPowerW/map.mainPacketPowerW-1);
    bandwidthChange=abs(map.halfBandwidthPeakWm2/map.spot.peakWm2-1);
    assert(leak<.005,'rp:Fig1Density','Density smoothing loses significant packet power.');
    core=map.spot.peakWm2;
    risPeaks=cellfun(@(b) b.spot.peakWm2,data.panels(j,2:4));
    % A sample maximum is bandwidth-sensitive. Report that uncertainty and
    % require the physical intensity ordering to survive halving bandwidth.
    assert(max(core,map.halfBandwidthPeakWm2)<min(risPeaks),'rp:Fig1Intensity', ...
        'Direct MC peak is not below every RIS peak under bandwidth refinement.');
    halves=cell(1,2);
    for part=1:2
        keep=mod(p.inputPacketIds,2)==part-1;
        subset=struct('landingXY',p.landingXY(keep,:),'powerW',2*p.powerW(keep), ...
            'ballistic',p.ballistic(keep));halves{part}=rp_fig1_mc_map(subset);
    end
    peakDifference=abs(halves{1}.spot.peakWm2/halves{2}.spot.peakWm2-1);
    widthDifference=abs(halves{1}.spot.equivalentRadiusM/halves{2}.spot.equivalentRadiusM-1);
    assert(peakDifference<.06 && widthDifference<.06,'rp:Fig1MCConvergence','Independent packet halves disagree.');
    assert(max(halves{1}.spot.peakWm2,halves{2}.spot.peakWm2)<min(risPeaks), ...
        'rp:Fig1HalfIntensity','Independent packet halves reverse the peak ordering.');
    rows(j,:)=[d.distanceM,p.incidentPowerW,p.receiverPowerW,balance, ...
        map.ballisticPowerFractionInMainView,map.totalPowerFractionInMainView, ...
        leak,bandwidthChange,peakDifference,widthDifference];
end
statistics=array2table(rows,'VariableNames',{'distanceM','incidentPowerW','allReceiverPowerW', ...
    'energyResidualW','ballisticFractionInMain','totalFractionInMain', ...
    'densityBoundaryLossFraction','halfBandwidthPeakRelativeChange', ...
    'independentHalfPeakRelativeDifference','independentHalfWidthRelativeDifference'});
f=openfig(fullfile(c.output,'Fig1_irradiance.fig'),'invisible');clean=onCleanup(@() close(f));drawnow;
axesAll=findall(f,'Type','axes');assert(numel(axesAll)==15,'rp:Fig1Panels','Need 12 main panels and 3 complete-MC context insets.');
allLimits=zeros(12,7);n=0;
for k=1:numel(axesAll)
    ax=axesAll(k);assert(ax.FontSize==10,'rp:Fig1TickFont','Ticks must be 10 pt.');
    if startsWith(ax.Tag,'Fig1_MC_full_halo_'),continue;end
    info=ax.UserData;map=data.panels{info.row,info.column};w=map.wave;
    edgeX=100*([w.x(1),w.x(end)]+[-1 1]*diff(w.x(1:2))/2-map.spot.centerXY(1));
    edgeY=100*([w.y(1),w.y(end)]+[-1 1]*diff(w.y(1:2))/2-map.spot.centerXY(2));
    assert(ax.XLim(1)<=edgeX(1) && ax.XLim(2)>=edgeX(2) && ...
        ax.YLim(1)<=edgeY(1) && ax.YLim(2)>=edgeY(2),'rp:Fig1DisplayCrop','Image pixels are cropped.');
    line=findobj(ax,'Tag','waist_ellipse');assert(isscalar(line) && strcmp(line.LineStyle,'--'));
    assert(all(line.XData>=ax.XLim(1) & line.XData<=ax.XLim(2)) && ...
        all(line.YData>=ax.YLim(1) & line.YData<=ax.YLim(2)),'rp:Fig1EllipseCrop','Waist ellipse is cropped.');
    e=rp_spot_ellipse(map.spot);quad=sum((e.offsetXY/map.spot.covarianceM2).*e.offsetXY,2);
    assert(max(abs(quad-4))<1e-10,'rp:Fig1Ellipse','Contour is not the covariance ellipse.');
    n=n+1;allLimits(n,:)=[info.row,info.column,ax.XLim,ax.YLim,max(abs(quad-4))];
end
assert(numel(findall(f,'Tag','waist_axis_1'))==6 && numel(findall(f,'Tag','waist_axis_2'))==6);
legends=findall(f,'Type','legend');assert(numel(legends)==12);
for k=1:numel(legends),assert(legends(k).UserData.overlapSamples==0,'rp:Fig1Legend','Legend covers the ellipse or axes.');end
limits=array2table(allLimits,'VariableNames',{'row','column','xMinCM','xMaxCM','yMinCM','yMaxCM','ellipseEquationError'});
verification=struct('status','passed','scope','Fig1 revision only','statistics',statistics,'limits',limits);
folder=fullfile(c.root,'audit','fig1_revision');
writetable(statistics,fullfile(folder,'mc_verification.csv'));writetable(limits,fullfile(folder,'display_verification.csv'));
save(fullfile(folder,'verification.mat'),'verification');
fid=fopen(fullfile(folder,'verification.json'),'w');guard=onCleanup(@() fclose(fid));
fprintf(fid,'%s',jsonencode(verification,PrettyPrint=true));
end
