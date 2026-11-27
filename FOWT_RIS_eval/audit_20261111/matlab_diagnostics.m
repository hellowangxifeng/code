function matlab_diagnostics()
%MATLAB_DIAGNOSTICS Independent probes. Writes only this audit's evidence.
root=fileparts(fileparts(fileparts(mfilename('fullpath'))));
oldpath=path; cleanup=onCleanup(@()path(oldpath)); %#ok<NASGU>
addpath(fullfile(root,'RIS_phase'),fullfile(root,'FOWT_RIS_eval','src'),fullfile(root,'sea_wave'));
rp_setup();
S=load(fullfile(root,'RIS_phase','results','validated_optics','ris_spots.mat'));
s=S.s; q=s.cases{1}; c=s.config; g=q.geometry; a=q.illumination;
I=load(fullfile(root,'FOWT_RIS_eval','results','improved_curvature.mat'));
n=numel(a.activeIds); P=a.unitDirectPowerW(a.activeIds); P=P/sum(P);
e=struct();
e.matlabVersion=version;
e.nominalGeometry=struct('rotation',g.RRis,'pRis',g.pRis,'rxDesired',g.rxDesired,'rxPlaneZ',g.rxPlaneZ,'incidenceDeg',g.incidenceDeg);
e.config=struct('lambda',c.source.wavelength,'sourcePowerW',c.source.power,'nWater',c.optics.nWater,'nGlass',c.optics.nGlass,'nAir',c.optics.nAir,'apertureM',c.array.apertureM,'device',c.device);
for which=1:2
    target=struct('mode','direction','u',I.p.targetDirectionWorld,'curvature','cancel');
    if which==2, target.u=[0 0 1]; end
    p=fx_phase_design(c,g,a,target);
    asym=abs(p.phaseCurvature(1,2,:)-p.phaseCurvature(2,1,:));
    bad=false(n,1); maxOver=0;
    for corner=[-1 -1;1 -1;1 1;-1 1].'
        delta=a.unitCenters(a.activeIds,:)-p.referenceXY+corner.'.*c.array.apertureM/2;
        profile=rp_phase_profile(c,p,(1:n).',delta);
        over=max(abs(profile.gradientOverK0)-c.optics.nAir*sind(c.device.maxAngleDeg),[],2);
        maxOver=max(maxOver,max(over)); bad=bad | over>1e-12;
    end
    tag=sprintf('curvature%d',which);
    e.(tag)=struct('maxHAsymmetry',max(asym,[],'all'),'gradientBudgetViolationUnits',nnz(bad),'violationPowerFraction',sum(P(bad)),'maxExcessGradient',maxOver,'reportedControlledPowerFraction',p.controlledPowerFraction);
end
t1=struct('mode','direction','u',I.p.targetDirectionWorld);
t2=t1; t2.u=2*t2.u;
p1=fx_phase_design(c,g,a,t1); p2=fx_phase_design(c,g,a,t2);
e.directionScale=struct('gradientDifference',max(abs(p1.gradientOverK0-p2.gradientOverK0),[],'all'),'maxPistonPhasorDifference',max(abs(exp(1i*p1.phaseAtReferenceRad)-exp(1i*p2.phaseAtReferenceRad))));
target=struct('mode','point','p',[I.f.centerXY,g.rxPlaneZ]);
p=fx_phase_design(c,g,a,target);
L=p.targetsXY; rx=target.p(1:2);
ideal=rp_paraxial(c,g,p.referenceXY,repmat(rx,n,1));
B=c.optics.nAir*sind(floor(c.device.maxAngleDeg/c.device.angleStepDeg)*c.device.angleStepDeg);
reachable=all(abs(ideal.q-p.qInAtReference)<=B+1e-12,2) & ideal.valid;
e.point=struct('requestedXY',rx,'actualCentroidXY',sum(P.*L,1),'centroidErrorM',norm(sum(P.*L,1)-rx),'exactPointReachablePower',sum(P(reachable)),'reportedControlledPowerFraction',p.controlledPowerFraction);
idealOriginal=rp_paraxial(c,g,a.referenceXY,repmat(g.rxDesired(1:2).',n,1));
originalReachable=all(abs(idealOriginal.q-q.phase.qInAtReference)<=B+1e-12,2) & idealOriginal.valid;
e.originalReceiver=struct('reachablePowerFraction',sum(P(originalReachable)),'distanceFromProductionSpotM',norm(q.footprint.centerXY-g.rxDesired(1:2).'),'distanceFromImprovedSpotM',norm(I.f.centerXY-g.rxDesired(1:2).'));
% Missing files / field references in probes, and model mode indexing.
e.motionProbePathExists=isfile(fullfile(root,'浮式平台六自由度实验数据','FOWT Platform Testing','FOWTmotion_11_3_45mmp8Hz.mat'));
M=load(fullfile(root,'FOWT_RIS_eval','results','measured_wave_ris_8.mat'));
o=M.o;
e.hydro=struct('massMatrixSize',size(o.massMatrix),'stiffnessDiagonal',diag(o.stiffnessMatrix),'radiationDampingSize',size(o.radiationDamping),'addedMass',o.radiationModel.infiniteAddedMass);
e.cachedFootprintFields=fieldnames(q.footprint);
e.savedImproved=struct('target',I.p.targetDirectionWorld,'center',I.f.centerXY,'mapPowerW',I.wave.mapPowerW,'maxHAsymmetry',max(abs(I.p.phaseCurvature(1,2,:)-I.p.phaseCurvature(2,1,:)),[],'all'));
% FFT bridge: default range can include mirrored negative-frequency bins.
t=(0:1/4:7.75).'; eta=0.01*cos(2*pi*1.25*t);
w=fx_wave_from_gauge(t,eta,minAmplitude=0);
syn=real(exp(1i*t*w.omega.')*w.coefficient);
e.fftLowSampleRate=struct('inputAmplitude',.01,'reconstructedAmplitude',(max(syn)-min(syn))/2,'rmsError',sqrt(mean((syn-eta).^2)));
files=dir(fullfile(root,'FOWT_RIS_eval','**','*.m'));
st=struct('file',{},'issues',{});
for j=1:numel(files)
    if contains(files(j).folder,'audit_20261111'), continue; end
    file=fullfile(files(j).folder,files(j).name);
    st(end+1)=struct('file',file,'issues',checkcode(file,'-struct')); %#ok<AGROW>
end
e.staticAnalysis=st;
file=fullfile(root,'FOWT_RIS_eval','audit_20261111','matlab_evidence.json');
fid=fopen(file,'w','n','UTF-8'); fprintf(fid,'%s',jsonencode(e,PrettyPrint=true)); fclose(fid);
e=rmfield(e,{'staticAnalysis','cachedFootprintFields'});
fprintf('%s\n',jsonencode(e,PrettyPrint=true));
end
