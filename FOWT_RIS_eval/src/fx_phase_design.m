function p = fx_phase_design(c,g,a,target)
%FX_PHASE_DESIGN Feasibility-aware per-unit linear phase for a common direction or a point.
%   target.mode = 'direction': target.u is a world unit vector (collimated beam).
%   target.mode = 'point'    : target.p is a world [x y z] receiver point (focused beam).
% The command is snapped to the real device grid and checked against the
% per-axis air-calibrated limit.  Unlike rp_vertical_command, all units share
% one target, so the array stays collimated even when +z is unreachable.
arguments
    c (1,1) struct
    g (1,1) struct
    a (1,1) struct
    target (1,1) struct
end
ids=a.activeIds; n=numel(ids); ref=a.referenceXY;
incoming=rp_layer_path(ref-g.sourceLocal(1:2).',[-g.sourceLocal(3),g.hIn], ...
    [c.optics.nWater,c.optics.nGlass],c.numerics.bisectionIterations);
qIn=incoming.q;
% First guess of the upper-face exit points for a point target.
glassKz=sqrt(max(realmin,c.optics.nGlass^2-sum(qIn.^2,2)));
topBody=ref+g.hOut*qIn./glassKz;
topWorld=[topBody,repmat(g.zTop,n,1)]*g.RRis.'+g.pRis.';
if strcmpi(target.mode,'direction')
    dTarget=repmat(target.u(:).',n,1);
    dTarget=dTarget./vecnorm(dTarget,2,2);
elseif strcmpi(target.mode,'point')
    dTarget=target.p(:).'-topWorld;
    dTarget=dTarget./vecnorm(dTarget,2,2);
else
    error('fx:Target','Unknown target mode.');
end
qTarget=c.optics.nAir*(dTarget*g.RRis); qTarget=qTarget(:,1:2);
bRequested=qTarget-qIn;
maxStep=c.device.angleStepDeg;
maxCmd=floor(c.device.maxAngleDeg/maxStep+1e-12)*maxStep;
nAir=c.optics.nAir;
command=round(asind(max(-1,min(1,bRequested/nAir)))/maxStep)*maxStep;
command=max(-maxCmd,min(maxCmd,command));
qOut=qIn+nAir*sind(command);
[command,qOut,valid]=repairCommands(c,g,qIn,qTarget,command,nAir,maxCmd,maxStep);
glassKz=sqrt(max(realmin,c.optics.nGlass^2-sum(qOut.^2,2)));
topBody=ref+g.hOut*qOut./glassKz;
topWorld=[topBody,repmat(g.zTop,n,1)]*g.RRis.'+g.pRis.';
directionWorld=[qOut,sqrt(max(0,nAir^2-sum(qOut.^2,2)))]/nAir*g.RRis.';
bActual=qOut-qIn;
[targets,reaches]=deal(nan(n,2),false(n,1));
for m=1:n
    [targets(m,:),reaches(m)]=rp_carrier_landing(c,g,ref(m,:),qOut(m,:));
end
valid=valid & reaches;
targets(~valid,:)=repmat(g.pRis(1:2).',sum(~valid),1);
out=rp_paraxial(c,g,ref,targets); out.valid=out.valid & valid;
% Incident curvature (same construction as rp_phase.m).
q=reshape(qIn.',2,1,n); identity=repmat(eye(2),1,1,n);
kw=incoming.kz(:,1); kg=incoming.kz(:,2);
J=identity.*reshape(-g.sourceLocal(3)./kw+g.hIn./kg,1,1,n)+ ...
    pagemtimes(q,pagetranspose(q)).*reshape(-g.sourceLocal(3)./kw.^3+g.hIn./kg.^3,1,1,n);
incomingH=pagemldivide(J,identity);
incident=rp_incident_field(c,g,ref);
% Optional curvature cancellation within each unit, budget-limited.
curvature=zeros(2,2,n);
if isfield(target,'curvature') && strcmpi(target.curvature,'cancel')
    scale=1;
    if isfield(target,'curvatureScale'), scale=target.curvatureScale; end
    H=-scale*incomingH;
    extent=c.array.apertureM/2+abs(a.unitCenters(ids,:)-ref);
    h11=reshape(H(1,1,:),n,1); h12=reshape(H(1,2,:),n,1);
    h21=reshape(H(2,1,:),n,1); h22=reshape(H(2,2,:),n,1);
    worst=abs(bActual)+[abs(h11).*extent(:,1)+abs(h12).*extent(:,2), ...
        abs(h21).*extent(:,1)+abs(h22).*extent(:,2)];
    over=worst>c.optics.nAir*sind(c.device.maxAngleDeg)+1e-13;
    if any(over(:))
        bmax=c.optics.nAir*sind(c.device.maxAngleDeg);
        denom=max(worst-abs(bActual),realmin);
        shrink=min(ones(n,2),(bmax-abs(bActual))./denom);
        shrink=min(shrink,1);
        S=zeros(2,2,n);
        S(1,1,:)=reshape(shrink(:,1),1,1,n);
        S(2,2,:)=reshape(shrink(:,2),1,1,n);
        H=pagemtimes(S,H);
    end
    curvature=H;
end
glassOpl=c.optics.nGlass^2*g.hOut./glassKz;
if strcmpi(target.mode,'direction')
    coord=topWorld*target.u(:);
    phase0=mod(2*pi/c.source.wavelength*(nAir*coord-glassOpl)-incident.phaseRad,2*pi);
else
    airDist=vecnorm(target.p(:).'-topWorld,2,2);
    phase0=mod(-2*pi/c.source.wavelength*(airDist+glassOpl)-incident.phaseRad,2*pi);
end
weights=a.unitDirectPowerW(ids);
if sum(weights)<=0, weights=a.unitPowerW(ids); end
angle=atan2d(vecnorm(directionWorld(:,1:2),2,2),directionWorld(:,3));
p=struct('unitIds',ids,'referenceXY',ref,'targetsXY',targets, ...
    'qInAtReference',qIn,'incomingH',incomingH,'commandDeg',command, ...
    'gradientOverK0',bActual,'phaseCurvature',curvature, ...
    'phaseAtReferenceRad',phase0,'outgoing',out,'controlled',valid, ...
    'controlledPowerFraction',sum(weights(valid))/max(sum(weights),realmin), ...
    'referenceTransmittedPowerW',sum(weights.*rp_fresnel(c.optics.nGlass,c.optics.nAir,qOut)), ...
    'phaseModel','linear','requestedPhaseCurvature',zeros(2,2,n), ...
    'fullRequestedPhaseCurvature',zeros(2,2,n),'quadraticWorstGradientOverK0',abs(bActual), ...
    'incidentPhaseAtReferenceRad',incident.phaseRad,'referenceDirectionWorld',directionWorld, ...
    'referenceAngleToWorldZDeg',angle,'verticalCarrierValid',valid);
if strcmpi(target.mode,'direction')
    p.objective='common world direction on the discrete per-axis command grid';
    p.targetDirectionWorld=target.u(:).';
else
    p.objective='common receiver point on the discrete per-axis command grid';
    p.targetPointWorld=target.p(:).';
end
end

function [command,qOut,valid]=repairCommands(c,g,qIn,qTarget,command,nAir,maxCmd,maxStep)
% Choose, in a small grid neighbourhood, the legal command closest to target.
ids=1:size(qIn,1);
qOut=qIn+nAir*sind(command);
bad=find(~(sum(qOut.^2,2)<1-1e-12));
offsets=(-2:2)*maxStep;
for j=bad(:).'
    bestScore=-Inf; bestCmd=command(j,:); bestQ=qOut(j,:);
    for dx=offsets
        for dy=offsets
            cmd=[command(j,1)+dx,command(j,2)+dy];
            if any(abs(cmd)>maxCmd+1e-12), continue; end
            qq=qIn(j,:)+nAir*sind(cmd);
            if sum(qq.^2)>=1-1e-12, continue; end
            d=[qq,sqrt(1-sum(qq.^2))];
            score=dot(d,[qTarget(j,:),sqrt(max(0,1-sum(qTarget(j,:).^2)))]);
            if score>bestScore, bestScore=score; bestCmd=cmd; bestQ=qq; end
        end
    end
    command(j,:)=bestCmd; qOut(j,:)=bestQ;
end
valid=sum(qOut.^2,2)<1-1e-12;
end
