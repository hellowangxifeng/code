function p = rp_phase(c,g,a)
%RP_PHASE Steer each tile as close to world +z as its discrete commands allow.
% Linear steering is the production objective. Optional quadratic/tanh models
% compare cancellation of incident curvature, without circular expansion.
arguments
    c (1,1) struct
    g (1,1) struct
    a (1,1) struct
end
ids=a.activeIds; n=numel(ids); ref=a.referenceXY;
incoming=rp_layer_path(ref-g.sourceLocal(1:2).',[-g.sourceLocal(3),g.hIn], ...
    [c.optics.nWater,c.optics.nGlass],c.numerics.bisectionIterations);
[b0,command,direction,valid]=rp_vertical_command(c,g,incoming.q);
qActual=incoming.q+b0;
targets=repmat(g.pRis(1:2).',n,1);
for m=find(valid).'
    [target,reaches]=rp_carrier_landing(c,g,ref(m,:),qActual(m,:));
    valid(m)=reaches;
    if reaches, targets(m,:)=target; end
end
% Expand diffraction about each ACTUAL carrier, including saturated commands.
out=rp_paraxial(c,g,ref,targets); out.valid=out.valid & valid;
q=reshape(incoming.q.',2,1,n); identity=repmat(eye(2),1,1,n);
kw=incoming.kz(:,1); kg=incoming.kz(:,2);
J=identity.*reshape(-g.sourceLocal(3)./kw+g.hIn./kg,1,1,n)+ ...
    pagemtimes(q,pagetranspose(q)).*reshape(-g.sourceLocal(3)./kw.^3+g.hIn./kg.^3,1,1,n);
incomingH=pagemldivide(J,identity);
H=-incomingH; fullRequestedH=H;
if c.phase.separable, H(1,2,:)=0; H(2,1,:)=0; end
if strcmp(c.phase.model,'linear'), H(:)=0; end
requestedH=H;
extent=c.array.apertureM/2+abs(a.unitCenters(ids,:)-ref);
h11=reshape(H(1,1,:),n,1); h12=reshape(H(1,2,:),n,1);
h21=reshape(H(2,1,:),n,1); h22=reshape(H(2,2,:),n,1);
worst=abs(b0)+[abs(h11).*extent(:,1)+abs(h12).*extent(:,2), ...
    abs(h21).*extent(:,1)+abs(h22).*extent(:,2)];
bound=c.optics.nAir*sind(c.device.maxAngleDeg);
controlled=all(worst<=bound+1e-13,2) & out.valid;
if ismember(c.phase.model,{'tanh','quadratic_raw'}), controlled=out.valid; end
curvature=H.*reshape(controlled,1,1,n);
incident=rp_incident_field(c,g,ref);
% Common world-horizontal phase at each actual upper-face exit point.
% RX horizontal motion does not enter the steering objective.
glassKz=sqrt(max(realmin,c.optics.nGlass^2-sum(qActual.^2,2)));
top=[ref+g.hOut*qActual./glassKz,repmat(g.zTop,n,1)]*g.RRis.'+g.pRis.';
glassOpl=c.optics.nGlass^2*g.hOut./glassKz;
phase0=mod(2*pi/c.source.wavelength*(c.optics.nAir*top(:,3)-glassOpl)-incident.phaseRad,2*pi);
weights=a.unitDirectPowerW(ids);
if sum(weights)<=0, weights=a.unitPowerW(ids); end
angle=atan2d(vecnorm(direction(:,1:2),2,2),direction(:,3));
p=struct('unitIds',ids,'referenceXY',ref,'targetsXY',targets, ...
    'qInAtReference',incoming.q,'incomingH',incomingH,'commandDeg',command, ...
    'gradientOverK0',b0,'phaseCurvature',curvature,'phaseAtReferenceRad',phase0, ...
    'outgoing',out,'controlled',controlled, ...
    'controlledPowerFraction',sum(weights(controlled))/max(sum(weights),realmin), ...
    'referenceTransmittedPowerW',sum(weights.*rp_fresnel(c.optics.nGlass,c.optics.nAir,qActual)), ...
    'phaseModel',c.phase.model,'requestedPhaseCurvature',requestedH, ...
    'fullRequestedPhaseCurvature',fullRequestedH,'quadraticWorstGradientOverK0',worst, ...
    'incidentPhaseAtReferenceRad',incident.phaseRad,'referenceDirectionWorld',direction, ...
    'referenceAngleToWorldZDeg',angle,'verticalCarrierValid',valid);
p.objective='minimum angle to world +z on the discrete per-axis command grid';
end
