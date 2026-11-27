function p=rp_finite_phase(c,g,a,family,width,margin,baseline)
%RP_FINITE_PHASE Quantized affine tiles, finite target transport and piston.
% All propagated designs use legal reference tilts; no unbounded lens phase.
if nargin<7, baseline=rp_phase(c,g,a); end
p=baseline;
if strcmp(family,'steering'), p.family='steering'; return; end
if strcmp(family,'unmodulated'), p=rp_zero_phase(c,g,a); return; end
assert(ismember(family,{'gaussian','flat_top','robust'}),'rp:Family','Unknown shaping family.');
n=numel(p.unitIds); weights=a.unitDirectPowerW(p.unitIds);
weights=weights/max(sum(weights),realmin);
ref=p.referenceXY; meanXY=sum(ref.*weights,1);
D=ref-meanXY; covariance=(D.*weights).'*D;
[V,L]=eig(covariance); W=V*diag(1./sqrt(max(diag(L),1e-12)))*V.';
white=D*W; radius=vecnorm(white,2,2);
[~,order]=sort(radius); cumulative=cumsum(weights(order));
u=zeros(n,1); u(order)=cumulative-weights(order)/2;
u=max(1e-7,min(1-1e-7,u));
angle=atan2(white(:,2),white(:,1));
if strcmp(family,'flat_top')
    radial=width*sqrt(u);
elseif strcmp(family,'robust')
    % Generalized Gaussian transport with a smooth shoulder/tail. Width and
    % carrier are selected against actual finite-aperture EGG/GG ABER/OP.
    exponent=c.study.robustExponent;
    radial=width*(-log1p(-u)).^(1/exponent)/sqrt(2*gamma(1+2/exponent));
else
    radial=width*sqrt(-log1p(-u)/2);
end
offset=radial.*[cos(angle),sin(angle)];
if strcmp(family,'gaussian')
    % Ordinary astigmatic Gaussian: keep incident anisotropy at equal RMS width.
    shape=V*diag(sqrt(max(diag(L),1e-12)))*V.';
    shape=shape/sqrt(trace(shape*shape.')/2);
    offset=offset*shape;
end
bound=c.optics.nAir*sind(floor(c.device.maxAngleDeg/c.device.angleStepDeg)*c.device.angleStepDeg);
q0=zeros(1,2);
tail=c.study.maximumTailPower/4;
for d=1:2
    lo=weightedQuantile(p.qInAtReference(:,d)-bound,weights,1-tail);
    hi=weightedQuantile(p.qInAtReference(:,d)+bound,weights,tail);
    assert(lo<=hi,'rp:CommonDirection','No common carrier for the retained power.');
    q0(d)=min(hi,max(lo,0));
end
% Move inside the feasible direction box to leave room for finite-width shaping.
q0=q0+margin*sign(q0);
assert(sum(q0.^2)<c.optics.nAir^2,'rp:Direction','Common direction is evanescent.');
[center,valid]=rp_carrier_landing(c,g,[0 0],q0);
assert(valid,'rp:Direction','Chosen common carrier misses the forward plane.');
targets=center+offset;
requested=rp_paraxial(c,g,ref,targets);
[b,command,saturated]=rp_command(c,requested.q-p.qInAtReference);
q=p.qInAtReference+b;
valid=sum(q.^2,2)<c.optics.nAir^2;
actual=nan(n,2);
for m=find(valid).'
    [actual(m,:),valid(m)]=rp_carrier_landing(c,g,ref(m,:),q(m,:));
end
actual(~valid,:)=targets(~valid,:);
p.gradientOverK0=b; p.commandDeg=command; p.phaseCurvature(:)=0;
p.targetsXY=actual; p.outgoing=rp_paraxial(c,g,ref,actual);
p.outgoing.valid=p.outgoing.valid & valid;
input=rp_incident_field(c,g,ref);
% Equal path phase at EACH assigned target; all tile pistons are wrapped.
p.phaseAtReferenceRad=mod(-input.phaseRad-2*pi/c.source.wavelength*p.outgoing.opl,2*pi);
p.phaseModel='linear'; p.family=family; p.targetWidthM=width;
p.targetCenterXY=center; p.targetMargin=margin;
p.requestedTargetsXY=targets; p.saturated=saturated;
p.saturatedPowerFraction=sum(weights(saturated));
p.referenceDirectionWorld=[q,sqrt(max(0,c.optics.nAir^2-sum(q.^2,2)))]/c.optics.nAir*g.RRis.';
p.objective='finite width transport with quantized affine tile gradients';
end

function v=weightedQuantile(x,w,f)
[x,idx]=sort(x); t=cumsum(w(idx)); v=x(find(t>=f*t(end),1));
end
