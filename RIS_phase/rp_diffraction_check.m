function audit = rp_diffraction_check(c,panel)
%RP_DIFFRACTION_CHECK Compare the outgoing quadratic OPL against exact layered OPL.
% Independent direct aperture quadrature for one bright transmitting tile.
% No FFT, frequency deposition or pixel-coherence approximation is used here.
arguments
    c (1,1) struct
    panel (1,1) struct
end
g=panel.geometry; a=panel.illumination; p=panel.phase;
weight=a.unitDirectPowerW(p.unitIds).*p.outgoing.valid;
if isfield(panel,'wave') && isfield(panel.wave,'unitMapPowerW') && any(panel.wave.unitMapPowerW>0)
    weight=panel.wave.unitMapPowerW.*p.outgoing.valid;
end
[~,m]=max(weight);
ref=p.referenceXY(m,:); center=a.unitCenters(p.unitIds(m),:);
% Test the geometric carrier and several interior points of the transmitted field.
M=p.outgoing.A(:,:,m)+p.incomingH(:,:,m)+p.phaseCurvature(:,:,m);
map=p.outgoing.B(:,:,m)\M;
first=0.4*(map*[c.array.apertureM(1)/2;0]).';
second=0.4*(map*[0;c.array.apertureM(2)/2]).';
first=first*min(1,0.25/max(norm(first),realmin));
second=second*min(1,0.25/max(norm(second),realmin));
eta=[0 0;first;-first;second;-second];
receiver=p.targetsXY(m,:)+eta;
rows=zeros(size(receiver,1),5); localRows=zeros(size(receiver,1),4);
warpedChange=zeros(size(receiver,1),1);
for j=1:size(receiver,1)
    [exactLow,quadraticLow,~,warpedLow]=integralAt(1,receiver(j,:),m,c,g,a,p,ref,center);
    [exactHigh,quadraticHigh,localHigh,warpedHigh]=integralAt(2,receiver(j,:),m,c,g,a,p,ref,center);
    fprintf('Direct check point %d/5, base refinement change %.3g / %.3g.\n',j,abs(exactHigh-exactLow)/max(exactHigh,realmin),abs(quadraticHigh-quadraticLow)/max(quadraticHigh,realmin));
    for refinement=[4 8 16]
        changes=[abs(exactHigh-exactLow)/max(exactHigh,realmin), ...
            abs(quadraticHigh-quadraticLow)/max(quadraticHigh,realmin), ...
            abs(warpedHigh-warpedLow)/max(warpedHigh,realmin)];
        if max(changes)<0.01, break; end
        exactLow=exactHigh; quadraticLow=quadraticHigh; warpedLow=warpedHigh;
        [exactHigh,quadraticHigh,localHigh,warpedHigh]=integralAt(refinement,receiver(j,:),m,c,g,a,p,ref,center);
        fprintf('  refinement %d change %.3g / %.3g.\n',refinement, ...
            abs(exactHigh-exactLow)/max(exactHigh,realmin),abs(quadraticHigh-quadraticLow)/max(quadraticHigh,realmin));
    end
    rows(j,:)=[exactHigh,quadraticHigh,abs(exactHigh-quadraticHigh)/max(exactHigh,realmin), ...
        abs(exactHigh-exactLow)/max(exactHigh,realmin), ...
        abs(quadraticHigh-quadraticLow)/max(quadraticHigh,realmin)];
    localRows(j,:)=[localHigh,abs(exactHigh-localHigh)/max(exactHigh,realmin), ...
        warpedHigh,abs(exactHigh-warpedHigh)/max(exactHigh,realmin)];
    warpedChange(j)=abs(warpedHigh-warpedLow)/max(warpedHigh,realmin);
end
audit=struct('unitId',p.unitIds(m),'receiverXY',receiver, ...
    'comparison',array2table(rows,'VariableNames', ...
    {'exactIntensity','quadraticIntensity','relativeKernelError','exactQuadratureChange','quadraticQuadratureChange'}));
audit.receiverLocalComparison=array2table(localRows,'VariableNames', ...
    {'localQuadraticIntensity','relativeKernelError','warpedIntensity','warpedRelativeKernelError'});
audit.fixedCenterComparison=audit.comparison;
audit.comparison.quadraticIntensity=localRows(:,3);
audit.comparison.relativeKernelError=localRows(:,4);
audit.comparison.quadraticQuadratureChange=warpedChange;
end

function [exact,quadratic,receiverLocal,warped]=integralAt(N,receiver,m,c,g,a,p,ref,center)
lambda=c.source.wavelength; k0=2*pi/lambda;
toReceiver=rp_paraxial(c,g,ref,receiver);
M=toReceiver.A+p.incomingH(:,:,m)+p.phaseCurvature(:,:,m);
[rotation,eigenvalues]=eig((M+M.')/2);
[~,order]=sort(abs(diag(eigenvalues)),'descend');
rotation=rotation(:,order); eigenvalues=eigenvalues(order,order);
extent=abs(rotation.')*c.array.apertureM.';
required=1.3*abs(diag(eigenvalues)).*extent.^2/lambda;
linear=(p.qInAtReference(m,:)+p.gradientOverK0(m,:)-toReceiver.q)*rotation;
required=required+2*extent.*abs(linear.')/lambda;
counts=N*2.^nextpow2(max(64,required));
fprintf('  Reference grid %d x %d (strong axis first).\n',counts(1),counts(2));
% Integrate the exact rotated rectangle, without a staircased binary mask.
% Split the first coordinate at polygon vertices; its second-coordinate
% limits are linear within each segment. Chunking bounds memory at any N.
corners=[-1 -1;-1 1;1 -1;1 1].*c.array.apertureM/2;
vertices=sort(unique(corners*rotation(:,1)));
u=[]; du=[];
for segment=1:numel(vertices)-1
    width=vertices(segment+1)-vertices(segment);
    cells=max(1,ceil(counts(1)*width/extent(1)));
    u=[u,vertices(segment)+((1:cells)-0.5)*width/cells]; %#ok<AGROW>
    du=[du,repmat(width/cells,1,cells)]; %#ok<AGROW>
end
low=-Inf(1,numel(u)); high=Inf(1,numel(u));
for dimension=1:2
    coefficient=rotation(dimension,2);
    if abs(coefficient)<1e-12, continue; end
    first=(-c.array.apertureM(dimension)/2-rotation(dimension,1)*u)/coefficient;
    second=(c.array.apertureM(dimension)/2-rotation(dimension,1)*u)/coefficient;
    low=max(low,min(first,second)); high=min(high,max(first,second));
end
vNodes=((1:counts(2)).'-0.5)/counts(2);
exactField=0; quadraticField=0; localField=0; warpedField=0;
eta=receiver-p.targetsXY(m,:);
for firstColumn=1:64:numel(u)
    columns=firstColumn:min(firstColumn+63,numel(u));
    uu=repmat(u(columns),counts(2),1);
    vv=low(columns)+vNodes.*(high(columns)-low(columns));
    weights=repmat(du(columns).*(high(columns)-low(columns))/counts(2),counts(2),1);
    xy=[uu(:),vv(:)]*rotation.'+center; delta=xy-ref;
    incident=rp_incident_field(c,g,xy);
    profile=rp_phase_profile(c,p,m,delta);
    device=p.phaseAtReferenceRad(m)+k0*profile.phaseOverK0;
    local=([receiver,g.rxPlaneZ]-g.pRis.')*g.RRis;
    path=rp_layer_path(local(1:2)-xy,[g.hOut,local(3)-g.zTop], ...
        [c.optics.nGlass,c.optics.nAir],c.numerics.bisectionIterations);
    transmission=rp_fresnel(c.optics.nGlass,c.optics.nAir,path.q);
    amplitude=sqrt(incident.irradianceWm2*a.ballisticModelScale).*weights(:);
    exactField=exactField+sum(amplitude.*sqrt(transmission).* ...
        exp(1i*(incident.phaseRad+device+k0*path.opl)));
    phaseOut=p.outgoing.opl(m)-delta*p.outgoing.q(m,:).'+ ...
        eta*p.outgoing.slope(m,:).'+sum((delta*p.outgoing.A(:,:,m)).*delta,2)/2- ...
        delta*p.outgoing.B(:,:,m)*eta.'+eta*p.outgoing.C(:,:,m)*eta.'/2;
    qLoaded=p.qInAtReference(m,:)+profile.gradientOverK0+delta*p.incomingH(:,:,m).';
    transmission=rp_fresnel(c.optics.nGlass,c.optics.nAir,qLoaded);
    quadraticField=quadraticField+sum(amplitude.*sqrt(transmission).* ...
        exp(1i*(incident.phaseRad+device+k0*phaseOut)));
    localPhase=toReceiver.opl-delta*toReceiver.q.'+sum((delta*toReceiver.A).*delta,2)/2;
    localField=localField+sum(amplitude.*sqrt(transmission).* ...
        exp(1i*(incident.phaseRad+device+k0*localPhase)));
    warpedPhase=toReceiver.opl-delta*toReceiver.q.'+sum((delta*p.outgoing.A(:,:,m)).*delta,2)/2;
    warpedField=warpedField+sum(amplitude.*sqrt(transmission).* ...
        exp(1i*(incident.phaseRad+device+k0*warpedPhase)));
end
exact=abs(det(toReceiver.B))/lambda^2*abs(exactField)^2;
quadratic=abs(det(p.outgoing.B(:,:,m)))/lambda^2*abs(quadraticField)^2;
receiverLocal=abs(det(toReceiver.B))/lambda^2*abs(localField)^2;
warped=abs(det(toReceiver.B))/lambda^2*abs(warpedField)^2;
end
