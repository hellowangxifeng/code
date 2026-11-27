function w = rp_wave(c,g,a,mc,p,r)
%RP_WAVE Per-unit canonical Fresnel integrals and detector-pixel coherent averaging.
arguments
    c (1,1) struct
    g (1,1) struct
    a (1,1) struct
    mc (1,1) struct
    p (1,1) struct
    r (1,1) struct
end
half=max(p.receiverHalfWidthM,c.wave.minimumHalfWidthM);
nx=c.spot.mapBins(1); ny=c.spot.mapBins(2);
xe=linspace(p.centerXY(1)-half(1),p.centerXY(1)+half(1),nx+1);
ye=linspace(p.centerXY(2)-half(2),p.centerXY(2)+half(2),ny+1);
x=(xe(1:end-1)+xe(2:end))/2; y=(ye(1:end-1)+ye(2:end))/2;
[xx,yy]=meshgrid(x,y); points=[xx(:),yy(:)]; pixelSize=[xe(2)-xe(1),ye(2)-ye(1)];
n=numel(p.unitIds);
unitPower=a.unitDirectPowerW(p.unitIds);
coherentUnits=find(unitPower>0 & p.outgoing.valid);
[~,order]=sort(unitPower(coherentUnits));
weak=coherentUnits(order(cumsum(unitPower(coherentUnits(order)))<= ...
    c.wave.maximumOmittedPowerFraction*a.activeDirectPowerW));
coherentUnits=setdiff(coherentUnits,weak,'stable');
omittedPower=sum(unitPower(weak));
fields=complex(zeros(numel(xx),numel(coherentUnits)));
slopeX=zeros(size(fields)); slopeY=slopeX; pupilPower=zeros(n,1); gridSamples=zeros(n,1);
paddingUsed=zeros(2,n); quadratureCapped=false(n,1); fourierSteps=zeros(2,n); gridShape=zeros(2,n);
for column=1:numel(coherentUnits)
    m=coherentUnits(column);
    power=a.unitDirectPowerW(p.unitIds(m));
    if power<=0 || ~p.outgoing.valid(m), continue; end
    [field,dxPhase,dyPhase,integral,N,padding,capped,step,shape]=unitField(c,g,a,p,m,points,power,x,y);
    fields(:,column)=field; slopeX(:,column)=dxPhase; slopeY(:,column)=dyPhase;
    pupilPower(m)=integral; gridSamples(m)=N;
    paddingUsed(:,m)=padding; quadratureCapped(m)=capped; fourierSteps(:,m)=step; gridShape(:,m)=shape;
    if mod(column,100)==0 || column==numel(coherentUnits)
        fprintf('    propagated %d/%d coherent units\n',column,numel(coherentUnits));
    end
end
incoherent=sum(abs(fields).^2,2);
unitMapPower=zeros(n,1);
unitMapPower(coherentUnits)=sum(abs(fields).^2,1).'*prod(pixelSize);
if c.wave.pixelCoherence
    % At most 0.02% of the irradiance threshold at ANY pixel, by a sinc
    % envelope bound, rather than an unquantified truncation of interference.
    tolerance=2e-4*c.spot.minIrradianceWm2;
    [coherent,coherenceError]=rp_pixel_intensity(fields,slopeX,slopeY,pixelSize,tolerance);
else
    coherent=abs(sum(fields,2)).^2; coherenceError=0;
end
% Scattered packets have no measured optical path phase: add their MC radiance.
use=r.valid & ~a.ballistic;
ix=discretize(r.landing(use,1),xe); iy=discretize(r.landing(use,2),ye);
power=r.powerW(use); inside=isfinite(ix) & isfinite(iy);
halo=accumarray([iy(inside),ix(inside)],power(inside),[ny nx],@sum,0)/prod(pixelSize);
intensity=reshape(coherent,ny,nx)+halo;
w=struct('x',x,'y',y,'irradianceWm2',intensity,'ballisticPixelMeanWm2',reshape(coherent,ny,nx), ...
    'ballisticIncoherentWm2',reshape(incoherent,ny,nx),'scatteredIrradianceWm2',halo, ...
    'cellAreaM2',prod(pixelSize),'pupilPowerW',sum(pupilPower), ...
    'inputDirectPowerW',a.activeDirectPowerW, ...
    'interfaceLossW',a.activeDirectPowerW-sum(pupilPower)-omittedPower, ...
    'mapPowerW',sum(intensity,'all')*prod(pixelSize), ...
    'diagonalMapPowerW',sum(incoherent)*prod(pixelSize), ...
    'scatteredMapPowerW',sum(halo,'all')*prod(pixelSize), ...
    'apertureSamples',gridSamples,'pixelAveraged',c.wave.pixelCoherence);
w.unitMapPowerW=unitMapPower;
w.ballisticPointField=reshape(sum(fields,2),ny,nx);
w.fieldConvention='flux phasor sqrt(W/m^2); point-field magnitude is not a pixel mean';
w.scatterFieldModel='zero-mean incoherent ensemble; HG transport fixes power, not a unique phasor';
w.omittedBallisticPowerW=omittedPower;
w.coherenceErrorBoundWm2=coherenceError;
w.ballisticBoardMCPowerW=sum(mc.powerW(mc.ballistic));
w.fftPaddingUsed=paddingUsed;
w.fourierStepPixels=fourierSteps;
w.apertureGridSamples=gridShape;
w.fourierMethod='windowed CZT, exact sampled Fourier integral';
w.propagationModel='exact loaded phase profile; frozen outgoing aperture curvature; exact receiver endpoints';
w.unresolvedQuadraturePowerFraction=sum(unitMapPower(quadratureCapped))/max(sum(unitMapPower),realmin);
end

function [field,slopeX,slopeY,pupilPower,N,padding,capped,fourierStep,sampleCounts]=unitField(c,g,a,p,m,points,power,pixelX,pixelY)
lambda=c.source.wavelength; k0=2*pi/lambda;
A=p.outgoing.A(:,:,m); B=p.outgoing.B(:,:,m);
H=p.phaseCurvature(:,:,m); M=A+p.incomingH(:,:,m)+H;
[rotation,~,~]=svd(B);
rotatedCurvature=rotation.'*M*rotation;
extent=abs(rotation.')*c.array.apertureM.';
pixelSize=[pixelX(2)-pixelX(1),pixelY(2)-pixelY(1)];
tileCenter=a.unitCenters(p.unitIds(m),:); ref=p.referenceXY(m,:);
receiver=rp_paraxial(c,g,ref,points);
corners=[pixelX(1)-pixelSize(1),pixelY(1)-pixelSize(2); ...
    pixelX(end)+pixelSize(1),pixelY(1)-pixelSize(2); ...
    pixelX(1)-pixelSize(1),pixelY(end)+pixelSize(2); ...
    pixelX(end)+pixelSize(1),pixelY(end)+pixelSize(2)];
boundary=rp_paraxial(c,g,ref,corners);
jacobians=cat(3,receiver.B(:,:,receiver.valid),boundary.B(:,:,boundary.valid));
if isempty(jacobians), jacobians=B; end
inverseMap=pagemldivide(jacobians,repmat(rotation,1,1,size(jacobians,3)));
scaled=abs(lambda*inverseMap)./reshape(pixelSize,2,1,1);
stretch=reshape(max(scaled,[],[1 3]),2,1);
requiredPadding=2.^nextpow2(c.wave.receiverSamplingMargin*stretch./extent);
padding=c.wave.fftPadding;
if ~isscalar(padding), padding=padding(:,m); end
capped=false;
padding=max(padding,requiredPadding);
fourierStep=stretch./(padding.*extent);
required=1.1*extent.*(abs(rotatedCurvature)*extent)/lambda;
if isfield(p,'phaseModel') && strcmp(p.phaseModel,'tanh')
    offsets=[-1 -1;1 -1;1 1;-1 1].*c.array.apertureM/2;
    profile=rp_phase_profile(c,p,m,offsets+tileCenter-ref);
    centerGradient=p.gradientOverK0(m,:)+(tileCenter-ref)*H;
    low=min(profile.gradientOverK0,[],1)-centerGradient;
    high=max(profile.gradientOverK0,[],1)-centerGradient;
    smoothBound=abs(A+p.incomingH(:,:,m))*c.array.apertureM.'/2;
    gradientBound=abs(rotation.')*(max(abs([low;high]),[],1).'+smoothBound);
    required=max(required,2.2*extent.*gradientBound/lambda);
end
oversampling=1;
if isfield(c.wave,'apertureOversampling'), oversampling=c.wave.apertureOversampling; end
sampleCounts=oversampling*2.^nextpow2(max(c.wave.apertureSamples,required)); N=max(sampleCounts);
assert(N<=c.wave.maxApertureSamples,'rp:WaveSampling', ...
    'Phase curvature needs %d aperture samples; increase wave.maxApertureSamples.',N);
dx=extent(1)/sampleCounts(1); dy=extent(2)/sampleCounts(2);
u=((1:sampleCounts(1))-(sampleCounts(1)+1)/2)*dx;
v=((1:sampleCounts(2))-(sampleCounts(2)+1)/2)*dy;
[ux,uy]=meshgrid(u,v);
localX=rotation(1,1)*ux+rotation(1,2)*uy;
localY=rotation(2,1)*ux+rotation(2,2)*uy;
xi=localX+tileCenter(1)-ref(1); yi=localY+tileCenter(2)-ref(2);
qin=p.qInAtReference(m,:);
% Interpolate a smooth envelope and residual phase, handling the carrier
% analytically. An optical-wavelength grid across the whole panel is unnecessary.
sourceAxis=linspace(-0.5,0.5,25);
[sx,sy]=meshgrid(sourceAxis*c.array.apertureM(1),sourceAxis*c.array.apertureM(2));
input=rp_incident_field(c,g,[sx(:),sy(:)]+tileCenter);
atReference=rp_incident_field(c,g,ref);
delta=[sx(:),sy(:)]+tileCenter-ref;
localPhase=atReference.phaseRad+k0*(delta*qin.'+sum((delta*p.incomingH(:,:,m)).*delta,2)/2);
phaseResidual=reshape(input.phaseRad-localPhase,size(sx));
amplitude=interp2(sx,sy,reshape(sqrt(input.irradianceWm2),size(sx)),localX,localY,'cubic',0);
residualPhase=interp2(sx,sy,phaseResidual,localX,localY,'cubic',0);
illuminated=abs(localX)<=c.array.apertureM(1)/2 & abs(localY)<=c.array.apertureM(2)/2;
amplitude(~illuminated)=0;
if ~any(amplitude,'all'), amplitude=double(illuminated); end
integral=sum(amplitude.^2,'all')*dx*dy;
if integral<=0
    field=complex(zeros(size(points,1),1)); slopeX=zeros(size(field)); slopeY=slopeX;
    pupilPower=0; return
end
amplitude=amplitude*sqrt(power/integral);
profile=rp_phase_profile(c,p,m,[xi(:),yi(:)]);
incomingGradient=[xi(:),yi(:)]*p.incomingH(:,:,m);
outX=reshape(qin(1)+incomingGradient(:,1)+profile.gradientOverK0(:,1),size(xi));
outY=reshape(qin(2)+incomingGradient(:,2)+profile.gradientOverK0(:,2),size(yi));
transmission=reshape(rp_fresnel(c.optics.nGlass,c.optics.nAir,[outX(:),outY(:)]),size(ux));
amplitude=amplitude.*sqrt(transmission);
pupilPower=sum(amplitude.^2,'all')*dx*dy;
if ~any(receiver.valid)
    field=complex(zeros(size(points,1),1)); slopeX=zeros(size(field)); slopeY=slopeX;
    return
end
% Handle the reference-point offset as an exact frequency shift, keeping the
% sampled quadratic chirp centered on the physical tile (no extra carrier alias).
chirp=amplitude.*exp(1i*(k0/2*(M(1,1)*localX.^2+2*M(1,2)*localX.*localY+M(2,2)*localY.^2)+residualPhase));
if isfield(p,'phaseModel') && strcmp(p.phaseModel,'tanh')
    quadratic=[xi(:),yi(:)]*p.gradientOverK0(m,:).'+ ...
        sum(([xi(:),yi(:)]*H).*[xi(:),yi(:)],2)/2;
    chirp=chirp.*exp(1i*k0*reshape(profile.phaseOverK0-quadratic,size(xi)));
end
offsetReference=tileCenter-ref;
qCenter=qin+p.gradientOverK0(m,:)+offsetReference*M;
frequency=((receiver.q-qCenter)/lambda)*rotation;
% The singular-vector basis of B avoids a long, thin, rotated frequency
% footprint becoming a huge axis-aligned rectangle at near-grazing emergence.
% Exact endpoint directions remove the large receiver-window Taylor error.
% Only the aperture-coordinate curvature is frozen at the unit reference.
frequencyBoundary=((boundary.q(boundary.valid,:)-qCenter)/lambda)*rotation;
frequencyRange=[frequency(receiver.valid,:);frequencyBoundary];
guardQ=max(abs(jacobians),[],3)*(2*pixelSize.');
guardFrequency=abs(rotation.')*guardQ/lambda;
maximumPoints=32000000;
if isfield(c.wave,'maximumFrequencySamples'), maximumPoints=c.wave.maximumFrequencySamples; end
[F,frequencyX,frequencyY,df]=rp_fourier_window(chirp,[dx dy],padding.', ...
    [min(frequencyRange,[],1)-guardFrequency.';max(frequencyRange,[],1)+guardFrequency.'],maximumPoints);
if isempty(F)
    field=complex(zeros(size(points,1),1)); slopeX=zeros(size(field)); slopeY=slopeX;
    return
end
[fgx,fgy]=meshgrid(frequencyX,frequencyY);
spectralPower=abs(F).^2*prod(df);
pixelIntensity=rp_deposit_power([fgx(:),fgy(:)]*rotation.',spectralPower(:), ...
    c,g,ref,qCenter,lambda,pixelX,pixelY);
% Remove the known stationary-phase chirp before interpolation; otherwise a
% coarse frequency grid attenuates a rapidly rotating complex Fresnel field.
[phaseBasis,phaseValues]=eig((rotatedCurvature+rotatedCurvature.')/2);
phaseValues=diag(phaseValues);
phaseExtent=abs(phaseBasis.')*extent;
strong=abs(phaseValues).*phaseExtent.^2/lambda>1;
inverseValues=zeros(2,1); inverseValues(strong)=1./phaseValues(strong);
inverse=phaseBasis*diag(inverseValues)*phaseBasis.';
% A near-zero curvature has a smooth Fourier phase. Dividing by that tiny
% curvature invents a fast chirp and aliases relative phase between tiles.
if any(strong)
    F=F.*exp(1i*pi*lambda*(inverse(1,1)*fgx.^2+ ...
        2*inverse(1,2)*fgx.*fgy+inverse(2,2)*fgy.^2));
end
value=interp2(frequencyX,frequencyY,F,frequency(:,1),frequency(:,2),'cubic',0);
dfx=(frequencyX(2)-frequencyX(1))/4; dfy=(frequencyY(2)-frequencyY(1))/4;
plusX=interp2(frequencyX,frequencyY,F,frequency(:,1)+dfx,frequency(:,2),'cubic',0);
minusX=interp2(frequencyX,frequencyY,F,frequency(:,1)-dfx,frequency(:,2),'cubic',0);
plusY=interp2(frequencyX,frequencyY,F,frequency(:,1),frequency(:,2)+dfy,'cubic',0);
minusY=interp2(frequencyX,frequencyY,F,frequency(:,1),frequency(:,2)-dfy,'cubic',0);
frequencySlope=[angle(plusX.*conj(minusX))/(2*dfx),angle(plusY.*conj(minusY))/(2*dfy)] ...
    -2*pi*lambda*frequency*inverse;
remodulation=exp(-1i*pi*lambda*sum((frequency*inverse).*frequency,2));
value=value.*remodulation;
carrier=k0*(offsetReference*M*offsetReference.'/2+ ...
    (qin+p.gradientOverK0(m,:)-receiver.q)*offsetReference.');
% Source phase + programmed tile piston + propagation phase must all enter.
carrier=carrier+mod(atReference.phaseRad+p.phaseAtReferenceRad(m)+k0*receiver.opl,2*pi);
% Pixel powers come from a conservative Fourier-area pushforward (Parseval),
% while complex interpolation supplies phase and local interference carriers.
field=sqrt(pixelIntensity).*exp(1i*(angle(value)+carrier));
field(~receiver.valid)=0;
transport=frequencySlope*rotation.'/lambda-k0*offsetReference;
correction=pagemtimes(reshape(transport.',1,2,[]),receiver.B);
phaseSlope=k0*receiver.slope+reshape(correction,2,[]).';
slopeX=phaseSlope(:,1); slopeY=phaseSlope(:,2);
floorValue=max(abs(field).^2)*c.wave.supportRelativeFloor;
negligible=abs(field).^2<floorValue;
field(negligible)=0; slopeX(negligible)=0; slopeY(negligible)=0;
end

