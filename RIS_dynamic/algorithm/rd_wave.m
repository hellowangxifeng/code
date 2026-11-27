function w = rd_wave(c,g,a,mc,p,r)
%RD_WAVE Receiver-plane irradiance: ballistic complex field + scattered packets.
%   w = rd_wave(c,g,a,mc,p,r) is RIS_phase/rp_wave (MODEL_AND_FORMULAS eqs.
%   54-66) for an arbitrary frozen board pose g, with two changes that make
%   30-50 m receiver planes affordable; neither changes the field model:
%
%   1. Beamlet spectral window. Tile m's Fresnel/CZT spectrum F_m(f) (eq. 56)
%      is concentrated where its stationary rays land, f = M_m l / lambda for
%      l inside the 14.4 mm aperture. The CZT window is that box widened by
%      c.wave.beamletGuardLobes diffraction lobes (1/d) on each side; outside
%      it lies only the sinc^2 edge-diffraction tail (<~2/(pi^2 N) per axis).
%      Pixels outside the landed window get exactly zero from this tile. The
%      dropped spectral power is booked in w.beamletTruncatedPowerW.
%   2. Neighbour coherence. The pixel-averaged |sum_m U_m|^2 (eq. 64) keeps
%      the exact sinc-weighted cross terms of every tile pair up to
%      c.wave.coherentNeighbours rings apart (rd_pixel_intensity). Tiles r
%      pitches apart arrive with phase-slope difference ~ k0 r p / L, so the
%      pixel sinc suppresses the omitted terms ~ 1/(pi r p dx/(lambda L)).
%      w.coherenceRingPowerW lists the map power of each kept ring so the
%      truncation is visible; aperture integration averages them further.
%
%   Scattered packets (non-ballistic, through the LC masks, refracted by the
%   loaded local gradient) add incoherently, exactly as in rp_wave.
half = max(p.receiverHalfWidthM,c.wave.minimumHalfWidthM);
nx = c.spot.mapBins(1);
ny = c.spot.mapBins(2);
xe = linspace(p.centerXY(1)-half(1),p.centerXY(1)+half(1),nx+1);
ye = linspace(p.centerXY(2)-half(2),p.centerXY(2)+half(2),ny+1);
x = (xe(1:end-1)+xe(2:end))/2;
y = (ye(1:end-1)+ye(2:end))/2;
[xx,yy] = meshgrid(x,y);
points = [xx(:),yy(:)];
pixelSize = [xe(2)-xe(1),ye(2)-ye(1)];
n = numel(p.unitIds);
unitPower = a.unitDirectPowerW(p.unitIds);
coherentUnits = find(unitPower>0 & p.outgoing.valid);
[~,order] = sort(unitPower(coherentUnits));
weak = coherentUnits(order(cumsum(unitPower(coherentUnits(order)))<= ...
    c.wave.maximumOmittedPowerFraction*a.activeDirectPowerW));
coherentUnits = setdiff(coherentUnits,weak,'stable');
omittedPower = sum(unitPower(weak));
K = numel(coherentUnits);
support = cell(1,K); values = cell(1,K); slopeX = cell(1,K); slopeY = cell(1,K);
pupilPower = zeros(n,1); windowPower = zeros(n,1); gridSamples = zeros(n,1);
incoherent = zeros(numel(xx),1);
for column = 1:K
    m = coherentUnits(column);
    [idx,field,sx,sy,integral,captured,N] = unitField(c,g,a,p,m,points,unitPower(m),x,y);
    support{column} = idx; values{column} = field; slopeX{column} = sx; slopeY{column} = sy;
    incoherent(idx) = incoherent(idx)+abs(field).^2;
    pupilPower(m) = integral; windowPower(m) = captured; gridSamples(m) = N;
    if mod(column,400)==0 || column==K
        fprintf('    propagated %d/%d coherent tiles\n',column,K);
    end
end
id = p.unitIds(coherentUnits);
tile = [mod(id-1,c.array.count(1))+1,floor((id-1)/c.array.count(1))+1];
[coherent,ringIrradiance] = rd_pixel_intensity(support,values,slopeX,slopeY,numel(xx), ...
    pixelSize,tile,c.wave.coherentNeighbours);
use = r.valid & ~a.ballistic;
ix = discretize(r.landing(use,1),xe);
iy = discretize(r.landing(use,2),ye);
power = r.powerW(use);
inside = isfinite(ix) & isfinite(iy);
halo = accumarray([iy(inside),ix(inside)],power(inside),[ny nx],@sum,0)/prod(pixelSize);
intensity = reshape(coherent,ny,nx)+halo;
w = struct('x',x,'y',y,'irradianceWm2',intensity,'ballisticPixelMeanWm2',reshape(coherent,ny,nx), ...
    'ballisticIncoherentWm2',reshape(incoherent,ny,nx),'scatteredIrradianceWm2',halo, ...
    'cellAreaM2',prod(pixelSize),'pupilPowerW',sum(pupilPower), ...
    'inputDirectPowerW',a.activeDirectPowerW, ...
    'interfaceLossW',a.activeDirectPowerW-sum(pupilPower)-omittedPower, ...
    'mapPowerW',sum(intensity,'all')*prod(pixelSize), ...
    'diagonalMapPowerW',sum(incoherent)*prod(pixelSize), ...
    'scatteredMapPowerW',sum(halo,'all')*prod(pixelSize), ...
    'omittedBallisticPowerW',omittedPower, ...
    'beamletTruncatedPowerW',sum(pupilPower-windowPower), ...
    'coherenceRingPowerW',ringIrradiance*prod(pixelSize), ...
    'coherentNeighbours',c.wave.coherentNeighbours,'apertureSamples',gridSamples, ...
    'pixelAveraged',true,'unresolvedQuadraturePowerFraction',0);
w.fieldConvention = 'flux phasor sqrt(W/m^2); pixel-averaged coherent sum';
w.scatterFieldModel = 'zero-mean incoherent ensemble; HG transport fixes power, not a unique phasor';
w.propagationModel = ['RIS_phase rp_wave unit integrals (frozen aperture curvature, exact ' ...
    'receiver endpoints), beamlet CZT windows, neighbour-ring pixel coherence'];
w.ballisticBoardMCPowerW = sum(mc.powerW(mc.ballistic));
end

function [idx,field,slopeX,slopeY,pupilPower,windowPower,N] = unitField(c,g,a,p,m,points,power,pixelX,pixelY)
% rp_wave/unitField for the linear phase model, restricted to the beamlet.
lambda = c.source.wavelength;
k0 = 2*pi/lambda;
idx = zeros(0,1); field = complex(zeros(0,1)); slopeX = zeros(0,1); slopeY = zeros(0,1);
pupilPower = 0; windowPower = 0;
A = p.outgoing.A(:,:,m); B = p.outgoing.B(:,:,m);
H = p.phaseCurvature(:,:,m); M = A+p.incomingH(:,:,m)+H;
[rotation,~,~] = svd(B);
rotatedCurvature = rotation.'*M*rotation;
extent = abs(rotation.')*c.array.apertureM.';
pixelSize = [pixelX(2)-pixelX(1),pixelY(2)-pixelY(1)];
tileCenter = a.unitCenters(p.unitIds(m),:);
ref = p.referenceXY(m,:);
qin = p.qInAtReference(m,:);
offsetReference = tileCenter-ref;
qCenter = qin+p.gradientOverK0(m,:)+offsetReference*M;
% Stationary-ray frequency box of the tile, widened by the guard lobes.
corners = [-1 -1;1 -1;1 1;-1 1].*c.array.apertureM/2;
cornerFrequency = (corners*M/lambda)*rotation;
guard = c.wave.beamletGuardLobes./extent.';
beam = [min(cornerFrequency,[],1)-guard;max(cornerFrequency,[],1)+guard];
% Land the box outline to find the receiver pixels this tile can reach.
[bx,by] = meshgrid(linspace(beam(1,1),beam(2,1),5),linspace(beam(1,2),beam(2,2),5));
outline = [bx(:),by(:)];
[landed,reach] = rp_carrier_landing(c,g,ref,lambda*outline*rotation.'+qCenter);
if ~any(reach), N = 0; return; end
low = min(landed(reach,:),[],1)-2*pixelSize;
high = max(landed(reach,:),[],1)+2*pixelSize;
selected = find(all(points>=low & points<=high,2));
if isempty(selected), N = 0; return; end
receiver = rp_paraxial(c,g,ref,points(selected,:));
box = [low;high(1) low(2);low(1) high(2);high];
boundary = rp_paraxial(c,g,ref,box);
jacobians = cat(3,receiver.B(:,:,receiver.valid),boundary.B(:,:,boundary.valid));
if isempty(jacobians), jacobians = B; end
inverseMap = pagemldivide(jacobians,repmat(rotation,1,1,size(jacobians,3)));
scaled = abs(lambda*inverseMap)./reshape(pixelSize,2,1,1);
stretch = reshape(max(scaled,[],[1 3]),2,1);
requiredPadding = 2.^nextpow2(c.wave.receiverSamplingMargin*stretch./extent);
padding = max(c.wave.fftPadding,requiredPadding);
required = 1.1*extent.*(abs(rotatedCurvature)*extent)/lambda;
oversampling = 1;
if isfield(c.wave,'apertureOversampling'), oversampling = c.wave.apertureOversampling; end
sampleCounts = oversampling*2.^nextpow2(max(c.wave.apertureSamples,required));
N = max(sampleCounts);
assert(N<=c.wave.maxApertureSamples,'rd:WaveSampling', ...
    'Phase curvature needs %d aperture samples; increase wave.maxApertureSamples.',N);
dx = extent(1)/sampleCounts(1); dy = extent(2)/sampleCounts(2);
u = ((1:sampleCounts(1))-(sampleCounts(1)+1)/2)*dx;
v = ((1:sampleCounts(2))-(sampleCounts(2)+1)/2)*dy;
[ux,uy] = meshgrid(u,v);
localX = rotation(1,1)*ux+rotation(1,2)*uy;
localY = rotation(2,1)*ux+rotation(2,2)*uy;
xi = localX+tileCenter(1)-ref(1); yi = localY+tileCenter(2)-ref(2);
% Smooth envelope and residual Gaussian phase on a 25 x 25 tile grid (rp_wave).
sourceAxis = linspace(-0.5,0.5,25);
[sx,sy] = meshgrid(sourceAxis*c.array.apertureM(1),sourceAxis*c.array.apertureM(2));
input = rp_incident_field(c,g,[sx(:),sy(:)]+tileCenter);
atReference = rp_incident_field(c,g,ref);
delta = [sx(:),sy(:)]+tileCenter-ref;
localPhase = atReference.phaseRad+k0*(delta*qin.'+sum((delta*p.incomingH(:,:,m)).*delta,2)/2);
phaseResidual = reshape(input.phaseRad-localPhase,size(sx));
amplitude = interp2(sx,sy,reshape(sqrt(input.irradianceWm2),size(sx)),localX,localY,'cubic',0);
residualPhase = interp2(sx,sy,phaseResidual,localX,localY,'cubic',0);
illuminated = abs(localX)<=c.array.apertureM(1)/2 & abs(localY)<=c.array.apertureM(2)/2;
amplitude(~illuminated) = 0;
if ~any(amplitude,'all'), amplitude = double(illuminated); end
integral = sum(amplitude.^2,'all')*dx*dy;
if integral<=0, return; end
amplitude = amplitude*sqrt(power/integral);
profile = rp_phase_profile(c,p,m,[xi(:),yi(:)]);
incomingGradient = [xi(:),yi(:)]*p.incomingH(:,:,m);
outX = reshape(qin(1)+incomingGradient(:,1)+profile.gradientOverK0(:,1),size(xi));
outY = reshape(qin(2)+incomingGradient(:,2)+profile.gradientOverK0(:,2),size(yi));
transmission = reshape(rp_fresnel(c.optics.nGlass,c.optics.nAir,[outX(:),outY(:)]),size(ux));
amplitude = amplitude.*sqrt(transmission);
pupilPower = sum(amplitude.^2,'all')*dx*dy;
if ~any(receiver.valid), return; end
chirp = amplitude.*exp(1i*(k0/2*(M(1,1)*localX.^2+2*M(1,2)*localX.*localY+M(2,2)*localY.^2)+residualPhase));
frequency = ((receiver.q-qCenter)/lambda)*rotation;
frequencyBoundary = ((boundary.q(boundary.valid,:)-qCenter)/lambda)*rotation;
frequencyRange = [frequency(receiver.valid,:);frequencyBoundary];
guardQ = max(abs(jacobians),[],3)*(2*pixelSize.');
guardFrequency = abs(rotation.')*guardQ/lambda;
window = [max(min(frequencyRange,[],1)-guardFrequency.',beam(1,:)); ...
    min(max(frequencyRange,[],1)+guardFrequency.',beam(2,:))];
if any(window(1,:)>window(2,:)), return; end
[F,frequencyX,frequencyY,df] = rp_fourier_window(chirp,[dx dy],padding.',window, ...
    c.wave.maximumFrequencySamples);
if isempty(F), return; end
[fgx,fgy] = meshgrid(frequencyX,frequencyY);
spectralPower = abs(F).^2*prod(df);
windowPower = min(pupilPower,sum(spectralPower,'all'));
deposited = rp_deposit_power([fgx(:),fgy(:)]*rotation.',spectralPower(:), ...
    c,g,ref,qCenter,lambda,pixelX,pixelY);
pixelIntensity = deposited(selected);
[phaseBasis,phaseValues] = eig((rotatedCurvature+rotatedCurvature.')/2);
phaseValues = diag(phaseValues);
phaseExtent = abs(phaseBasis.')*extent;
strong = abs(phaseValues).*phaseExtent.^2/lambda>1;
inverseValues = zeros(2,1); inverseValues(strong) = 1./phaseValues(strong);
inverse = phaseBasis*diag(inverseValues)*phaseBasis.';
if any(strong)
    F = F.*exp(1i*pi*lambda*(inverse(1,1)*fgx.^2+2*inverse(1,2)*fgx.*fgy+inverse(2,2)*fgy.^2));
end
value = interp2(frequencyX,frequencyY,F,frequency(:,1),frequency(:,2),'cubic',0);
dfx = (frequencyX(2)-frequencyX(1))/4; dfy = (frequencyY(2)-frequencyY(1))/4;
plusX = interp2(frequencyX,frequencyY,F,frequency(:,1)+dfx,frequency(:,2),'cubic',0);
minusX = interp2(frequencyX,frequencyY,F,frequency(:,1)-dfx,frequency(:,2),'cubic',0);
plusY = interp2(frequencyX,frequencyY,F,frequency(:,1),frequency(:,2)+dfy,'cubic',0);
minusY = interp2(frequencyX,frequencyY,F,frequency(:,1),frequency(:,2)-dfy,'cubic',0);
frequencySlope = [angle(plusX.*conj(minusX))/(2*dfx),angle(plusY.*conj(minusY))/(2*dfy)] ...
    -2*pi*lambda*frequency*inverse;
value = value.*exp(-1i*pi*lambda*sum((frequency*inverse).*frequency,2));
carrier = k0*(offsetReference*M*offsetReference.'/2+ ...
    (qin+p.gradientOverK0(m,:)-receiver.q)*offsetReference.');
% Source phase + programmed tile piston + propagation phase all enter (eq. 58).
carrier = carrier+mod(atReference.phaseRad+p.phaseAtReferenceRad(m)+k0*receiver.opl,2*pi);
local = sqrt(pixelIntensity).*exp(1i*(angle(value)+carrier));
local(~receiver.valid) = 0;
transport = frequencySlope*rotation.'/lambda-k0*offsetReference;
correction = pagemtimes(reshape(transport.',1,2,[]),receiver.B);
phaseSlope = k0*receiver.slope+reshape(correction,2,[]).';
keep = abs(local).^2>=max(abs(local).^2)*c.wave.supportRelativeFloor & abs(local)>0;
idx = selected(keep);
field = local(keep);
slopeX = phaseSlope(keep,1);
slopeY = phaseSlope(keep,2);
end
