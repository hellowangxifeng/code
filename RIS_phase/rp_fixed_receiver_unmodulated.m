function result=rp_fixed_receiver_unmodulated(c,d,target)
%RP_FIXED_RECEIVER_UNMODULATED Fixed RX scatter plus an optimistic tail envelope.
% A remote receiver is outside the native field map. No extrapolated zero is
% passed off as a physical diffraction calculation: bound the unresolved tail.
rp_setup(); folder=fullfile(c.output,'fixed_receiver'); file=fullfile(folder,'unmodulated.mat');
p=d.bank{1}.phase; half=.50; n=c.spot.mapBins(1);
digest=java.security.MessageDigest.getInstance('SHA-256');
digest.update(unicode2native([jsonencode(c.study),jsonencode(target), ...
    rp_stage2_optical_signature(c),fileread([mfilename('fullpath') '.m']), ...
    fileread(which('rp_nonstationary_rect_bound'))],'UTF-8'));
signature=lower(reshape(dec2hex(typecast(digest.digest(),'uint8'),2).',1,[]));
if isfile(file)
    old=load(file,'result','signature');
    if strcmp(old.signature,signature), result=old.result; return; end
end
xe=linspace(target(1)-half,target(1)+half,n+1);
ye=linspace(target(2)-half,target(2)+half,n+1);
x=(xe(1:end-1)+xe(2:end))/2; y=(ye(1:end-1)+ye(2:end))/2;
r=rp_trace(c,d.geometry,d.illumination,p);
use=r.valid & ~d.illumination.ballistic;
ix=discretize(r.landing(use,1),xe); iy=discretize(r.landing(use,2),ye);
power=r.powerW(use); inside=isfinite(ix) & isfinite(iy); area=(xe(2)-xe(1))*(ye(2)-ye(1));
scatter=accumarray([iy(inside),ix(inside)],power(inside),[n n],@sum,0)/area;
wave=struct('x',x,'y',y,'irradianceWm2',scatter,'cellAreaM2',area);
samples=rp_link_samples(c); channelLower=rp_aperture_channel(c,d.geometry,wave,target,samples);
grid=[17 33]; envelopes=zeros(1,2); minOscillations=zeros(1,2);
for level=1:2
    [envelopes(level),minOscillations(level)]=tailEnvelope(c,d,p,target,half,grid(level));
end
assert(abs(envelopes(2)/envelopes(1)-1)<.03,'rp:TailBoundGrid','Refine tail-envelope quadrature.');
assert(min(minOscillations)>100,'rp:FarTail','Nonstationary far-tail approximation is inappropriate.');
% The integration-by-parts inequality is exact for the scalar aperture kernel.
% Its derivatives are evaluated numerically, with two-grid convergence and
% 20 percent headroom. This is an optimistic baseline, not an estimated tail.
irradianceEnvelope=1.2*max(envelopes);
wave.irradianceWm2=scatter+irradianceEnvelope;
channelUpper=rp_aperture_channel(c,d.geometry,wave,target,samples);
result=struct('channelLower',channelLower,'channelUpper',channelUpper,'wave',wave, ...
    'tailIrradianceEnvelopeWm2',irradianceEnvelope,'tailGridValuesWm2',envelopes, ...
    'tailMinimumOscillations',minOscillations,'scatterPacketsInWindow',sum(inside), ...
    'scatterPowerInWindowW',sum(power(inside)),'receiverCenterXY',target, ...
    'meaning','Optimistic ABER/OP lower envelope: empirical MC scatter plus bounded coherent ballistic tail');
save(file,'result','signature','-v7.3');
end

function [intensity,minOscillation]=tailEnvelope(c,d,p,target,half,n)
g=d.geometry; a=d.illumination; k0=2*pi/c.source.wavelength;
x=linspace(-c.array.apertureM(1)/2,c.array.apertureM(1)/2,n);
y=linspace(-c.array.apertureM(2)/2,c.array.apertureM(2)/2,n);
[X,Y]=meshgrid(x,y); offsets=[X(:),Y(:)];
targets=target+half*[-1 -1;-1 1;1 -1;1 1;0 0];
fieldBounds=zeros(numel(p.unitIds),1); minOscillation=Inf;
for m=1:numel(p.unitIds)
    inputPower=a.unitDirectPowerW(p.unitIds(m)); if inputPower<=0, continue; end
    center=a.unitCenters(p.unitIds(m),:); points=center+offsets;
    incident=rp_incident_field(c,g,points); atCenter=rp_incident_field(c,g,center);
    inputAmplitude=reshape(sqrt(incident.irradianceWm2),n,n);
    normalization=sqrt(inputPower/trapz(y,trapz(x,inputAmplitude.^2,2)));
    for j=1:size(targets,1)
        out=rp_paraxial(c,g,points,targets(j,:)); ref=rp_paraxial(c,g,center,targets(j,:));
        residual=incident.phaseRad-atCenter.phaseRad-k0*(offsets*atCenter.q.')+ ...
            k0*(out.opl-ref.opl+offsets*ref.q.');
        [rx,ry]=gradient(reshape(residual,n,n),x(2)-x(1),y(2)-y(1));
        difference=k0*(atCenter.q-ref.q); gx=rx+difference(1); gy=ry+difference(2);
        minOscillation=min(minOscillation,min([min(abs(gx),[],'all')*range(x), ...
            min(abs(gy),[],'all')*range(y)])/(2*pi));
        determinant=squeeze(out.B(1,1,:).*out.B(2,2,:)-out.B(1,2,:).*out.B(2,1,:));
        transmission=rp_fresnel(c.optics.nGlass,c.optics.nAir,out.q);
        amplitude=normalization*inputAmplitude.*reshape(sqrt(abs(determinant).*transmission)/c.source.wavelength,n,n);
        fieldBounds(m)=max(fieldBounds(m),rp_nonstationary_rect_bound(amplitude,gx,gy,x,y));
    end
end
intensity=sum(fieldBounds)^2;
end
