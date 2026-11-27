function window = rp_receiver_window(c,g,a,p,r)
%RP_RECEIVER_WINDOW Deterministic power envelope of the loaded phase, in world XY.
% A finite main-beam window is independent of RX jitter and a 100 m search disk.
% Diffuse scatter is retained wherever it falls inside this explicitly saved window.
arguments
    c (1,1) struct
    g (1,1) struct
    a (1,1) struct
    p (1,1) struct
    r (1,1) struct
end
count=c.wave.windowSamples;
sampleAxis=((1:count)-0.5)/count-0.5;
[x,y]=meshgrid(sampleAxis*c.array.apertureM(1),sampleAxis*c.array.apertureM(2));
n=numel(p.unitIds); samples=numel(x);
xy=repelem(a.unitCenters(p.unitIds,:),samples,1)+repmat([x(:),y(:)],n,1);
incoming=rp_layer_path(xy-g.sourceLocal(1:2).',[-g.sourceLocal(3),g.hIn], ...
    [c.optics.nWater,c.optics.nGlass],c.numerics.bisectionIterations);
packets=struct('lcXY',xy,'qIn',incoming.q,'sampleUnitIndex',repelem((1:n).',samples), ...
    'powerW',rp_bsf_ballistic(c,g,xy,g.hIn)*prod(c.array.apertureM)/samples, ...
    'unitId',repelem(p.unitIds,samples));
traced=rp_trace(c,g,packets,p);
use=traced.valid & traced.powerW>0 & p.outgoing.valid(packets.sampleUnitIndex);
landing=traced.landing(use,1:2); weights=traced.powerW(use);
source='deterministic BSF ballistic aperture quadrature';
if isempty(weights)
    use=r.valid & r.powerW>0;
    landing=r.landing(use,1:2); weights=r.powerW(use);
    source='MC transmitted packets (no ballistic carrier)';
end
diffraction=[0 0];
bright=p.outgoing.valid & a.unitDirectPowerW(p.unitIds)> ...
    c.wave.maximumOmittedPowerFraction*a.activeDirectPowerW;
for m=find(bright).'
    % Map the tile's first Fourier zeros to world XY through the actual carrier.
    width=abs(p.outgoing.B(:,:,m)\diag(c.source.wavelength./c.array.apertureM));
    diffraction=max(diffraction,2*sum(width,2).');
end
if isempty(weights)
    center=g.pRis(1:2).'; half=repmat(c.wave.minimumHalfWidthM,1,2);
    captured=0;
else
    tail=(1-c.wave.windowPowerFraction)/4;
    low=[weightedQuantile(landing(:,1),weights,tail),weightedQuantile(landing(:,2),weights,tail)];
    high=[weightedQuantile(landing(:,1),weights,1-tail),weightedQuantile(landing(:,2),weights,1-tail)];
    center=(low+high)/2;
    half=max((high-low)/2*(1+c.wave.windowPadding)+diffraction,c.wave.minimumHalfWidthM);
    captured=sum(weights(all(abs(landing-center)<=half,2)))/sum(weights);
end
window=struct('centerXY',center,'halfWidthM',half,'rayCapturedPowerFraction',captured, ...
    'rayTransmittedPowerW',sum(weights),'source',source,'diffractionMarginM',diffraction, ...
    'boundsXY',[center-half;center+half]);
end

function value=weightedQuantile(x,weights,fraction)
[x,order]=sort(x); cumulative=cumsum(weights(order));
index=find(cumulative>=fraction*cumulative(end),1);
value=x(index);
end
